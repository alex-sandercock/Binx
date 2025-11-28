pub mod io;
pub mod math;
pub mod model;

use ndarray::Array1;
use io::LocusData;
use rayon::prelude::*;
use std::io::{BufWriter, Write};

// Re-export FitMode for CLI usage
pub use model::FitMode;

/// Represents the result of a genotyping run for a single locus
#[derive(Debug, Clone)]
pub struct GenotypeResult {
    pub locus_id: String,
    pub genotype_probs: Vec<Vec<f64>>, // [Sample][Genotype]
    pub best_genotypes: Vec<usize>,    // [Sample]
    pub bias: f64,
    pub seq_error: f64,
    pub overdispersion: f64,           // Added rho parameter
    pub model_mu: f64,
    pub model_sigma: f64,
    pub final_log_lik: f64,
}

/// Input source for dosage estimation.
pub enum InputSource {
    /// Legacy two-line CSV (alternating ref/total rows).
    TwoLineCsv(String),
    /// Separate ref/total matrices (markers in rows, samples in columns).
    RefTotalMatrices { ref_path: String, total_path: String },
    /// VCF input (gzipped or plain), optionally streamed in chunks.
    Vcf { path: String, chunk_size: Option<usize> },
}

pub fn run_norm_model(
    ref_counts: &Array1<u32>,
    total_counts: &Array1<u32>,
    ploidy: usize,
    mode: FitMode,
) -> anyhow::Result<GenotypeResult> {
    model::fit_norm_with_mode(ref_counts, total_counts, ploidy, mode)
}

/// Top-level entry point used by the CLI for genotype dosage estimation.
pub fn run_dosage(
    input: InputSource,
    ploidy: usize,
    mode: FitMode,
    verbose: bool,
    threads: Option<usize>,
) -> anyhow::Result<()> {
    // Configure rayon thread pool if requested.
    if let Some(n) = threads {
        if n > 0 {
            let _ = rayon::ThreadPoolBuilder::new().num_threads(n).build_global();
        }
    }

    if verbose {
        println!("Starting dosage estimation (mode: {:?})", mode);
    }
    println!("Locus\tBias\tRho\tMu\tSigma\tLogLik\tGenotypes");

    let process_locus = |locus_id: String, ref_counts: Array1<u32>, total_counts: Array1<u32>| {
        match run_norm_model(&ref_counts, &total_counts, ploidy, mode) {
            Ok(mut res) => {
                res.locus_id = locus_id;
                let geno_str: Vec<String> = res.best_genotypes.iter().map(|g| g.to_string()).collect();
                Some(format!("{}\t{:.3}\t{:.4}\t{:.3}\t{:.3}\t{:.2}\t{}",
                    res.locus_id,
                    res.bias,
                    res.overdispersion,
                    res.model_mu,
                    res.model_sigma,
                    res.final_log_lik,
                    geno_str.join(",")
                ))
            },
            Err(e) => {
                eprintln!("Error processing locus: {}", e);
                None
            },
        }
    };

    match input {
        InputSource::TwoLineCsv(csv_file) => {
            if verbose {
                println!("Parsing {}...", csv_file);
            }
            let loci = io::parse_two_line_csv(&csv_file)
                .map_err(|e| anyhow::anyhow!("Failed to parse CSV file: {}", e))?;
            if verbose {
                println!("Found {} loci from CSV", loci.len());
            }
            let mut out = BufWriter::new(std::io::stdout().lock());
            for locus in loci {
                if let Some(line) = process_locus(locus.id, locus.ref_counts, locus.total_counts) {
                    writeln!(out, "{}", line)?;
                }
            }
        }
        InputSource::RefTotalMatrices { ref_path, total_path } => {
            if verbose {
                println!("Parsing ref matrix {} and total matrix {}...", ref_path, total_path);
            }
            let matrices = io::parse_ref_total_matrices(&ref_path, &total_path)
                .map_err(|e| anyhow::anyhow!("Failed to parse ref/total matrices: {}", e))?;
            if verbose {
                println!("Found {} loci from matrices", matrices.marker_ids.len());
            }
            let mut out = BufWriter::new(std::io::stdout().lock());
            for (row_idx, locus_id) in matrices.marker_ids.iter().enumerate() {
                let ref_counts = matrices.ref_counts.row(row_idx).to_owned();
                let total_counts = matrices.total_counts.row(row_idx).to_owned();
                if let Some(line) = process_locus(locus_id.clone(), ref_counts, total_counts) {
                    writeln!(out, "{}", line)?;
                }
            }
        }
        InputSource::Vcf { path, chunk_size } => {
            if verbose {
                println!("Streaming VCF {}{}", path, chunk_size.map(|c| format!(" in chunks of {}", c)).unwrap_or_default());
            }
            let mut buffer: Vec<LocusData> = Vec::new();
            let chunk_size = chunk_size.unwrap_or(0);
            let mut out = BufWriter::new(std::io::stdout().lock());
            io::stream_vcf_records(&path, |rec| {
                if chunk_size == 0 {
                    if let Some(line) = process_locus(rec.id, rec.ref_counts, rec.total_counts) {
                        // Best-effort: ignore I/O errors here; they will be caught on flush.
                        let _ = writeln!(out, "{}", line);
                    }
                } else {
                    buffer.push(LocusData {
                        id: rec.id,
                        ref_counts: rec.ref_counts,
                        total_counts: rec.total_counts,
                    });
                    if buffer.len() >= chunk_size {
                        let chunk: Vec<LocusData> = buffer.drain(..).collect();
                        let results: Vec<String> = chunk
                            .into_par_iter()
                            .filter_map(|locus| {
                                process_locus(locus.id, locus.ref_counts, locus.total_counts)
                            })
                            .collect();
                        for line in results {
                            let _ = writeln!(out, "{}", line);
                        }
                    }
                }
            }).map_err(|e| anyhow::anyhow!("Failed to read VCF: {}", e))?;

            // Flush remaining buffered loci
            if !buffer.is_empty() {
                let chunk: Vec<LocusData> = buffer.drain(..).collect();
                let results: Vec<String> = chunk
                    .into_par_iter()
                    .filter_map(|locus| process_locus(locus.id, locus.ref_counts, locus.total_counts))
                    .collect();
                for line in results {
                    writeln!(out, "{}", line)?;
                }
            }
            out.flush()?;
        }
    }

    Ok(())
}
