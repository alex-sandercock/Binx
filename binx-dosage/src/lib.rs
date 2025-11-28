pub mod io;
pub mod math;
pub mod model;

use ndarray::Array1;
use io::{LocusData, MatrixData};

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
) -> anyhow::Result<()> {
    let loci: Vec<LocusData>;
    let matrix_input: Option<MatrixData>;

    match input {
        InputSource::TwoLineCsv(csv_file) => {
            if verbose {
                println!("Parsing {}...", csv_file);
            }
            loci = io::parse_two_line_csv(&csv_file)
                .map_err(|e| anyhow::anyhow!("Failed to parse CSV file: {}", e))?;
            matrix_input = None;
        }
        InputSource::RefTotalMatrices { ref_path, total_path } => {
            if verbose {
                println!("Parsing ref matrix {} and total matrix {}...", ref_path, total_path);
            }
            let matrices = io::parse_ref_total_matrices(&ref_path, &total_path)
                .map_err(|e| anyhow::anyhow!("Failed to parse ref/total matrices: {}", e))?;
            matrix_input = Some(matrices);
            loci = Vec::new(); // not used in this branch
        }
    }

    if verbose {
        let count = if let Some(mat) = &matrix_input { mat.marker_ids.len() } else { loci.len() };
        println!("Found {} loci. Running Norm model...", count);
        println!("Locus\tBias\tRho\tMu\tSigma\tLogLik\tGenotypes");
    } else {
        println!("Locus\tBias\tRho\tMu\tSigma\tLogLik\tGenotypes");
    }

    if let Some(mat) = matrix_input {
        for (row_idx, locus_id) in mat.marker_ids.iter().enumerate() {
            let ref_counts = mat.ref_counts.row(row_idx).to_owned();
            let total_counts = mat.total_counts.row(row_idx).to_owned();

            match run_norm_model(&ref_counts, &total_counts, ploidy, mode) {
                Ok(mut res) => {
                    res.locus_id = locus_id.clone();
                    let geno_str: Vec<String> = res.best_genotypes.iter().map(|g| g.to_string()).collect();
                    println!("{}\t{:.3}\t{:.4}\t{:.3}\t{:.3}\t{:.2}\t{}",
                        res.locus_id,
                        res.bias,
                        res.overdispersion,
                        res.model_mu,
                        res.model_sigma,
                        res.final_log_lik,
                        geno_str.join(",")
                    );
                },
                Err(e) => eprintln!("Error processing locus {}: {}", locus_id, e),
            }
        }
    } else {
        for locus in loci {
            match run_norm_model(
                &locus.ref_counts,
                &locus.total_counts,
                ploidy,
                mode,
            ) {
                Ok(mut res) => {
                    res.locus_id = locus.id;

                    let geno_str: Vec<String> = res.best_genotypes.iter().map(|g| g.to_string()).collect();

                    println!("{}\t{:.3}\t{:.4}\t{:.3}\t{:.3}\t{:.2}\t{}",
                        res.locus_id,
                        res.bias,
                        res.overdispersion,
                        res.model_mu,
                        res.model_sigma,
                        res.final_log_lik,
                        geno_str.join(",")
                    );
                },
                Err(e) => eprintln!("Error processing locus {}: {}", locus.id, e),
            }
        }
    }

    Ok(())
}
