pub mod io;
pub mod math;
pub mod model;

use ndarray::Array1;
use io::LocusData;
use rayon::prelude::*;
use std::cell::RefCell;
use std::io::Write;
use std::sync::Arc;
use std::time::Instant;

// Re-export FitMode for CLI usage
pub use model::FitMode;

/// Calculates adaptive chunk size based on available parallelism.
/// Conservative sizing to minimize overhead while maintaining parallelism.
fn adaptive_chunk_size() -> usize {
    let threads = rayon::current_num_threads();
    if threads <= 1 {
        512 // Single-threaded: use larger chunks
    } else if threads <= 16 {
        256 // Default for most common cases (2-16 threads)
    } else if threads <= 32 {
        192 // High parallelism: moderate reduction
    } else {
        128 // Very high parallelism (32+ threads): smaller chunks
    }
}

/// Formats seconds into a human-readable time string.
fn format_duration(seconds: f64) -> String {
    if seconds < 60.0 {
        format!("{:.1}s", seconds)
    } else if seconds < 3600.0 {
        let mins = (seconds / 60.0).floor();
        let secs = seconds % 60.0;
        format!("{}m {:.0}s", mins, secs)
    } else {
        let hours = (seconds / 3600.0).floor();
        let mins = ((seconds % 3600.0) / 60.0).floor();
        format!("{}h {}m", hours, mins)
    }
}

/// Creates a simple text-based progress bar.
fn progress_bar(percentage: f64, width: usize) -> String {
    let filled = ((percentage / 100.0) * width as f64).round() as usize;
    let empty = width.saturating_sub(filled);
    format!("[{}{}]", "█".repeat(filled), "░".repeat(empty))
}

fn maybe_report_progress(
    processed: usize,
    total: Option<usize>,
    start: &Instant,
    last_report: &mut Instant,
) {
    if last_report.elapsed().as_secs_f64() < 2.0 {
        return;
    }
    let elapsed = start.elapsed().as_secs_f64();
    let rate = processed as f64 / elapsed.max(1e-3);

    match total {
        Some(t) => {
            let pct = (processed as f64 / t.max(1) as f64) * 100.0;
            let remaining = t.saturating_sub(processed);
            let eta_secs = if rate > 0.0 {
                remaining as f64 / rate
            } else {
                0.0
            };

            let bar = progress_bar(pct, 20);
            eprint!(
                "\r{} {}/{} ({:.1}%) | Elapsed: {} | Rate: {:.0}/s | ETA: {}",
                bar,
                processed,
                t,
                pct,
                format_duration(elapsed),
                rate,
                format_duration(eta_secs)
            );
            let _ = std::io::stderr().flush();
        }
        None => {
            // Spinner for streaming mode
            let spinner = ["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"];
            let spinner_idx = (elapsed as usize / 2) % spinner.len();
            eprint!(
                "\r{} Processing {} loci | Elapsed: {} | Rate: {:.0}/s   ",
                spinner[spinner_idx],
                processed,
                format_duration(elapsed),
                rate
            );
            let _ = std::io::stderr().flush();
        }
    }
    *last_report = Instant::now();
}

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

#[derive(Clone, Copy, Debug)]
pub enum OutputFormat {
    Matrix,
    Stats,
    Beagle,
    Vcf,
    PlinkRaw,
    GwasPoly,
}

impl OutputFormat {
    /// Returns the standard file extension for this format (without the dot).
    pub fn extension(&self) -> &'static str {
        match self {
            OutputFormat::Matrix => "matrix.tsv",
            OutputFormat::Stats => "stats.tsv",
            OutputFormat::Beagle => "bgl",
            OutputFormat::Vcf => "vcf",
            OutputFormat::PlinkRaw => "raw",
            OutputFormat::GwasPoly => "gwaspoly.tsv",
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub enum CompressMode {
    None,
    Gzip,
}

/// Applies the correct file suffix to an output path based on format and compression.
/// Strips common format-related extensions first, then appends the correct one.
pub fn apply_output_suffix(
    path: &str,
    format: OutputFormat,
    compress: CompressMode,
) -> String {
    // Strip existing format-related extensions (order matters - longest first)
    let base = path
        .trim_end_matches(".gz")
        .trim_end_matches(".gwaspoly.tsv")
        .trim_end_matches(".matrix.tsv")
        .trim_end_matches(".stats.tsv")
        .trim_end_matches(".tsv")
        .trim_end_matches(".bgl")
        .trim_end_matches(".vcf")
        .trim_end_matches(".raw")
        .trim_end_matches(".txt")
        .to_string();

    // Build new path with correct extension
    let mut result = format!("{}.{}", base, format.extension());

    if matches!(compress, CompressMode::Gzip) {
        result.push_str(".gz");
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_apply_output_suffix() {
        // Basic suffix addition
        assert_eq!(
            apply_output_suffix("output", OutputFormat::Vcf, CompressMode::None),
            "output.vcf"
        );

        // With compression
        assert_eq!(
            apply_output_suffix("output", OutputFormat::Vcf, CompressMode::Gzip),
            "output.vcf.gz"
        );

        // Strip existing extension
        assert_eq!(
            apply_output_suffix("output.txt", OutputFormat::Vcf, CompressMode::None),
            "output.vcf"
        );

        // Keep correct extension
        assert_eq!(
            apply_output_suffix("output.vcf", OutputFormat::Vcf, CompressMode::None),
            "output.vcf"
        );

        // Replace extension and add compression
        assert_eq!(
            apply_output_suffix("output.tsv", OutputFormat::Beagle, CompressMode::Gzip),
            "output.bgl.gz"
        );

        // PLINK raw format
        assert_eq!(
            apply_output_suffix("genotypes", OutputFormat::PlinkRaw, CompressMode::None),
            "genotypes.raw"
        );

        // Stats format with specific suffix
        assert_eq!(
            apply_output_suffix("results", OutputFormat::Stats, CompressMode::Gzip),
            "results.stats.tsv.gz"
        );

        // Matrix format with specific suffix
        assert_eq!(
            apply_output_suffix("dosages", OutputFormat::Matrix, CompressMode::None),
            "dosages.matrix.tsv"
        );

        // Matrix format - strip existing .matrix.tsv
        assert_eq!(
            apply_output_suffix("dosages.matrix.tsv", OutputFormat::Matrix, CompressMode::None),
            "dosages.matrix.tsv"
        );

        // Stats format - strip existing .stats.tsv
        assert_eq!(
            apply_output_suffix("results.stats.tsv", OutputFormat::Stats, CompressMode::Gzip),
            "results.stats.tsv.gz"
        );

        // Matrix with compression
        assert_eq!(
            apply_output_suffix("data", OutputFormat::Matrix, CompressMode::Gzip),
            "data.matrix.tsv.gz"
        );
    }
}

// Re-export LocusOutput from output module
use io::output::LocusOutput;

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
    output: Option<&str>,
    format: OutputFormat,
    compress: CompressMode,
) -> anyhow::Result<()> {
    // Apply correct suffix to output path
    let output_with_suffix = output.map(|p| {
        let suffixed = apply_output_suffix(p, format, compress);
        // Inform user if we modified the path
        if suffixed != p {
            eprintln!("Output path adjusted to: {}", suffixed);
        }
        suffixed
    });
    let output_ref = output_with_suffix.as_deref();

    // PLINK .raw requires transposition, so we must collect all results
    if matches!(format, OutputFormat::PlinkRaw) {
        return run_dosage_collected(input, ploidy, mode, verbose, threads, output_ref, format, compress);
    }

    // All other formats can stream
    run_dosage_streaming(input, ploidy, mode, verbose, threads, output_ref, format, compress)
}

/// Streaming path: processes loci in chunks and writes immediately.
/// Used for matrix, stats, beagle, and vcf formats.
fn run_dosage_streaming(
    input: InputSource,
    ploidy: usize,
    mode: FitMode,
    verbose: bool,
    threads: Option<usize>,
    output: Option<&str>,
    format: OutputFormat,
    compress: CompressMode,
) -> anyhow::Result<()> {
    // Configure rayon thread pool if requested.
    if let Some(n) = threads {
        if n > 0 {
            if let Err(e) = rayon::ThreadPoolBuilder::new().num_threads(n).build_global() {
                eprintln!("Warning: Failed to set Rayon thread pool to {} threads: {}", n, e);
                eprintln!("         Rayon global pool may have already been initialized.");
            }
        }
    }

    // Always log actual thread count being used
    let actual_threads = rayon::current_num_threads();
    if let Some(requested) = threads {
        if actual_threads != requested {
            eprintln!("Note: Requested {} threads but Rayon is using {} threads", requested, actual_threads);
        }
    }

    let input_desc = match &input {
        InputSource::TwoLineCsv(path) => format!("csv={}", path),
        InputSource::RefTotalMatrices { ref_path, total_path } => {
            format!("ref={} total={}", ref_path, total_path)
        }
        InputSource::Vcf { path, chunk_size } => {
            format!("vcf={} chunk-size={}", path, chunk_size.unwrap_or_else(adaptive_chunk_size))
        }
    };

    let thread_desc = threads
        .map(|t| t.to_string())
        .unwrap_or_else(|| "auto".to_string());
    let output_desc = output.unwrap_or("stdout");
    let format_desc = match format {
        OutputFormat::Matrix => "matrix",
        OutputFormat::Stats => "stats",
        OutputFormat::Beagle => "beagle",
        OutputFormat::Vcf => "vcf",
        OutputFormat::PlinkRaw => "plink-raw",
        OutputFormat::GwasPoly => "gwaspoly",
    };
    let compress_desc = match compress {
        CompressMode::None => "none",
        CompressMode::Gzip => "gzip",
    };

    eprintln!(
        "{} v{} | mode={:?} | ploidy={} | {} | threads={} | output={} | format={} | compress={}",
        env!("CARGO_PKG_NAME"),
        env!("CARGO_PKG_VERSION"),
        mode,
        ploidy,
        input_desc,
        thread_desc,
        output_desc,
        format_desc,
        compress_desc
    );
    if verbose {
        let threads = rayon::current_num_threads();
        eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
        eprintln!("binx-dosage: Genotype Likelihood Estimation");
        eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
        eprintln!("  Ploidy: {}", ploidy);
        eprintln!("  Mode: {:?}", mode);
        eprintln!("  Threads: {}", threads);
        eprintln!("  Chunk size: {}", adaptive_chunk_size());
        eprintln!("  Input: {}", input_desc);
        eprintln!("  Output: {} ({})", format_desc, compress_desc);
        eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    }

    // Create streaming writer
    let mut writer = io::output::StreamingWriter::new(format, compress, output, ploidy)?;

    let start = Instant::now();
    let mut last_report = Instant::now();
    let mut processed: usize = 0;

    let process_locus = |locus_id: String,
                          ref_counts: Array1<u32>,
                          total_counts: Array1<u32>,
                          vcf_meta: Option<(Arc<String>, u64, Arc<String>, Arc<String>)>| {
        match run_norm_model(&ref_counts, &total_counts, ploidy, mode) {
            Ok(res) => {
                let (vcf_chrom, vcf_pos, vcf_ref, vcf_alt) = if let Some((c, p, r, a)) = vcf_meta {
                    (Some(c), Some(p), Some(r), Some(a))
                } else {
                    (None, None, None, None)
                };

                // Flatten probabilities: Vec<Vec<f64>> -> Vec<f64>
                let num_genotypes = if res.genotype_probs.is_empty() {
                    0
                } else {
                    res.genotype_probs[0].len()
                };
                let num_samples = res.genotype_probs.len();
                let mut flat_probs = Vec::with_capacity(num_samples * num_genotypes);
                for sample_probs in res.genotype_probs {
                    flat_probs.extend(sample_probs);
                }

                Some(LocusOutput {
                    id: locus_id,
                    best: res.best_genotypes,
                    probs: flat_probs,
                    num_genotypes,
                    bias: res.bias,
                    rho: res.overdispersion,
                    mu: res.model_mu,
                    sigma: res.model_sigma,
                    loglik: res.final_log_lik,
                    ref_counts: Arc::new(ref_counts),
                    total_counts: Arc::new(total_counts),
                    vcf_chrom,
                    vcf_pos,
                    vcf_ref,
                    vcf_alt,
                })
            },
            Err(e) => {
                eprintln!("Error processing locus {}: {}", locus_id, e);
                None
            },
        }
    };

    match input {
        InputSource::TwoLineCsv(csv_file) => {
            if verbose {
                eprintln!("Parsing {}...", csv_file);
            }
            let data = io::parse_two_line_csv(&csv_file)
                .map_err(|e| anyhow::anyhow!("Failed to parse CSV file: {}", e))?;
            if verbose {
                eprintln!("Found {} loci from CSV", data.loci.len());
            }

            writer.write_header(&data.sample_names)?;

            let chunk_size = adaptive_chunk_size();
            let total = data.loci.len();
            let mut iter = data.loci.into_iter();
            loop {
                let chunk: Vec<LocusData> = iter.by_ref().take(chunk_size).collect();
                if chunk.is_empty() {
                    break;
                }
                let chunk_len = chunk.len();
                let chunk_results: Vec<LocusOutput> = if matches!(mode, FitMode::Turbo | FitMode::TurboAuto) {
                    chunk
                        .into_iter()
                        .filter_map(|locus| {
                            let vcf_meta = if let (Some(c), Some(p), Some(r), Some(a)) =
                                (&locus.vcf_chrom, locus.vcf_pos, &locus.vcf_ref, &locus.vcf_alt) {
                                Some((c.clone(), p, r.clone(), a.clone()))
                            } else {
                                None
                            };
                            process_locus(locus.id, locus.ref_counts, locus.total_counts, vcf_meta)
                        })
                        .collect()
                } else {
                    chunk
                        .into_par_iter()
                        .filter_map(|locus| {
                            let vcf_meta = if let (Some(c), Some(p), Some(r), Some(a)) =
                                (&locus.vcf_chrom, locus.vcf_pos, &locus.vcf_ref, &locus.vcf_alt) {
                                Some((c.clone(), p, r.clone(), a.clone()))
                            } else {
                                None
                            };
                            process_locus(locus.id, locus.ref_counts, locus.total_counts, vcf_meta)
                        })
                        .collect()
                };

                writer.write_chunk(&chunk_results)?;
                processed += chunk_len;
                maybe_report_progress(processed, Some(total), &start, &mut last_report);
            }
        }
        InputSource::RefTotalMatrices { ref_path, total_path } => {
            if verbose {
                eprintln!("Parsing ref matrix {} and total matrix {}...", ref_path, total_path);
            }
            let matrices = io::parse_ref_total_matrices(&ref_path, &total_path)
                .map_err(|e| anyhow::anyhow!("Failed to parse ref/total matrices: {}", e))?;
            if verbose {
                eprintln!("Found {} loci from matrices", matrices.marker_ids.len());
            }

            writer.write_header(&matrices.sample_names)?;

            let chunk_size = adaptive_chunk_size();
            let total = matrices.marker_ids.len();
            let mut start_idx = 0;
            while start_idx < total {
                let end = (start_idx + chunk_size).min(total);
                let chunk: Vec<usize> = (start_idx..end).collect();
                let chunk_results: Vec<LocusOutput> = if matches!(mode, FitMode::Turbo | FitMode::TurboAuto) {
                    chunk
                        .into_iter()
                        .filter_map(|row_idx| {
                            let ref_counts = matrices.ref_counts.row(row_idx).to_owned();
                            let total_counts = matrices.total_counts.row(row_idx).to_owned();
                            process_locus(matrices.marker_ids[row_idx].clone(), ref_counts, total_counts, None)
                        })
                        .collect()
                } else {
                    chunk
                        .into_par_iter()
                        .filter_map(|row_idx| {
                            let ref_counts = matrices.ref_counts.row(row_idx).to_owned();
                            let total_counts = matrices.total_counts.row(row_idx).to_owned();
                            process_locus(matrices.marker_ids[row_idx].clone(), ref_counts, total_counts, None)
                        })
                        .collect()
                };

                writer.write_chunk(&chunk_results)?;
                processed += end - start_idx;
                maybe_report_progress(processed, Some(total), &start, &mut last_report);
                start_idx = end;
            }
        }
        InputSource::Vcf { path, chunk_size } => {
            let mut chunk_size = chunk_size.unwrap_or_else(adaptive_chunk_size);
            // Turbo modes benefit from smaller chunks to surface progress sooner and avoid long, silent batches
            if matches!(mode, FitMode::Turbo | FitMode::TurboAuto) {
                chunk_size = chunk_size.min(64).max(1);
            }

            // Pre-count VCF records for progress bar with ETA
            let total_loci = if verbose {
                eprintln!("Counting VCF records...");
                match io::count_vcf_records(&path) {
                    Ok(count) => {
                        eprintln!("Found {} loci", count);
                        Some(count)
                    }
                    Err(e) => {
                        eprintln!("Warning: Could not count records ({}), progress will show without ETA", e);
                        None
                    }
                }
            } else {
                match io::count_vcf_records(&path) {
                    Ok(count) => Some(count),
                    Err(_) => None,
                }
            };

            if verbose {
                let desc = if chunk_size == 0 {
                    "one by one".to_string()
                } else {
                    format!("chunks of {}", chunk_size)
                };
                eprintln!("Streaming VCF {} ({})", path, desc);
            }
            let mut buffer: Vec<LocusData> = Vec::with_capacity(chunk_size);
            let sample_names_opt = RefCell::new(None::<Vec<String>>);
            let header_written = RefCell::new(false);
            let writer_cell = RefCell::new(writer);
            let write_error = RefCell::new(None::<anyhow::Error>);

            io::stream_vcf_records(&path, |names| {
                *sample_names_opt.borrow_mut() = Some(names.to_vec());
            }, |rec| {
                // Write header on first record (after we have sample names)
                if !*header_written.borrow() {
                    if let Some(ref names) = *sample_names_opt.borrow() {
                        let _ = writer_cell.borrow_mut().write_header(names);
                        *header_written.borrow_mut() = true;
                    }
                }

                if chunk_size == 0 {
                    let vcf_meta = Some((rec.chrom.clone(), rec.pos, rec.ref_allele.clone(), rec.alt_allele.clone()));
                    if let Some(res) = process_locus(rec.id, rec.ref_counts, rec.total_counts, vcf_meta) {
                        if let Err(e) = writer_cell.borrow_mut().write_chunk(&[res]) {
                            *write_error.borrow_mut() = Some(e);
                            return;
                        }
                        processed += 1;
                        maybe_report_progress(processed, total_loci, &start, &mut last_report);
                    }
                } else {
                    buffer.push(LocusData {
                        id: rec.id,
                        ref_counts: rec.ref_counts,
                        total_counts: rec.total_counts,
                        vcf_chrom: Some(rec.chrom),
                        vcf_pos: Some(rec.pos),
                        vcf_ref: Some(rec.ref_allele),
                        vcf_alt: Some(rec.alt_allele),
                    });
                    if buffer.len() >= chunk_size {
                        let chunk: Vec<LocusData> = buffer.drain(..).collect();
                        let chunk_len = chunk.len();
                        // Conditional parallelism: Turbo modes process loci sequentially (inner parallelism),
                        // other modes process loci in parallel (outer parallelism)
                        let chunk_results: Vec<LocusOutput> = if matches!(mode, FitMode::Turbo | FitMode::TurboAuto) {
                            chunk
                                .into_iter()  // Sequential for Turbo modes
                                .filter_map(|locus| {
                                    let vcf_meta = if let (Some(c), Some(p), Some(r), Some(a)) =
                                        (&locus.vcf_chrom, locus.vcf_pos, &locus.vcf_ref, &locus.vcf_alt) {
                                        Some((c.clone(), p, r.clone(), a.clone()))
                                    } else {
                                        None
                                    };
                                    process_locus(locus.id, locus.ref_counts, locus.total_counts, vcf_meta)
                                })
                                .collect()
                        } else {
                            chunk
                                .into_par_iter()  // Parallel for other modes
                                .filter_map(|locus| {
                                    let vcf_meta = if let (Some(c), Some(p), Some(r), Some(a)) =
                                        (&locus.vcf_chrom, locus.vcf_pos, &locus.vcf_ref, &locus.vcf_alt) {
                                        Some((c.clone(), p, r.clone(), a.clone()))
                                    } else {
                                        None
                                    };
                                    process_locus(locus.id, locus.ref_counts, locus.total_counts, vcf_meta)
                                })
                                .collect()
                        };
                        processed += chunk_len;
                        if let Err(e) = writer_cell.borrow_mut().write_chunk(&chunk_results) {
                            *write_error.borrow_mut() = Some(e);
                            return;
                        }
                        maybe_report_progress(processed, total_loci, &start, &mut last_report);
                    }
                }
            }).map_err(|e| anyhow::anyhow!("Failed to read VCF: {}", e))?;

            let mut writer = writer_cell.into_inner();

            // Flush remaining buffered loci
            if !buffer.is_empty() {
                let chunk: Vec<LocusData> = buffer.drain(..).collect();
                let chunk_len = chunk.len();
                // Conditional parallelism for flush (same as main processing)
                let chunk_results: Vec<LocusOutput> = if matches!(mode, FitMode::Turbo | FitMode::TurboAuto) {
                    chunk
                        .into_iter()  // Sequential for Turbo modes
                        .filter_map(|locus| {
                            let vcf_meta = if let (Some(c), Some(p), Some(r), Some(a)) =
                                (&locus.vcf_chrom, locus.vcf_pos, &locus.vcf_ref, &locus.vcf_alt) {
                                Some((c.clone(), p, r.clone(), a.clone()))
                            } else {
                                None
                            };
                            process_locus(locus.id, locus.ref_counts, locus.total_counts, vcf_meta)
                        })
                        .collect()
                } else {
                    chunk
                        .into_par_iter()  // Parallel for other modes
                        .filter_map(|locus| {
                            let vcf_meta = if let (Some(c), Some(p), Some(r), Some(a)) =
                                (&locus.vcf_chrom, locus.vcf_pos, &locus.vcf_ref, &locus.vcf_alt) {
                                Some((c.clone(), p, r.clone(), a.clone()))
                            } else {
                                None
                            };
                            process_locus(locus.id, locus.ref_counts, locus.total_counts, vcf_meta)
                        })
                        .collect()
                };
                processed += chunk_len;
                writer.write_chunk(&chunk_results)?;
                maybe_report_progress(processed, None, &start, &mut last_report);
            }

            if sample_names_opt.borrow().is_none() {
                return Err(anyhow::anyhow!("VCF header with sample names was not found"));
            }

            if let Some(e) = write_error.into_inner() {
                return Err(anyhow::anyhow!("Failed to write output chunk: {}", e));
            }

            writer.finish()?;

            let elapsed = start.elapsed().as_secs_f64();
            let rate = processed as f64 / elapsed.max(1e-3);
            eprintln!("\n✓ Completed {} loci in {} ({:.0} loci/s)", processed, format_duration(elapsed), rate);

            return Ok(());
        }
    }

    writer.finish()?;

    let elapsed = start.elapsed().as_secs_f64();
    let rate = processed as f64 / elapsed.max(1e-3);
    eprintln!("\n✓ Completed {} loci in {} ({:.0} loci/s)", processed, format_duration(elapsed), rate);

    Ok(())
}

/// Collection path: collects all results in memory before writing.
/// Used only for PLINK .raw format which requires transposition.
fn run_dosage_collected(
    input: InputSource,
    ploidy: usize,
    mode: FitMode,
    verbose: bool,
    threads: Option<usize>,
    output: Option<&str>,
    format: OutputFormat,
    compress: CompressMode,
) -> anyhow::Result<()> {
    // Configure rayon thread pool if requested.
    if let Some(n) = threads {
        if n > 0 {
            if let Err(e) = rayon::ThreadPoolBuilder::new().num_threads(n).build_global() {
                eprintln!("Warning: Failed to set Rayon thread pool to {} threads: {}", n, e);
                eprintln!("         Rayon global pool may have already been initialized.");
            }
        }
    }

    // Always log actual thread count being used
    let actual_threads = rayon::current_num_threads();
    if let Some(requested) = threads {
        if actual_threads != requested {
            eprintln!("Note: Requested {} threads but Rayon is using {} threads", requested, actual_threads);
        }
    }

    let input_desc = match &input {
        InputSource::TwoLineCsv(path) => format!("csv={}", path),
        InputSource::RefTotalMatrices { ref_path, total_path } => {
            format!("ref={} total={}", ref_path, total_path)
        }
        InputSource::Vcf { path, chunk_size } => {
            format!("vcf={} chunk-size={}", path, chunk_size.unwrap_or_else(adaptive_chunk_size))
        }
    };

    let thread_desc = threads
        .map(|t| t.to_string())
        .unwrap_or_else(|| "auto".to_string());
    let output_desc = output.unwrap_or("stdout");
    let format_desc = match format {
        OutputFormat::Matrix => "matrix",
        OutputFormat::Stats => "stats",
        OutputFormat::Beagle => "beagle",
        OutputFormat::Vcf => "vcf",
        OutputFormat::PlinkRaw => "plink-raw",
        OutputFormat::GwasPoly => "gwaspoly",
    };
    let compress_desc = match compress {
        CompressMode::None => "none",
        CompressMode::Gzip => "gzip",
    };

    eprintln!(
        "{} v{} | mode={:?} | ploidy={} | {} | threads={} | output={} | format={} | compress={}",
        env!("CARGO_PKG_NAME"),
        env!("CARGO_PKG_VERSION"),
        mode,
        ploidy,
        input_desc,
        thread_desc,
        output_desc,
        format_desc,
        compress_desc
    );
    if verbose {
        let threads = rayon::current_num_threads();
        eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
        eprintln!("binx-dosage: Genotype Likelihood Estimation");
        eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
        eprintln!("  Ploidy: {}", ploidy);
        eprintln!("  Mode: {:?}", mode);
        eprintln!("  Threads: {}", threads);
        eprintln!("  Chunk size: {}", adaptive_chunk_size());
        eprintln!("  Input: {}", input_desc);
        eprintln!("  Output: {} ({})", format_desc, compress_desc);
        eprintln!("  ⚠ Collection mode (PLINK format)");
        eprintln!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");
    }

    let start = Instant::now();
    let mut last_report = Instant::now();
    let mut processed: usize = 0;

    let process_locus = |locus_id: String,
                          ref_counts: Array1<u32>,
                          total_counts: Array1<u32>,
                          vcf_meta: Option<(Arc<String>, u64, Arc<String>, Arc<String>)>| {
        match run_norm_model(&ref_counts, &total_counts, ploidy, mode) {
            Ok(res) => {
                let (vcf_chrom, vcf_pos, vcf_ref, vcf_alt) = if let Some((c, p, r, a)) = vcf_meta {
                    (Some(c), Some(p), Some(r), Some(a))
                } else {
                    (None, None, None, None)
                };

                // Flatten probabilities: Vec<Vec<f64>> -> Vec<f64>
                let num_genotypes = if res.genotype_probs.is_empty() {
                    0
                } else {
                    res.genotype_probs[0].len()
                };
                let num_samples = res.genotype_probs.len();
                let mut flat_probs = Vec::with_capacity(num_samples * num_genotypes);
                for sample_probs in res.genotype_probs {
                    flat_probs.extend(sample_probs);
                }

                Some(LocusOutput {
                    id: locus_id,
                    best: res.best_genotypes,
                    probs: flat_probs,
                    num_genotypes,
                    bias: res.bias,
                    rho: res.overdispersion,
                    mu: res.model_mu,
                    sigma: res.model_sigma,
                    loglik: res.final_log_lik,
                    ref_counts: Arc::new(ref_counts),
                    total_counts: Arc::new(total_counts),
                    vcf_chrom,
                    vcf_pos,
                    vcf_ref,
                    vcf_alt,
                })
            },
            Err(e) => {
                eprintln!("Error processing locus {}: {}", locus_id, e);
                None
            },
        }
    };

    let (results, sample_names) = match input {
        InputSource::TwoLineCsv(csv_file) => {
            if verbose {
                eprintln!("Parsing {}...", csv_file);
            }
            let data = io::parse_two_line_csv(&csv_file)
                .map_err(|e| anyhow::anyhow!("Failed to parse CSV file: {}", e))?;
            if verbose {
                eprintln!("Found {} loci from CSV", data.loci.len());
            }
            let sample_names = data.sample_names.clone();
            let chunk_size = adaptive_chunk_size();
            let total = data.loci.len();
            let mut iter = data.loci.into_iter();
            let mut all_results = Vec::with_capacity(total);
            loop {
                let chunk: Vec<LocusData> = iter.by_ref().take(chunk_size).collect();
                if chunk.is_empty() {
                    break;
                }
                let chunk_len = chunk.len();
                let chunk_results: Vec<LocusOutput> = if matches!(mode, FitMode::Turbo | FitMode::TurboAuto) {
                    chunk
                        .into_iter()
                        .filter_map(|locus| {
                            let vcf_meta = if let (Some(c), Some(p), Some(r), Some(a)) =
                                (&locus.vcf_chrom, locus.vcf_pos, &locus.vcf_ref, &locus.vcf_alt) {
                                Some((c.clone(), p, r.clone(), a.clone()))
                            } else {
                                None
                            };
                            process_locus(locus.id, locus.ref_counts, locus.total_counts, vcf_meta)
                        })
                        .collect()
                } else {
                    chunk
                        .into_par_iter()
                        .filter_map(|locus| {
                            let vcf_meta = if let (Some(c), Some(p), Some(r), Some(a)) =
                                (&locus.vcf_chrom, locus.vcf_pos, &locus.vcf_ref, &locus.vcf_alt) {
                                Some((c.clone(), p, r.clone(), a.clone()))
                            } else {
                                None
                            };
                            process_locus(locus.id, locus.ref_counts, locus.total_counts, vcf_meta)
                        })
                        .collect()
                };
                all_results.extend(chunk_results);
                processed += chunk_len;
                maybe_report_progress(processed, Some(total), &start, &mut last_report);
            }
            (all_results, sample_names)
        }
        InputSource::RefTotalMatrices { ref_path, total_path } => {
            if verbose {
                eprintln!("Parsing ref matrix {} and total matrix {}...", ref_path, total_path);
            }
            let matrices = io::parse_ref_total_matrices(&ref_path, &total_path)
                .map_err(|e| anyhow::anyhow!("Failed to parse ref/total matrices: {}", e))?;
            if verbose {
                eprintln!("Found {} loci from matrices", matrices.marker_ids.len());
            }
            let sample_names = matrices.sample_names.clone();
            let chunk_size = adaptive_chunk_size();
            let total = matrices.marker_ids.len();
            let mut start_idx = 0;
            let mut all_results = Vec::with_capacity(total);
            while start_idx < total {
                let end = (start_idx + chunk_size).min(total);
                let chunk: Vec<usize> = (start_idx..end).collect();
                let chunk_results: Vec<LocusOutput> = if matches!(mode, FitMode::Turbo | FitMode::TurboAuto) {
                    chunk
                        .into_iter()
                        .filter_map(|row_idx| {
                            let ref_counts = matrices.ref_counts.row(row_idx).to_owned();
                            let total_counts = matrices.total_counts.row(row_idx).to_owned();
                            process_locus(matrices.marker_ids[row_idx].clone(), ref_counts, total_counts, None)
                        })
                        .collect()
                } else {
                    chunk
                        .into_par_iter()
                        .filter_map(|row_idx| {
                            let ref_counts = matrices.ref_counts.row(row_idx).to_owned();
                            let total_counts = matrices.total_counts.row(row_idx).to_owned();
                            process_locus(matrices.marker_ids[row_idx].clone(), ref_counts, total_counts, None)
                        })
                        .collect()
                };
                all_results.extend(chunk_results);
                processed += end - start_idx;
                maybe_report_progress(processed, Some(total), &start, &mut last_report);
                start_idx = end;
            }
            (all_results, sample_names)
        }
        InputSource::Vcf { path, chunk_size } => {
            let mut chunk_size = chunk_size.unwrap_or_else(adaptive_chunk_size);
            if matches!(mode, FitMode::Turbo | FitMode::TurboAuto) {
                chunk_size = chunk_size.min(64).max(1);
            }

            // Pre-count VCF records for progress bar with ETA
            let total_loci = if verbose {
                eprintln!("Counting VCF records...");
                match io::count_vcf_records(&path) {
                    Ok(count) => {
                        eprintln!("Found {} loci", count);
                        Some(count)
                    }
                    Err(e) => {
                        eprintln!("Warning: Could not count records ({}), progress will show without ETA", e);
                        None
                    }
                }
            } else {
                match io::count_vcf_records(&path) {
                    Ok(count) => Some(count),
                    Err(_) => None,
                }
            };

            if verbose {
                let desc = if chunk_size == 0 {
                    "one by one".to_string()
                } else {
                    format!("chunks of {}", chunk_size)
                };
                eprintln!("Streaming VCF {} ({})", path, desc);
            }
            let mut buffer: Vec<LocusData> = Vec::with_capacity(chunk_size);
            let mut sample_names: Option<Vec<String>> = None;
            let mut all_results = Vec::with_capacity(1024); // Estimate: start with capacity for 1K loci

            io::stream_vcf_records(&path, |names| {
                sample_names = Some(names.to_vec());
            }, |rec| {
                if chunk_size == 0 {
                    let vcf_meta = Some((rec.chrom.clone(), rec.pos, rec.ref_allele.clone(), rec.alt_allele.clone()));
                    if let Some(res) = process_locus(rec.id, rec.ref_counts, rec.total_counts, vcf_meta) {
                        all_results.push(res);
                        processed += 1;
                        maybe_report_progress(processed, total_loci, &start, &mut last_report);
                    }
                } else {
                    buffer.push(LocusData {
                        id: rec.id,
                        ref_counts: rec.ref_counts,
                        total_counts: rec.total_counts,
                        vcf_chrom: Some(rec.chrom),
                        vcf_pos: Some(rec.pos),
                        vcf_ref: Some(rec.ref_allele),
                        vcf_alt: Some(rec.alt_allele),
                    });
                    if buffer.len() >= chunk_size {
                        let chunk: Vec<LocusData> = buffer.drain(..).collect();
                        let chunk_results: Vec<LocusOutput> = chunk
                            .into_par_iter()
                            .filter_map(|locus| {
                                let vcf_meta = if let (Some(c), Some(p), Some(r), Some(a)) =
                                    (&locus.vcf_chrom, locus.vcf_pos, &locus.vcf_ref, &locus.vcf_alt) {
                                    Some((c.clone(), p, r.clone(), a.clone()))
                                } else {
                                    None
                                };
                                process_locus(locus.id, locus.ref_counts, locus.total_counts, vcf_meta)
                            })
                            .collect();
                        processed += chunk_results.len();
                        all_results.extend(chunk_results);
                        maybe_report_progress(processed, total_loci, &start, &mut last_report);
                    }
                }
            }).map_err(|e| anyhow::anyhow!("Failed to read VCF: {}", e))?;

            // Flush remaining buffered loci
            if !buffer.is_empty() {
                let chunk: Vec<LocusData> = buffer.drain(..).collect();
                let chunk_len = chunk.len();
                // Conditional parallelism for flush (same as main processing)
                let chunk_results: Vec<LocusOutput> = if matches!(mode, FitMode::Turbo | FitMode::TurboAuto) {
                    chunk
                        .into_iter()  // Sequential for Turbo modes
                        .filter_map(|locus| {
                            let vcf_meta = if let (Some(c), Some(p), Some(r), Some(a)) =
                                (&locus.vcf_chrom, locus.vcf_pos, &locus.vcf_ref, &locus.vcf_alt) {
                                Some((c.clone(), p, r.clone(), a.clone()))
                            } else {
                                None
                            };
                            process_locus(locus.id, locus.ref_counts, locus.total_counts, vcf_meta)
                        })
                        .collect()
                } else {
                    chunk
                        .into_par_iter()  // Parallel for other modes
                        .filter_map(|locus| {
                            let vcf_meta = if let (Some(c), Some(p), Some(r), Some(a)) =
                                (&locus.vcf_chrom, locus.vcf_pos, &locus.vcf_ref, &locus.vcf_alt) {
                                Some((c.clone(), p, r.clone(), a.clone()))
                            } else {
                                None
                            };
                            process_locus(locus.id, locus.ref_counts, locus.total_counts, vcf_meta)
                        })
                        .collect()
                };
                processed += chunk_len;
                all_results.extend(chunk_results);
                maybe_report_progress(processed, None, &start, &mut last_report);
            }

            let sample_names = sample_names
                .ok_or_else(|| anyhow::anyhow!("VCF header with sample names was not found"))?;

            (all_results, sample_names)
        }
    };

    if verbose {
        eprintln!("Writing output...");
    }

    io::output::write_output(results, &sample_names, format, compress, output, ploidy)?;

    let elapsed = start.elapsed().as_secs_f64();
    let rate = processed as f64 / elapsed.max(1e-3);
    eprintln!("\n✓ Completed {} loci in {} ({:.0} loci/s)", processed, format_duration(elapsed), rate);

    Ok(())
}
