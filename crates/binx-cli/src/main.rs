use anyhow::Result;
use clap::{Parser, Subcommand};
use std::cell::RefCell;
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::str::FromStr;
use std::time::Instant;

/// binx: a Rust-based genomic analysis toolkit
#[derive(Parser)]
#[command(
    name = "binx",
    version,
    about = "binx: a Rust-based genomic analysis toolkit (GWAS, multiallelic GWAS, more to come)"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Multiallelic GWAS (stub for now)
    Multigwas {
        #[arg(long)]
        geno: String,
        #[arg(long)]
        pheno: String,
        #[arg(long)]
        trait_name: String,
        #[arg(long)]
        kinship: String,
        #[arg(long)]
        pcs: Option<String>,
        #[arg(long)]
        ploidy: u8,
        #[arg(long, default_value = "reference")]
        encoding: String,
        #[arg(long, default_value = "general")]
        model: String,
        #[arg(long)]
        out: String,
    },

    /// Run genome-wide association study
    #[command(after_help = "EXAMPLES:
    # Basic GWAS with additive model
    binx gwas --geno geno.tsv --pheno pheno.csv --trait yield --ploidy 4 --out results.csv

    # Multiple models with threshold calculation
    binx gwas --geno geno.tsv --pheno pheno.csv --trait yield --ploidy 4 \\
        --models additive,general --threshold m.eff --out results.csv

    # With kinship matrix and plots
    binx gwas --geno geno.tsv --pheno pheno.csv --trait yield --ploidy 4 \\
        --kinship kinship.tsv --plot both --out results.csv

    # Use a different method (when available)
    binx gwas --method gwaspoly --geno geno.tsv --pheno pheno.csv --trait yield --ploidy 4 --out results.csv

METHODS:
    gwaspoly        GWASpoly-style GWAS for polyploids (default)

MODELS (for gwaspoly method):
    additive        Linear effect of allele dosage
    general         Separate effect for each dosage class
    1-dom-ref       Single-dose dominance (reference)
    1-dom-alt       Single-dose dominance (alternate)
    2-dom-ref       Double-dose dominance (reference)
    2-dom-alt       Double-dose dominance (alternate)
    diplo-general   Diploidized general
    diplo-additive  Diploidized additive")]
    Gwas {
        // === Method ===
        /// GWAS method to use
        #[arg(long, default_value = "gwaspoly", help_heading = "Method")]
        method: String,

        // === Input/Output ===
        /// Genotype dosage file (TSV: marker, chr, pos, samples...)
        #[arg(long, help_heading = "Input/Output")]
        geno: String,

        /// Phenotype file (CSV: sample_id, traits...)
        #[arg(long, help_heading = "Input/Output")]
        pheno: String,

        /// Output results CSV
        #[arg(long, help_heading = "Input/Output")]
        out: String,

        /// Optional kinship matrix TSV
        #[arg(long, help_heading = "Input/Output")]
        kinship: Option<String>,

        // === Analysis ===
        /// Trait name to analyze
        #[arg(long, help_heading = "Analysis")]
        r#trait: String,

        /// Ploidy level (e.g., 2, 4, 6)
        #[arg(long, help_heading = "Analysis")]
        ploidy: u8,

        /// Gene action models (comma-separated)
        #[arg(long, default_value = "additive,general", help_heading = "Analysis")]
        models: String,

        /// Use LOCO (Leave-One-Chromosome-Out)
        #[arg(long, default_value_t = false, help_heading = "Analysis")]
        loco: bool,

        /// Covariates from phenotype file (comma-separated)
        #[arg(long, help_heading = "Analysis")]
        covariates: Option<String>,

        /// Number of principal components to include as fixed effects (P+K model)
        #[arg(long, default_value = "0", help_heading = "Analysis")]
        n_pc: usize,

        // === QC Filters ===
        /// Minimum minor allele frequency (0.0-0.5)
        #[arg(long, default_value = "0.0", help_heading = "QC Filters")]
        min_maf: f64,

        /// Maximum genotype frequency (0.0-1.0, 0=auto)
        #[arg(long, default_value = "0.0", help_heading = "QC Filters")]
        max_geno_freq: f64,

        /// Allow samples in geno but not in pheno
        #[arg(long, default_value_t = false, help_heading = "QC Filters")]
        allow_missing_samples: bool,

        // === Threshold ===
        /// Method: bonferroni, m.eff, or fdr
        #[arg(long, help_heading = "Threshold")]
        threshold: Option<String>,

        /// Significance level (default: 0.05)
        #[arg(long, default_value = "0.05", help_heading = "Threshold")]
        alpha: f64,

        // === Output Options ===
        /// Generate plots: manhattan, qq, or both
        #[arg(long, help_heading = "Output Options")]
        plot: Option<String>,

        /// Custom path for plot files
        #[arg(long, help_heading = "Output Options")]
        plot_output: Option<String>,

        /// Use parallel marker testing
        #[arg(long, default_value_t = false, help_heading = "Output Options")]
        parallel: bool,
    },

    /// Compute kinship matrix from biallelic dosages
    Kinship {
        /// Genotype dosage file (biallelic; markers x samples)
        #[arg(long)]
        geno: String,

        /// Ploidy (e.g., 2, 4, 6)
        #[arg(long)]
        ploidy: u8,

        /// Kinship method: vanraden (default) or gwaspoly
        #[arg(long, default_value = "vanraden")]
        method: String,

        /// Output TSV path for kinship matrix
        #[arg(long)]
        out: String,
    },

    /// Calculate significance thresholds for GWAS results
    #[command(after_help = "EXAMPLES:
    # Bonferroni correction
    binx threshold --results gwas.csv --method bonferroni

    # M.eff (effective number of tests) - requires genotype data
    binx threshold --results gwas.csv --method m.eff --geno geno.tsv --ploidy 4

    # FDR (Benjamini-Hochberg)
    binx threshold --results gwas.csv --method fdr --alpha 0.1

METHODS:
    bonferroni  Simple Bonferroni: -log10(alpha / n_markers)
    m.eff       Effective tests (Moskvina & Schmidt 2008), accounts for LD
    fdr         False Discovery Rate (Benjamini-Hochberg)")]
    Threshold {
        /// GWAS results CSV file
        #[arg(long)]
        results: String,

        /// Threshold method
        #[arg(long, value_parser = ["bonferroni", "m.eff", "fdr"])]
        method: String,

        /// Significance level
        #[arg(long, default_value = "0.05")]
        alpha: f64,

        /// Genotype file (required for m.eff)
        #[arg(long)]
        geno: Option<String>,

        /// Ploidy (required for m.eff)
        #[arg(long)]
        ploidy: Option<u8>,
    },

    /// Genotype dosage estimation from sequencing read counts
    Dosage {
        /// CSV file with alternating lines of Ref and Total counts per locus
        #[arg(long)]
        csv: Option<String>,

        /// Use matrix mode (requires --ref and --total)
        #[arg(long, default_value_t = false)]
        counts: bool,

        /// Ref count matrix (markers in rows, samples in columns; first column marker ID)
        #[arg(long)]
        ref_path: Option<String>,

        /// Total count matrix (markers in rows, samples in columns; first column marker ID)
        #[arg(long)]
        total_path: Option<String>,

        /// VCF file (plain or gzipped) with FORMAT/AD allele depths
        #[arg(long)]
        vcf: Option<String>,

        /// Chunk size for streaming VCF markers (0 or omit = stream one by one)
        #[arg(long)]
        chunk_size: Option<usize>,

        /// Number of threads to use for parallel dosage (default: num_cpus)
        #[arg(long)]
        threads: Option<usize>,

        /// Output path (defaults to stdout)
        #[arg(long)]
        output: Option<String>,

        /// Ploidy (e.g., 2, 4, 6)
        #[arg(long)]
        ploidy: usize,

        /// Multi-start optimization mode (auto, updog, updog-fast, updog-exact, fast, turbo, turboauto, turboauto-safe)
        #[arg(long, default_value = "auto")]
        mode: String,

        /// Output format (matrix, stats, beagle, vcf, plink, gwaspoly)
        #[arg(long, default_value = "matrix")]
        format: String,

        /// Compression mode (none, gzip)
        #[arg(long, default_value = "none")]
        compress: String,

        /// Enable verbose output
        #[arg(long, default_value_t = false)]
        verbose: bool,
    },

    /// Convert VCF to other formats
    Convert {
        /// Input VCF file (plain or gzipped)
        #[arg(long)]
        vcf: String,

        /// Output path
        #[arg(long)]
        output: String,

        /// Output format: csv (two-line ref/total counts from AD) or gwaspoly (dosages from GT)
        #[arg(long, default_value = "csv")]
        format: String,

        /// Enable verbose progress output
        #[arg(long, default_value_t = false)]
        verbose: bool,
    },

    /// Identify significant QTLs from GWAS results
    #[command(after_help = "EXAMPLES:
    # Extract QTLs from GWAS results
    binx qtl --input gwas_results.csv --output qtls.csv

    # Prune nearby signals within 1 Mb window
    binx qtl --input gwas_results.csv --bp-window 1000000

    # Pipe from gwaspoly (requires threshold in results)
    binx gwaspoly ... --threshold m.eff --out /dev/stdout 2>/dev/null | binx qtl --bp-window 1000000

NOTE: Input must have 'threshold' column (use --threshold with gwaspoly).
      Only markers with score >= threshold are reported as QTLs.")]
    Qtl {
        /// Input GWAS results file (stdin if omitted)
        #[arg(long)]
        input: Option<String>,

        /// Prune signals within this window (bp)
        #[arg(long)]
        bp_window: Option<u64>,

        /// Output file (stdout if omitted)
        #[arg(long)]
        output: Option<String>,
    },

    /// Generate GWAS plots (Manhattan, QQ, LD) from results
    #[command(after_help = "EXAMPLES:
    # Manhattan plot from GWAS results
    binx plot --input gwas.csv --output manhattan.svg

    # QQ plot for a specific model
    binx plot --input gwas.csv --output qq.svg --plot-type qq --model additive

    # LD decay plot with threshold annotation
    binx plot --input geno.tsv --output ld.svg --plot-type ld --ploidy 4 --r2-threshold 0.2

    # LD plot for specific chromosomes
    binx plot --input geno.tsv --output ld.svg --plot-type ld --ploidy 4 --chromosomes chr05,chr09")]
    Plot {
        // === Input/Output ===
        /// Input file: GWAS results CSV (manhattan/qq) or genotype TSV (ld)
        #[arg(long, help_heading = "Input/Output")]
        input: String,

        /// Output file path (.svg or .png)
        #[arg(long, help_heading = "Input/Output")]
        output: String,

        /// Plot type
        #[arg(long, default_value = "manhattan", value_parser = ["manhattan", "qq", "ld"], help_heading = "Input/Output")]
        plot_type: String,

        // === Appearance ===
        /// Plot title
        #[arg(long, help_heading = "Appearance")]
        title: Option<String>,

        /// Color theme
        #[arg(long, default_value = "classic", value_parser = ["classic", "nature", "colorful", "dark", "high_contrast"], help_heading = "Appearance")]
        theme: String,

        /// Plot width in pixels
        #[arg(long, default_value = "1200", help_heading = "Appearance")]
        width: u32,

        /// Plot height in pixels
        #[arg(long, default_value = "600", help_heading = "Appearance")]
        height: u32,

        // === Manhattan/QQ Options ===
        /// Filter to specific model (e.g., additive, general)
        #[arg(long, help_heading = "Manhattan/QQ Options")]
        model: Option<String>,

        /// Significance threshold as -log10(p)
        #[arg(long, default_value = "5.0", help_heading = "Manhattan/QQ Options")]
        threshold: f64,

        /// Suggestive threshold as -log10(p) (0 to disable)
        #[arg(long, default_value = "3.0", help_heading = "Manhattan/QQ Options")]
        suggestive: f64,

        /// Filter to specific chromosomes (comma-separated)
        #[arg(long, help_heading = "Manhattan/QQ Options")]
        chromosomes: Option<String>,

        // === LD Plot Options ===
        /// Ploidy level (required for ld plot)
        #[arg(long, help_heading = "LD Plot Options")]
        ploidy: Option<u8>,

        /// R² threshold to mark on plot (draws line, reports distance)
        #[arg(long, help_heading = "LD Plot Options")]
        r2_threshold: Option<f64>,

        /// Maximum marker pairs to sample
        #[arg(long, default_value = "10000", help_heading = "LD Plot Options")]
        max_pairs: usize,

        /// Maximum markers per chromosome
        #[arg(long, help_heading = "LD Plot Options")]
        max_loci: Option<usize>,

        /// Number of distance bins for smoothing
        #[arg(long, default_value = "50", help_heading = "LD Plot Options")]
        n_bins: usize,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Multigwas {
            geno,
            pheno,
            trait_name,
            kinship,
            pcs,
            ploidy,
            encoding,
            model,
            out,
        } => {
            binx_multigwas::run_multigwas(
                &geno,
                &pheno,
                &trait_name,
                &kinship,
                pcs.as_ref(),
                ploidy,
                &encoding,
                &model,
                &out,
            )?;
        }
        Commands::Kinship { geno, ploidy, method, out } => {
            let kin_method = binx_kinship::KinshipMethod::from_str(&method).unwrap_or_else(|e| {
                eprintln!("{}", e);
                std::process::exit(1);
            });
            binx_kinship::run_kinship(&geno, ploidy, kin_method, &out)?;
        }
        Commands::Threshold { results, method, alpha, geno, ploidy } => {
            // Parse threshold method
            let threshold_method = gwaspoly_rs::ThresholdMethod::from_str(&method).unwrap_or_else(|e| {
                eprintln!("{}", e);
                std::process::exit(1);
            });

            // Load GWAS results
            eprintln!("Loading GWAS results from {}...", results);
            let gwas_results = gwaspoly_rs::load_gwas_results_for_threshold(&results)?;
            eprintln!("Loaded {} results", gwas_results.len());

            // Calculate thresholds
            let thresholds = if threshold_method == gwaspoly_rs::ThresholdMethod::Meff {
                // M.eff requires genotype data
                let geno_path = geno.as_ref().ok_or_else(|| {
                    anyhow::anyhow!("M.eff method requires --geno argument")
                })?;
                let ploidy_val = ploidy.ok_or_else(|| {
                    anyhow::anyhow!("M.eff method requires --ploidy argument")
                })?;

                eprintln!("Loading genotype data for M.eff calculation...");
                let geno_data = gwaspoly_rs::load_genotypes(geno_path, ploidy_val)?;
                gwaspoly_rs::set_threshold(&gwas_results, &geno_data, threshold_method, alpha)?
            } else {
                // Bonferroni and FDR don't need genotype data
                gwaspoly_rs::set_threshold_simple(&gwas_results, threshold_method, alpha)?
            };

            // Print results
            gwaspoly_rs::print_thresholds(&thresholds);
        }
        Commands::Gwas {
            method,
            geno,
            pheno,
            r#trait,
            covariates,
            n_pc,
            kinship,
            allow_missing_samples,
            ploidy,
            models,
            loco,
            min_maf,
            max_geno_freq,
            out,
            parallel,
            plot,
            plot_output,
            threshold,
            alpha,
        } => {
            // Parse GWAS method
            let gwas_method = binx_gwas::GwasMethod::from_str(&method).unwrap_or_else(|e| {
                eprintln!("{}", e);
                std::process::exit(1);
            });

            let covariate_list = covariates.as_deref().map(parse_csv_list);
            let model_list: Vec<binx_gwas::GeneActionModel> = models
                .split(',')
                .map(|s| s.trim())
                .filter(|s| !s.is_empty())
                .map(|s| {
                    binx_gwas::GeneActionModel::from_str(s).unwrap_or_else(|e| {
                        eprintln!("Invalid model '{}': {}", s, e);
                        std::process::exit(1);
                    })
                })
                .collect();

            if model_list.is_empty() {
                eprintln!("No valid models specified");
                std::process::exit(1);
            }

            // Parse threshold method if provided
            let threshold_method = threshold.as_ref().map(|t| {
                binx_gwas::ThresholdMethod::from_str(t).unwrap_or_else(|e| {
                    eprintln!("{}", e);
                    std::process::exit(1);
                })
            });

            if parallel {
                eprintln!("Using parallel marker testing...");
            }
            if n_pc > 0 {
                eprintln!("Including {} principal components as fixed effects (P+K model)", n_pc);
            }
            eprintln!("Running GWAS with method: {}", gwas_method);
            binx_gwas::run_gwas(
                &geno,
                &pheno,
                &r#trait,
                covariate_list.as_deref(),
                kinship.as_deref(),
                allow_missing_samples,
                ploidy,
                &model_list,
                loco,
                min_maf,
                max_geno_freq,
                &out,
                parallel,
                threshold_method,
                alpha,
                gwas_method,
                n_pc,
            )?;

            // Generate plots if requested
            if let Some(plot_types) = plot {
                use binx_plotting::{manhattan_plot, qq_plot, PlotConfig, load_gwas_results};

                eprintln!("Loading results for plotting...");
                let results = load_gwas_results(&out)?;

                let base_path = plot_output.as_deref().unwrap_or(&out);
                let config = PlotConfig::default();

                for plot_type in plot_types.split(',').map(|s| s.trim()) {
                    match plot_type {
                        "manhattan" => {
                            let path = format!("{}.manhattan.svg", base_path.trim_end_matches(".csv"));
                            eprintln!("Generating Manhattan plot: {}", path);
                            manhattan_plot(&results, &path, config.clone())?;
                        }
                        "qq" => {
                            let path = format!("{}.qq.svg", base_path.trim_end_matches(".csv"));
                            eprintln!("Generating QQ plot: {}", path);
                            qq_plot(&results, &path, config.clone())?;
                        }
                        "both" => {
                            let manhattan_path = format!("{}.manhattan.svg", base_path.trim_end_matches(".csv"));
                            let qq_path = format!("{}.qq.svg", base_path.trim_end_matches(".csv"));
                            eprintln!("Generating Manhattan plot: {}", manhattan_path);
                            manhattan_plot(&results, &manhattan_path, config.clone())?;
                            eprintln!("Generating QQ plot: {}", qq_path);
                            qq_plot(&results, &qq_path, config.clone())?;
                        }
                        _ => {
                            eprintln!("Unknown plot type '{}', skipping. Use: manhattan, qq, or both", plot_type);
                        }
                    }
                }
            }
        }
        Commands::Dosage { csv, counts, ref_path, total_path, vcf, chunk_size, threads, output, ploidy, mode, format, compress, verbose } => {
            let fit_mode = match mode.as_str() {
                "auto" => binx_dosage::FitMode::Auto,
                "updog" => binx_dosage::FitMode::Updog,
                "updog-fast" => binx_dosage::FitMode::UpdogFast,
                "updog-exact" => binx_dosage::FitMode::UpdogExact,
                "fast" => binx_dosage::FitMode::Fast,
                "turbo" => binx_dosage::FitMode::Turbo,
                "turboauto" => binx_dosage::FitMode::TurboAuto,
                "turboauto-safe" => binx_dosage::FitMode::TurboAutoSafe,
                _ => {
                    eprintln!("Invalid mode: {}. Use 'auto', 'updog', 'updog-fast', 'updog-exact', 'fast', 'turbo', 'turboauto', or 'turboauto-safe'", mode);
                    std::process::exit(1);
                }
            };
            let output_format = match format.as_str() {
                "matrix" => binx_dosage::OutputFormat::Matrix,
                "stats" => binx_dosage::OutputFormat::Stats,
                "beagle" => binx_dosage::OutputFormat::Beagle,
                "vcf" => binx_dosage::OutputFormat::Vcf,
                "plink" => binx_dosage::OutputFormat::PlinkRaw,
                "gwaspoly" => binx_dosage::OutputFormat::GwasPoly,
                _ => {
                    eprintln!("Invalid format: {}. Use 'matrix', 'stats', 'beagle', 'vcf', 'plink', or 'gwaspoly'", format);
                    std::process::exit(1);
                }
            };
            let compress_mode = match compress.as_str() {
                "none" => binx_dosage::CompressMode::None,
                "gzip" => binx_dosage::CompressMode::Gzip,
                _ => {
                    eprintln!("Invalid compress: {}. Use 'none' or 'gzip'", compress);
                    std::process::exit(1);
                }
            };
            let modes_selected = counts as u8 + (vcf.is_some() as u8) + (csv.is_some() as u8);
            if modes_selected == 0 {
                eprintln!("Provide one of: --csv, --counts (with --ref/--total), or --vcf");
                std::process::exit(1);
            }
            if modes_selected > 1 {
                eprintln!("Inputs are mutually exclusive: choose only one of --csv, --counts, or --vcf");
                std::process::exit(1);
            }

            let input = if let Some(vcf_path) = vcf {
                binx_dosage::InputSource::Vcf { path: vcf_path, chunk_size }
            } else if counts {
                let ref_path = ref_path.unwrap_or_else(|| {
                    eprintln!("--counts requires --ref-path <path>");
                    std::process::exit(1);
                });
                let total_path = total_path.unwrap_or_else(|| {
                    eprintln!("--counts requires --total-path <path>");
                    std::process::exit(1);
                });
                if csv.is_some() {
                    eprintln!("Do not provide --csv when using --counts");
                    std::process::exit(1);
                }
                binx_dosage::InputSource::RefTotalMatrices { ref_path, total_path }
            } else {
                let csv_path = csv.unwrap_or_else(|| {
                    eprintln!("--csv is required unless --counts or --vcf is set");
                    std::process::exit(1);
                });
                binx_dosage::InputSource::TwoLineCsv(csv_path)
            };
            binx_dosage::run_dosage(
                input,
                ploidy,
                fit_mode,
                verbose,
                threads,
                output.as_deref(),
                output_format,
                compress_mode,
            )?;
        }
        Commands::Convert { vcf, output, format, verbose } => {
            match format.as_str() {
                "csv" => convert_vcf_to_csv(&vcf, &output, verbose)?,
                "gwaspoly" => convert_vcf_to_gwaspoly(&vcf, &output, verbose)?,
                _ => {
                    eprintln!("Invalid format: {}. Use 'csv' or 'gwaspoly'", format);
                    std::process::exit(1);
                }
            }
        }
        Commands::Qtl { input, bp_window, output } => {
            run_qtl(input.as_deref(), bp_window, output.as_deref())?;
        }
        Commands::Plot {
            input,
            output,
            plot_type,
            model,
            threshold,
            suggestive,
            title,
            theme,
            width,
            height,
            chromosomes,
            ploidy,
            max_pairs,
            max_loci,
            n_bins,
            r2_threshold,
        } => {
            run_plot(
                &input,
                &output,
                &plot_type,
                model.as_deref(),
                threshold,
                suggestive,
                title,
                &theme,
                width,
                height,
                chromosomes.as_deref(),
                ploidy,
                max_pairs,
                max_loci,
                n_bins,
                r2_threshold,
            )?;
        }
    }

    Ok(())
}

fn parse_csv_list(s: &str) -> Vec<String> {
    s.split(',')
        .map(|v| v.trim())
        .filter(|v| !v.is_empty())
        .map(|v| v.to_string())
        .collect()
}

fn convert_vcf_to_csv(vcf_path: &str, output_path: &str, verbose: bool) -> Result<()> {
    use binx_dosage::io::{stream_vcf_records, VcfRecordCounts};

    let writer = RefCell::new(BufWriter::with_capacity(64 * 1024, File::create(output_path)?));
    let header_written = RefCell::new(false);
    let write_error = RefCell::new(None::<anyhow::Error>);
    let processed = RefCell::new(0usize);
    let last_report = RefCell::new(Instant::now());

    stream_vcf_records(
        vcf_path,
        |names| {
            let mut w = writer.borrow_mut();
            let header = std::iter::once("locus".to_string())
                .chain(names.iter().cloned())
                .collect::<Vec<_>>()
                .join(",");
            if let Err(e) = writeln!(w, "{}", header) {
                *write_error.borrow_mut() = Some(anyhow::anyhow!(e));
            } else {
                *header_written.borrow_mut() = true;
            }
        },
        |rec: VcfRecordCounts| {
            if write_error.borrow().is_some() {
                return;
            }
            if !*header_written.borrow() {
                *write_error.borrow_mut() =
                    Some(anyhow::anyhow!("VCF header with sample names was not found"));
                return;
            }

            let mut w = writer.borrow_mut();
            // Ref line
            if let Err(e) = write!(w, "{}", rec.id) {
                *write_error.borrow_mut() = Some(anyhow::anyhow!(e));
                return;
            }
            for &c in rec.ref_counts.iter() {
                if let Err(e) = write!(w, ",{}", c) {
                    *write_error.borrow_mut() = Some(anyhow::anyhow!(e));
                    return;
                }
            }
            if let Err(e) = writeln!(w) {
                *write_error.borrow_mut() = Some(anyhow::anyhow!(e));
                return;
            }

            // Total line
            if let Err(e) = write!(w, "{}", rec.id) {
                *write_error.borrow_mut() = Some(anyhow::anyhow!(e));
                return;
            }
            for &c in rec.total_counts.iter() {
                if let Err(e) = write!(w, ",{}", c) {
                    *write_error.borrow_mut() = Some(anyhow::anyhow!(e));
                    return;
                }
            }
            if let Err(e) = writeln!(w) {
                *write_error.borrow_mut() = Some(anyhow::anyhow!(e));
                return;
            }

            if verbose && last_report.borrow().elapsed().as_secs_f64() >= 2.0 {
                let mut lr = last_report.borrow_mut();
                *lr = Instant::now();
                let count = {
                    let mut p = processed.borrow_mut();
                    *p += 1;
                    *p
                };
                eprint!("\rConverted {} loci...", count);
                let _ = std::io::stderr().flush();
            } else {
                let mut p = processed.borrow_mut();
                *p += 1;
            }
        },
    )
    .map_err(|e| anyhow::anyhow!("Failed to read VCF: {}", e))?;

    if let Some(e) = write_error.into_inner() {
        return Err(e);
    }

    writer.into_inner().flush()?;
    if verbose {
        let count = processed.into_inner();
        eprintln!("\nWrote {} loci to {}", count, output_path);
    }
    Ok(())
}

fn convert_vcf_to_gwaspoly(vcf_path: &str, output_path: &str, verbose: bool) -> Result<()> {
    use binx_dosage::io::{stream_vcf_gt_dosages, VcfRecordDosage};

    let writer = RefCell::new(BufWriter::with_capacity(64 * 1024, File::create(output_path)?));
    let header_written = RefCell::new(false);
    let write_error = RefCell::new(None::<anyhow::Error>);
    let processed = RefCell::new(0usize);
    let last_report = RefCell::new(Instant::now());

    stream_vcf_gt_dosages(
        vcf_path,
        |names| {
            let mut w = writer.borrow_mut();
            // GWASpoly format header: Marker, Chrom, Position, Sample1, Sample2, ...
            if let Err(e) = write!(w, "Marker\tChrom\tPosition") {
                *write_error.borrow_mut() = Some(anyhow::anyhow!(e));
                return;
            }
            for name in names {
                if let Err(e) = write!(w, "\t{}", name) {
                    *write_error.borrow_mut() = Some(anyhow::anyhow!(e));
                    return;
                }
            }
            if let Err(e) = writeln!(w) {
                *write_error.borrow_mut() = Some(anyhow::anyhow!(e));
            } else {
                *header_written.borrow_mut() = true;
            }
        },
        |rec: VcfRecordDosage| {
            if write_error.borrow().is_some() {
                return;
            }
            if !*header_written.borrow() {
                *write_error.borrow_mut() =
                    Some(anyhow::anyhow!("VCF header with sample names was not found"));
                return;
            }

            let mut w = writer.borrow_mut();
            // Write marker ID, chromosome, position
            if let Err(e) = write!(w, "{}\t{}\t{}", rec.id, rec.chrom, rec.pos) {
                *write_error.borrow_mut() = Some(anyhow::anyhow!(e));
                return;
            }

            // Write dosages (NA for missing)
            for dosage in &rec.dosages {
                let val = match dosage {
                    Some(d) => d.to_string(),
                    None => "NA".to_string(),
                };
                if let Err(e) = write!(w, "\t{}", val) {
                    *write_error.borrow_mut() = Some(anyhow::anyhow!(e));
                    return;
                }
            }
            if let Err(e) = writeln!(w) {
                *write_error.borrow_mut() = Some(anyhow::anyhow!(e));
                return;
            }

            if verbose && last_report.borrow().elapsed().as_secs_f64() >= 2.0 {
                let mut lr = last_report.borrow_mut();
                *lr = Instant::now();
                let count = {
                    let mut p = processed.borrow_mut();
                    *p += 1;
                    *p
                };
                eprint!("\rConverted {} loci...", count);
                let _ = std::io::stderr().flush();
            } else {
                let mut p = processed.borrow_mut();
                *p += 1;
            }
        },
    )
    .map_err(|e| anyhow::anyhow!("Failed to read VCF: {}", e))?;

    if let Some(e) = write_error.into_inner() {
        return Err(e);
    }

    writer.into_inner().flush()?;
    if verbose {
        let count = processed.into_inner();
        eprintln!("\nWrote {} loci to {}", count, output_path);
    }
    Ok(())
}

fn run_qtl(input: Option<&str>, bp_window: Option<u64>, output: Option<&str>) -> Result<()> {
    use std::io::{stdin, BufReader};

    // Load results from file or stdin
    let results = if let Some(path) = input {
        eprintln!("Loading GWAS results from {}...", path);
        gwaspoly_rs::get_qtl_from_file(path, bp_window)?
    } else {
        // Read from stdin
        eprintln!("Reading GWAS results from stdin...");
        let stdin = stdin();
        let reader = BufReader::new(stdin.lock());
        let results = load_marker_results_from_reader(reader)?;
        eprintln!("Loaded {} results", results.len());
        gwaspoly_rs::get_qtl(&results, bp_window)?
    };

    eprintln!("Found {} significant QTLs", results.len());

    // Write output
    if let Some(path) = output {
        gwaspoly_rs::write_qtl_results(&results, path)?;
        eprintln!("QTL results written to {}", path);
    } else {
        // Write to stdout
        println!("marker_id,chrom,pos,model,score,effect,threshold");
        for qtl in &results {
            let effect_str = qtl.effect.map(|e| format!("{:.6}", e)).unwrap_or_else(|| "NA".to_string());
            println!(
                "{},{},{:.0},{},{:.4},{},{}",
                qtl.marker_id,
                qtl.chrom,
                qtl.pos,
                qtl.model,
                qtl.score,
                effect_str,
                qtl.threshold
            );
        }
    }

    Ok(())
}

fn load_marker_results_from_reader<R: BufRead>(reader: R) -> Result<Vec<gwaspoly_rs::MarkerResult>> {
    let mut results = Vec::new();
    let mut lines = reader.lines();

    // Parse header
    let header = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("Empty input"))??;
    let columns: Vec<&str> = header.split(',').map(|s| s.trim()).collect();

    // Find column indices
    let find_col = |name: &str| -> Option<usize> {
        columns.iter().position(|&c| c == name)
    };

    let idx_marker = find_col("marker_id").ok_or_else(|| anyhow::anyhow!("Missing marker_id column"))?;
    let idx_chrom = find_col("chrom");
    let idx_pos = find_col("pos");
    let idx_model = find_col("model").ok_or_else(|| anyhow::anyhow!("Missing model column"))?;
    let idx_score = find_col("score").ok_or_else(|| anyhow::anyhow!("Missing score column"))?;
    let idx_pvalue = find_col("p_value");
    let idx_effect = find_col("effect");
    let idx_nobs = find_col("n_obs");
    let idx_threshold = find_col("threshold");

    for line in lines {
        let line = line?;
        let fields: Vec<&str> = line.split(',').map(|s| s.trim()).collect();

        if fields.len() <= idx_marker {
            continue;
        }

        let parse_optional_f64 = |idx: Option<usize>| -> Option<f64> {
            idx.and_then(|i| fields.get(i))
                .and_then(|s| {
                    if *s == "NA" || s.is_empty() {
                        None
                    } else {
                        s.parse().ok()
                    }
                })
        };

        let parse_optional_str = |idx: Option<usize>| -> Option<String> {
            idx.and_then(|i| fields.get(i))
                .and_then(|s| {
                    if *s == "NA" || s.is_empty() {
                        None
                    } else {
                        Some(s.to_string())
                    }
                })
        };

        let score = fields.get(idx_score)
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.0);

        let p_value = parse_optional_f64(idx_pvalue).unwrap_or_else(|| {
            10_f64.powf(-score)
        });

        let n_obs = idx_nobs
            .and_then(|i| fields.get(i))
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        results.push(gwaspoly_rs::MarkerResult {
            marker_id: fields[idx_marker].to_string(),
            chrom: parse_optional_str(idx_chrom),
            pos: parse_optional_f64(idx_pos),
            model: fields[idx_model].to_string(),
            score,
            p_value,
            effect: parse_optional_f64(idx_effect),
            n_obs,
            threshold: parse_optional_f64(idx_threshold),
        });
    }

    Ok(results)
}

fn run_plot(
    input: &str,
    output: &str,
    plot_type: &str,
    model: Option<&str>,
    threshold: f64,
    suggestive: f64,
    title: Option<String>,
    theme_name: &str,
    width: u32,
    height: u32,
    chromosomes: Option<&str>,
    ploidy: Option<u8>,
    max_pairs: usize,
    max_loci: Option<usize>,
    n_bins: usize,
    r2_threshold: Option<f64>,
) -> Result<()> {
    use binx_plotting::{manhattan_plot, qq_plot, PlotConfig, load_gwas_results, filter_by_model, themes::Theme};
    use binx_plotting::{ld_plot_from_file, LDPlotConfig};

    // Parse theme
    let theme = match theme_name {
        "classic" => Theme::classic(),
        "nature" => Theme::nature(),
        "colorful" => Theme::colorful(),
        "dark" => Theme::dark(),
        "high_contrast" => Theme::high_contrast(),
        _ => {
            eprintln!("Unknown theme '{}', using classic", theme_name);
            Theme::classic()
        }
    };

    // Parse chromosome filter
    let chrom_filter: Option<Vec<String>> = chromosomes.map(|c| {
        c.split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect()
    });

    // Handle LD plot separately (different input format)
    if plot_type == "ld" {
        let ploidy_val = ploidy.ok_or_else(|| {
            anyhow::anyhow!("LD plot requires --ploidy argument")
        })?;

        let ld_config = LDPlotConfig {
            max_pairs,
            max_loci_per_chrom: max_loci,
            n_bins,
            show_points: false,
            r2_threshold,
            chromosomes: chrom_filter,
            plot_config: PlotConfig {
                width,
                height,
                significance_threshold: threshold,
                suggestive_threshold: if suggestive > 0.0 { Some(suggestive) } else { None },
                title,
                theme,
                point_size: 3,
                show_chrom_labels: true,
                chromosomes: None,
            },
        };

        ld_plot_from_file(input, output, ploidy_val, ld_config)?;
        eprintln!("LD plot saved to {}", output);
        return Ok(());
    }

    // For manhattan/qq plots, load GWAS results
    eprintln!("Loading GWAS results from {}...", input);
    let mut results = load_gwas_results(input)?;

    // Filter by model if specified
    if let Some(m) = model {
        let before = results.len();
        results = filter_by_model(&results, m);
        eprintln!("Filtered to model '{}': {} -> {} results", m, before, results.len());
    }

    if results.is_empty() {
        anyhow::bail!("No results to plot after filtering");
    }

    if let Some(ref chroms) = chrom_filter {
        eprintln!("Filtering to chromosomes: {:?}", chroms);
    }

    // Build config
    let config = PlotConfig {
        width,
        height,
        significance_threshold: threshold,
        suggestive_threshold: if suggestive > 0.0 { Some(suggestive) } else { None },
        title,
        theme,
        point_size: 3,
        show_chrom_labels: true,
        chromosomes: chrom_filter,
    };

    // Generate plot
    match plot_type {
        "manhattan" => {
            eprintln!("Generating Manhattan plot...");
            manhattan_plot(&results, output, config)?;
        }
        "qq" => {
            eprintln!("Generating QQ plot...");
            qq_plot(&results, output, config)?;
        }
        _ => {
            anyhow::bail!("Unknown plot type '{}'. Use 'manhattan', 'qq', or 'ld'", plot_type);
        }
    }

    eprintln!("Plot saved to {}", output);
    Ok(())
}
