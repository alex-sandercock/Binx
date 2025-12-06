use anyhow::Result;
use clap::{Parser, Subcommand};
use std::cell::RefCell;
use std::fs::File;
use std::io::{BufWriter, Write};
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

    /// GWASpoly-style GWAS for polyploids with multiple genetic models
    Gwaspoly {
        /// Genotype dosage file (biallelic; markers x samples with chr/pos)
        #[arg(long)]
        geno: String,

        /// Phenotype file (sample_id + traits)
        #[arg(long)]
        pheno: String,

        /// Trait name to analyze
        #[arg(long)]
        trait_name: String,

        /// Comma-separated covariate names from phenotype file
        #[arg(long)]
        covariates: Option<String>,

        /// Optional kinship matrix TSV (sample_id + samples)
        #[arg(long)]
        kinship: Option<String>,

        /// Allow dropping genotype samples that lack phenotypes
        #[arg(long, default_value_t = false)]
        allow_missing_samples: bool,

        /// Optional environment column name to filter phenotype rows
        #[arg(long)]
        env_column: Option<String>,

        /// Environment value to keep (used with --env-column)
        #[arg(long)]
        env_value: Option<String>,

        /// Ploidy (e.g., 2, 4, 6)
        #[arg(long)]
        ploidy: u8,

        /// Comma-separated gene action models (additive,general,1-dom-ref,1-dom-alt,2-dom-ref,2-dom-alt,diplo-general,diplo-additive)
        #[arg(long, default_value = "additive,general")]
        models: String,

        /// Use LOCO (Leave-One-Chromosome-Out) kinship matrices
        #[arg(long, default_value_t = false)]
        loco: bool,

        /// Minimum minor allele frequency (0.0-0.5)
        #[arg(long, default_value = "0.0")]
        min_maf: f64,

        /// Maximum genotype frequency for QC (0.0-1.0, 0 for auto)
        #[arg(long, default_value = "0.0")]
        max_geno_freq: f64,

        /// Output results CSV
        #[arg(long)]
        out: String,

        /// Use parallel marker testing (experimental)
        #[arg(long, default_value_t = false)]
        parallel: bool,

        /// Generate plots after GWAS (manhattan, qq, or both)
        #[arg(long)]
        plot: Option<String>,

        /// Output path for plots (default: <out>.manhattan.svg, <out>.qq.svg)
        #[arg(long)]
        plot_output: Option<String>,
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

    /// Generate GWAS plots (Manhattan, QQ) from results
    Plot {
        /// Input GWAS results CSV file
        #[arg(long)]
        input: String,

        /// Output file path (extension determines format: .svg or .png)
        #[arg(long)]
        output: String,

        /// Plot type: manhattan or qq
        #[arg(long, default_value = "manhattan")]
        plot_type: String,

        /// Filter to specific model (e.g., additive, general)
        #[arg(long)]
        model: Option<String>,

        /// Significance threshold as -log10(p), default 5.0
        #[arg(long, default_value = "5.0")]
        threshold: f64,

        /// Suggestive threshold as -log10(p), default 3.0 (use 0 to disable)
        #[arg(long, default_value = "3.0")]
        suggestive: f64,

        /// Plot title
        #[arg(long)]
        title: Option<String>,

        /// Color theme: classic, nature, colorful, dark, high_contrast
        #[arg(long, default_value = "classic")]
        theme: String,

        /// Plot width in pixels
        #[arg(long, default_value = "1200")]
        width: u32,

        /// Plot height in pixels
        #[arg(long, default_value = "600")]
        height: u32,
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
        Commands::Gwaspoly {
            geno,
            pheno,
            trait_name,
            covariates,
            kinship,
            allow_missing_samples,
            env_column,
            env_value,
            ploidy,
            models,
            loco,
            min_maf,
            max_geno_freq,
            out,
            parallel,
            plot,
            plot_output,
        } => {
            let covariate_list = covariates.as_deref().map(parse_csv_list);
            let model_list: Vec<gwaspoly_rs::GeneActionModel> = models
                .split(',')
                .map(|s| s.trim())
                .filter(|s| !s.is_empty())
                .map(|s| {
                    gwaspoly_rs::GeneActionModel::from_str(s).unwrap_or_else(|e| {
                        eprintln!("Invalid model '{}': {}", s, e);
                        std::process::exit(1);
                    })
                })
                .collect();

            if model_list.is_empty() {
                eprintln!("No valid models specified");
                std::process::exit(1);
            }

            if parallel {
                eprintln!("Using parallel marker testing...");
            }
            gwaspoly_rs::run_gwaspoly(
                &geno,
                &pheno,
                &trait_name,
                covariate_list.as_deref(),
                kinship.as_deref(),
                allow_missing_samples,
                env_column.as_deref(),
                env_value.as_deref(),
                ploidy,
                &model_list,
                loco,
                min_maf,
                max_geno_freq,
                &out,
                parallel,
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
        } => {
            run_plot(&input, &output, &plot_type, model.as_deref(), threshold, suggestive, title, &theme, width, height)?;
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
) -> Result<()> {
    use binx_plotting::{manhattan_plot, qq_plot, PlotConfig, load_gwas_results, filter_by_model, themes::Theme};

    // Load results
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
            anyhow::bail!("Unknown plot type '{}'. Use 'manhattan' or 'qq'", plot_type);
        }
    }

    eprintln!("Plot saved to {}", output);
    Ok(())
}
