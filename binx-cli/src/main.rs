use anyhow::Result;
use clap::{Parser, Subcommand};

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
    /// Biallelic GWASpoly-style GWAS (supports LM and LMM; use --kinship for LMM)
    Gwas {
        /// Genotype dosage file (biallelic; markers x samples)
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

        /// Optional TSV with PCs (sample_id + PC1..PCn)
        #[arg(long)]
        pcs: Option<String>,

        /// Optional kinship matrix TSV (sample_id + samples)
        #[arg(long)]
        kinship: Option<String>,

        /// Allow dropping genotype samples that lack phenotypes
        #[arg(long, default_value_t = false)]
        allow_missing_samples: bool,

        /// Optional environment column name to filter phenotype rows (e.g., env)
        #[arg(long)]
        env_column: Option<String>,

        /// Environment value to keep (used with --env-column)
        #[arg(long)]
        env_value: Option<String>,

        /// Ploidy (e.g., 2, 4, 6)
        #[arg(long)]
        ploidy: u8,

        /// Gene action model (additive, general, etc.)
        #[arg(long, default_value = "additive")]
        model: String,

        /// Output results TSV
        #[arg(long)]
        out: String,
    },

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
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Gwas {
        geno,
        pheno,
        trait_name,
        covariates,
        pcs,
        kinship,
        allow_missing_samples,
        env_column,
        env_value,
        ploidy,
        model,
        out,
    } => {
        let covariate_list = covariates.as_deref().map(parse_csv_list);
        binx_gwas::run_gwas(
            &geno,
            &pheno,
            &trait_name,
            covariate_list.as_deref(),
            pcs.as_deref(),
            kinship.as_deref(),
            allow_missing_samples,
            env_column.as_deref(),
            env_value.as_deref(),
            ploidy,
            &model,
            &out,
        )?;
    }
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
            )?;
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
