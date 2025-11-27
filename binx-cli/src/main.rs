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

    /// Compute kinship matrix (VanRaden additive) from biallelic dosages
    Kinship {
        /// Genotype dosage file (biallelic; markers x samples)
        #[arg(long)]
        geno: String,

        /// Ploidy (e.g., 2, 4, 6)
        #[arg(long)]
        ploidy: u8,

        /// Output TSV path for kinship matrix
        #[arg(long)]
        out: String,
    },

    /// Genotype dosage estimation from sequencing read counts
    Dosage {
        /// CSV file with alternating lines of Ref and Total counts per locus
        #[arg(long)]
        csv: String,

        /// Ploidy (e.g., 2, 4, 6)
        #[arg(long)]
        ploidy: usize,

        /// Multi-start optimization mode (auto, updog, fast)
        #[arg(long, default_value = "auto")]
        mode: String,

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
        Commands::Kinship { geno, ploidy, out } => {
            binx_kinship::run_kinship(&geno, ploidy, &out)?;
        }
        Commands::Dosage { csv, ploidy, mode, verbose } => {
            let fit_mode = match mode.as_str() {
                "auto" => binx_dosage::FitMode::Auto,
                "updog" => binx_dosage::FitMode::Updog,
                "fast" => binx_dosage::FitMode::Fast,
                _ => {
                    eprintln!("Invalid mode: {}. Use 'auto', 'updog', or 'fast'", mode);
                    std::process::exit(1);
                }
            };
            binx_dosage::run_dosage(&csv, ploidy, fit_mode, verbose)?;
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
