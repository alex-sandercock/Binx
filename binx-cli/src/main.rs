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
    /// Biallelic GWASpoly-style GWAS (Phase 2: LM only, no K)
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
