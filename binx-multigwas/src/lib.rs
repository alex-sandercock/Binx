use anyhow::Result;
use std::path::Path;

/// Placeholder for future multiallelic GWAS implementation.
pub fn run_multigwas<P: AsRef<Path>>(
    geno_path: P,
    pheno_path: P,
    trait_name: &str,
    kinship_path: P,
    pcs_path: Option<P>,
    ploidy: u8,
    encoding: &str,
    model: &str,
    out_path: P,
) -> Result<()> {
    let _ = (
        geno_path,
        pheno_path,
        trait_name,
        kinship_path,
        pcs_path,
        ploidy,
        encoding,
        model,
        out_path,
    );
    Err(anyhow::anyhow!(
        "binx multigwas not yet implemented (stub placeholder)"
    ))
}
