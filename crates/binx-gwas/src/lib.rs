//! binx-gwas: GWAS orchestration and workflow management for binx.
//!
//! This crate provides the orchestration layer for running GWAS analyses:
//! - Sample and marker alignment between genotype and phenotype data
//! - Multi-method GWAS comparison (gwaspoly, rrblup-assoc, etc.)
//! - LOCO (Leave-One-Chromosome-Out) coordination
//! - Workflow management and result aggregation
//!
//! This crate serves as the integration layer between standalone method crates
//! (gwaspoly-rs, rrblup-rs) and the binx ecosystem.

use anyhow::Result;
use std::str::FromStr;

// Re-export types from gwaspoly-rs that CLI needs
pub use gwaspoly_rs::{GeneActionModel, ThresholdMethod};

/// GWAS method to use for analysis.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum GwasMethod {
    /// GWASpoly method - polyploid GWAS with multiple genetic models
    #[default]
    Gwaspoly,
    // Future: FastGwas, etc.
}

impl FromStr for GwasMethod {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "gwaspoly" => Ok(GwasMethod::Gwaspoly),
            _ => Err(anyhow::anyhow!(
                "Unknown GWAS method '{}'. Available: gwaspoly",
                s
            )),
        }
    }
}

impl std::fmt::Display for GwasMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GwasMethod::Gwaspoly => write!(f, "gwaspoly"),
        }
    }
}

/// Run GWAS analysis using the specified method.
///
/// This is the main orchestration function that routes to the appropriate
/// GWAS implementation based on the selected method.
#[allow(clippy::too_many_arguments)]
pub fn run_gwas(
    geno_path: &str,
    pheno_path: &str,
    trait_name: &str,
    covariate_names: Option<&[String]>,
    kinship_path: Option<&str>,
    allow_missing_samples: bool,
    ploidy: u8,
    models: &[GeneActionModel],
    loco: bool,
    min_maf: f64,
    max_geno_freq: f64,
    out_path: &str,
    use_parallel: bool,
    threshold_method: Option<ThresholdMethod>,
    alpha: f64,
    method: GwasMethod,
    n_pc: usize,
) -> Result<()> {
    match method {
        GwasMethod::Gwaspoly => {
            gwaspoly_rs::gwaspoly(
                geno_path,
                pheno_path,
                trait_name,
                covariate_names,
                kinship_path,
                allow_missing_samples,
                None, // env_column (deprecated)
                None, // env_value (deprecated)
                ploidy,
                models,
                loco,
                min_maf,
                max_geno_freq,
                out_path,
                use_parallel,
                threshold_method,
                alpha,
                n_pc,
            )
        }
    }
}

/// Placeholder for GWAS workflow configuration
#[derive(Debug, Clone, Default)]
pub struct GwasConfig {
    pub parallel: bool,
    pub loco: bool,
    pub min_maf: f64,
    pub max_geno_freq: f64,
}

/// Placeholder for aligned dataset ready for GWAS
#[derive(Debug, Clone)]
pub struct AlignedData {
    pub n_samples: usize,
    pub n_markers: usize,
    pub sample_ids: Vec<String>,
}

/// Align samples between genotype and phenotype data
///
/// This function will handle:
/// - Finding common samples between geno and pheno
/// - Reordering matrices to match
/// - Handling missing data
pub fn align_samples(
    geno_ids: &[String],
    pheno_ids: &[String],
    allow_missing: bool,
) -> Result<Vec<usize>> {
    let mut indices = Vec::new();

    for pheno_id in pheno_ids {
        if let Some(idx) = geno_ids.iter().position(|g| g == pheno_id) {
            indices.push(idx);
        } else if !allow_missing {
            anyhow::bail!("Sample '{}' in phenotype not found in genotype", pheno_id);
        }
    }

    Ok(indices)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_align_samples_exact_match() {
        let geno = vec!["A".into(), "B".into(), "C".into()];
        let pheno = vec!["A".into(), "B".into(), "C".into()];

        let indices = align_samples(&geno, &pheno, false).unwrap();
        assert_eq!(indices, vec![0, 1, 2]);
    }

    #[test]
    fn test_align_samples_reorder() {
        let geno = vec!["A".into(), "B".into(), "C".into()];
        let pheno = vec!["C".into(), "A".into(), "B".into()];

        let indices = align_samples(&geno, &pheno, false).unwrap();
        assert_eq!(indices, vec![2, 0, 1]);
    }

    #[test]
    fn test_align_samples_missing_allowed() {
        let geno = vec!["A".into(), "B".into()];
        let pheno = vec!["A".into(), "B".into(), "C".into()];

        let indices = align_samples(&geno, &pheno, true).unwrap();
        assert_eq!(indices, vec![0, 1]); // C is skipped
    }
}
