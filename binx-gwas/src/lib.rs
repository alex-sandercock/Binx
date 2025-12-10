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
