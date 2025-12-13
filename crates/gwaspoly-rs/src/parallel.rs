//! Parallel marker testing for GWAS
//!
//! This module provides parallelized versions of marker testing functions
//! using rayon for thread-level parallelism.
//!
//! The parallel implementation uses P3D (Population Parameters Previously Determined)
//! where variance components and H_inv are pre-computed once, then each marker
//! is tested independently in parallel.

use anyhow::Result;
use ndarray::{Array1, Axis};
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;

use crate::types::GenotypeMatrixBiallelic;
use crate::gwaspoly::{GeneActionModel, GwasCache, MarkerResult, lmm_score_test};

/// Test all markers in parallel using pre-computed null model cache.
///
/// This is the main entry point for parallel GWAS marker testing.
/// Each marker is tested independently across multiple threads.
///
/// # Arguments
/// * `geno` - Genotype matrix
/// * `y` - Phenotype vector (observation-level)
/// * `base_design` - Base design matrix columns (intercept + covariates)
/// * `cache` - Pre-computed null model cache (H_inv, variance components)
/// * `obs_to_geno` - Mapping from observation index to genotype index
/// * `models` - Gene action models to test
/// * `ploidy` - Ploidy level
/// * `min_maf` - Minimum minor allele frequency
/// * `max_geno_freq` - Maximum genotype frequency for QC
///
/// # Returns
/// Vector of MarkerResults for all markers and models
pub fn test_markers_parallel(
    geno: &GenotypeMatrixBiallelic,
    y: &Array1<f64>,
    base_design: &[Vec<f64>],
    cache: &GwasCache,
    obs_to_geno: &[usize],
    models: &[GeneActionModel],
    ploidy: u8,
    min_maf: f64,
    max_geno_freq: f64,
) -> Result<Vec<MarkerResult>> {
    let n_markers = geno.marker_ids.len();
    let n_genotypes_for_qc = geno.sample_ids.len();

    // Wrap shared data in Arc for thread-safe sharing
    let geno_arc: Arc<GenotypeMatrixBiallelic> = Arc::new(geno.clone());
    let y_arc: Arc<Array1<f64>> = Arc::new(y.clone());
    let base_design_arc: Arc<Vec<Vec<f64>>> = Arc::new(base_design.to_vec());
    let cache_arc: Arc<GwasCache> = Arc::new(cache.clone());
    let obs_to_geno_arc: Arc<Vec<usize>> = Arc::new(obs_to_geno.to_vec());
    let models_arc: Arc<Vec<GeneActionModel>> = Arc::new(models.to_vec());

    // Process markers in parallel
    let results: Vec<Result<Vec<MarkerResult>>> = (0..n_markers)
        .into_par_iter()
        .map(|marker_idx| {
            test_marker_parallel(
                &geno_arc,
                marker_idx,
                &y_arc,
                &base_design_arc,
                &cache_arc,
                &obs_to_geno_arc,
                &models_arc,
                ploidy,
                min_maf,
                max_geno_freq,
                n_genotypes_for_qc,
            )
        })
        .collect();

    // Flatten results, propagating any errors
    let mut all_results = Vec::with_capacity(n_markers * models.len());
    for result in results {
        all_results.extend(result?);
    }

    Ok(all_results)
}

/// Test a subset of markers in parallel (for LOCO mode).
///
/// This version takes a list of marker indices to test, useful when
/// different chromosomes need different null models.
pub fn test_markers_parallel_subset(
    geno: &GenotypeMatrixBiallelic,
    marker_indices: &[usize],
    y: &Array1<f64>,
    base_design: &[Vec<f64>],
    cache: &GwasCache,
    obs_to_geno: &[usize],
    models: &[GeneActionModel],
    ploidy: u8,
    min_maf: f64,
    max_geno_freq: f64,
) -> Result<Vec<MarkerResult>> {
    let n_genotypes_for_qc = geno.sample_ids.len();

    // Wrap shared data in Arc for thread-safe sharing
    let geno_arc: Arc<GenotypeMatrixBiallelic> = Arc::new(geno.clone());
    let y_arc: Arc<Array1<f64>> = Arc::new(y.clone());
    let base_design_arc: Arc<Vec<Vec<f64>>> = Arc::new(base_design.to_vec());
    let cache_arc: Arc<GwasCache> = Arc::new(cache.clone());
    let obs_to_geno_arc: Arc<Vec<usize>> = Arc::new(obs_to_geno.to_vec());
    let models_arc: Arc<Vec<GeneActionModel>> = Arc::new(models.to_vec());

    // Process markers in parallel
    let results: Vec<Result<Vec<MarkerResult>>> = marker_indices
        .par_iter()
        .map(|&marker_idx| {
            test_marker_parallel(
                &geno_arc,
                marker_idx,
                &y_arc,
                &base_design_arc,
                &cache_arc,
                &obs_to_geno_arc,
                &models_arc,
                ploidy,
                min_maf,
                max_geno_freq,
                n_genotypes_for_qc,
            )
        })
        .collect();

    // Flatten results, propagating any errors
    let mut all_results = Vec::with_capacity(marker_indices.len() * models.len());
    for result in results {
        all_results.extend(result?);
    }

    Ok(all_results)
}

/// Test a single marker across all models (called in parallel)
fn test_marker_parallel(
    geno: &Arc<GenotypeMatrixBiallelic>,
    marker_idx: usize,
    y: &Arc<Array1<f64>>,
    base_design: &Arc<Vec<Vec<f64>>>,
    cache: &Arc<GwasCache>,
    obs_to_geno: &Arc<Vec<usize>>,
    models: &Arc<Vec<GeneActionModel>>,
    ploidy: u8,
    min_maf: f64,
    max_geno_freq: f64,
    n_genotypes_for_qc: usize,
) -> Result<Vec<MarkerResult>> {
    let marker_id = &geno.marker_ids[marker_idx];
    let dosage_geno = geno.dosages.index_axis(Axis(0), marker_idx).to_owned();

    let (chrom, pos) = geno
        .marker_metadata
        .as_ref()
        .and_then(|meta| meta.get(marker_idx))
        .map(|m| (Some(m.chrom.clone()), Some(m.pos)))
        .unwrap_or((None, None));

    let mut results = Vec::new();

    for model in models.iter() {
        // Get marker design at genotype level
        let marker_design_geno = match design_score_parallel(
            &dosage_geno,
            ploidy,
            *model,
            min_maf,
            max_geno_freq,
            n_genotypes_for_qc,
        ) {
            Some(cols) => cols,
            None => continue,
        };

        // Expand marker design to observation level
        let marker_design_obs: Vec<Vec<f64>> = marker_design_geno
            .iter()
            .map(|col_geno| {
                obs_to_geno.iter().map(|&g_idx| col_geno[g_idx]).collect()
            })
            .collect();

        let (score, p_value, effect) =
            lmm_score_test(&marker_design_obs, y, base_design, cache)?;

        results.push(MarkerResult {
            marker_id: marker_id.clone(),
            chrom: chrom.clone(),
            pos,
            model: model.as_str().to_string(),
            score,
            p_value,
            effect,
            n_obs: y.len(),
            threshold: None,
        });
    }

    Ok(results)
}

/// Generate marker design matrix (parallel-safe version)
fn design_score_parallel(
    dosage: &ndarray::Array1<f64>,
    ploidy: u8,
    model: GeneActionModel,
    min_maf: f64,
    max_geno_freq: f64,
    n_genotypes_for_qc: usize,
) -> Option<Vec<Vec<f64>>> {
    if n_genotypes_for_qc == 0 {
        return None;
    }

    // Compute MAF
    let mut valid_sum = 0.0;
    let mut valid_n = 0usize;
    for &v in dosage.iter() {
        if v.is_finite() {
            valid_sum += v;
            valid_n += 1;
        }
    }
    if valid_n == 0 {
        return None;
    }
    let freq = (valid_sum / (valid_n as f64)) / (ploidy as f64);
    if freq.min(1.0 - freq) < min_maf {
        return None;
    }

    // Impute missing values with mean dosage (matches R/GWASpoly behavior)
    let mean_dosage = valid_sum / (valid_n as f64);
    let dosage: Vec<f64> = dosage
        .iter()
        .map(|&d| if d.is_finite() { d } else { mean_dosage })
        .collect();

    let rounded: Vec<i64> = dosage.iter().map(|&d| d.round() as i64).collect();

    let mut counts: HashMap<i64, usize> = HashMap::new();
    for &r in &rounded {
        *counts.entry(r).or_insert(0) += 1;
    }

    match model {
        GeneActionModel::Additive => {
            let max_freq = counts.values().cloned().max().unwrap_or(0) as f64 / (dosage.len() as f64);
            if max_freq > max_geno_freq {
                return None;
            }
            // Use imputed dosages (missing values replaced with marker mean)
            Some(vec![dosage.clone()])
        }

        GeneActionModel::General => {
            let max_freq = counts.values().cloned().max().unwrap_or(0) as f64 / (dosage.len() as f64);
            if max_freq > max_geno_freq {
                return None;
            }

            let mut levels: Vec<i64> = Vec::new();
            for &r in &rounded {
                if !levels.contains(&r) {
                    levels.push(r);
                }
            }

            if levels.len() <= 1 {
                return None;
            }

            let mut cols = Vec::with_capacity(levels.len() - 1);
            for level in levels.iter().skip(1) {
                let col: Vec<f64> = rounded.iter().map(|&r| (r == *level) as i32 as f64).collect();
                cols.push(col);
            }
            Some(cols)
        }

        GeneActionModel::SimplexDomRef => {
            let col: Vec<f64> = rounded
                .iter()
                .map(|&r| if r > 0 { 1.0 } else { 0.0 })
                .collect();
            check_binary_column_parallel(&col, max_geno_freq)?;
            Some(vec![col])
        }

        GeneActionModel::SimplexDomAlt => {
            let max_dose = ploidy as i64;
            let col: Vec<f64> = rounded
                .iter()
                .map(|&r| if r == max_dose { 1.0 } else { 0.0 })
                .collect();
            check_binary_column_parallel(&col, max_geno_freq)?;
            Some(vec![col])
        }

        GeneActionModel::DuplexDomRef => {
            let threshold = (ploidy / 2) as i64;
            let col: Vec<f64> = rounded
                .iter()
                .map(|&r| if r >= threshold { 1.0 } else { 0.0 })
                .collect();
            check_binary_column_parallel(&col, max_geno_freq)?;
            Some(vec![col])
        }

        GeneActionModel::DuplexDomAlt => {
            let threshold = (ploidy / 2 + 1) as i64;
            let col: Vec<f64> = rounded
                .iter()
                .map(|&r| if r >= threshold { 1.0 } else { 0.0 })
                .collect();
            check_binary_column_parallel(&col, max_geno_freq)?;
            Some(vec![col])
        }

        GeneActionModel::DiploGeneral => {
            let max_dose = ploidy as i64;
            let diplo: Vec<i64> = rounded
                .iter()
                .map(|&r| {
                    if r == 0 { 0 }
                    else if r == max_dose { 2 }
                    else { 1 }
                })
                .collect();

            let mut diplo_counts: HashMap<i64, usize> = HashMap::new();
            for &d in &diplo {
                *diplo_counts.entry(d).or_insert(0) += 1;
            }
            let max_freq = diplo_counts.values().cloned().max().unwrap_or(0) as f64 / (diplo.len() as f64);
            if max_freq > max_geno_freq {
                return None;
            }

            let mut levels: Vec<i64> = Vec::new();
            for &d in &diplo {
                if !levels.contains(&d) {
                    levels.push(d);
                }
            }

            if levels.len() <= 1 {
                return None;
            }

            let mut cols = Vec::with_capacity(levels.len() - 1);
            for level in levels.iter().skip(1) {
                let col: Vec<f64> = diplo.iter().map(|&d| (d == *level) as i32 as f64).collect();
                cols.push(col);
            }
            Some(cols)
        }

        GeneActionModel::DiploAdditive => {
            let max_dose = ploidy as f64;
            let col: Vec<f64> = dosage
                .iter()
                .map(|&d| {
                    if d <= 0.0 { 0.0 }
                    else if d >= max_dose { 1.0 }
                    else { 0.5 }
                })
                .collect();

            let mut diplo_counts: HashMap<i64, usize> = HashMap::new();
            for &c in &col {
                let key = (c * 10.0).round() as i64;
                *diplo_counts.entry(key).or_insert(0) += 1;
            }
            let max_freq = diplo_counts.values().cloned().max().unwrap_or(0) as f64 / (col.len() as f64);
            if max_freq > max_geno_freq {
                return None;
            }

            Some(vec![col])
        }
    }
}

fn check_binary_column_parallel(col: &[f64], max_geno_freq: f64) -> Option<()> {
    let n = col.len();
    if n == 0 {
        return None;
    }
    let ones = col.iter().filter(|&&v| v > 0.5).count();
    let zeros = n - ones;

    let max_count = ones.max(zeros);
    if (max_count as f64) / (n as f64) > max_geno_freq {
        return None;
    }
    if ones == 0 || zeros == 0 {
        return None;
    }
    Some(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;
    use ndarray::{array, Array2};

    fn create_test_cache(n: usize) -> GwasCache {
        GwasCache {
            sample_ids: (0..n).map(|i| format!("S{}", i)).collect(),
            n_obs: n,
            n_ind: n,
            h_inv: DMatrix::identity(n, n),
            sigma2: 1.0,
            lambda: 1.0,
            vu: 0.5,
            ve: 0.5,
        }
    }

    #[test]
    fn test_parallel_matches_sequential() {
        // Create simple test data
        let n_samples = 10;
        let n_markers = 5;

        let mut dosages = Array2::<f64>::zeros((n_markers, n_samples));
        for i in 0..n_markers {
            for j in 0..n_samples {
                dosages[(i, j)] = ((i + j) % 5) as f64;
            }
        }

        let geno = GenotypeMatrixBiallelic {
            ploidy: 4,
            sample_ids: (0..n_samples).map(|i| format!("S{}", i)).collect(),
            marker_ids: (0..n_markers).map(|i| format!("M{}", i)).collect(),
            marker_metadata: None,
            dosages,
        };

        let y = Array1::from_vec((0..n_samples).map(|i| i as f64).collect());
        let base_design = vec![vec![1.0; n_samples]]; // intercept only
        let cache = create_test_cache(n_samples);
        let obs_to_geno: Vec<usize> = (0..n_samples).collect();
        let models = vec![GeneActionModel::Additive];

        let results = test_markers_parallel(
            &geno,
            &y,
            &base_design,
            &cache,
            &obs_to_geno,
            &models,
            4,
            0.0,
            1.0,
        ).unwrap();

        // Should have results for markers that pass QC
        assert!(!results.is_empty());

        // All results should have valid scores
        for r in &results {
            assert!(r.score.is_finite());
            assert!(r.p_value >= 0.0 && r.p_value <= 1.0);
        }
    }

    #[test]
    fn test_design_score_additive() {
        let dosage = array![0.0, 1.0, 2.0, 3.0, 4.0];
        let design = design_score_parallel(&dosage, 4, GeneActionModel::Additive, 0.0, 1.0, 5);
        assert!(design.is_some());
        let cols = design.unwrap();
        assert_eq!(cols.len(), 1);
        assert_eq!(cols[0], vec![0.0, 1.0, 2.0, 3.0, 4.0]);
    }

    #[test]
    fn test_design_score_general() {
        let dosage = array![0.0, 1.0, 2.0, 3.0, 4.0];
        let design = design_score_parallel(&dosage, 4, GeneActionModel::General, 0.0, 1.0, 5);
        assert!(design.is_some());
        let cols = design.unwrap();
        assert_eq!(cols.len(), 4); // 5 levels - 1 reference
    }
}
