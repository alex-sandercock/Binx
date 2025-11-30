//! gwaspoly.rs: Main GWASpoly function implementation
//!
//! This module provides a faithful Rust implementation of R/GWASpoly's GWASpoly() function
//! for genome-wide association studies in polyploid organisms.
//!
//! ## Features
//! - Eight genetic models (additive, general, 1-dom-ref, 1-dom-alt, 2-dom-ref, 2-dom-alt, diplo-general, diplo-additive)
//! - LOCO (Leave-One-Chromosome-Out) support
//! - Multiple model testing per marker
//! - Q+K mixed model framework
//! - Support for repeated observations per genotype (multi-environment trials)

use anyhow::{anyhow, Result};
use binx_core::{
    load_genotypes_biallelic_from_tsv, load_kinship_from_tsv,
    load_phenotypes_from_tsv, GenotypeMatrixBiallelic, KinshipMatrix,
};
use binx_kinship::compute_kinship_gwaspoly;
use nalgebra::{DMatrix, DVector};
use ndarray::{Array1, Array2, Axis};
use rrblup_rs::{mixed_solve_new, MixedSolveOptions};

use crate::set_k::compute_loco_kinship_gwaspoly;
use statrs::distribution::{ContinuousCDF, FisherSnedecor};
use std::collections::HashMap;

use crate::parallel;

/// Cached quantities from fitting the null mixed model.
/// Used for efficient marker testing via P3D (population parameters previously determined).
#[derive(Debug, Clone)]
pub struct GwasCache {
    pub sample_ids: Vec<String>,
    pub n_obs: usize,
    pub n_ind: usize,
    pub h_inv: DMatrix<f64>,
    pub sigma2: f64,
    pub lambda: f64,
    pub vu: f64,
    pub ve: f64,
}

/// Gene action models supported by GWASpoly
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum GeneActionModel {
    /// Additive: effect proportional to allele dosage
    Additive,
    /// General: no constraints, separate effect per dosage class
    General,
    /// Simplex dominant (reference): 0 vs (1,2,...,ploidy)
    SimplexDomRef,
    /// Simplex dominant (alternate): (0,1,...,ploidy-1) vs ploidy
    SimplexDomAlt,
    /// Duplex dominant (reference): (0,1) vs (2,...,ploidy) for tetraploid
    DuplexDomRef,
    /// Duplex dominant (alternate): (0,...,ploidy-2) vs (ploidy-1,ploidy) for tetraploid
    DuplexDomAlt,
    /// Diploidized general: heterozygotes collapsed, homozygotes separate
    DiploGeneral,
    /// Diploidized additive: heterozygotes halfway between homozygotes
    DiploAdditive,
}

impl GeneActionModel {
    /// Parse model name from string
    pub fn from_str(s: &str) -> Result<Self> {
        match s.to_lowercase().as_str() {
            "additive" => Ok(GeneActionModel::Additive),
            "general" => Ok(GeneActionModel::General),
            "1-dom-ref" | "simplex-dom-ref" | "1-dom" => Ok(GeneActionModel::SimplexDomRef),
            "1-dom-alt" | "simplex-dom-alt" => Ok(GeneActionModel::SimplexDomAlt),
            "2-dom-ref" | "duplex-dom-ref" | "2-dom" => Ok(GeneActionModel::DuplexDomRef),
            "2-dom-alt" | "duplex-dom-alt" => Ok(GeneActionModel::DuplexDomAlt),
            "diplo-general" => Ok(GeneActionModel::DiploGeneral),
            "diplo-additive" => Ok(GeneActionModel::DiploAdditive),
            other => Err(anyhow!("Unknown gene action model: {}", other)),
        }
    }

    /// Get string representation
    pub fn as_str(&self) -> &'static str {
        match self {
            GeneActionModel::Additive => "additive",
            GeneActionModel::General => "general",
            GeneActionModel::SimplexDomRef => "1-dom-ref",
            GeneActionModel::SimplexDomAlt => "1-dom-alt",
            GeneActionModel::DuplexDomRef => "2-dom-ref",
            GeneActionModel::DuplexDomAlt => "2-dom-alt",
            GeneActionModel::DiploGeneral => "diplo-general",
            GeneActionModel::DiploAdditive => "diplo-additive",
        }
    }

    /// Get all models for multi-model testing
    pub fn all_models() -> Vec<GeneActionModel> {
        vec![
            GeneActionModel::Additive,
            GeneActionModel::General,
            GeneActionModel::SimplexDomRef,
            GeneActionModel::SimplexDomAlt,
            GeneActionModel::DuplexDomRef,
            GeneActionModel::DuplexDomAlt,
            GeneActionModel::DiploGeneral,
            GeneActionModel::DiploAdditive,
        ]
    }
}

/// Result for a single marker test
#[derive(Debug, Clone)]
pub struct MarkerResult {
    pub marker_id: String,
    pub chrom: Option<String>,
    pub pos: Option<f64>,
    pub model: String,
    pub score: f64,
    pub p_value: f64,
    pub effect: Option<f64>,
    pub n_obs: usize,
}

/// Fit null mixed model using the validated rrblup-rs mixed_solve implementation.
///
/// This wraps rrblup_rs::mixed_solve_new (faithful R/rrBLUP::mixed.solve implementation)
/// for use in GWAS marker testing.
pub(crate) fn fit_null_model(
    y: &Array1<f64>,
    x0: &Array2<f64>,
    kinship: &KinshipMatrix,
    z: Option<&Array2<f64>>,
    obs_ids: Option<&[String]>,
) -> Result<GwasCache> {
    let n = y.len();
    let m = kinship.matrix.nrows();

    // Convert y to slice
    let y_slice: Vec<f64> = y.iter().cloned().collect();

    // Convert X to DMatrix
    let x_dmat = DMatrix::from_fn(x0.nrows(), x0.ncols(), |i, j| x0[(i, j)]);

    // Convert K to DMatrix
    let k_dmat = DMatrix::from_fn(m, m, |i, j| kinship.matrix[(i, j)]);

    // Convert Z to DMatrix if provided
    let z_dmat = z.map(|z_arr| {
        DMatrix::from_fn(z_arr.nrows(), z_arr.ncols(), |i, j| z_arr[(i, j)])
    });

    // Call validated mixed_solve_new with return_hinv=true
    let opts = MixedSolveOptions {
        return_hinv: true,
        ..Default::default()
    };

    let result = mixed_solve_new(
        &y_slice,
        z_dmat.as_ref(),
        Some(&k_dmat),
        Some(&x_dmat),
        Some(opts),
    )?;

    let h_inv = result.hinv
        .ok_or_else(|| anyhow!("H_inv not computed by mixed_solve"))?;

    let ids = obs_ids.map(|s| s.to_vec())
        .unwrap_or_else(|| kinship.sample_ids.clone());

    Ok(GwasCache {
        sample_ids: ids,
        n_obs: n,
        n_ind: m,
        h_inv,
        sigma2: result.ve,
        lambda: result.ve / result.vu.max(1e-10),
        vu: result.vu,
        ve: result.ve,
    })
}

/// Main GWASpoly analysis function
///
/// This function faithfully implements R/GWASpoly's approach:
/// - Supports repeated observations per genotype (multi-environment trials)
/// - Uses Z incidence matrix to map observations to genotypes
/// - Fixed effects (like environment) are included in the model, not filtered
pub fn run_gwaspoly(
    geno_path: &str,
    pheno_path: &str,
    trait_name: &str,
    covariate_names: Option<&[String]>,
    kinship_path: Option<&str>,
    allow_missing_samples: bool,
    _env_column: Option<&str>,  // Deprecated: use covariate_names instead
    _env_value: Option<&str>,   // Deprecated: include all data, use fixed effects
    ploidy: u8,
    models: &[GeneActionModel],
    loco: bool,
    min_maf: f64,
    max_geno_freq: f64,
    out_path: &str,
    use_parallel: bool,
) -> Result<()> {
    // Load data - don't filter by environment, include all observations
    let geno_full = load_genotypes_biallelic_from_tsv(geno_path, ploidy)?;
    let pheno = load_phenotypes_from_tsv(pheno_path)?;

    // Get unique genotype IDs that have phenotype data (following R/GWASpoly)
    // R/GWASpoly subsets genotypes to only those with phenotype observations
    let pheno_gid_set: std::collections::HashSet<&str> = pheno.sample_ids
        .iter()
        .map(|s| s.as_str())
        .collect();

    // Find which genotype samples have phenotype data and their indices
    let mut geno_keep_indices: Vec<usize> = Vec::new();
    let mut geno_keep_ids: Vec<String> = Vec::new();
    for (i, gid) in geno_full.sample_ids.iter().enumerate() {
        if pheno_gid_set.contains(gid.as_str()) {
            geno_keep_indices.push(i);
            geno_keep_ids.push(gid.clone());
        }
    }

    if geno_keep_indices.is_empty() {
        return Err(anyhow!("No overlapping sample_ids between genotype and phenotype files"));
    }

    // Subset genotype matrix to only samples with phenotype data
    let n_markers = geno_full.marker_ids.len();
    let n_ind = geno_keep_indices.len();
    let mut dosages_subset = Array2::<f64>::zeros((n_markers, n_ind));
    for (new_j, &old_j) in geno_keep_indices.iter().enumerate() {
        for i in 0..n_markers {
            dosages_subset[(i, new_j)] = geno_full.dosages[(i, old_j)];
        }
    }

    let geno = GenotypeMatrixBiallelic {
        ploidy: geno_full.ploidy,
        sample_ids: geno_keep_ids,
        marker_ids: geno_full.marker_ids.clone(),
        marker_metadata: geno_full.marker_metadata.clone(),
        dosages: dosages_subset,
    };

    eprintln!("DEBUG: Subsetted genotypes from {} to {} samples (matching phenotype)",
              geno_full.sample_ids.len(), n_ind);

    // Build mapping from phenotype observations to genotypes
    // This follows R/GWASpoly: Z[cbind(1:n, match(pheno.gid, geno.gid))] <- 1
    let mut obs_indices: Vec<usize> = Vec::new();      // Which pheno rows to keep
    let mut obs_to_geno: Vec<usize> = Vec::new();      // Maps obs -> geno index (in subsetted geno)
    let mut missing_ids: Vec<String> = Vec::new();

    // Build geno ID -> index map (using subsetted geno)
    let geno_id_to_idx: HashMap<&str, usize> = geno.sample_ids
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();

    for (obs_idx, pheno_id) in pheno.sample_ids.iter().enumerate() {
        if let Some(&geno_idx) = geno_id_to_idx.get(pheno_id.as_str()) {
            obs_indices.push(obs_idx);
            obs_to_geno.push(geno_idx);
        } else {
            missing_ids.push(pheno_id.clone());
        }
    }

    if obs_indices.is_empty() {
        return Err(anyhow!("No overlapping sample_ids between genotype and phenotype files"));
    }

    if !allow_missing_samples && !missing_ids.is_empty() {
        // Only warn if there are many missing - some missing is expected
        if missing_ids.len() > pheno.sample_ids.len() / 2 {
            eprintln!("Warning: {} phenotype observations have no matching genotype", missing_ids.len());
        }
    }

    let n_obs = obs_indices.len();

    // Extract trait vector y for observations with matching genotypes
    let y_full = pheno
        .traits
        .get(trait_name)
        .ok_or_else(|| anyhow!("Trait '{}' not found in phenotype file", trait_name))?;

    let y: Array1<f64> = Array1::from_iter(obs_indices.iter().map(|&i| y_full[i]));

    // Build observation-level sample IDs (for reporting)
    let obs_sample_ids: Vec<String> = obs_indices
        .iter()
        .map(|&i| pheno.sample_ids[i].clone())
        .collect();

    // Build Z incidence matrix (n_obs × n_ind)
    // Z[i, j] = 1 if observation i belongs to genotype j
    let mut z_mat = Array2::<f64>::zeros((n_obs, n_ind));
    for (obs_row, &geno_col) in obs_to_geno.iter().enumerate() {
        z_mat[(obs_row, geno_col)] = 1.0;
    }

    // Build base design matrix X (n_obs × p) with intercept and covariates
    let base_design = build_base_design_obs_level(
        &pheno,
        covariate_names,
        &obs_indices,
        n_obs,
    )?;
    let x0 = base_design_to_array2(&base_design)?;

    let n_genotypes_for_qc = n_ind;
    let max_gf = if max_geno_freq <= 0.0 {
        (1.0 - 5.0 / (n_ind.max(1) as f64)).clamp(0.01, 0.99)
    } else {
        max_geno_freq
    };

    let mut all_results: Vec<MarkerResult> = Vec::new();

    if loco {
        // LOCO mode: compute chromosome-specific kinships using GWASpoly method
        // This matches R/GWASpoly's set.K(LOCO=TRUE) + makeLOCO approach:
        // 1. Compute per-chromosome K from only that chromosome's markers
        // 2. For testing chrX, average K's from all other chromosomes
        let loco_kinships = compute_loco_kinship_gwaspoly(&geno)?;

        // Group markers by chromosome
        let meta = geno
            .marker_metadata
            .as_ref()
            .ok_or_else(|| anyhow!("Marker metadata required for LOCO"))?;

        let mut chrom_markers: HashMap<String, Vec<usize>> = HashMap::new();
        for (idx, m) in meta.iter().enumerate() {
            chrom_markers.entry(m.chrom.clone()).or_default().push(idx);
        }

        // Process each chromosome
        for (chrom, marker_indices) in chrom_markers {
            let kin = loco_kinships
                .get(&chrom)
                .ok_or_else(|| anyhow!("Missing LOCO kinship for chromosome {}", chrom))?;

            let kin_aligned = align_kinship_to_genotypes_clone(kin, &geno.sample_ids)?;

            // Fit null model with Z matrix for repeated observations
            let cache = fit_null_model(
                &y,
                &x0,
                &kin_aligned,
                Some(&z_mat),
                Some(&obs_sample_ids)
            )?;

            // Test markers on this chromosome (parallel or sequential)
            if use_parallel {
                // Filter geno to only markers on this chromosome for parallel testing
                let chrom_results = parallel::test_markers_parallel_subset(
                    &geno,
                    &marker_indices,
                    &y,
                    &base_design,
                    &cache,
                    &obs_to_geno,
                    models,
                    ploidy,
                    min_maf,
                    max_gf,
                )?;
                all_results.extend(chrom_results);
            } else {
                for &marker_idx in &marker_indices {
                    let marker_results = test_marker_all_models_with_z(
                        &geno,
                        marker_idx,
                        &y,
                        &base_design,
                        &cache,
                        &z_mat,
                        &obs_to_geno,
                        models,
                        ploidy,
                        min_maf,
                        max_gf,
                        n_genotypes_for_qc,
                    )?;
                    all_results.extend(marker_results);
                }
            }
        }
    } else {
        // Standard mode: single kinship matrix
        let kin = if let Some(k_path) = kinship_path {
            load_kinship_from_tsv(k_path)?
        } else {
            compute_kinship_gwaspoly(&geno)?
        };
        let kin_aligned = align_kinship_to_genotypes_owned(kin, &geno.sample_ids)?;

        // Debug: print dimensions and sample IDs
        eprintln!("DEBUG: n_obs={}, n_ind={}, y.len()={}, x0.shape={:?}, z_mat.shape={:?}, kin.shape={:?}",
            n_obs, n_ind, y.len(), x0.shape(), z_mat.shape(), kin_aligned.matrix.shape());

        // Check kinship diagonal
        let kin_diag_mean = (0..n_ind).map(|i| kin_aligned.matrix[(i, i)]).sum::<f64>() / n_ind as f64;
        let kin_diag_min = (0..n_ind).map(|i| kin_aligned.matrix[(i, i)]).fold(f64::INFINITY, f64::min);
        let kin_diag_max = (0..n_ind).map(|i| kin_aligned.matrix[(i, i)]).fold(f64::NEG_INFINITY, f64::max);
        eprintln!("DEBUG: Kinship diagonal: mean={:.4}, min={:.4}, max={:.4}", kin_diag_mean, kin_diag_min, kin_diag_max);

        // Check Z matrix row sums (should all be 1.0)
        let z_row_sum_first5: Vec<f64> = (0..5.min(n_obs)).map(|i| {
            (0..n_ind).map(|j| z_mat[(i, j)]).sum()
        }).collect();
        eprintln!("DEBUG: Z matrix row sums (first 5): {:?}", z_row_sum_first5);

        // Check y variance
        let y_mean = y.iter().sum::<f64>() / n_obs as f64;
        let y_var = y.iter().map(|v| (v - y_mean).powi(2)).sum::<f64>() / (n_obs - 1) as f64;
        eprintln!("DEBUG: y mean={:.4}, var={:.4}", y_mean, y_var);

        // Fit null model with Z matrix for repeated observations
        // Note: We use the validated mixed_solve_new for accuracy.
        // The speedup comes from parallel marker testing (rayon).
        let cache = fit_null_model(
            &y,
            &x0,
            &kin_aligned,
            Some(&z_mat),
            Some(&obs_sample_ids)
        )?;

        // Debug: print variance components
        eprintln!("DEBUG: lambda={}, vu={}, ve={}, sigma2={}", cache.lambda, cache.vu, cache.ve, cache.sigma2);

        // Debug: print H_inv corner values for comparison with R
        eprintln!("DEBUG: H_inv[0:3,0:3]:");
        for i in 0..3.min(cache.h_inv.nrows()) {
            let row: Vec<f64> = (0..3.min(cache.h_inv.ncols()))
                .map(|j| cache.h_inv[(i, j)])
                .collect();
            eprintln!("  {:?}", row);
        }

        // Test all markers (parallel or sequential)
        if use_parallel {
            let parallel_results = parallel::test_markers_parallel(
                &geno,
                &y,
                &base_design,
                &cache,
                &obs_to_geno,
                models,
                ploidy,
                min_maf,
                max_gf,
            )?;
            all_results.extend(parallel_results);
        } else {
            for marker_idx in 0..geno.marker_ids.len() {
                let marker_results = test_marker_all_models_with_z(
                    &geno,
                    marker_idx,
                    &y,
                    &base_design,
                    &cache,
                    &z_mat,
                    &obs_to_geno,
                    models,
                    ploidy,
                    min_maf,
                    max_gf,
                    n_genotypes_for_qc,
                )?;
                all_results.extend(marker_results);
            }
        }
    }

    // Write results
    write_gwaspoly_results(out_path, &all_results)?;

    Ok(())
}

/// Test a single marker across all specified models with Z incidence matrix
pub(crate) fn test_marker_all_models_with_z(
    geno: &GenotypeMatrixBiallelic,
    marker_idx: usize,
    y: &Array1<f64>,
    base_design: &[Vec<f64>],
    cache: &GwasCache,
    _z_mat: &Array2<f64>,  // Z matrix (kept for API consistency, expansion done via obs_to_geno)
    obs_to_geno: &[usize],
    models: &[GeneActionModel],
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

    for model in models {
        // Get marker design at genotype level
        let marker_design_geno = match design_score(
            &dosage_geno,
            ploidy,
            *model,
            min_maf,
            max_geno_freq,
            n_genotypes_for_qc
        ) {
            Some(cols) => cols,
            None => continue, // Skip this model if QC fails
        };

        // Expand marker design to observation level via Z: Z %*% S
        // Each column of S (genotype-level) becomes a column at observation level
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
        });
    }

    Ok(results)
}

/// Generate marker design matrix based on gene action model (at genotype level)
pub(crate) fn design_score(
    dosage: &Array1<f64>,
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

    // Round dosages to integers for categorical models
    let rounded: Vec<i64> = dosage.iter().map(|&d| d.round() as i64).collect();

    // Count genotype frequencies
    let mut counts: HashMap<i64, usize> = HashMap::new();
    for &r in &rounded {
        *counts.entry(r).or_insert(0) += 1;
    }

    match model {
        GeneActionModel::Additive => {
            // Check max genotype frequency
            let max_freq = counts.values().cloned().max().unwrap_or(0) as f64 / (dosage.len() as f64);
            if max_freq > max_geno_freq {
                return None;
            }
            // Use original (possibly fractional) dosages
            Some(vec![dosage.to_vec()])
        }

        GeneActionModel::General => {
            // Check max genotype frequency
            let max_freq = counts.values().cloned().max().unwrap_or(0) as f64 / (dosage.len() as f64);
            if max_freq > max_geno_freq {
                return None;
            }

            // Get unique levels in order of appearance
            let mut levels: Vec<i64> = Vec::new();
            for &r in &rounded {
                if !levels.contains(&r) {
                    levels.push(r);
                }
            }

            if levels.len() <= 1 {
                return None;
            }

            // Dummy coding with first level as reference
            let mut cols = Vec::with_capacity(levels.len() - 1);
            for level in levels.iter().skip(1) {
                let col: Vec<f64> = rounded.iter().map(|&r| (r == *level) as i32 as f64).collect();
                cols.push(col);
            }
            Some(cols)
        }

        GeneActionModel::SimplexDomRef => {
            // 0 vs (1, 2, ..., ploidy): reference allele is dominant
            let col: Vec<f64> = rounded
                .iter()
                .map(|&r| if r > 0 { 1.0 } else { 0.0 })
                .collect();
            check_binary_column(&col, max_geno_freq)?;
            Some(vec![col])
        }

        GeneActionModel::SimplexDomAlt => {
            // (0, 1, ..., ploidy-1) vs ploidy: alternate allele is dominant
            let max_dose = ploidy as i64;
            let col: Vec<f64> = rounded
                .iter()
                .map(|&r| if r == max_dose { 1.0 } else { 0.0 })
                .collect();
            check_binary_column(&col, max_geno_freq)?;
            Some(vec![col])
        }

        GeneActionModel::DuplexDomRef => {
            // For tetraploid: (0, 1) vs (2, 3, 4)
            let threshold = (ploidy / 2) as i64;
            let col: Vec<f64> = rounded
                .iter()
                .map(|&r| if r >= threshold { 1.0 } else { 0.0 })
                .collect();
            check_binary_column(&col, max_geno_freq)?;
            Some(vec![col])
        }

        GeneActionModel::DuplexDomAlt => {
            // For tetraploid: (0, 1, 2) vs (3, 4)
            let threshold = (ploidy / 2 + 1) as i64;
            let col: Vec<f64> = rounded
                .iter()
                .map(|&r| if r >= threshold { 1.0 } else { 0.0 })
                .collect();
            check_binary_column(&col, max_geno_freq)?;
            Some(vec![col])
        }

        GeneActionModel::DiploGeneral => {
            // Diploidized: collapse all heterozygotes
            let max_dose = ploidy as i64;
            let diplo: Vec<i64> = rounded
                .iter()
                .map(|&r| {
                    if r == 0 { 0 }
                    else if r == max_dose { 2 }
                    else { 1 }
                })
                .collect();

            // Check frequencies
            let mut diplo_counts: HashMap<i64, usize> = HashMap::new();
            for &d in &diplo {
                *diplo_counts.entry(d).or_insert(0) += 1;
            }
            let max_freq = diplo_counts.values().cloned().max().unwrap_or(0) as f64 / (diplo.len() as f64);
            if max_freq > max_geno_freq {
                return None;
            }

            // Get levels
            let mut levels: Vec<i64> = Vec::new();
            for &d in &diplo {
                if !levels.contains(&d) {
                    levels.push(d);
                }
            }

            if levels.len() <= 1 {
                return None;
            }

            // Dummy coding
            let mut cols = Vec::with_capacity(levels.len() - 1);
            for level in levels.iter().skip(1) {
                let col: Vec<f64> = diplo.iter().map(|&d| (d == *level) as i32 as f64).collect();
                cols.push(col);
            }
            Some(cols)
        }

        GeneActionModel::DiploAdditive => {
            // Diploidized additive: 0, 0.5, 1 coding
            let max_dose = ploidy as f64;
            let col: Vec<f64> = dosage
                .iter()
                .map(|&d| {
                    if d <= 0.0 { 0.0 }
                    else if d >= max_dose { 1.0 }
                    else { 0.5 }
                })
                .collect();

            // Check frequencies
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

/// Check if a binary column passes QC
fn check_binary_column(col: &[f64], max_geno_freq: f64) -> Option<()> {
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

/// LMM score test for marker effect
pub(crate) fn lmm_score_test(
    marker_design: &[Vec<f64>],
    y: &Array1<f64>,
    base_design: &[Vec<f64>],
    cache: &GwasCache,
) -> Result<(f64, f64, Option<f64>)> {
    let n = y.len();
    let p0 = base_design.len();
    let p_marker = marker_design.len();

    if p_marker == 0 {
        return Err(anyhow!("Empty marker design"));
    }

    // Build full design X = [base_design | marker_design]
    let p = p0 + p_marker;
    let mut data = Vec::with_capacity(n * p);
    for row in 0..n {
        for col in base_design {
            data.push(col[row]);
        }
        for col in marker_design {
            data.push(col[row]);
        }
    }
    let x = DMatrix::from_row_slice(n, p, &data);

    let y_vec = DVector::from_row_slice(
        y.as_slice().ok_or_else(|| anyhow!("y not contiguous"))?,
    );

    // W = X' H^-1 X
    let w = &x.transpose() * &cache.h_inv * &x;
    let w_inv = match w.clone().try_inverse() {
        Some(inv) => inv,
        None => {
            let mut w_eps = w.clone();
            for i in 0..p {
                w_eps[(i, i)] += 1e-8;
            }
            w_eps.try_inverse().ok_or_else(|| anyhow!("W not invertible"))?
        }
    };

    // beta = W^-1 X' H^-1 y
    let beta = &w_inv * (&x.transpose() * &cache.h_inv * &y_vec);

    // Residual variance
    let resid = &y_vec - &x * &beta;
    let v2 = (n as f64) - (p as f64);
    if v2 <= 0.0 {
        return Err(anyhow!("Degrees of freedom <= 0"));
    }
    let s2 = (resid.transpose() * &cache.h_inv * &resid)[(0, 0)] / v2;

    // Extract marker portion
    let marker_beta = beta.rows(p0, p_marker).into_owned();
    let w_marker = w_inv.slice((p0, p0), (p_marker, p_marker)).into_owned();

    let inv_w_marker = match w_marker.clone().try_inverse() {
        Some(m) => m,
        None => {
            let mut jittered = w_marker.clone();
            for i in 0..p_marker {
                jittered[(i, i)] += 1e-8;
            }
            jittered
                .try_inverse()
                .ok_or_else(|| anyhow!("W_marker not invertible"))?
        }
    };

    // F-statistic
    let denom = s2 * (p_marker as f64);
    let f_stat = if denom > 0.0 {
        (marker_beta.transpose() * &inv_w_marker * &marker_beta)[(0, 0)] / denom
    } else {
        0.0
    };

    // P-value from F distribution
    let p_value = if f_stat.is_finite() && f_stat >= 0.0 {
        let f_dist = FisherSnedecor::new(p_marker as f64, v2)?;
        (1.0 - f_dist.cdf(f_stat)).max(0.0)
    } else {
        1.0
    };

    // Score = -log10(p)
    // Cap at a practical minimum p-value to avoid numerical overflow differences with R
    // R's pf() function has different precision limits than statrs FisherSnedecor
    // A p-value of 1e-300 (score=300) is already infinitely significant
    let p_clamped = p_value.max(1e-300);
    let score = -p_clamped.log10();

    // Effect estimate (only for single-column design)
    let effect = if p_marker == 1 {
        Some(marker_beta[0])
    } else {
        None
    };

    Ok((score, p_value, effect))
}

// Helper functions

/// Build base design matrix at observation level (n_obs × p)
fn build_base_design_obs_level(
    pheno: &binx_core::PhenotypeTable,
    covariate_names: Option<&[String]>,
    obs_indices: &[usize],
    n_obs: usize,
) -> Result<Vec<Vec<f64>>> {
    let mut cols = Vec::new();

    // Intercept
    cols.push(vec![1.0f64; n_obs]);

    // Covariates
    if let Some(names) = covariate_names {
        for name in names {
            if let Some(cov) = pheno.traits.get(name).or_else(|| pheno.covariates.get(name)) {
                // Numeric covariate
                let col: Vec<f64> = obs_indices.iter().map(|&i| cov[i]).collect();
                cols.push(col);
            } else if let Some(fvals) = pheno.factor_covariates.get(name) {
                // Factor covariate: dummy-code with reference level dropped
                let values: Vec<String> = obs_indices.iter().map(|&i| fvals[i].clone()).collect();

                // Get unique levels in order of appearance
                let mut levels: Vec<String> = Vec::new();
                for v in &values {
                    if !levels.contains(v) {
                        levels.push(v.clone());
                    }
                }

                if levels.len() > 1 {
                    // Dummy code: skip first level (reference)
                    for level in levels.iter().skip(1) {
                        let col: Vec<f64> = values
                            .iter()
                            .map(|v| (v == level) as i32 as f64)
                            .collect();
                        cols.push(col);
                    }
                }
            } else {
                return Err(anyhow!("Covariate '{}' not found", name));
            }
        }
    }

    Ok(cols)
}

pub(crate) fn base_design_to_array2(cols: &[Vec<f64>]) -> Result<Array2<f64>> {
    if cols.is_empty() {
        return Err(anyhow!("Base design must have at least one column"));
    }
    let n_samples = cols[0].len();
    let p = cols.len();
    let mut mat = Array2::<f64>::zeros((n_samples, p));
    for (j, col) in cols.iter().enumerate() {
        for (i, val) in col.iter().enumerate() {
            mat[(i, j)] = *val;
        }
    }
    Ok(mat)
}

pub(crate) fn align_kinship_to_genotypes_clone(
    kin: &KinshipMatrix,
    geno_ids: &[String],
) -> Result<KinshipMatrix> {
    let mut map = HashMap::new();
    for (i, sid) in kin.sample_ids.iter().enumerate() {
        map.insert(sid.as_str(), i);
    }

    let mut idx = Vec::with_capacity(geno_ids.len());
    for sid in geno_ids {
        let &i = map
            .get(sid.as_str())
            .ok_or_else(|| anyhow!("Kinship missing sample {}", sid))?;
        idx.push(i);
    }

    let n = geno_ids.len();
    let mut matrix = Array2::<f64>::zeros((n, n));
    for (i_new, &i_old) in idx.iter().enumerate() {
        for (j_new, &j_old) in idx.iter().enumerate() {
            matrix[(i_new, j_new)] = kin.matrix[(i_old, j_old)];
        }
    }

    Ok(KinshipMatrix {
        sample_ids: geno_ids.to_vec(),
        matrix,
    })
}

fn align_kinship_to_genotypes_owned(
    kin: KinshipMatrix,
    geno_ids: &[String],
) -> Result<KinshipMatrix> {
    let mut map = HashMap::new();
    for (i, sid) in kin.sample_ids.iter().enumerate() {
        map.insert(sid.as_str(), i);
    }

    let mut idx = Vec::with_capacity(geno_ids.len());
    for sid in geno_ids {
        let &i = map
            .get(sid.as_str())
            .ok_or_else(|| anyhow!("Kinship missing sample {}", sid))?;
        idx.push(i);
    }

    let n = geno_ids.len();
    let mut matrix = Array2::<f64>::zeros((n, n));
    for (i_new, &i_old) in idx.iter().enumerate() {
        for (j_new, &j_old) in idx.iter().enumerate() {
            matrix[(i_new, j_new)] = kin.matrix[(i_old, j_old)];
        }
    }

    Ok(KinshipMatrix {
        sample_ids: geno_ids.to_vec(),
        matrix,
    })
}

fn write_gwaspoly_results(path: &str, results: &[MarkerResult]) -> Result<()> {
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b',')
        .from_path(path)?;

    wtr.write_record(&[
        "marker_id",
        "chrom",
        "pos",
        "model",
        "score",
        "p_value",
        "effect",
        "n_obs",
    ])?;

    for r in results {
        let chrom = r.chrom.as_deref().unwrap_or("");
        let pos = r.pos.map(|p| p.to_string()).unwrap_or_default();
        let effect = r
            .effect
            .map(|e| e.to_string())
            .unwrap_or_else(|| "NA".to_string());
        wtr.write_record(&[
            &r.marker_id,
            chrom,
            &pos,
            &r.model,
            &r.score.to_string(),
            &r.p_value.to_string(),
            &effect,
            &r.n_obs.to_string(),
        ])?;
    }

    wtr.flush()?;
    Ok(())
}

/// Write results in GWASpoly-style wide format (markers as rows, models as columns)
pub fn write_gwaspoly_scores_wide(
    path: &str,
    results: &[MarkerResult],
    models: &[GeneActionModel],
) -> Result<()> {
    // Group results by marker
    let mut marker_scores: HashMap<String, HashMap<String, f64>> = HashMap::new();

    for r in results {
        marker_scores
            .entry(r.marker_id.clone())
            .or_default()
            .insert(r.model.clone(), r.score);
    }

    // Get marker order (preserve original order)
    let mut markers: Vec<String> = Vec::new();
    for r in results {
        if !markers.contains(&r.marker_id) {
            markers.push(r.marker_id.clone());
        }
    }

    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b',')
        .from_path(path)?;

    // Header: marker_id, model1, model2, ...
    let mut header = vec!["".to_string()];
    for model in models {
        header.push(model.as_str().to_string());
    }
    wtr.write_record(&header)?;

    // Data rows
    for marker in &markers {
        let mut row = vec![marker.clone()];
        if let Some(scores) = marker_scores.get(marker) {
            for model in models {
                let score = scores.get(model.as_str()).copied().unwrap_or(f64::NAN);
                row.push(score.to_string());
            }
        } else {
            for _ in models {
                row.push("NA".to_string());
            }
        }
        wtr.write_record(&row)?;
    }

    wtr.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn test_design_additive() {
        let dosage = array![0.0, 1.0, 2.0, 3.0, 4.0];
        let design = design_score(&dosage, 4, GeneActionModel::Additive, 0.0, 1.0, 5);
        assert!(design.is_some());
        let cols = design.unwrap();
        assert_eq!(cols.len(), 1);
        assert_eq!(cols[0], vec![0.0, 1.0, 2.0, 3.0, 4.0]);
    }

    #[test]
    fn test_design_general() {
        let dosage = array![0.0, 1.0, 2.0, 3.0, 4.0];
        let design = design_score(&dosage, 4, GeneActionModel::General, 0.0, 1.0, 5);
        assert!(design.is_some());
        let cols = design.unwrap();
        // Should have ploidy columns (5 levels - 1 reference = 4)
        assert_eq!(cols.len(), 4);
    }

    #[test]
    fn test_design_simplex_dom_ref() {
        let dosage = array![0.0, 1.0, 2.0, 3.0, 4.0];
        let design = design_score(&dosage, 4, GeneActionModel::SimplexDomRef, 0.0, 1.0, 5);
        assert!(design.is_some());
        let cols = design.unwrap();
        assert_eq!(cols.len(), 1);
        // 0 -> 0, everything else -> 1
        assert_eq!(cols[0], vec![0.0, 1.0, 1.0, 1.0, 1.0]);
    }

    #[test]
    fn test_design_diplo_additive() {
        let dosage = array![0.0, 1.0, 2.0, 3.0, 4.0];
        let design = design_score(&dosage, 4, GeneActionModel::DiploAdditive, 0.0, 1.0, 5);
        assert!(design.is_some());
        let cols = design.unwrap();
        assert_eq!(cols.len(), 1);
        // 0 -> 0, heterozygotes -> 0.5, ploidy -> 1
        assert_eq!(cols[0], vec![0.0, 0.5, 0.5, 0.5, 1.0]);
    }

    #[test]
    fn test_model_parsing() {
        assert_eq!(
            GeneActionModel::from_str("additive").unwrap(),
            GeneActionModel::Additive
        );
        assert_eq!(
            GeneActionModel::from_str("1-dom-ref").unwrap(),
            GeneActionModel::SimplexDomRef
        );
    }
}
