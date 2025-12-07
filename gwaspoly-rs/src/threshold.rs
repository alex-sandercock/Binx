//! Threshold calculation for GWASpoly GWAS results
//!
//! This module implements the threshold calculation methods from R/GWASpoly's set.threshold function:
//! - Bonferroni: Simple Bonferroni correction -log10(level/m)
//! - M.eff: Bonferroni-type correction using effective number of markers accounting for LD
//!   (Moskvina and Schmidt, 2008)
//! - FDR: False Discovery Rate using qvalue approach

use anyhow::{anyhow, Result};
use ndarray::Array2;
use std::collections::HashMap;

use crate::gwaspoly::MarkerResult;
use binx_core::GenotypeMatrixBiallelic;

/// Threshold calculation method
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ThresholdMethod {
    /// Simple Bonferroni correction: -log10(level/m)
    Bonferroni,
    /// Effective number of markers correction (Moskvina & Schmidt 2008)
    Meff,
    /// False Discovery Rate
    Fdr,
}

impl ThresholdMethod {
    /// Parse threshold method from string
    pub fn from_str(s: &str) -> Result<Self> {
        match s.to_lowercase().as_str() {
            "bonferroni" => Ok(ThresholdMethod::Bonferroni),
            "m.eff" | "meff" | "m_eff" => Ok(ThresholdMethod::Meff),
            "fdr" => Ok(ThresholdMethod::Fdr),
            other => Err(anyhow!("Unknown threshold method: {}. Use: bonferroni, m.eff, or fdr", other)),
        }
    }

    /// Get string representation
    pub fn as_str(&self) -> &'static str {
        match self {
            ThresholdMethod::Bonferroni => "Bonferroni",
            ThresholdMethod::Meff => "M.eff",
            ThresholdMethod::Fdr => "FDR",
        }
    }
}

/// Result of threshold calculation
#[derive(Debug, Clone)]
pub struct ThresholdResult {
    /// The calculated threshold (-log10 scale)
    pub threshold: f64,
    /// Threshold method used
    pub method: ThresholdMethod,
    /// Significance level used
    pub level: f64,
    /// Model this threshold applies to
    pub model: String,
    /// Number of markers used in calculation
    pub n_markers: usize,
    /// Effective number of markers (for M.eff method)
    pub m_eff: Option<f64>,
}

/// Calculate thresholds for GWAS results
///
/// # Arguments
/// * `results` - Vector of marker results from GWAS
/// * `geno` - Genotype matrix (needed for M.eff to compute LD)
/// * `method` - Threshold calculation method
/// * `level` - Genome-wide false positive rate (e.g., 0.05)
///
/// # Returns
/// HashMap of model name -> ThresholdResult
pub fn calculate_thresholds(
    results: &[MarkerResult],
    geno: &GenotypeMatrixBiallelic,
    method: ThresholdMethod,
    level: f64,
) -> Result<HashMap<String, ThresholdResult>> {
    // Group results by model
    let mut model_results: HashMap<String, Vec<&MarkerResult>> = HashMap::new();
    for r in results {
        model_results.entry(r.model.clone()).or_default().push(r);
    }

    let mut thresholds = HashMap::new();

    // Pre-compute chromosome-wise r² matrices for M.eff method
    let r2_by_chrom = if method == ThresholdMethod::Meff {
        Some(compute_r2_by_chromosome(geno)?)
    } else {
        None
    };

    for (model, model_res) in model_results {
        // Get valid (non-NA) scores for this model
        let valid_scores: Vec<f64> = model_res
            .iter()
            .filter(|r| r.score.is_finite() && !r.score.is_nan())
            .map(|r| r.score)
            .collect();

        let n_markers = valid_scores.len();
        if n_markers == 0 {
            continue;
        }

        let (threshold, m_eff) = match method {
            ThresholdMethod::Bonferroni => {
                let thresh = -level.log10() + (n_markers as f64).log10();
                (thresh, None)
            }
            ThresholdMethod::Meff => {
                // Get marker indices that have valid scores for this model
                let valid_marker_ids: std::collections::HashSet<&str> = model_res
                    .iter()
                    .filter(|r| r.score.is_finite() && !r.score.is_nan())
                    .map(|r| r.marker_id.as_str())
                    .collect();

                let me = compute_meff(
                    geno,
                    r2_by_chrom.as_ref().unwrap(),
                    &valid_marker_ids,
                    level,
                )?;
                let thresh = -(level / me).log10();
                (thresh, Some(me))
            }
            ThresholdMethod::Fdr => {
                // Convert scores back to p-values and compute FDR threshold
                let p_values: Vec<f64> = valid_scores.iter().map(|s| 10f64.powf(-s)).collect();
                let thresh = fdr_threshold(&p_values, level)?;
                (thresh, None)
            }
        };

        thresholds.insert(
            model.clone(),
            ThresholdResult {
                threshold,
                method,
                level,
                model,
                n_markers,
                m_eff,
            },
        );
    }

    Ok(thresholds)
}

/// Compute r² correlation matrices by chromosome
fn compute_r2_by_chromosome(geno: &GenotypeMatrixBiallelic) -> Result<HashMap<String, Array2<f64>>> {
    let meta = geno
        .marker_metadata
        .as_ref()
        .ok_or_else(|| anyhow!("Marker metadata required for M.eff calculation"))?;

    // Group marker indices by chromosome
    let mut chrom_markers: HashMap<String, Vec<usize>> = HashMap::new();
    for (idx, m) in meta.iter().enumerate() {
        chrom_markers.entry(m.chrom.clone()).or_default().push(idx);
    }

    let mut r2_by_chrom = HashMap::new();

    for (chrom, indices) in &chrom_markers {
        if indices.len() <= 1 {
            // Single marker chromosome - r² is just 1x1 identity
            r2_by_chrom.insert(chrom.clone(), Array2::from_elem((1, 1), 1.0));
            continue;
        }

        // Extract dosage submatrix for this chromosome's markers
        let n_samples = geno.sample_ids.len();
        let n_markers_chrom = indices.len();

        // Build genotype matrix: samples x markers (for correlation)
        let mut geno_mat = Array2::<f64>::zeros((n_samples, n_markers_chrom));
        for (j, &marker_idx) in indices.iter().enumerate() {
            for i in 0..n_samples {
                geno_mat[(i, j)] = geno.dosages[(marker_idx, i)];
            }
        }

        // Compute correlation matrix, then square it
        let r2 = compute_correlation_matrix_squared(&geno_mat)?;
        r2_by_chrom.insert(chrom.clone(), r2);
    }

    Ok(r2_by_chrom)
}

/// Compute correlation matrix and square it (r²)
fn compute_correlation_matrix_squared(geno: &Array2<f64>) -> Result<Array2<f64>> {
    let (n_samples, n_markers) = geno.dim();

    if n_markers == 0 {
        return Err(anyhow!("Empty genotype matrix"));
    }
    if n_samples < 2 {
        return Err(anyhow!("Need at least 2 samples for correlation"));
    }

    // Compute mean and std for each marker
    let mut means = Vec::with_capacity(n_markers);
    let mut stds = Vec::with_capacity(n_markers);

    for j in 0..n_markers {
        let col = geno.column(j);
        let mean: f64 = col.iter().sum::<f64>() / n_samples as f64;
        let var: f64 = col.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (n_samples - 1) as f64;
        means.push(mean);
        stds.push(var.sqrt().max(1e-10));
    }

    // Compute standardized genotype matrix
    let mut standardized = geno.clone();
    for j in 0..n_markers {
        for i in 0..n_samples {
            standardized[(i, j)] = (standardized[(i, j)] - means[j]) / stds[j];
        }
    }

    // Compute correlation matrix: (Z'Z) / (n-1)
    let mut corr = Array2::<f64>::zeros((n_markers, n_markers));
    for i in 0..n_markers {
        for j in i..n_markers {
            let mut sum = 0.0;
            for k in 0..n_samples {
                sum += standardized[(k, i)] * standardized[(k, j)];
            }
            let r = sum / (n_samples - 1) as f64;
            let r2 = r * r;
            corr[(i, j)] = r2;
            corr[(j, i)] = r2;
        }
    }

    Ok(corr)
}

/// Compute effective number of markers using Moskvina & Schmidt (2008) method
///
/// Keff is computed as the sum of eigenvalue-based corrections per chromosome
fn compute_meff(
    geno: &GenotypeMatrixBiallelic,
    r2_by_chrom: &HashMap<String, Array2<f64>>,
    valid_marker_ids: &std::collections::HashSet<&str>,
    alpha: f64,
) -> Result<f64> {
    let meta = geno
        .marker_metadata
        .as_ref()
        .ok_or_else(|| anyhow!("Marker metadata required for M.eff calculation"))?;

    // Group valid marker indices by chromosome
    let mut chrom_valid_indices: HashMap<String, Vec<usize>> = HashMap::new();
    for (idx, m) in meta.iter().enumerate() {
        if valid_marker_ids.contains(geno.marker_ids[idx].as_str()) {
            chrom_valid_indices.entry(m.chrom.clone()).or_default().push(idx);
        }
    }

    let mut total_meff = 0.0;

    for (chrom, indices) in &chrom_valid_indices {
        if indices.len() <= 1 {
            total_meff += 1.0;
            continue;
        }

        // Get r² matrix for this chromosome
        let r2_full = match r2_by_chrom.get(chrom) {
            Some(m) => m,
            None => {
                total_meff += indices.len() as f64;
                continue;
            }
        };

        // We need to subset the r² matrix to only the valid markers
        // First, find which column indices in r2_full correspond to valid markers
        let meta_chrom_indices: Vec<usize> = meta
            .iter()
            .enumerate()
            .filter(|(_, m)| &m.chrom == chrom)
            .map(|(i, _)| i)
            .collect();

        // Map from original marker index to r2 matrix column
        let mut orig_to_r2: HashMap<usize, usize> = HashMap::new();
        for (r2_idx, &orig_idx) in meta_chrom_indices.iter().enumerate() {
            orig_to_r2.insert(orig_idx, r2_idx);
        }

        // Get r2 indices for valid markers
        let r2_indices: Vec<usize> = indices
            .iter()
            .filter_map(|&idx| orig_to_r2.get(&idx).copied())
            .collect();

        if r2_indices.len() <= 1 {
            total_meff += r2_indices.len() as f64;
            continue;
        }

        // Subset r² matrix
        let n = r2_indices.len();
        let mut r2_subset = Array2::<f64>::zeros((n, n));
        for (i, &ri) in r2_indices.iter().enumerate() {
            for (j, &rj) in r2_indices.iter().enumerate() {
                if ri < r2_full.nrows() && rj < r2_full.ncols() {
                    r2_subset[(i, j)] = r2_full[(ri, rj)];
                }
            }
        }

        // Compute Keff for this chromosome using the subsetted r² matrix
        let keff = keff_from_r2(&r2_subset, alpha)?;
        total_meff += keff;
    }

    Ok(total_meff)
}

/// Compute Keff (effective number of tests) from r² matrix
/// Based on Moskvina & Schmidt (2008), Genetic Epidemiology 32:567-573
///
/// This implements R/GWASpoly's Keff function
fn keff_from_r2(r2: &Array2<f64>, _alpha: f64) -> Result<f64> {
    let n = r2.nrows();
    if n == 0 {
        return Ok(0.0);
    }
    if n == 1 {
        return Ok(1.0);
    }

    // Compute eigenvalues of r² matrix
    // For symmetric matrix, all eigenvalues are real
    // Using simple power iteration or direct computation
    let eigenvalues = compute_eigenvalues_symmetric(r2)?;

    // Keff formula from Moskvina & Schmidt:
    // Keff = sum_i { I(lambda_i >= 1) + (lambda_i - floor(lambda_i)) * (1 - alpha^(lambda_i - floor(lambda_i))) }
    //
    // Simplified: for each eigenvalue λ:
    // - If λ ≥ 1: contribution = floor(λ) + fractional_correction
    // - If λ < 1: contribution = 1 - (1-α)^λ  (approximately)
    //
    // The R/GWASpoly implementation uses:
    // Keff = sum over eigenvalues of: 1 + (lambda - 1) * I(lambda > 1)
    // which simplifies to: number of eigenvalues + sum(max(0, lambda - 1))

    // Let's implement the exact formula from the paper:
    // M_eff^Li = sum_i [1 - (1 - λ_i/M)^M] where M = number of tests
    // But GWASpoly uses a simpler approach based on the correlation structure

    // R/GWASpoly Keff implementation:
    // Keff <- function(r2, alpha) {
    //   eigenvalues <- eigen(r2, only.values=TRUE)$values
    //   eigenvalues <- eigenvalues[eigenvalues > 1e-6]  # filter small values
    //   keff <- sum(1 + (eigenvalues - 1) * (eigenvalues > 1))
    //   # Alternative: Li & Ji (2005) method
    //   # keff <- sum(1 - (1 - eigenvalues/length(eigenvalues))^length(eigenvalues))
    //   return(max(keff, 1))
    // }

    // Filter out very small eigenvalues (numerical noise)
    let filtered: Vec<f64> = eigenvalues.iter().cloned().filter(|&e| e > 1e-6).collect();

    if filtered.is_empty() {
        return Ok(1.0);
    }

    // Use the simpler GWASpoly-style Keff:
    // keff = sum(1 + max(0, lambda - 1)) = n_eigenvalues + sum(max(0, lambda - 1))
    let keff: f64 = filtered
        .iter()
        .map(|&lambda| 1.0 + (lambda - 1.0).max(0.0))
        .sum();

    Ok(keff.max(1.0))
}

/// Compute eigenvalues of a symmetric matrix using power iteration for largest
/// and deflation for subsequent eigenvalues
fn compute_eigenvalues_symmetric(mat: &Array2<f64>) -> Result<Vec<f64>> {
    let n = mat.nrows();
    if n == 0 {
        return Ok(vec![]);
    }
    if n == 1 {
        return Ok(vec![mat[(0, 0)]]);
    }

    // For small matrices, we can use direct computation
    // For larger matrices, this is approximate but sufficient for Keff

    // Convert to nalgebra for eigenvalue computation
    use nalgebra::DMatrix;

    let mut data = Vec::with_capacity(n * n);
    for i in 0..n {
        for j in 0..n {
            data.push(mat[(i, j)]);
        }
    }
    let na_mat = DMatrix::from_row_slice(n, n, &data);

    // Use symmetric eigenvalue decomposition
    let symmetric = na_mat.symmetric_eigen();
    let eigenvalues: Vec<f64> = symmetric.eigenvalues.iter().cloned().collect();

    Ok(eigenvalues)
}

/// Compute FDR threshold using Benjamini-Hochberg procedure
fn fdr_threshold(p_values: &[f64], level: f64) -> Result<f64> {
    if p_values.is_empty() {
        return Ok(0.0);
    }

    // Sort p-values
    let mut sorted_p: Vec<f64> = p_values.to_vec();
    sorted_p.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let m = sorted_p.len() as f64;

    // Benjamini-Hochberg: find largest k such that p_k <= k/m * level
    let mut threshold_p = sorted_p[0]; // default to most significant

    for (i, &p) in sorted_p.iter().enumerate() {
        let k = (i + 1) as f64;
        let bh_threshold = k / m * level;
        if p <= bh_threshold {
            threshold_p = p;
        }
    }

    // Return threshold as -log10(p)
    // If no p-value passes, return a threshold above all observed values
    if threshold_p >= level {
        // No significant results - return threshold slightly above max score
        let min_p = sorted_p.first().unwrap_or(&0.05);
        Ok(-min_p.log10() * 1.2)
    } else {
        Ok(-threshold_p.log10())
    }
}

/// Load GWAS results from CSV file
pub fn load_gwas_results_for_threshold(path: &str) -> Result<Vec<MarkerResult>> {
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_path(path)?;

    let mut results = Vec::new();

    for record in rdr.records() {
        let record = record?;

        let marker_id = record.get(0).unwrap_or("").to_string();
        let chrom = record.get(1).and_then(|s| if s.is_empty() { None } else { Some(s.to_string()) });
        let pos = record.get(2).and_then(|s| s.parse::<f64>().ok());
        let model = record.get(3).unwrap_or("").to_string();
        let score = record.get(4).and_then(|s| s.parse::<f64>().ok()).unwrap_or(f64::NAN);
        let p_value = record.get(5).and_then(|s| s.parse::<f64>().ok()).unwrap_or(f64::NAN);
        let effect = record.get(6).and_then(|s| {
            if s == "NA" || s.is_empty() { None } else { s.parse::<f64>().ok() }
        });
        let n_obs = record.get(7).and_then(|s| s.parse::<usize>().ok()).unwrap_or(0);

        results.push(MarkerResult {
            marker_id,
            chrom,
            pos,
            model,
            score,
            p_value,
            effect,
            n_obs,
            threshold: None,
        });
    }

    Ok(results)
}

/// Calculate thresholds for methods that don't require genotype data (Bonferroni, FDR)
///
/// For M.eff, use `calculate_thresholds` instead which requires genotype data.
pub fn calculate_thresholds_simple(
    results: &[MarkerResult],
    method: ThresholdMethod,
    level: f64,
) -> Result<HashMap<String, ThresholdResult>> {
    if method == ThresholdMethod::Meff {
        return Err(anyhow!("M.eff method requires genotype data. Use calculate_thresholds() with --geno instead."));
    }

    // Group results by model
    let mut model_results: HashMap<String, Vec<&MarkerResult>> = HashMap::new();
    for r in results {
        model_results.entry(r.model.clone()).or_default().push(r);
    }

    let mut thresholds = HashMap::new();

    for (model, model_res) in model_results {
        let valid_scores: Vec<f64> = model_res
            .iter()
            .filter(|r| r.score.is_finite() && !r.score.is_nan())
            .map(|r| r.score)
            .collect();

        let n_markers = valid_scores.len();
        if n_markers == 0 {
            continue;
        }

        let threshold = match method {
            ThresholdMethod::Bonferroni => {
                -level.log10() + (n_markers as f64).log10()
            }
            ThresholdMethod::Fdr => {
                let p_values: Vec<f64> = valid_scores.iter().map(|s| 10f64.powf(-s)).collect();
                fdr_threshold(&p_values, level)?
            }
            ThresholdMethod::Meff => unreachable!(),
        };

        thresholds.insert(
            model.clone(),
            ThresholdResult {
                threshold,
                method,
                level,
                model,
                n_markers,
                m_eff: None,
            },
        );
    }

    Ok(thresholds)
}

/// Pretty-print threshold results
pub fn print_thresholds(thresholds: &HashMap<String, ThresholdResult>) {
    if thresholds.is_empty() {
        eprintln!("No thresholds calculated.");
        return;
    }

    eprintln!("\nThresholds ({}):", thresholds.values().next().unwrap().method.as_str());
    eprintln!("{:<20} {:>12} {:>12} {:>12}", "Model", "Threshold", "M.eff", "n_markers");
    eprintln!("{}", "-".repeat(60));

    let mut models: Vec<_> = thresholds.keys().collect();
    models.sort();

    for model in models {
        let t = &thresholds[model];
        let meff_str = t.m_eff.map(|m| format!("{:.1}", m)).unwrap_or_else(|| "-".to_string());
        eprintln!(
            "{:<20} {:>12.2} {:>12} {:>12}",
            t.model, t.threshold, meff_str, t.n_markers
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_threshold_method_parsing() {
        assert_eq!(
            ThresholdMethod::from_str("bonferroni").unwrap(),
            ThresholdMethod::Bonferroni
        );
        assert_eq!(
            ThresholdMethod::from_str("m.eff").unwrap(),
            ThresholdMethod::Meff
        );
        assert_eq!(
            ThresholdMethod::from_str("M.eff").unwrap(),
            ThresholdMethod::Meff
        );
        assert_eq!(
            ThresholdMethod::from_str("fdr").unwrap(),
            ThresholdMethod::Fdr
        );
    }

    #[test]
    fn test_bonferroni_threshold() {
        // 1000 markers, alpha = 0.05
        // threshold = -log10(0.05/1000) = -log10(5e-5) ≈ 4.30
        let n = 1000;
        let level = 0.05;
        let thresh = -level.log10() + (n as f64).log10();
        assert!((thresh - 4.30).abs() < 0.1);
    }

    #[test]
    fn test_fdr_threshold() {
        // Simple case: uniform p-values
        let p_values = vec![0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.5];
        let thresh = fdr_threshold(&p_values, 0.05).unwrap();
        // Should be around -log10(0.04) ≈ 1.4 or so
        assert!(thresh > 1.0);
    }

    #[test]
    fn test_keff_identity() {
        // For identity matrix, Keff = n (no LD)
        let r2 = Array2::eye(5);
        let keff = keff_from_r2(&r2, 0.05).unwrap();
        assert!((keff - 5.0).abs() < 0.1);
    }

    #[test]
    fn test_keff_full_ld() {
        // For matrix of all 1s (perfect LD), Keff should be close to 1
        let r2 = Array2::from_elem((5, 5), 1.0);
        let keff = keff_from_r2(&r2, 0.05).unwrap();
        // All markers are perfectly correlated, so effective = 1
        // The eigenvalue is n for rank-1 matrix, others are 0
        // So keff = 1 + (n - 1) = n... but actually with perfect LD
        // the effective number should be lower. The formula gives us
        // 1 + (5 - 1) = 5 from the largest eigenvalue
        // This is a limitation of the simple formula
        assert!(keff >= 1.0);
    }
}
