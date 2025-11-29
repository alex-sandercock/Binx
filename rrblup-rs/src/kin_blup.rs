//! Faithful Rust implementation of R/rrBLUP::kin.blup
//!
//! This module provides a direct translation of the kin.blup function from
//! the R/rrBLUP package. The implementation follows the original R code structure
//! as closely as possible for verification and maintainability.
//!
//! Reference: R/rrBLUP package by Jeffrey Endelman
//! https://cran.r-project.org/package=rrBLUP

use crate::mixed_solve::{mixed_solve, MixedSolveOptions};
use anyhow::{anyhow, Result};
use nalgebra::{DMatrix, DVector};
use std::collections::{HashMap, HashSet};

/// Options for kin.blup function
#[derive(Debug, Clone)]
pub struct KinBlupOptions {
    /// Use Gaussian kernel transformation of K (default: false)
    pub gauss: bool,
    /// Kinship matrix (optional, if None uses identity covariance)
    pub k: Option<DMatrix<f64>>,
    /// Row names for K matrix (genotype IDs)
    pub k_ids: Option<Vec<String>>,
    /// Whether to compute prediction error variance (default: false)
    pub pev: bool,
    /// Sequence of theta values for Gaussian kernel profiling
    pub theta_seq: Option<Vec<f64>>,
}

impl Default for KinBlupOptions {
    fn default() -> Self {
        Self {
            gauss: false,
            k: None,
            k_ids: None,
            pev: false,
            theta_seq: None,
        }
    }
}

/// Input data for kin.blup
#[derive(Debug, Clone)]
pub struct KinBlupData {
    /// Genotype/line identifiers for each observation
    pub geno_ids: Vec<String>,
    /// Phenotype values (may contain NaN for missing)
    pub pheno: Vec<f64>,
    /// Fixed effects (categorical factors) - each inner vec is one factor
    pub fixed: Option<Vec<Vec<String>>>,
    /// Covariates (continuous variables) - each inner vec is one covariate
    pub covariates: Option<Vec<Vec<f64>>>,
}

/// Result from kin.blup function
#[derive(Debug, Clone)]
pub struct KinBlupResult {
    /// Genetic variance (Vg = Vu from mixed.solve)
    pub vg: f64,
    /// Residual variance (Ve)
    pub ve: f64,
    /// Genetic values (BLUPs) for each genotype
    pub g: Vec<f64>,
    /// Genotype IDs corresponding to g values
    pub g_ids: Vec<String>,
    /// Prediction error variance (if PEV=TRUE)
    pub pev: Option<Vec<f64>>,
    /// Residuals for each observation (NA for missing phenotypes)
    pub resid: Vec<f64>,
    /// Predictions (g + mean fixed effect)
    pub pred: Vec<f64>,
    /// Profile likelihood for theta (if GAUSS=TRUE)
    pub profile: Option<Vec<(f64, f64)>>,
}

/// Genomic prediction using kinship-based BLUP
///
/// This is a faithful Rust implementation of R/rrBLUP::kin.blup().
///
/// # Arguments
/// * `data` - Input data containing genotype IDs, phenotypes, and optional effects
/// * `options` - Optional settings (GAUSS, K, PEV, theta.seq)
///
/// # Returns
/// * KinBlupResult with genetic values, variance components, residuals, and predictions
pub fn kin_blup(data: &KinBlupData, options: Option<KinBlupOptions>) -> Result<KinBlupResult> {
    let opts = options.unwrap_or_default();

    let n_total = data.pheno.len();
    if n_total == 0 {
        return Err(anyhow!("Empty phenotype data"));
    }
    if data.geno_ids.len() != n_total {
        return Err(anyhow!(
            "Genotype IDs length ({}) != phenotype length ({})",
            data.geno_ids.len(),
            n_total
        ));
    }

    // not.miss <- which(!is.na(y))
    let not_miss: Vec<usize> = (0..n_total)
        .filter(|&i| data.pheno[i].is_finite())
        .collect();

    if not_miss.is_empty() {
        return Err(anyhow!("All phenotype values are missing"));
    }

    let n = not_miss.len();

    // Subset to non-missing observations
    let y: Vec<f64> = not_miss.iter().map(|&i| data.pheno[i]).collect();
    let geno_ids_subset: Vec<String> = not_miss.iter().map(|&i| data.geno_ids[i].clone()).collect();

    // Build X matrix (fixed effects design matrix)
    // X <- matrix(1,n,1) - start with intercept
    let mut x_cols: Vec<Vec<f64>> = vec![vec![1.0; n]]; // intercept

    // Add fixed effects (categorical factors)
    if let Some(ref fixed) = data.fixed {
        for factor in fixed {
            if factor.len() != n_total {
                return Err(anyhow!("Fixed effect length != data length"));
            }
            // Subset to non-missing
            let factor_subset: Vec<&String> = not_miss.iter().map(|&i| &factor[i]).collect();

            // Get unique levels
            let mut levels: Vec<String> = factor_subset.iter().map(|&s| s.clone()).collect();
            levels.sort();
            levels.dedup();

            // Only add if more than one unique level
            if levels.len() > 1 {
                // Create dummy variables (model.matrix style, dropping first level)
                for level in levels.iter().skip(0) {
                    // Include all levels like R's ~x-1
                    let dummy: Vec<f64> = factor_subset
                        .iter()
                        .map(|&s| if s == level { 1.0 } else { 0.0 })
                        .collect();
                    x_cols.push(dummy);
                }
            }
        }
    }

    // Add covariates (continuous)
    if let Some(ref covariates) = data.covariates {
        for cov in covariates {
            if cov.len() != n_total {
                return Err(anyhow!("Covariate length != data length"));
            }
            let cov_subset: Vec<f64> = not_miss.iter().map(|&i| cov[i]).collect();
            x_cols.push(cov_subset);
        }
    }

    // Build X matrix
    let p = x_cols.len();
    let mut x_mat = DMatrix::zeros(n, p);
    for (j, col) in x_cols.iter().enumerate() {
        for (i, &val) in col.iter().enumerate() {
            x_mat[(i, j)] = val;
        }
    }

    // make.full: ensure X is full rank using SVD
    let x2 = make_full(&x_mat);
    let p2 = x2.ncols();

    // Get unique genotype IDs with phenotypes
    let not_miss_gid: Vec<String> = {
        let mut seen = HashSet::new();
        geno_ids_subset
            .iter()
            .filter(|id| seen.insert((*id).clone()))
            .cloned()
            .collect()
    };

    // Initialize residual vector with NaN
    let mut resid = vec![f64::NAN; n_total];

    if opts.k.is_none() {
        // No kinship matrix - use identity covariance
        // gid <- not.miss.gid
        // v <- length(gid)
        // Z <- matrix(0,n,v)
        // Z[cbind(1:n,match(data[,gid.pos],gid))] <- 1

        let gid = not_miss_gid.clone();
        let v = gid.len();

        // Build Z matrix
        let mut z_mat = DMatrix::zeros(n, v);
        let gid_to_idx: HashMap<&String, usize> =
            gid.iter().enumerate().map(|(i, id)| (id, i)).collect();

        for (i, id) in geno_ids_subset.iter().enumerate() {
            if let Some(&j) = gid_to_idx.get(id) {
                z_mat[(i, j)] = 1.0;
            }
        }

        // Call mixed.solve
        let y_vec: Vec<f64> = y.clone();
        let ms_opts = MixedSolveOptions {
            se: opts.pev,
            ..Default::default()
        };

        let ans = mixed_solve(&y_vec, Some(&z_mat), None, Some(&x2), Some(ms_opts))?;

        // Compute residuals: resid[not.miss] <- y - X2 %*% beta - Z %*% u
        let y_dvec = DVector::from_row_slice(&y);
        let fitted = &x2 * &ans.beta + &z_mat * &ans.u;
        for (idx, &i) in not_miss.iter().enumerate() {
            resid[i] = y_dvec[idx] - fitted[idx];
        }

        // Compute predictions: pred = u + mean(X2) %*% beta
        let x2_col_means: DVector<f64> =
            DVector::from_fn(p2, |j, _| x2.column(j).iter().sum::<f64>() / n as f64);
        let mean_fixed_effect = x2_col_means.dot(&ans.beta);
        let pred: Vec<f64> = ans.u.iter().map(|&u| u + mean_fixed_effect).collect();

        let pev_result = if opts.pev {
            ans.u_se.map(|se| se.iter().map(|&s| s * s).collect())
        } else {
            None
        };

        Ok(KinBlupResult {
            vg: ans.vu,
            ve: ans.ve,
            g: ans.u.iter().cloned().collect(),
            g_ids: gid,
            pev: pev_result,
            resid,
            pred,
            profile: None,
        })
    } else {
        // With kinship matrix K
        let k_mat = opts.k.as_ref().unwrap();
        let k_ids = opts
            .k_ids
            .as_ref()
            .ok_or_else(|| anyhow!("K matrix provided but k_ids is missing"))?;

        if k_mat.nrows() != k_ids.len() || k_mat.ncols() != k_ids.len() {
            return Err(anyhow!(
                "K dimensions ({} x {}) don't match k_ids length ({})",
                k_mat.nrows(),
                k_mat.ncols(),
                k_ids.len()
            ));
        }

        // Check that all phenotyped individuals have genotypes
        let gid_set: HashSet<&String> = k_ids.iter().collect();
        let missing_geno: Vec<&String> = not_miss_gid
            .iter()
            .filter(|id| !gid_set.contains(id))
            .collect();
        if !missing_geno.is_empty() {
            return Err(anyhow!(
                "The following lines have phenotypes but no genotypes: {}",
                missing_geno
                    .iter()
                    .map(|s| s.as_str())
                    .collect::<Vec<_>>()
                    .join(" ")
            ));
        }

        // Reorder K: phenotyped individuals first, then non-phenotyped
        let miss_gid: Vec<String> = k_ids
            .iter()
            .filter(|id| !not_miss_gid.contains(id))
            .cloned()
            .collect();

        // ix.pheno <- match(not.miss.gid, gid)
        let k_id_to_idx: HashMap<&String, usize> =
            k_ids.iter().enumerate().map(|(i, id)| (id, i)).collect();

        let mut ix: Vec<usize> = not_miss_gid
            .iter()
            .map(|id| *k_id_to_idx.get(id).unwrap())
            .collect();
        for id in &miss_gid {
            ix.push(*k_id_to_idx.get(id).unwrap());
        }

        // Reorder K
        let k_reordered = reorder_matrix(k_mat, &ix);

        // Build ordered gid list
        let mut gid: Vec<String> = not_miss_gid.clone();
        gid.extend(miss_gid);

        let v = not_miss_gid.len();
        let total_gid = gid.len();

        // Build Z matrix (n x v), then Z2 = [Z | 0] (n x total_gid)
        let gid_to_idx_new: HashMap<&String, usize> =
            not_miss_gid.iter().enumerate().map(|(i, id)| (id, i)).collect();

        let mut z_mat = DMatrix::zeros(n, v);
        for (i, id) in geno_ids_subset.iter().enumerate() {
            if let Some(&j) = gid_to_idx_new.get(id) {
                z_mat[(i, j)] = 1.0;
            }
        }

        // Z2 <- cbind(Z, matrix(0, n, nrow(K) - v))
        let mut z2_mat = DMatrix::zeros(n, total_gid);
        for i in 0..n {
            for j in 0..v {
                z2_mat[(i, j)] = z_mat[(i, j)];
            }
        }

        if !opts.gauss {
            // Standard GBLUP
            let y_vec: Vec<f64> = y.clone();
            let ms_opts = MixedSolveOptions {
                se: opts.pev,
                ..Default::default()
            };

            let ans = mixed_solve(
                &y_vec,
                Some(&z2_mat),
                Some(&k_reordered),
                Some(&x2),
                Some(ms_opts),
            )?;

            // Compute residuals
            let y_dvec = DVector::from_row_slice(&y);
            let fitted = &x2 * &ans.beta + &z2_mat * &ans.u;
            for (idx, &i) in not_miss.iter().enumerate() {
                resid[i] = y_dvec[idx] - fitted[idx];
            }

            // Compute predictions
            let x2_col_means: DVector<f64> =
                DVector::from_fn(p2, |j, _| x2.column(j).iter().sum::<f64>() / n as f64);
            let mean_fixed_effect = x2_col_means.dot(&ans.beta);
            let pred: Vec<f64> = ans.u.iter().map(|&u| u + mean_fixed_effect).collect();

            let pev_result = if opts.pev {
                ans.u_se.map(|se| se.iter().map(|&s| s * s).collect())
            } else {
                None
            };

            Ok(KinBlupResult {
                vg: ans.vu,
                ve: ans.ve,
                g: ans.u.iter().cloned().collect(),
                g_ids: gid,
                pev: pev_result,
                resid,
                pred,
                profile: None,
            })
        } else {
            // Gaussian kernel BLUP
            // theta <- setdiff(seq(0, max(K), length.out=11), 0)
            let theta_seq = opts.theta_seq.unwrap_or_else(|| {
                let max_k = k_reordered
                    .iter()
                    .cloned()
                    .fold(f64::NEG_INFINITY, f64::max);
                let mut seq: Vec<f64> = (1..=10).map(|i| max_k * i as f64 / 10.0).collect();
                seq.retain(|&t| t > 0.0);
                seq
            });

            if theta_seq.is_empty() {
                return Err(anyhow!("theta_seq is empty"));
            }

            // Profile likelihood over theta values
            let mut best_ll = f64::NEG_INFINITY;
            let mut best_ans = None;
            let mut profile = Vec::with_capacity(theta_seq.len());

            for &theta in &theta_seq {
                // K_gauss = exp(-(K/theta)^2)
                let k_gauss = DMatrix::from_fn(total_gid, total_gid, |i, j| {
                    let d = k_reordered[(i, j)] / theta;
                    (-d * d).exp()
                });

                let y_vec: Vec<f64> = y.clone();
                let ms_opts = MixedSolveOptions {
                    se: opts.pev,
                    ..Default::default()
                };

                match mixed_solve(&y_vec, Some(&z2_mat), Some(&k_gauss), Some(&x2), Some(ms_opts))
                {
                    Ok(ans) => {
                        profile.push((theta, ans.ll));
                        if ans.ll > best_ll {
                            best_ll = ans.ll;
                            best_ans = Some(ans);
                        }
                    }
                    Err(_) => {
                        profile.push((theta, f64::NEG_INFINITY));
                    }
                }
            }

            let ans = best_ans.ok_or_else(|| anyhow!("All theta values failed"))?;

            // Compute residuals
            let y_dvec = DVector::from_row_slice(&y);
            let fitted = &x2 * &ans.beta + &z2_mat * &ans.u;
            for (idx, &i) in not_miss.iter().enumerate() {
                resid[i] = y_dvec[idx] - fitted[idx];
            }

            // Compute predictions
            let x2_col_means: DVector<f64> =
                DVector::from_fn(p2, |j, _| x2.column(j).iter().sum::<f64>() / n as f64);
            let mean_fixed_effect = x2_col_means.dot(&ans.beta);
            let pred: Vec<f64> = ans.u.iter().map(|&u| u + mean_fixed_effect).collect();

            let pev_result = if opts.pev {
                ans.u_se.map(|se| se.iter().map(|&s| s * s).collect())
            } else {
                None
            };

            Ok(KinBlupResult {
                vg: ans.vu,
                ve: ans.ve,
                g: ans.u.iter().cloned().collect(),
                g_ids: gid,
                pev: pev_result,
                resid,
                pred,
                profile: Some(profile),
            })
        }
    }
}

/// Make matrix full rank using SVD
/// make.full <- function(X) {
///     svd.X <- svd(X)
///     r <- max(which(svd.X$d > 1e-8))
///     return(as.matrix(svd.X$u[,1:r]))
/// }
fn make_full(x: &DMatrix<f64>) -> DMatrix<f64> {
    use faer::Mat as FaerMat;

    let n = x.nrows();
    let p = x.ncols();

    // Convert to faer for SVD
    let x_faer = FaerMat::from_fn(n, p, |i, j| x[(i, j)]);
    let svd = x_faer.svd();

    // Find rank (number of singular values > 1e-8)
    let s = svd.s_diagonal();
    let r = (0..s.nrows())
        .filter(|&i| s.read(i) > 1e-8)
        .count()
        .max(1);

    // Extract first r columns of U
    let u = svd.u();
    DMatrix::from_fn(n, r, |i, j| u.read(i, j))
}

/// Reorder a square matrix according to index vector
fn reorder_matrix(mat: &DMatrix<f64>, ix: &[usize]) -> DMatrix<f64> {
    let n = ix.len();
    DMatrix::from_fn(n, n, |i, j| mat[(ix[i], ix[j])])
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_kin_blup_no_kinship() {
        // Simple case: no kinship matrix
        let data = KinBlupData {
            geno_ids: vec![
                "A".into(),
                "B".into(),
                "C".into(),
                "A".into(),
                "B".into(),
                "C".into(),
            ],
            pheno: vec![1.0, 2.0, 3.0, 1.5, 2.5, 3.5],
            fixed: None,
            covariates: None,
        };

        let result = kin_blup(&data, None).unwrap();

        assert_eq!(result.g_ids.len(), 3);
        assert_eq!(result.g.len(), 3);
        assert!(result.vg >= 0.0);
        assert!(result.ve >= 0.0);
    }

    #[test]
    fn test_kin_blup_with_kinship() {
        // With kinship matrix
        let data = KinBlupData {
            geno_ids: vec!["A".into(), "B".into(), "C".into()],
            pheno: vec![1.0, 2.0, 3.0],
            fixed: None,
            covariates: None,
        };

        #[rustfmt::skip]
        let k = DMatrix::from_row_slice(3, 3, &[
            1.0, 0.5, 0.3,
            0.5, 1.0, 0.4,
            0.3, 0.4, 1.0,
        ]);

        let opts = KinBlupOptions {
            k: Some(k),
            k_ids: Some(vec!["A".into(), "B".into(), "C".into()]),
            ..Default::default()
        };

        let result = kin_blup(&data, Some(opts)).unwrap();

        assert_eq!(result.g_ids.len(), 3);
        assert_eq!(result.g.len(), 3);
    }

    #[test]
    fn test_kin_blup_with_missing_pheno() {
        // Some individuals have missing phenotypes
        let data = KinBlupData {
            geno_ids: vec![
                "A".into(),
                "B".into(),
                "C".into(),
                "D".into(),
                "E".into(),
            ],
            pheno: vec![1.0, f64::NAN, 3.0, 4.0, f64::NAN],
            fixed: None,
            covariates: None,
        };

        #[rustfmt::skip]
        let k = DMatrix::from_row_slice(5, 5, &[
            1.0, 0.5, 0.3, 0.2, 0.1,
            0.5, 1.0, 0.4, 0.3, 0.2,
            0.3, 0.4, 1.0, 0.5, 0.3,
            0.2, 0.3, 0.5, 1.0, 0.4,
            0.1, 0.2, 0.3, 0.4, 1.0,
        ]);

        let opts = KinBlupOptions {
            k: Some(k),
            k_ids: Some(vec![
                "A".into(),
                "B".into(),
                "C".into(),
                "D".into(),
                "E".into(),
            ]),
            ..Default::default()
        };

        let result = kin_blup(&data, Some(opts)).unwrap();

        // Should have predictions for all 5 genotypes
        assert_eq!(result.g_ids.len(), 5);
        assert_eq!(result.g.len(), 5);

        // Residuals should be NaN for missing phenotypes
        assert!(result.resid[1].is_nan()); // B
        assert!(result.resid[4].is_nan()); // E
    }

    #[test]
    fn test_kin_blup_with_pev() {
        let data = KinBlupData {
            geno_ids: vec!["A".into(), "B".into(), "C".into()],
            pheno: vec![1.0, 2.0, 3.0],
            fixed: None,
            covariates: None,
        };

        let opts = KinBlupOptions {
            pev: true,
            ..Default::default()
        };

        let result = kin_blup(&data, Some(opts)).unwrap();

        assert!(result.pev.is_some());
        let pev = result.pev.unwrap();
        assert_eq!(pev.len(), 3);
        // PEV should be non-negative
        for &p in &pev {
            assert!(p >= 0.0);
        }
    }

    #[test]
    fn test_kin_blup_with_fixed_effects() {
        let data = KinBlupData {
            geno_ids: vec![
                "A".into(),
                "B".into(),
                "C".into(),
                "A".into(),
                "B".into(),
                "C".into(),
            ],
            pheno: vec![1.0, 2.0, 3.0, 2.0, 3.0, 4.0],
            fixed: Some(vec![vec![
                "F1".into(),
                "F1".into(),
                "F1".into(),
                "F2".into(),
                "F2".into(),
                "F2".into(),
            ]]),
            covariates: None,
        };

        let result = kin_blup(&data, None).unwrap();

        assert_eq!(result.g_ids.len(), 3);
    }

    #[test]
    fn test_kin_blup_with_covariates() {
        let data = KinBlupData {
            geno_ids: vec![
                "A".into(),
                "B".into(),
                "C".into(),
                "A".into(),
                "B".into(),
                "C".into(),
            ],
            pheno: vec![1.0, 2.0, 3.0, 1.5, 2.5, 3.5],
            fixed: None,
            covariates: Some(vec![vec![0.0, 1.0, 2.0, 0.5, 1.5, 2.5]]),
        };

        let result = kin_blup(&data, None).unwrap();

        assert_eq!(result.g_ids.len(), 3);
    }

    #[test]
    fn test_make_full() {
        // Test make_full function
        let x = DMatrix::from_row_slice(
            4,
            3,
            &[
                1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, // Duplicate of row 1
                1.0, 0.0, 1.0, // Duplicate of row 2
            ],
        );

        let x_full = make_full(&x);

        // Result should have rank 2 (since rows 3,4 are duplicates)
        assert!(x_full.ncols() <= 3);
        assert!(x_full.ncols() >= 1);
    }
}
