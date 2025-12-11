//! rrblup-rs: Rust implementation of R/rrBLUP package
//!
//! This crate provides a Rust implementation of the core functionality from the
//! R/rrBLUP package for mixed model analysis, including:
//!
//! - `a_mat()`: Additive relationship matrix from markers (A.mat in R)
//! - `mixed_solve()`: REML-based mixed model solver (mixed.solve in R)
//! - `kin_blup()`: Genomic prediction with kinship matrix (kin.blup in R)
//!
//! # Example
//!
//! ```ignore
//! use rrblup_rs::mixed_solve::{mixed_solve, MixedSolveOptions};
//! use rrblup_rs::a_mat::{a_mat, AMatOptions};
//!
//! // Compute relationship matrix
//! let geno = DMatrix::from_row_slice(3, 5, &[...]); // {-1, 0, 1} coded
//! let a_result = a_mat(&geno, None)?;
//!
//! // Fit mixed model
//! let y = vec![1.0, 2.0, 3.0];
//! let result = mixed_solve(&y, None, Some(&a_result.a), None, None)?;
//! println!("Vu = {}, Ve = {}", result.vu, result.ve);
//! ```

// Module for mixed.solve function
pub mod mixed_solve;

// Module for A.mat function
pub mod a_mat;

// Module for kin.blup function
pub mod kin_blup;

// Re-export main types and functions from mixed_solve
pub use mixed_solve::{
    mixed_solve as mixed_solve_reml, Method, MixedSolveOptions, MixedSolveResult,
};

// Re-export main types and functions from a_mat
pub use a_mat::{a_mat, AMatOptions, AMatResult, ImputeMethod, ShrinkConfig, ShrinkMethod};

// Re-export main types and functions from kin_blup
pub use kin_blup::{kin_blup, KinBlupData, KinBlupOptions, KinBlupResult};

use anyhow::{anyhow, Result};
use faer::Mat as FaerMat;
use nalgebra::{DMatrix, DVector, SymmetricEigen};
use ndarray::{Array1, Array2};
use statrs::distribution::{ContinuousCDF, FisherSnedecor};

/// Cached quantities from fitting the null mixed model.
/// Used for efficient marker testing via P3D (population parameters previously determined).
#[derive(Debug, Clone)]
pub struct MixedModelCache {
    pub sample_ids: Vec<String>,
    pub n_obs: usize,
    pub n_ind: usize,
    pub h_inv: DMatrix<f64>,
    pub sigma2: f64,
    pub lambda: f64,
    pub vu: f64,
    pub ve: f64,
}

/// Solve mixed model y = Xb + Zu + e using REML spectral decomposition.
///
/// This is a Rust implementation of R/rrBLUP::mixed.solve().
///
/// # Arguments
/// * `y` - Response vector (n x 1)
/// * `x` - Fixed effects design matrix (n x p), defaults to intercept if None
/// * `z` - Random effects design matrix (n x m), defaults to identity if None
/// * `k` - Covariance matrix for random effects (m x m), defaults to identity if None
/// * `return_hinv` - Whether to compute and return H inverse
///
/// # Returns
/// * MixedSolveResult with variance components and BLUPs
///
/// # Note
/// This is a legacy interface. Prefer using `mixed_solve::mixed_solve` directly
/// which provides more options and matches R/rrBLUP more closely.
pub fn mixed_solve_legacy(
    y: &Array1<f64>,
    x: Option<&Array2<f64>>,
    z: Option<&Array2<f64>>,
    k: Option<&Array2<f64>>,
    return_hinv: bool,
) -> Result<LegacyMixedSolveResult> {
    let n_full = y.len();

    // Handle missing values (NA in R)
    let not_na: Vec<usize> = (0..n_full).filter(|&i| y[i].is_finite()).collect();
    let n = not_na.len();

    if n == 0 {
        return Err(anyhow!("All y values are NA"));
    }

    // Subset y to non-NA
    let y_vec: DVector<f64> = DVector::from_iterator(n, not_na.iter().map(|&i| y[i]));

    // Default X to intercept
    let x_default = Array2::from_shape_fn((n_full, 1), |_| 1.0);
    let x_full = x.unwrap_or(&x_default);
    let p = x_full.ncols();

    // Subset X to non-NA rows
    let x_mat = DMatrix::from_fn(n, p, |i, j| x_full[(not_na[i], j)]);

    // Default Z to identity
    let m = match z {
        Some(z_arr) => z_arr.ncols(),
        None => n,
    };

    // Subset Z to non-NA rows
    let z_mat = match z {
        Some(z_arr) => DMatrix::from_fn(n, m, |i, j| z_arr[(not_na[i], j)]),
        None => DMatrix::identity(n, m),
    };

    // K matrix (m x m)
    let k_mat: Option<DMatrix<f64>> = k.map(|k_arr| {
        DMatrix::from_fn(k_arr.nrows(), k_arr.ncols(), |i, j| k_arr[(i, j)])
    });

    // Check dimensions
    if z_mat.nrows() != n {
        return Err(anyhow!("Z rows ({}) != n ({})", z_mat.nrows(), n));
    }
    if x_mat.nrows() != n {
        return Err(anyhow!("X rows ({}) != n ({})", x_mat.nrows(), n));
    }
    if let Some(ref k) = k_mat {
        if k.nrows() != m || k.ncols() != m {
            return Err(anyhow!("K must be {} x {}", m, m));
        }
    }

    // XtX and rank check
    let xtx = x_mat.transpose() * &x_mat;
    let xtx_inv = xtx
        .clone()
        .try_inverse()
        .ok_or_else(|| anyhow!("X not full rank"))?;

    // S = I - X(X'X)^-1 X' (projection onto space orthogonal to X)
    let s_mat = DMatrix::identity(n, n) - &x_mat * &xtx_inv * x_mat.transpose();

    // Determine spectral method
    let use_cholesky = n > m + p;

    let (phi, theta, u_mat, omega_sq): (Vec<f64>, Vec<f64>, DMatrix<f64>, Vec<f64>);

    if use_cholesky {
        // CHOLESKY BRANCH
        let zbt = if let Some(ref k) = k_mat {
            let mut k_jittered = k.clone();
            for i in 0..m {
                k_jittered[(i, i)] += 1e-6;
            }
            let chol = k_jittered
                .cholesky()
                .ok_or_else(|| anyhow!("K not positive semi-definite"))?;
            let b = chol.l();
            &z_mat * b.transpose()
        } else {
            z_mat.clone()
        };

        let zbt_faer = nalgebra_to_faer(&zbt);
        let svd_zbt = zbt_faer.svd();
        let u_full = faer_to_nalgebra(&svd_zbt.u());
        let d_vals = svd_zbt.s_diagonal();

        u_mat = u_full;

        phi = (0..n)
            .map(|i| {
                if i < d_vals.nrows() {
                    let d = d_vals.read(i);
                    d * d
                } else {
                    0.0
                }
            })
            .collect();

        let szbt = &s_mat * &zbt;
        let szbt_faer = nalgebra_to_faer(&szbt);
        let svd_szbt = szbt_faer.thin_svd();
        let u_szbt = faer_to_nalgebra(&svd_szbt.u());
        let d_szbt = svd_szbt.s_diagonal();

        let n_u_cols = u_szbt.ncols();
        let mut combined = DMatrix::zeros(n, p + n_u_cols);
        for i in 0..n {
            for j in 0..p {
                combined[(i, j)] = x_mat[(i, j)];
            }
            for j in 0..n_u_cols {
                combined[(i, p + j)] = u_szbt[(i, j)];
            }
        }

        let combined_faer = nalgebra_to_faer(&combined);
        let qr = combined_faer.qr();
        let q_full = faer_to_nalgebra(&qr.compute_q().as_ref());
        let q_complement = q_full.columns(p, n - p).into_owned();
        let r_faer = qr.compute_r();
        let r_mat = faer_mat_to_nalgebra(&r_faer);

        let r22_rows = m.min(r_mat.nrows().saturating_sub(p));
        let r22_cols = m.min(r_mat.ncols().saturating_sub(p));
        let r22_size = r22_rows.min(r22_cols);

        if r22_size > 0 && d_szbt.nrows() > 0 {
            let mut r22_sq = DMatrix::zeros(r22_size, r22_size);
            for i in 0..r22_size {
                for j in 0..r22_size {
                    let val = r_mat[(p + i, p + j)];
                    r22_sq[(i, j)] = val * val;
                }
            }

            let t_r22_sq = r22_sq.transpose();
            let d_sq_len = r22_size.min(d_szbt.nrows());
            let d_sq: Vec<f64> = (0..d_sq_len)
                .map(|i| {
                    let d = d_szbt.read(i);
                    d * d
                })
                .collect();
            let d_sq_vec = DVector::from_row_slice(&d_sq);

            let theta_sub = t_r22_sq.clone().try_inverse().map(|inv| inv * &d_sq_vec);

            match theta_sub {
                Some(ts) => {
                    let n_theta = n - p;
                    theta = (0..n_theta)
                        .map(|i| {
                            if i < ts.len() {
                                ts[i].max(0.0)
                            } else {
                                0.0
                            }
                        })
                        .collect();
                }
                None => {
                    theta = vec![0.0; n - p];
                }
            }
        } else {
            theta = vec![0.0; n - p];
        }

        let omega_vec = q_complement.transpose() * &y_vec;
        omega_sq = omega_vec.iter().map(|v| v * v).collect();
    } else {
        // EIGEN BRANCH
        let offset = (n as f64).sqrt();

        let hb = if let Some(ref k) = k_mat {
            let zk = &z_mat * k;
            let zkzt = &zk * z_mat.transpose();
            let mut hb = zkzt;
            for i in 0..n {
                hb[(i, i)] += offset;
            }
            hb
        } else {
            let zzt = &z_mat * z_mat.transpose();
            let mut hb = zzt;
            for i in 0..n {
                hb[(i, i)] += offset;
            }
            hb
        };

        let hb_eig = SymmetricEigen::new(hb.clone());
        phi = hb_eig.eigenvalues.iter().map(|v| v - offset).collect();

        if phi.iter().cloned().fold(f64::INFINITY, f64::min) < -1e-6 {
            return Err(anyhow!("K not positive semi-definite (phi < 0)"));
        }

        u_mat = hb_eig.eigenvectors.clone();

        let shbs = &s_mat * &hb * &s_mat;
        let shbs_eig = SymmetricEigen::new(shbs);

        let n_theta = n - p;
        theta = shbs_eig
            .eigenvalues
            .iter()
            .take(n_theta)
            .map(|v| (v - offset).max(0.0))
            .collect();

        let q = shbs_eig.eigenvectors.columns(0, n_theta).into_owned();
        let omega_vec = q.transpose() * &y_vec;
        omega_sq = omega_vec.iter().map(|v| v * v).collect();
    }

    // REML optimization
    let df = n - p;
    let bounds = (1e-9f64, 1e9f64);

    let reml_obj = |lambda: f64| -> f64 {
        if lambda <= 0.0 {
            return f64::INFINITY;
        }
        let sum_ratio: f64 = omega_sq
            .iter()
            .zip(theta.iter())
            .map(|(o, t)| o / (t + lambda))
            .sum();
        if sum_ratio <= 0.0 {
            return f64::INFINITY;
        }
        let sum_log: f64 = theta.iter().map(|t| (t + lambda).ln()).sum();
        (df as f64) * sum_ratio.ln() + sum_log
    };

    // Golden section search
    let (mut a, mut b) = bounds;
    let gr = 0.5 * (1.0 + 5f64.sqrt());
    let mut c = b - (b - a) / gr;
    let mut d = a + (b - a) / gr;
    let mut fc = reml_obj(c);
    let mut fd = reml_obj(d);

    for _ in 0..100 {
        if (b - a).abs() < 1e-8 {
            break;
        }
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - (b - a) / gr;
            fc = reml_obj(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + (b - a) / gr;
            fd = reml_obj(d);
        }
    }

    let lambda_opt = if fc < fd { c } else { d };

    let vu_opt: f64 = omega_sq
        .iter()
        .zip(theta.iter())
        .map(|(o, t)| o / (t + lambda_opt))
        .sum::<f64>()
        / (df as f64);

    let ve_opt = lambda_opt * vu_opt;

    let h_inv = if return_hinv {
        let mut hinv = DMatrix::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                let mut sum = 0.0;
                for k in 0..n {
                    sum += u_mat[(i, k)] * u_mat[(j, k)] / (phi[k] + lambda_opt);
                }
                hinv[(i, j)] = sum;
            }
        }
        Some(hinv)
    } else {
        None
    };

    let h_inv_for_beta = if let Some(ref hinv) = h_inv {
        hinv.clone()
    } else {
        let mut hinv = DMatrix::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                let mut sum = 0.0;
                for k in 0..n {
                    sum += u_mat[(i, k)] * u_mat[(j, k)] / (phi[k] + lambda_opt);
                }
                hinv[(i, j)] = sum;
            }
        }
        hinv
    };

    let hinv_x = &h_inv_for_beta * &x_mat;
    let w = x_mat.transpose() * &hinv_x;
    let w_inv = w
        .clone()
        .try_inverse()
        .ok_or_else(|| anyhow!("W not invertible"))?;
    let hinv_y = &h_inv_for_beta * &y_vec;
    let beta = &w_inv * (x_mat.transpose() * &hinv_y);

    let resid = &y_vec - &x_mat * &beta;
    let hinv_resid = &h_inv_for_beta * &resid;

    let u_blup = if let Some(ref k) = k_mat {
        let kzt = k * z_mat.transpose();
        vu_opt * &kzt * &hinv_resid
    } else {
        vu_opt * z_mat.transpose() * &hinv_resid
    };

    let obj_val = reml_obj(lambda_opt);
    let ll = -0.5
        * (obj_val
            + (df as f64)
            + (df as f64) * (2.0 * std::f64::consts::PI / (df as f64)).ln());

    Ok(LegacyMixedSolveResult {
        vu: vu_opt,
        ve: ve_opt,
        beta: beta.iter().cloned().collect(),
        u: u_blup.iter().cloned().collect(),
        ll,
        lambda: lambda_opt,
        h_inv,
    })
}

/// Result from legacy mixed model solving interface
#[derive(Debug, Clone)]
pub struct LegacyMixedSolveResult {
    /// Variance of random effects (Vu)
    pub vu: f64,
    /// Residual variance (Ve)
    pub ve: f64,
    /// Fixed effects coefficients (beta)
    pub beta: Vec<f64>,
    /// Random effects BLUPs (u)
    pub u: Vec<f64>,
    /// Log-likelihood (LL)
    pub ll: f64,
    /// Ridge parameter (lambda = Ve/Vu)
    pub lambda: f64,
    /// Inverse of H matrix (optional, for GWAS)
    pub h_inv: Option<DMatrix<f64>>,
}

/// Score test for a single marker using pre-computed cache.
pub fn score_test_marker(
    marker_design: &[Vec<f64>],
    y: &Array1<f64>,
    x0: &Array2<f64>,
    cache: &MixedModelCache,
) -> Result<(f64, f64, Option<f64>)> {
    let n = y.len();
    let p0 = x0.ncols();
    let p_marker = marker_design.len();

    if p_marker == 0 {
        return Err(anyhow!("Empty marker design"));
    }

    // Build full design X = [X0 | marker_design]
    let p = p0 + p_marker;
    let mut data = Vec::with_capacity(n * p);
    for row in 0..n {
        for col in 0..p0 {
            data.push(x0[(row, col)]);
        }
        for col in marker_design {
            data.push(col[row]);
        }
    }
    let x = DMatrix::from_row_slice(n, p, &data);

    let y_vec = DVector::from_row_slice(
        y.as_slice()
            .ok_or_else(|| anyhow!("y not contiguous"))?,
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
            w_eps
                .try_inverse()
                .ok_or_else(|| anyhow!("W not invertible"))?
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
    let cov_beta_marker = s2 * w_inv.slice((p0, p0), (p_marker, p_marker)).into_owned();

    let f_stat = if p_marker == 1 {
        let cov = cov_beta_marker[(0, 0)];
        if cov > 0.0 {
            marker_beta[0] * marker_beta[0] / cov
        } else {
            0.0
        }
    } else {
        let cov_inv = cov_beta_marker.clone().try_inverse().unwrap_or_else(|| {
            let mut jittered = cov_beta_marker.clone();
            for i in 0..p_marker {
                jittered[(i, i)] += 1e-8;
            }
            jittered
                .try_inverse()
                .unwrap_or(DMatrix::identity(p_marker, p_marker))
        });
        (marker_beta.transpose() * &cov_inv * &marker_beta)[(0, 0)] / (p_marker as f64)
    };

    let p_value = if f_stat.is_finite() && f_stat >= 0.0 {
        let f_dist = FisherSnedecor::new(p_marker as f64, v2)?;
        (1.0 - f_dist.cdf(f_stat)).max(0.0)
    } else {
        1.0
    };

    let effect = if p_marker == 1 {
        Some(marker_beta[0])
    } else {
        None
    };

    Ok((f_stat, p_value, effect))
}

// Helper functions for matrix conversions

fn nalgebra_to_faer(m: &DMatrix<f64>) -> FaerMat<f64> {
    let nrows = m.nrows();
    let ncols = m.ncols();
    FaerMat::from_fn(nrows, ncols, |i, j| m[(i, j)])
}

fn faer_to_nalgebra(m: &faer::MatRef<f64>) -> DMatrix<f64> {
    let nrows = m.nrows();
    let ncols = m.ncols();
    DMatrix::from_fn(nrows, ncols, |i, j| m.read(i, j))
}

fn faer_mat_to_nalgebra(m: &FaerMat<f64>) -> DMatrix<f64> {
    let nrows = m.nrows();
    let ncols = m.ncols();
    DMatrix::from_fn(nrows, ncols, |i, j| m.read(i, j))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::array;

    #[test]
    fn test_mixed_solve_basic() {
        let y = array![1.0, 2.0, 3.0, 4.0, 5.0];
        let result = mixed_solve_legacy(&y, None, None, None, false).unwrap();

        assert!(result.vu >= 0.0);
        assert!(result.ve > 0.0);
        assert_relative_eq!(result.beta[0], 3.0, epsilon = 1.0);
    }
}
