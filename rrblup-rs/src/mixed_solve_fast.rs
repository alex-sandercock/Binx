//! Fast mixed model solver using faer throughout
//!
//! This module provides an optimized version of mixed.solve that:
//! 1. Uses faer matrices throughout (no nalgebra conversion overhead)
//! 2. Uses faer's SIMD-optimized eigendecomposition
//! 3. Leverages faer's cache-friendly matrix operations
//!
//! API is compatible with mixed_solve.rs for drop-in replacement.

use anyhow::{anyhow, Result};
use faer::prelude::*;
use faer::solvers::SolverCore;
use faer::{Col, Mat};
use nalgebra::{DMatrix, DVector};

// Re-use types from mixed_solve
pub use crate::mixed_solve::{Method, MixedSolveOptions, MixedSolveResult};

/// Fast mixed.solve using faer matrices throughout
///
/// This is functionally identical to mixed_solve() but uses faer's
/// optimized linear algebra operations throughout.
pub fn mixed_solve_fast(
    y: &[f64],
    z: Option<&DMatrix<f64>>,
    k: Option<&DMatrix<f64>>,
    x: Option<&DMatrix<f64>>,
    options: Option<MixedSolveOptions>,
) -> Result<MixedSolveResult> {
    let opts = options.unwrap_or_default();
    let pi = std::f64::consts::PI;

    let n_full = y.len();

    // Handle missing values
    let not_na: Vec<usize> = (0..n_full).filter(|&i| y[i].is_finite()).collect();
    if not_na.is_empty() {
        return Err(anyhow!("All y values are NA"));
    }
    let n = not_na.len();

    // Convert inputs to faer matrices (only conversion point)
    let y_vec = Col::<f64>::from_fn(n, |i| y[not_na[i]]);

    let p = x.map(|m| m.ncols()).unwrap_or(1);
    let x_mat: Mat<f64> = match x {
        Some(x_na) => Mat::from_fn(n, p, |i, j| x_na[(not_na[i], j)]),
        None => Mat::from_fn(n, 1, |_, _| 1.0),
    };

    let m = z.map(|m| m.ncols()).unwrap_or(n);
    let z_mat: Mat<f64> = match z {
        Some(z_na) => Mat::from_fn(n, m, |i, j| z_na[(not_na[i], j)]),
        None => Mat::identity(n, n),
    };

    // K matrix (optional)
    let k_mat: Option<Mat<f64>> = k.map(|k_na| {
        Mat::from_fn(k_na.nrows(), k_na.ncols(), |i, j| k_na[(i, j)])
    });

    // Dimension checks
    if let Some(ref k_m) = k_mat {
        if k_m.nrows() != m || k_m.ncols() != m {
            return Err(anyhow!("K must be {} x {}", m, m));
        }
    }

    // XtX and inverse
    let xtx = x_mat.transpose() * &x_mat;
    let xtx_lu = xtx.partial_piv_lu();
    let xtx_inv = xtx_lu.inverse();

    // S = I - X @ (X'X)^-1 @ X'
    let x_xtxinv = &x_mat * &xtx_inv;
    let s_mat = Mat::<f64>::identity(n, n) - &x_xtxinv * x_mat.transpose();

    // Determine spectral method
    let use_eigen = n <= m + p;

    let phi: Vec<f64>;
    let theta: Vec<f64>;
    let u_mat: Mat<f64>;
    let omega_sq: Vec<f64>;

    if use_eigen {
        // ============================================================
        // EIGEN BRANCH: faer eigendecomposition
        // ============================================================
        let offset = (n as f64).sqrt();

        // Hb = Z @ K @ Z' + offset * I
        let mut hb: Mat<f64> = if let Some(ref k_m) = k_mat {
            let zk = &z_mat * k_m;
            &zk * z_mat.transpose()
        } else {
            &z_mat * z_mat.transpose()
        };
        for i in 0..n {
            hb.write(i, i, hb.read(i, i) + offset);
        }

        // Eigendecomposition of Hb
        let hb_eig = hb.selfadjoint_eigendecomposition(faer::Side::Lower);
        let hb_vals = hb_eig.s().column_vector();
        let hb_vecs = hb_eig.u();

        // faer returns ascending order, need descending
        let mut idx: Vec<usize> = (0..n).collect();
        idx.sort_by(|&a, &b| {
            hb_vals.read(b).partial_cmp(&hb_vals.read(a)).unwrap()
        });

        phi = idx.iter().map(|&i| hb_vals.read(i) - offset).collect();

        let min_phi = phi.iter().cloned().fold(f64::INFINITY, f64::min);
        if min_phi < -1e-6 {
            return Err(anyhow!("K not positive semi-definite"));
        }

        u_mat = Mat::from_fn(n, n, |i, j| hb_vecs.read(i, idx[j]));

        // SHbS eigendecomposition
        let shbs = &s_mat * &hb * &s_mat;
        let shbs_eig = shbs.selfadjoint_eigendecomposition(faer::Side::Lower);
        let shbs_vals = shbs_eig.s().column_vector();
        let shbs_vecs = shbs_eig.u();

        let mut shbs_idx: Vec<usize> = (0..n).collect();
        shbs_idx.sort_by(|&a, &b| {
            shbs_vals.read(b).partial_cmp(&shbs_vals.read(a)).unwrap()
        });

        let n_theta = n - p;
        theta = shbs_idx.iter()
            .take(n_theta)
            .map(|&i| shbs_vals.read(i) - offset)
            .collect();

        // Q = columns of SHbS eigenvectors
        let q_mat = Mat::from_fn(n, n_theta, |i, j| shbs_vecs.read(i, shbs_idx[j]));

        // omega = Q' @ y
        let omega = q_mat.transpose() * y_vec.as_ref();
        omega_sq = (0..n_theta).map(|i| {
            let v = omega.read(i);
            v * v
        }).collect();

    } else {
        // ============================================================
        // CHOLESKY BRANCH: SVD-based
        // ============================================================
        let zbt: Mat<f64> = if let Some(ref k_m) = k_mat {
            let mut k_jit = k_m.clone();
            for i in 0..m {
                k_jit.write(i, i, k_jit.read(i, i) + 1e-6);
            }
            match k_jit.cholesky(faer::Side::Lower) {
                Ok(chol) => {
                    let l = chol.compute_l();
                    &z_mat * l.transpose()
                }
                Err(_) => return Err(anyhow!("K not positive semi-definite")),
            }
        } else {
            z_mat.clone()
        };

        let svd_zbt = zbt.svd();
        u_mat = svd_zbt.u().to_owned();
        let d_vals = svd_zbt.s_diagonal();

        phi = (0..n).map(|i| {
            if i < d_vals.nrows() {
                let d = d_vals.read(i);
                d * d
            } else {
                0.0
            }
        }).collect();

        // SZBt and its SVD
        let szbt = &s_mat * &zbt;
        let svd_szbt = szbt.thin_svd();
        let u_szbt = svd_szbt.u();
        let d_szbt = svd_szbt.s_diagonal();

        // QR of [X | U_szbt]
        let n_u = u_szbt.ncols();
        let mut combined = Mat::<f64>::zeros(n, p + n_u);
        for i in 0..n {
            for j in 0..p {
                combined.write(i, j, x_mat.read(i, j));
            }
            for j in 0..n_u {
                combined.write(i, p + j, u_szbt.read(i, j));
            }
        }

        let qr = combined.qr();
        let q_full = qr.compute_q();
        let r_mat = qr.compute_r();

        // theta computation
        let r22_size = m.min(r_mat.nrows().saturating_sub(p)).min(r_mat.ncols().saturating_sub(p));
        let n_theta = n - p;

        theta = if r22_size > 0 && d_szbt.nrows() > 0 {
            let mut r22_sq = Mat::<f64>::zeros(r22_size, r22_size);
            for i in 0..r22_size {
                for j in 0..r22_size {
                    let v = r_mat.read(p + i, p + j);
                    r22_sq.write(i, j, v * v);
                }
            }

            let d_len = r22_size.min(d_szbt.nrows());
            let mut d_sq = Col::<f64>::zeros(d_len);
            for i in 0..d_len {
                let d = d_szbt.read(i);
                d_sq.write(i, d * d);
            }

            let t_r22_sq = r22_sq.transpose().to_owned();
            let ans = t_r22_sq.partial_piv_lu().solve(&d_sq);

            (0..n_theta).map(|i| {
                if i < ans.nrows() { ans.read(i) } else { 0.0 }
            }).collect()
        } else {
            vec![0.0; n_theta]
        };

        // omega = Q[p:]' @ y
        let q_compl = Mat::from_fn(n, n_theta, |i, j| q_full.read(i, p + j));
        let omega = q_compl.transpose() * y_vec.as_ref();
        omega_sq = (0..n_theta).map(|i| {
            let v = omega.read(i);
            v * v
        }).collect();
    }

    // ============================================================
    // REML/ML Optimization
    // ============================================================
    let df: usize;
    let lambda_opt: f64;
    let obj_val: f64;

    if opts.method == Method::ML {
        let f_ml = |lambda: f64| -> f64 {
            if lambda <= 0.0 { return f64::INFINITY; }
            let sum_ratio: f64 = omega_sq.iter().zip(&theta)
                .map(|(o, t)| o / (t + lambda)).sum();
            if sum_ratio <= 0.0 { return f64::INFINITY; }
            let sum_log: f64 = phi.iter().map(|p| (p + lambda).ln()).sum();
            (n as f64) * sum_ratio.ln() + sum_log
        };
        let (lo, ov) = golden_section_minimize(f_ml, opts.bounds.0, opts.bounds.1);
        lambda_opt = lo;
        obj_val = ov;
        df = n;
    } else {
        let n_p = n - p;
        let f_reml = |lambda: f64| -> f64 {
            if lambda <= 0.0 { return f64::INFINITY; }
            let sum_ratio: f64 = omega_sq.iter().zip(&theta)
                .map(|(o, t)| o / (t + lambda)).sum();
            if sum_ratio <= 0.0 { return f64::INFINITY; }
            let sum_log: f64 = theta.iter().map(|t| (t + lambda).ln()).sum();
            (n_p as f64) * sum_ratio.ln() + sum_log
        };
        let (lo, ov) = golden_section_minimize(f_reml, opts.bounds.0, opts.bounds.1);
        lambda_opt = lo;
        obj_val = ov;
        df = n_p;
    }

    // Variance components
    let vu_opt: f64 = omega_sq.iter().zip(&theta)
        .map(|(o, t)| o / (t + lambda_opt))
        .sum::<f64>() / (df as f64);
    let ve_opt = lambda_opt * vu_opt;

    // H_inv = U @ diag(1/(phi+lambda)) @ U'
    let mut hinv = Mat::<f64>::zeros(n, n);
    for i in 0..n {
        for j in 0..=i {
            let mut sum = 0.0;
            for kk in 0..n {
                sum += u_mat.read(i, kk) * u_mat.read(j, kk) / (phi[kk] + lambda_opt);
            }
            hinv.write(i, j, sum);
            hinv.write(j, i, sum);
        }
    }

    // beta = (X' H^-1 X)^-1 X' H^-1 y
    let hinv_x = &hinv * &x_mat;
    let w = x_mat.transpose() * &hinv_x;
    let w_inv = w.partial_piv_lu().inverse();
    let hinv_y = &hinv * y_vec.as_ref();
    let beta_faer = &w_inv * (x_mat.transpose() * &hinv_y);

    // u = K Z' H^-1 (y - X beta)
    let resid = y_vec.as_ref() - &x_mat * &beta_faer;
    let hinv_resid = &hinv * &resid;

    let u_blup_faer: Col<f64> = if let Some(ref k_m) = k_mat {
        let kzt = k_m * z_mat.transpose();
        &kzt * &hinv_resid
    } else {
        z_mat.transpose() * &hinv_resid
    };

    // Log-likelihood
    let ll = -0.5 * (obj_val + (df as f64) + (df as f64) * (2.0 * pi / (df as f64)).ln());

    // Convert output to nalgebra (only conversion at output)
    let beta = DVector::from_fn(p, |i, _| beta_faer.read(i));
    let u_blup = DVector::from_fn(u_blup_faer.nrows(), |i, _| u_blup_faer.read(i));

    // Standard errors
    let (beta_se, u_se) = if opts.se {
        let beta_se_vec = DVector::from_fn(p, |i, _| (vu_opt * w_inv.read(i, i)).sqrt());

        let kzt: Mat<f64> = if let Some(ref k_m) = k_mat {
            k_m * z_mat.transpose()
        } else {
            z_mat.transpose().to_owned()
        };
        let kzt_hinv = &kzt * &hinv;
        let ww = &kzt_hinv * kzt.transpose();
        let www = &kzt_hinv * &x_mat;
        let www_winv = &www * &w_inv;
        let www_term = &www_winv * www.transpose();

        let u_se_vec = if k_mat.is_none() {
            DVector::from_fn(m, |i, _| {
                let val = vu_opt * (1.0 - ww.read(i, i) + www_term.read(i, i));
                if val > 0.0 { val.sqrt() } else { 0.0 }
            })
        } else {
            let k_m = k_mat.as_ref().unwrap();
            DVector::from_fn(m, |i, _| {
                let val = vu_opt * (k_m.read(i, i) - ww.read(i, i) + www_term.read(i, i));
                if val > 0.0 { val.sqrt() } else { 0.0 }
            })
        };

        (Some(beta_se_vec), Some(u_se_vec))
    } else {
        (None, None)
    };

    // Convert H_inv if needed
    let hinv_out = if opts.return_hinv {
        Some(DMatrix::from_fn(n, n, |i, j| hinv.read(i, j)))
    } else {
        None
    };

    Ok(MixedSolveResult {
        vu: vu_opt,
        ve: ve_opt,
        beta,
        beta_se,
        u: u_blup,
        u_se,
        ll,
        hinv: hinv_out,
    })
}

/// Golden section search for minimization
fn golden_section_minimize<F>(f: F, mut a: f64, mut b: f64) -> (f64, f64)
where
    F: Fn(f64) -> f64,
{
    let gr = 0.5 * (1.0 + 5f64.sqrt());
    let tol = 1e-8;

    let mut c = b - (b - a) / gr;
    let mut d = a + (b - a) / gr;
    let mut fc = f(c);
    let mut fd = f(d);

    for _ in 0..100 {
        if (b - a).abs() < tol { break; }
        if fc < fd {
            b = d; d = c; fd = fc;
            c = b - (b - a) / gr;
            fc = f(c);
        } else {
            a = c; c = d; fc = fd;
            d = a + (b - a) / gr;
            fd = f(d);
        }
    }

    if fc < fd { (c, fc) } else { (d, fd) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mixed_solve::mixed_solve;
    use approx::assert_relative_eq;

    #[test]
    fn test_mixed_solve_fast_matches_original() {
        let y = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let result_orig = mixed_solve(&y, None, None, None, None).unwrap();
        let result_fast = mixed_solve_fast(&y, None, None, None, None).unwrap();

        assert_relative_eq!(result_orig.vu, result_fast.vu, epsilon = 1e-6);
        assert_relative_eq!(result_orig.ve, result_fast.ve, epsilon = 1e-6);
        assert_relative_eq!(result_orig.beta[0], result_fast.beta[0], epsilon = 1e-6);
        assert_relative_eq!(result_orig.ll, result_fast.ll, epsilon = 1e-6);
    }

    #[test]
    fn test_mixed_solve_fast_with_kinship() {
        let y = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let n = y.len();
        let mut k = DMatrix::identity(n, n);
        for i in 0..n {
            for j in 0..n {
                if i != j { k[(i, j)] = 0.1; }
            }
        }

        let result_orig = mixed_solve(&y, None, Some(&k), None, None).unwrap();
        let result_fast = mixed_solve_fast(&y, None, Some(&k), None, None).unwrap();

        assert_relative_eq!(result_orig.vu, result_fast.vu, epsilon = 1e-5);
        assert_relative_eq!(result_orig.ve, result_fast.ve, epsilon = 1e-5);
    }

    #[test]
    fn test_mixed_solve_fast_with_fixed_effects() {
        let y = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let x = DMatrix::from_row_slice(6, 2, &[
            1.0, 0.0, 1.0, 1.0, 1.0, 2.0,
            1.0, 3.0, 1.0, 4.0, 1.0, 5.0,
        ]);

        let result_orig = mixed_solve(&y, None, None, Some(&x), None).unwrap();
        let result_fast = mixed_solve_fast(&y, None, None, Some(&x), None).unwrap();

        for i in 0..result_orig.beta.len() {
            assert_relative_eq!(result_orig.beta[i], result_fast.beta[i], epsilon = 1e-5);
        }
    }

    #[test]
    fn test_mixed_solve_fast_with_na() {
        let y = vec![1.0, f64::NAN, 3.0, f64::NAN, 5.0];
        let result_orig = mixed_solve(&y, None, None, None, None).unwrap();
        let result_fast = mixed_solve_fast(&y, None, None, None, None).unwrap();
        assert_relative_eq!(result_orig.beta[0], result_fast.beta[0], epsilon = 1e-6);
    }

    #[test]
    fn test_mixed_solve_fast_se() {
        let y = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let opts = MixedSolveOptions { se: true, ..Default::default() };

        let result_fast = mixed_solve_fast(&y, None, None, None, Some(opts)).unwrap();
        assert!(result_fast.beta_se.is_some());
        assert!(result_fast.u_se.is_some());
    }

    #[test]
    fn test_mixed_solve_fast_return_hinv() {
        let y = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let opts = MixedSolveOptions { return_hinv: true, ..Default::default() };

        let result = mixed_solve_fast(&y, None, None, None, Some(opts)).unwrap();
        assert!(result.hinv.is_some());
        let hinv = result.hinv.unwrap();
        assert_eq!(hinv.nrows(), 5);
    }
}
