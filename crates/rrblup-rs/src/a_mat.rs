//! Rust implementation of R/rrBLUP::A.mat
//!
//! This module computes the additive relationship matrix from marker data,
//! following the VanRaden (2008) method as implemented in R/rrBLUP.
//!
//! # Example
//!
//! ```
//! use rrblup_rs::a_mat::{a_mat, AMatOptions};
//! use nalgebra::DMatrix;
//!
//! // Genotype matrix: 3 individuals x 4 markers, coded as {-1, 0, 1}
//! let geno = DMatrix::from_row_slice(3, 4, &[
//!     -1.0,  0.0,  1.0,  0.0,
//!      0.0,  1.0, -1.0,  1.0,
//!      1.0, -1.0,  0.0, -1.0,
//! ]);
//!
//! let result = a_mat(&geno, None).unwrap();
//! assert_eq!(result.a.nrows(), 3);
//! assert_eq!(result.a.ncols(), 3);
//! ```
//!
//! # Reference
//!
//! R/rrBLUP package by Jeffrey Endelman:
//! <https://cran.r-project.org/package=rrBLUP>

use anyhow::{anyhow, Result};
use nalgebra::{DMatrix, DVector};

/// Imputation method for missing values
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ImputeMethod {
    /// Impute with column mean (default)
    #[default]
    Mean,
    /// EM algorithm for imputation
    EM,
}

/// Shrinkage method for A matrix estimation
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShrinkMethod {
    /// Endelman-Jannink shrinkage (default when shrink=TRUE)
    EJ,
    /// Regression-based shrinkage
    REG,
}

/// Shrinkage configuration
#[derive(Debug, Clone)]
pub struct ShrinkConfig {
    /// Shrinkage method
    pub method: ShrinkMethod,
    /// Number of iterations (for REG method)
    pub n_iter: usize,
    /// Number of QTL to sample (for REG method)
    pub n_qtl: usize,
}

impl Default for ShrinkConfig {
    fn default() -> Self {
        Self {
            method: ShrinkMethod::EJ,
            n_iter: 100,
            n_qtl: 10,
        }
    }
}

/// Options for A.mat function
#[derive(Debug, Clone)]
pub struct AMatOptions {
    /// Minimum minor allele frequency (default: 1/(2*n))
    pub min_maf: Option<f64>,
    /// Maximum fraction of missing values per marker (default: 1 - 1/(2*n))
    pub max_missing: Option<f64>,
    /// Imputation method for missing values
    pub impute_method: ImputeMethod,
    /// Tolerance for EM convergence (default: 0.02)
    pub tol: f64,
    /// Whether to use shrinkage estimation
    pub shrink: Option<ShrinkConfig>,
    /// Whether to return the imputed genotype matrix
    pub return_imputed: bool,
}

impl Default for AMatOptions {
    fn default() -> Self {
        Self {
            min_maf: None,
            max_missing: None,
            impute_method: ImputeMethod::Mean,
            tol: 0.02,
            shrink: None,
            return_imputed: false,
        }
    }
}

/// Result from A.mat function
#[derive(Debug, Clone)]
pub struct AMatResult {
    /// Additive relationship matrix (n x n)
    pub a: DMatrix<f64>,
    /// Imputed genotype matrix (optional, n x m_filtered)
    pub imputed: Option<DMatrix<f64>>,
    /// Indices of markers that passed filtering
    pub marker_indices: Vec<usize>,
    /// Shrinkage intensity (if shrinkage was applied)
    pub shrink_intensity: Option<f64>,
}

/// Compute additive relationship matrix from marker data
///
/// This is a Rust implementation of R/rrBLUP::A.mat().
///
/// # Arguments
///
/// * `x` - Genotype matrix (n individuals x m markers) with values in {-1, 0, 1}
///         where -1 = homozygous for allele 1, 0 = heterozygous, 1 = homozygous for allele 2.
///         Missing values should be `NaN`.
/// * `options` - Optional settings (min.MAF, max.missing, impute.method, etc.)
///
/// # Returns
///
/// [`AMatResult`] containing the relationship matrix and optionally the imputed genotypes.
///
/// # Errors
///
/// Returns an error if:
/// - Genotype matrix is empty (zero rows or columns)
/// - No markers pass MAF/missing filters
/// - Genetic variance (`var.A`) is zero (all markers monomorphic)
///
/// # Example
///
/// ```
/// use rrblup_rs::a_mat::{a_mat, AMatOptions};
/// use nalgebra::DMatrix;
///
/// let geno = DMatrix::from_row_slice(3, 4, &[
///     -1.0,  0.0,  1.0,  0.0,
///      0.0,  1.0, -1.0,  1.0,
///      1.0, -1.0,  0.0, -1.0,
/// ]);
///
/// let result = a_mat(&geno, None).unwrap();
/// println!("Relationship matrix:\n{}", result.a);
/// ```
///
/// # Notes
///
/// The genotype coding follows R/rrBLUP convention:
/// - X values are {-1, 0, 1} representing genotype dosages
/// - Allele frequency is computed as (X + 1) / 2, so freq = mean(X + 1) / 2
pub fn a_mat(x: &DMatrix<f64>, options: Option<AMatOptions>) -> Result<AMatResult> {
    let opts = options.unwrap_or_default();
    let n = x.nrows(); // number of individuals
    let m_total = x.ncols(); // total number of markers

    if n == 0 || m_total == 0 {
        return Err(anyhow!("Empty genotype matrix"));
    }

    // Compute fraction of missing values per marker
    // R: frac.missing <- apply(X,2,function(x){length(which(is.na(x)))/n})
    let frac_missing: Vec<f64> = (0..m_total)
        .map(|j| {
            let n_missing = (0..n).filter(|&i| !x[(i, j)].is_finite()).count();
            n_missing as f64 / n as f64
        })
        .collect();

    let has_missing = frac_missing.iter().any(|&f| f > 0.0);

    // Compute allele frequencies
    // R: freq <- apply(X + 1, 2, function(x) {mean(x, na.rm = missing)})/2
    let freq: Vec<f64> = (0..m_total)
        .map(|j| {
            let mut sum = 0.0;
            let mut count = 0;
            for i in 0..n {
                let val = x[(i, j)];
                if val.is_finite() {
                    sum += val + 1.0; // X + 1
                    count += 1;
                }
            }
            if count > 0 {
                sum / (2.0 * count as f64) // divide by 2 to get frequency
            } else {
                0.5 // default for all missing
            }
        })
        .collect();

    // Compute MAF
    // R: MAF <- apply(rbind(freq,1-freq),2,min)
    let maf: Vec<f64> = freq.iter().map(|&p| p.min(1.0 - p)).collect();

    // Default thresholds
    // if (is.null(min.MAF)) {min.MAF <- 1/(2*n)}
    // if (is.null(max.missing)) {max.missing <- 1 - 1/(2*n)}
    let min_maf = opts.min_maf.unwrap_or(1.0 / (2.0 * n as f64));
    let max_missing = opts.max_missing.unwrap_or(1.0 - 1.0 / (2.0 * n as f64));

    // Filter markers
    // R: markers <- which((MAF >= min.MAF)&(frac.missing <= max.missing))
    let markers: Vec<usize> = (0..m_total)
        .filter(|&j| maf[j] >= min_maf && frac_missing[j] <= max_missing)
        .collect();

    let m = markers.len();
    if m == 0 {
        return Err(anyhow!(
            "No markers passed filtering (min.MAF={}, max.missing={})",
            min_maf,
            max_missing
        ));
    }

    // R: var.A <- 2 * mean(freq[markers] * (1 - freq[markers]))
    let var_a: f64 = 2.0
        * markers
            .iter()
            .map(|&j| freq[j] * (1.0 - freq[j]))
            .sum::<f64>()
        / m as f64;

    if var_a == 0.0 {
        return Err(anyhow!("var.A is zero; check allele frequencies"));
    }

    // Handle monomorphic markers by setting to expected value
    // R: mono <- which(freq*(1-freq)==0)
    // X[,mono] <- 2*tcrossprod(one,matrix(freq[mono],length(mono),1))-1
    // (This is handled implicitly by using filtered markers)

    // Build centered W matrix
    // R: freq.mat <- tcrossprod(one, matrix(freq[markers], m, 1))
    // R: W <- X[, markers] + 1 - 2 * freq.mat
    let freq_markers: Vec<f64> = markers.iter().map(|&j| freq[j]).collect();

    let mut w = DMatrix::zeros(n, m);
    for (col_idx, &marker_j) in markers.iter().enumerate() {
        let p = freq_markers[col_idx];
        for i in 0..n {
            let val = x[(i, marker_j)];
            if val.is_finite() {
                // W = X + 1 - 2*freq
                w[(i, col_idx)] = val + 1.0 - 2.0 * p;
            } else {
                // Mark as NA for now (will be imputed)
                w[(i, col_idx)] = f64::NAN;
            }
        }
    }

    // Track shrinkage intensity
    let mut shrink_intensity: Option<f64> = None;

    let (a_mat, imputed_w) = if !has_missing {
        // No missing values - straightforward computation
        let (a, delta) = compute_a_no_missing(&w, var_a, m, n, &opts, &freq_markers)?;
        shrink_intensity = delta;
        (a, if opts.return_imputed { Some(w.clone()) } else { None })
    } else {
        // Has missing values - need imputation
        // First, impute with zeros (mean-centered imputation)
        // R: isna <- which(is.na(W))
        // W[isna] <- 0
        for j in 0..m {
            for i in 0..n {
                if !w[(i, j)].is_finite() {
                    w[(i, j)] = 0.0;
                }
            }
        }

        match opts.impute_method {
            ImputeMethod::EM => {
                // EM algorithm for imputation
                if m < n {
                    // Linear dependency - fall back to mean imputation
                    let (a, delta) =
                        compute_a_no_missing(&w, var_a, m, n, &opts, &freq_markers)?;
                    shrink_intensity = delta;
                    (a, if opts.return_imputed { Some(w.clone()) } else { None })
                } else {
                    // Check for linear dependency
                    let cov_mat = compute_row_covariance(&w);
                    let rank = estimate_rank(&cov_mat);

                    if rank < n - 1 {
                        // Linear dependency - fall back to mean imputation
                        let (a, delta) =
                            compute_a_no_missing(&w, var_a, m, n, &opts, &freq_markers)?;
                        shrink_intensity = delta;
                        (a, if opts.return_imputed { Some(w.clone()) } else { None })
                    } else {
                        // Do EM algorithm
                        // Restore NA values for EM
                        for (col_idx, &marker_j) in markers.iter().enumerate() {
                            for i in 0..n {
                                if !x[(i, marker_j)].is_finite() {
                                    w[(i, col_idx)] = f64::NAN;
                                }
                            }
                        }

                        let (a, w_imputed) = em_imputation(&w, var_a, opts.tol)?;
                        (
                            a,
                            if opts.return_imputed {
                                Some(w_imputed)
                            } else {
                                None
                            },
                        )
                    }
                }
            }
            ImputeMethod::Mean => {
                // Mean imputation (already done by setting NaN to 0)
                let (a, delta) = compute_a_no_missing(&w, var_a, m, n, &opts, &freq_markers)?;
                shrink_intensity = delta;
                (a, if opts.return_imputed { Some(w.clone()) } else { None })
            }
        }
    };

    // Convert imputed W back to X scale if requested
    // R: Ximp <- W - 1 + 2*freq.mat
    let imputed_x = imputed_w.map(|w_imp| {
        let mut x_imp = DMatrix::zeros(n, m);
        for j in 0..m {
            let p = freq_markers[j];
            for i in 0..n {
                x_imp[(i, j)] = w_imp[(i, j)] - 1.0 + 2.0 * p;
            }
        }
        x_imp
    });

    Ok(AMatResult {
        a: a_mat,
        imputed: imputed_x,
        marker_indices: markers,
        shrink_intensity,
    })
}

/// Compute A matrix when there are no missing values
fn compute_a_no_missing(
    w: &DMatrix<f64>,
    var_a: f64,
    m: usize,
    n: usize,
    opts: &AMatOptions,
    freq_markers: &[f64],
) -> Result<(DMatrix<f64>, Option<f64>)> {
    match &opts.shrink {
        Some(shrink_config) => match shrink_config.method {
            ShrinkMethod::EJ => {
                // Endelman-Jannink shrinkage
                // R: W.mean <- rowMeans(W)
                // R: cov.W <- cov.W.shrink(W)
                // R: A <- (cov.W+tcrossprod(W.mean))/var.A
                let w_mean = compute_row_means(w);
                let (cov_w, delta) = cov_w_shrink(w);
                let w_mean_outer = &w_mean * w_mean.transpose();
                let a = (&cov_w + &w_mean_outer) / var_a;
                Ok((a, Some(delta)))
            }
            ShrinkMethod::REG => {
                // Regression-based shrinkage
                let mut deltas = Vec::with_capacity(shrink_config.n_iter);
                for _ in 0..shrink_config.n_iter {
                    if let Some(d) = shrink_coeff(w, shrink_config.n_qtl, freq_markers) {
                        deltas.push(d);
                    }
                }
                let delta = if deltas.is_empty() {
                    0.0
                } else {
                    deltas.iter().sum::<f64>() / deltas.len() as f64
                };

                // R: A <- tcrossprod(W)/var.A/m
                // R: A <- (1-delta)*A + delta*mean(diag(A))*diag(n)
                let a_raw = (w * w.transpose()) / (var_a * m as f64);
                let diag_mean: f64 = (0..n).map(|i| a_raw[(i, i)]).sum::<f64>() / n as f64;
                let a = (1.0 - delta) * &a_raw + delta * diag_mean * DMatrix::identity(n, n);
                Ok((a, Some(delta)))
            }
        },
        None => {
            // No shrinkage
            // R: A <- tcrossprod(W)/var.A/m
            let a = (w * w.transpose()) / (var_a * m as f64);
            Ok((a, None))
        }
    }
}

/// Compute row means of a matrix
fn compute_row_means(w: &DMatrix<f64>) -> DVector<f64> {
    let n = w.nrows();
    let m = w.ncols();
    DVector::from_fn(n, |i, _| {
        let sum: f64 = (0..m).map(|j| w[(i, j)]).sum();
        sum / m as f64
    })
}

/// Compute row covariance matrix (covariance between individuals across markers)
fn compute_row_covariance(w: &DMatrix<f64>) -> DMatrix<f64> {
    let n = w.nrows();
    let m = w.ncols();

    // Center rows
    let row_means = compute_row_means(w);
    let mut centered = w.clone();
    for i in 0..n {
        for j in 0..m {
            centered[(i, j)] -= row_means[i];
        }
    }

    // Covariance = (centered @ centered') / (m - 1)
    (&centered * centered.transpose()) / (m - 1) as f64
}

/// Estimate rank of a matrix (count eigenvalues > threshold)
fn estimate_rank(mat: &DMatrix<f64>) -> usize {
    use nalgebra::SymmetricEigen;
    let eig = SymmetricEigen::new(mat.clone());
    let threshold = 1e-10 * eig.eigenvalues.iter().cloned().fold(0.0, f64::max);
    eig.eigenvalues.iter().filter(|&&v| v.abs() > threshold).count()
}

/// Shrinkage estimator for covariance (Endelman-Jannink method)
/// R: cov.W.shrink <- function(W) { ... }
fn cov_w_shrink(w: &DMatrix<f64>) -> (DMatrix<f64>, f64) {
    let m = w.ncols() as f64;
    let n = w.nrows();

    // R: Z <- t(scale(t(W),scale=FALSE)) - center rows
    let row_means = compute_row_means(w);
    let mut z = w.clone();
    for i in 0..n {
        for j in 0..w.ncols() {
            z[(i, j)] -= row_means[i];
        }
    }

    // Z2 <- Z^2
    let z2 = z.map(|v| v * v);

    // R: S <- tcrossprod(Z)/m
    let s = (&z * z.transpose()) / m;

    // R: target <- mean(diag(S))*diag(n)
    let diag_mean: f64 = (0..n).map(|i| s[(i, i)]).sum::<f64>() / n as f64;
    let target = diag_mean * DMatrix::identity(n, n);

    // R: var.S <- tcrossprod(Z2)/m^2-S^2/m
    let z2_outer = (&z2 * z2.transpose()) / (m * m);
    let s_sq = s.map(|v| v * v) / m;
    let var_s = &z2_outer - &s_sq;

    // b2 <- sum(var.S)
    let b2: f64 = var_s.iter().sum();

    // d2 <- sum((S-target)^2)
    let diff = &s - &target;
    let d2: f64 = diff.iter().map(|v| v * v).sum();

    // R: delta <- max(0,min(1,b2/d2))
    let delta = if d2 > 0.0 { (b2 / d2).clamp(0.0, 1.0) } else { 0.0 };

    // return(target*delta + (1-delta)*S)
    let cov_shrink = delta * &target + (1.0 - delta) * &s;
    (cov_shrink, delta)
}

/// Compute shrinkage coefficient (for REG method)
/// R: shrink.coeff <- function(i,W,n.qtl,p) { ... }
fn shrink_coeff(w: &DMatrix<f64>, n_qtl: usize, p: &[f64]) -> Option<f64> {
    use rand::seq::SliceRandom;
    use rand::thread_rng;

    let m = w.ncols();
    let n = w.nrows();

    if n_qtl >= m {
        return None;
    }

    // R: qtl <- sample(1:m, n.qtl)
    let mut rng = thread_rng();
    let mut indices: Vec<usize> = (0..m).collect();
    indices.shuffle(&mut rng);
    let qtl: Vec<usize> = indices[..n_qtl].to_vec();
    let not_qtl: Vec<usize> = indices[n_qtl..].to_vec();

    // R: A.mark <- tcrossprod(W[,-qtl])/sum(2*p[-qtl]*(1-p[-qtl]))
    let denom_mark: f64 = not_qtl.iter().map(|&j| 2.0 * p[j] * (1.0 - p[j])).sum();
    if denom_mark == 0.0 {
        return None;
    }
    let w_mark = select_columns(w, &not_qtl);
    let a_mark = (&w_mark * w_mark.transpose()) / denom_mark;

    // R: A.qtl <- tcrossprod(W[,qtl])/sum(2*p[qtl]*(1-p[qtl]))
    let denom_qtl: f64 = qtl.iter().map(|&j| 2.0 * p[j] * (1.0 - p[j])).sum();
    if denom_qtl == 0.0 {
        return None;
    }
    let w_qtl = select_columns(w, &qtl);
    let a_qtl = (&w_qtl * w_qtl.transpose()) / denom_qtl;

    // R: x <- as.vector(A.mark - mean(diag(A.mark))*diag(n))
    let diag_mean_mark: f64 = (0..n).map(|i| a_mark[(i, i)]).sum::<f64>() / n as f64;
    let x_mat = &a_mark - diag_mean_mark * DMatrix::identity(n, n);
    let x: Vec<f64> = x_mat.iter().cloned().collect();

    // R: y <- as.vector(A.qtl - mean(diag(A.qtl))*diag(n))
    let diag_mean_qtl: f64 = (0..n).map(|i| a_qtl[(i, i)]).sum::<f64>() / n as f64;
    let y_mat = &a_qtl - diag_mean_qtl * DMatrix::identity(n, n);
    let y: Vec<f64> = y_mat.iter().cloned().collect();

    // return(1-cov(y,x)/var(x))
    let cov_xy = covariance(&y, &x);
    let var_x = variance(&x);

    if var_x == 0.0 {
        None
    } else {
        Some(1.0 - cov_xy / var_x)
    }
}

/// Select specific columns from a matrix
fn select_columns(mat: &DMatrix<f64>, cols: &[usize]) -> DMatrix<f64> {
    let n = mat.nrows();
    let m = cols.len();
    DMatrix::from_fn(n, m, |i, j| mat[(i, cols[j])])
}

/// Compute covariance between two vectors
fn covariance(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mean_x: f64 = x.iter().sum::<f64>() / n;
    let mean_y: f64 = y.iter().sum::<f64>() / n;
    let cov: f64 = x
        .iter()
        .zip(y.iter())
        .map(|(&xi, &yi)| (xi - mean_x) * (yi - mean_y))
        .sum();
    cov / (n - 1.0)
}

/// Compute variance of a vector
fn variance(x: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mean: f64 = x.iter().sum::<f64>() / n;
    let var: f64 = x.iter().map(|&xi| (xi - mean).powi(2)).sum();
    var / (n - 1.0)
}

/// EM algorithm for imputation
/// Returns (A matrix, imputed W matrix)
fn em_imputation(w: &DMatrix<f64>, var_a: f64, tol: f64) -> Result<(DMatrix<f64>, DMatrix<f64>)> {
    let n = w.nrows();
    let m = w.ncols();

    // Initialize with mean imputation
    let mut w_imp = w.clone();
    for j in 0..m {
        let mut sum = 0.0;
        let mut count = 0;
        for i in 0..n {
            if w[(i, j)].is_finite() {
                sum += w[(i, j)];
                count += 1;
            }
        }
        let mean = if count > 0 { sum / count as f64 } else { 0.0 };
        for i in 0..n {
            if !w_imp[(i, j)].is_finite() {
                w_imp[(i, j)] = mean;
            }
        }
    }

    // R: mean.vec.new <- matrix(rowMeans(W),n,1)
    let mut mean_vec = compute_row_means(&w_imp);

    // R: cov.mat.new <- cov(t(W))
    let mut cov_mat = compute_row_covariance(&w_imp);

    // R: A.new <- (cov.mat.new + tcrossprod(mean.vec.new))/var.A
    let mean_outer = &mean_vec * mean_vec.transpose();
    let mut a_new = (&cov_mat + &mean_outer) / var_a;

    let max_iter = 100;
    for _ in 0..max_iter {
        let a_old = a_new.clone();
        let cov_mat_old = cov_mat.clone();
        let mean_vec_old = mean_vec.clone();

        // Impute using EM
        let (s, w_new) = impute_em_step(&w, &cov_mat_old, &mean_vec_old)?;

        w_imp = w_new;
        mean_vec = compute_row_means(&w_imp);

        // R: cov.mat.new <- (S-tcrossprod(mean.vec.new)*m)/(m-1)
        let mean_outer_new = &mean_vec * mean_vec.transpose();
        cov_mat = (&s - &mean_outer_new * m as f64) / (m - 1) as f64;

        // R: A.new <- (cov.mat.new + tcrossprod(mean.vec.new))/var.A
        a_new = (&cov_mat + &mean_outer_new) / var_a;

        // R: err <- norm(A.old-A.new,type="F")/n
        let diff = &a_old - &a_new;
        let err = diff.iter().map(|v| v * v).sum::<f64>().sqrt() / n as f64;

        if err < tol {
            break;
        }
    }

    Ok((a_new, w_imp))
}

/// Single step of EM imputation
/// R: impute.EM <- function(W, cov.mat, mean.vec) { ... }
fn impute_em_step(
    w: &DMatrix<f64>,
    cov_mat: &DMatrix<f64>,
    mean_vec: &DVector<f64>,
) -> Result<(DMatrix<f64>, DMatrix<f64>)> {
    let n = w.nrows();
    let m = w.ncols();

    let mut s = DMatrix::zeros(n, n);
    let mut w_imp = w.clone();

    for j in 0..m {
        // R: Wi <- matrix(W[,i],n,1)
        let mut wi: Vec<f64> = (0..n).map(|i| w[(i, j)]).collect();

        // R: missing <- which(is.na(Wi))
        let missing: Vec<usize> = (0..n).filter(|&i| !wi[i].is_finite()).collect();

        let d = if !missing.is_empty() {
            let not_na: Vec<usize> = (0..n).filter(|&i| wi[i].is_finite()).collect();

            if not_na.is_empty() {
                // All missing - use mean
                for &i in &missing {
                    wi[i] = mean_vec[i];
                }
                // D = tcrossprod(Wi)
                let wi_vec = DVector::from_row_slice(&wi);
                &wi_vec * wi_vec.transpose()
            } else {
                // R: Bt <- solve(cov.mat[not.NA,not.NA], cov.mat[not.NA,missing])
                let cov_obs_obs = select_submatrix(cov_mat, &not_na, &not_na);
                let cov_obs_miss = select_submatrix(cov_mat, &not_na, &missing);

                let bt = match cov_obs_obs.clone().try_inverse() {
                    Some(inv) => inv * &cov_obs_miss,
                    None => {
                        // Fall back to mean imputation
                        for &i in &missing {
                            wi[i] = mean_vec[i];
                        }
                        let wi_vec = DVector::from_row_slice(&wi);
                        &wi_vec * wi_vec.transpose()
                    }
                };

                // Wi[missing] <- mean.vec[missing] + crossprod(Bt, Wi[not.NA]-mean.vec[not.NA])
                let wi_obs: DVector<f64> =
                    DVector::from_fn(not_na.len(), |i, _| wi[not_na[i]] - mean_vec[not_na[i]]);
                let imputed = DVector::from_fn(missing.len(), |i, _| mean_vec[missing[i]])
                    + bt.transpose() * &wi_obs;

                for (idx, &i) in missing.iter().enumerate() {
                    wi[i] = imputed[idx];
                }

                // R: C <- cov.mat[missing,missing] - crossprod(cov.mat[not.NA,missing],Bt)
                let cov_miss_miss = select_submatrix(cov_mat, &missing, &missing);
                let c_mat = &cov_miss_miss - cov_obs_miss.transpose() * &bt;

                // R: D <- tcrossprod(Wi)
                // D[missing,missing] <- D[missing,missing] + C
                let wi_vec = DVector::from_row_slice(&wi);
                let mut d = &wi_vec * wi_vec.transpose();
                for (i_idx, &i) in missing.iter().enumerate() {
                    for (j_idx, &j_val) in missing.iter().enumerate() {
                        d[(i, j_val)] += c_mat[(i_idx, j_idx)];
                    }
                }
                d
            }
        } else {
            // No missing values
            let wi_vec = DVector::from_row_slice(&wi);
            &wi_vec * wi_vec.transpose()
        };

        // R: S <- S + D
        s += d;

        // Update W.imp
        for i in 0..n {
            w_imp[(i, j)] = wi[i];
        }
    }

    Ok((s, w_imp))
}

/// Select a submatrix given row and column indices
fn select_submatrix(mat: &DMatrix<f64>, rows: &[usize], cols: &[usize]) -> DMatrix<f64> {
    DMatrix::from_fn(rows.len(), cols.len(), |i, j| mat[(rows[i], cols[j])])
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_a_mat_simple() {
        // Simple 3x4 genotype matrix
        #[rustfmt::skip]
        let x = DMatrix::from_row_slice(3, 4, &[
            -1.0,  0.0,  1.0,  0.0,
             0.0,  1.0, -1.0,  1.0,
             1.0, -1.0,  0.0, -1.0,
        ]);

        let result = a_mat(&x, None).unwrap();

        // A should be 3x3 (n x n)
        assert_eq!(result.a.nrows(), 3);
        assert_eq!(result.a.ncols(), 3);

        // A should be symmetric
        for i in 0..3 {
            for j in 0..3 {
                assert_relative_eq!(result.a[(i, j)], result.a[(j, i)], epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn test_a_mat_with_missing() {
        // Matrix with missing values
        #[rustfmt::skip]
        let x = DMatrix::from_row_slice(3, 4, &[
            -1.0,      0.0,  1.0,  0.0,
             0.0,  f64::NAN, -1.0,  1.0,
             1.0,     -1.0,  0.0, -1.0,
        ]);

        let opts = AMatOptions {
            impute_method: ImputeMethod::Mean,
            ..Default::default()
        };

        let result = a_mat(&x, Some(opts)).unwrap();
        assert_eq!(result.a.nrows(), 3);
        assert_eq!(result.a.ncols(), 3);
    }

    #[test]
    fn test_a_mat_maf_filter() {
        // Matrix where one marker is monomorphic
        #[rustfmt::skip]
        let x = DMatrix::from_row_slice(4, 3, &[
            -1.0, -1.0, 1.0,  // marker 1: all -1 (monomorphic)
            -1.0,  0.0, 0.0,
            -1.0,  1.0, -1.0,
            -1.0, -1.0, 0.0,
        ]);

        let opts = AMatOptions {
            min_maf: Some(0.1), // Should filter out monomorphic marker
            ..Default::default()
        };

        let result = a_mat(&x, Some(opts)).unwrap();

        // First marker should be filtered out
        assert!(!result.marker_indices.contains(&0));
    }

    #[test]
    fn test_a_mat_shrink_ej() {
        #[rustfmt::skip]
        let x = DMatrix::from_row_slice(5, 10, &[
            -1.0,  0.0,  1.0,  0.0, -1.0,  1.0,  0.0, -1.0,  1.0,  0.0,
             0.0,  1.0, -1.0,  1.0,  0.0, -1.0,  1.0,  0.0, -1.0,  1.0,
             1.0, -1.0,  0.0, -1.0,  1.0,  0.0, -1.0,  1.0,  0.0, -1.0,
            -1.0,  1.0,  1.0,  0.0,  0.0,  1.0, -1.0,  0.0,  1.0, -1.0,
             0.0,  0.0, -1.0,  1.0, -1.0, -1.0,  0.0,  1.0, -1.0,  0.0,
        ]);

        let opts = AMatOptions {
            shrink: Some(ShrinkConfig {
                method: ShrinkMethod::EJ,
                ..Default::default()
            }),
            ..Default::default()
        };

        let result = a_mat(&x, Some(opts)).unwrap();

        // Should have shrinkage intensity
        assert!(result.shrink_intensity.is_some());
        let delta = result.shrink_intensity.unwrap();
        assert!(delta >= 0.0 && delta <= 1.0);
    }

    #[test]
    fn test_a_mat_return_imputed() {
        #[rustfmt::skip]
        let x = DMatrix::from_row_slice(3, 4, &[
            -1.0,  0.0,  1.0,  0.0,
             0.0,  1.0, -1.0,  1.0,
             1.0, -1.0,  0.0, -1.0,
        ]);

        let opts = AMatOptions {
            return_imputed: true,
            ..Default::default()
        };

        let result = a_mat(&x, Some(opts)).unwrap();
        assert!(result.imputed.is_some());
    }
}
