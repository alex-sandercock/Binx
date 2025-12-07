//! binx-core: shared data structures and I/O for the binx toolkit.

use anyhow::{anyhow, Result};
use csv;
use nalgebra::{DMatrix, DVector, SymmetricEigen};
use ndarray::{Array1, Array2};
use std::collections::HashMap;
use std::path::Path;

// faer for full SVD and QR with complete Q matrix
use faer::prelude::*;
use faer::Mat as FaerMat;

pub type SampleId = String;
pub type MarkerId = String;

/// Optional marker metadata for genomic coordinates.
#[derive(Clone, Debug)]
pub struct MarkerMetadata {
    pub chrom: String,
    pub pos: f64,
}

/// Principal components matrix.
pub struct PcMatrix {
    pub sample_ids: Vec<SampleId>,
    /// shape: (n_samples, n_pcs)
    pub pcs: Array2<f64>,
}

/// Simple phenotype table: samples × (traits + covariates).
pub struct PhenotypeTable {
    pub sample_ids: Vec<SampleId>,
    pub traits: HashMap<String, Array1<f64>>,       // numeric trait_name -> y
    pub covariates: HashMap<String, Array1<f64>>,   // numeric cov_name -> x
    pub factor_covariates: HashMap<String, Vec<String>>, // factor cov_name -> levels per observation
}

/// Biallelic dosage matrix: markers × samples, entries 0..ploidy.
#[derive(Clone, Debug)]
pub struct GenotypeMatrixBiallelic {
    pub ploidy: u8,
    pub sample_ids: Vec<SampleId>,
    pub marker_ids: Vec<MarkerId>,
    pub marker_metadata: Option<Vec<MarkerMetadata>>,
    pub dosages: Array2<f64>,
}

/// Kinship matrix with sample labels.
#[derive(Clone, Debug)]
pub struct KinshipMatrix {
    pub sample_ids: Vec<SampleId>,
    pub matrix: Array2<f64>,
}

/// Cached quantities from fitting the null mixed model.
/// Used for efficient marker testing via P3D (population parameters previously determined).
pub struct MixedModelCache {
    pub sample_ids: Vec<SampleId>,
    pub n_obs: usize,
    pub n_ind: usize,
    pub u: DMatrix<f64>,
    pub d_inv_sqrt: DVector<f64>,
    pub h_inv: DMatrix<f64>,
    pub y_star: DVector<f64>,
    pub x0_star: DMatrix<f64>,
    pub xtx_inv: DMatrix<f64>,
    pub sigma2: f64,
    pub lambda: f64,
    pub vu: f64,
    pub ve: f64,
}

/// Fit null mixed model using rrBLUP::mixed.solve algorithm.
///
/// This is a Rust implementation of R/rrBLUP's mixed.solve function.
/// Uses the cholesky branch when n > m + p, otherwise the eigen branch.
/// Uses faer for full SVD (Mat::svd gives full U matrix).
pub fn fit_null_mixed_model(
    y: &Array1<f64>,
    x0: &Array2<f64>,
    kinship: &KinshipMatrix,
    z: Option<&Array2<f64>>,
    obs_ids: Option<&[SampleId]>,
) -> Result<MixedModelCache> {
    let n = y.len();  // n_obs
    let p = x0.ncols();

    if x0.nrows() != n {
        return Err(anyhow!("X0 rows ({}) do not match y length ({})", x0.nrows(), n));
    }
    if p == 0 {
        return Err(anyhow!("X0 must have at least one column (intercept)"));
    }
    if n <= p {
        return Err(anyhow!("Not enough samples (n={}) for fixed effects (p={})", n, p));
    }

    let y_vec = DVector::from_row_slice(
        y.as_slice().ok_or_else(|| anyhow!("Trait array is not contiguous"))?,
    );
    let x_mat = DMatrix::from_row_slice(
        n, p, x0.as_slice().ok_or_else(|| anyhow!("Design matrix is not contiguous"))?,
    );

    // Build Z matrix
    let m = kinship.matrix.nrows();  // n_ind
    if kinship.matrix.ncols() != m {
        return Err(anyhow!("Kinship matrix must be square"));
    }

    let z_mat = if let Some(z) = z {
        if z.ncols() != m {
            return Err(anyhow!("Z columns ({}) do not match kinship size ({})", z.ncols(), m));
        }
        if z.nrows() != n {
            return Err(anyhow!("Z rows ({}) do not match observation count ({})", z.nrows(), n));
        }
        DMatrix::from_row_slice(z.nrows(), z.ncols(), z.as_slice().unwrap())
    } else {
        if m != n {
            return Err(anyhow!("Kinship size ({}) does not match sample count ({})", m, n));
        }
        DMatrix::identity(n, m)
    };

    // S = I - X(X'X)^-1 X' (projection onto space orthogonal to X)
    let xtx = x_mat.transpose() * &x_mat;
    let xtx_inv = xtx.clone().try_inverse()
        .ok_or_else(|| anyhow!("X not full rank"))?;
    let s_mat = DMatrix::identity(n, n) - &x_mat * &xtx_inv * x_mat.transpose();

    let k_mat = DMatrix::from_row_slice(m, m, kinship.matrix.as_slice().unwrap());

    // Determine spectral method (R's mixed.solve logic)
    let use_cholesky = n > m + p;

    let (phi, theta, u, omega_sq): (Vec<f64>, Vec<f64>, DMatrix<f64>, Vec<f64>);
    let offset: f64;

    if use_cholesky {
        // ============================================================
        // CHOLESKY BRANCH (n > m + p): R's spectral.method == "cholesky"
        // Uses faer for full SVD (equivalent to R's svd(A, nu=n))
        // ============================================================
        offset = 0.0;

        // B = chol(K + 1e-6*I) using nalgebra
        let mut k_jittered = k_mat.clone();
        for i in 0..m {
            k_jittered[(i, i)] += 1e-6;
        }
        let chol = k_jittered.clone().cholesky()
            .ok_or_else(|| anyhow!("K not positive semi-definite for Cholesky"))?;
        let b_mat = chol.l();  // Lower triangular

        // ZBt = Z @ B' (tcrossprod(Z, B) in R)
        let zbt_nalg = &z_mat * b_mat.transpose();

        // Convert to faer matrix for full SVD
        let zbt_faer = nalgebra_to_faer(&zbt_nalg);

        // Full SVD using faer: svd() gives full U (n x n), equivalent to R's svd(A, nu=n)
        let svd_zbt = zbt_faer.svd();
        let u_full_faer = svd_zbt.u();  // Full U matrix (n x n)
        let s_vals = svd_zbt.s_diagonal();

        // Convert U back to nalgebra
        u = faer_to_nalgebra(&u_full_faer);

        // phi = c(d^2, rep(0, n-m)) - singular values squared, padded
        phi = (0..n).map(|i| {
            if i < s_vals.nrows() {
                let s = s_vals.read(i);  // singular values are a column vector
                s * s
            } else {
                0.0
            }
        }).collect();

        // SZBt = S @ ZBt
        let szbt_nalg = &s_mat * &zbt_nalg;
        let szbt_faer = nalgebra_to_faer(&szbt_nalg);

        // SVD of SZBt (thin SVD is sufficient here)
        let svd_szbt = szbt_faer.thin_svd();
        let u_szbt_faer = svd_szbt.u();
        let d_szbt = svd_szbt.s_diagonal();

        // Convert to nalgebra for QR
        let u_szbt = faer_to_nalgebra(&u_szbt_faer);

        // QR decomposition of [X | U_szbt] using faer for complete Q
        let mut combined = DMatrix::zeros(n, p + u_szbt.ncols());
        for i in 0..n {
            for j in 0..p {
                combined[(i, j)] = x_mat[(i, j)];
            }
            for j in 0..u_szbt.ncols() {
                combined[(i, p + j)] = u_szbt[(i, j)];
            }
        }

        // Use faer QR for complete Q matrix
        let combined_faer = nalgebra_to_faer(&combined);
        let qr = combined_faer.qr();

        // Get complete Q matrix (n x n) using faer's compute_q()
        let q_full_faer = qr.compute_q();
        let q_full = faer_to_nalgebra(&q_full_faer.as_ref());

        // Q = Q_full[, (p+1):n] - columns (p+1) to n (0-indexed: p to n-1)
        let q_complement = q_full.columns(p, n - p).into_owned();

        // Get R matrix from QR
        let r_faer = qr.compute_r();
        let r_mat = faer_mat_to_nalgebra(&r_faer);

        // R22 = R[(p+1):(p+m), (p+1):(p+m)] in R (1-indexed)
        // In 0-indexed: R[p:p+m, p:p+m]
        let r_size = m.min(r_mat.nrows().saturating_sub(p)).min(r_mat.ncols().saturating_sub(p));

        // theta = c(solve(t(R^2), d^2), rep(0, n-m-p))
        // R's code: ans <- solve(t(R^2), svd.SZBt$d^2)
        // This solves: t(R22^2) @ x = d^2
        if r_size > 0 && d_szbt.nrows() > 0 {
            // Extract R22 and compute R22^2 (element-wise square)
            let mut r22_sq = DMatrix::zeros(r_size, r_size);
            for i in 0..r_size {
                for j in 0..r_size {
                    let val = r_mat[(p + i, p + j)];
                    r22_sq[(i, j)] = val * val;
                }
            }

            // t(R22^2) - transpose of element-wise squared R22
            let t_r22_sq = r22_sq.transpose();

            // d^2 from SVD
            let d_sq: Vec<f64> = (0..r_size.min(d_szbt.nrows()))
                .map(|i| {
                    let s = d_szbt.read(i);  // singular values are a column vector
                    s * s
                })
                .collect();
            let d_sq_vec = DVector::from_row_slice(&d_sq);

            // solve(t(R^2), d^2) = inv(t(R22^2)) @ d^2
            eprintln!("DEBUG theta: r_size={}, d_sq.len={}, r22_sq[0,0]={:.4}, t_r22_sq[0,0]={:.4}",
                r_size, d_sq.len(), r22_sq[(0,0)], t_r22_sq[(0,0)]);
            let theta_sub = t_r22_sq.clone().try_inverse()
                .map(|inv| {
                    eprintln!("DEBUG theta: inv succeeded, inv[0,0]={:.6}", inv[(0,0)]);
                    inv * &d_sq_vec
                })
                .unwrap_or_else(|| {
                    eprintln!("DEBUG theta: inv failed, using regularized fallback");
                    // Regularized fallback
                    let mut reg = t_r22_sq.clone();
                    for i in 0..reg.nrows().min(reg.ncols()) {
                        reg[(i, i)] += 1e-10;
                    }
                    reg.try_inverse().map(|inv| inv * &d_sq_vec).unwrap_or(d_sq_vec.clone())
                });

            // Pad with zeros: c(theta_sub, rep(0, n-m-p))
            let n_theta = n - p;
            theta = (0..n_theta).map(|i| {
                if i < theta_sub.len() {
                    theta_sub[i].max(0.0)
                } else {
                    0.0
                }
            }).collect();
        } else {
            theta = vec![0.0; n - p];
        }

        // omega = crossprod(Q, y) = Q' @ y
        let omega_vec = q_complement.transpose() * &y_vec;
        omega_sq = omega_vec.iter().map(|v| v * v).collect();

    } else {
        // ============================================================
        // EIGEN BRANCH (n <= m + p): R's spectral.method == "eigen"
        // ============================================================
        offset = (n as f64).sqrt();

        // Hb = Z @ K @ Z' + offset * I
        let zk = &z_mat * &k_mat;
        let mut hb = &zk * z_mat.transpose();
        for i in 0..n {
            hb[(i, i)] += offset;
        }

        let hb_eig = SymmetricEigen::new(hb.clone());
        phi = hb_eig.eigenvalues.iter().map(|v| *v - offset).collect();

        if phi.iter().cloned().fold(f64::INFINITY, f64::min) < -1e-6 {
            return Err(anyhow!("K not positive semi-definite (phi < 0)"));
        }
        u = hb_eig.eigenvectors.clone();

        // SHbS = S @ Hb @ S
        let shbs = &s_mat * &hb * &s_mat;
        let shbs_eig = SymmetricEigen::new(shbs);

        let n_theta = n - p;
        theta = shbs_eig.eigenvalues.iter()
            .take(n_theta)
            .map(|v| *v - offset)
            .collect();

        let q = shbs_eig.eigenvectors.columns(0, n_theta).into_owned();
        let omega_vec = q.transpose() * &y_vec;
        omega_sq = omega_vec.iter().map(|v| v * v).collect();
    }

    // REML optimization (same for both branches)
    let df = n - p;

    let reml_obj = |lambda: f64| -> Option<f64> {
        if lambda <= 0.0 {
            return None;
        }
        let denom: f64 = omega_sq.iter()
            .zip(theta.iter())
            .map(|(o, t)| o / (t + lambda))
            .sum();
        if denom <= 0.0 {
            return None;
        }
        let term: f64 = theta.iter().map(|t| (t + lambda).ln()).sum();
        Some((df as f64) * denom.ln() + term)
    };

    // Golden section search for optimal lambda
    let (mut a, mut b) = (1e-9f64, 1e9f64);
    let gr = 0.5 * (1.0 + 5f64.sqrt());
    let mut c = b - (b - a) / gr;
    let mut d = a + (b - a) / gr;
    let mut fc = reml_obj(c).unwrap_or(f64::INFINITY);
    let mut fd = reml_obj(d).unwrap_or(f64::INFINITY);

    for _ in 0..120 {
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - (b - a) / gr;
            fc = reml_obj(c).unwrap_or(f64::INFINITY);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + (b - a) / gr;
            fd = reml_obj(d).unwrap_or(f64::INFINITY);
        }
    }

    let lambda_opt = if fc < fd { c } else { d };
    let vu_opt = omega_sq.iter()
        .zip(theta.iter())
        .map(|(o, t)| o / (t + lambda_opt))
        .sum::<f64>() / (df as f64);
    let ve_opt = lambda_opt * vu_opt;

    // Compute H^-1 = U @ diag(1/(phi + lambda)) @ U'
    let d_inv_vals: Vec<f64> = phi.iter().map(|v| 1.0 / (v + lambda_opt)).collect();
    let d_inv_sqrt = DVector::from_iterator(
        phi.len(),
        phi.iter().map(|v| 1.0 / (v + lambda_opt).sqrt()),
    );
    let d_inv = DMatrix::from_diagonal(&DVector::from_row_slice(&d_inv_vals));
    let h_inv = &u * d_inv * u.transpose();

    // Compute transformed quantities for efficient marker testing
    let y_t = u.transpose() * &y_vec;
    let mut y_star = y_t.clone();
    elementwise_scale_vec(&mut y_star, &d_inv_sqrt);

    let x_t = u.transpose() * &x_mat;
    let mut x0_star = x_t.clone();
    elementwise_scale_rows(&mut x0_star, &d_inv_sqrt);

    let xtx_star = x0_star.transpose() * &x0_star;
    let xtx_star_inv = xtx_star.try_inverse()
        .ok_or_else(|| anyhow!("Failed to invert X0*'X0*"))?;

    let cache_sample_ids = if let Some(ids) = obs_ids {
        ids.to_vec()
    } else {
        kinship.sample_ids.clone()
    };

    // Debug output
    let theta_sum: f64 = theta.iter().sum();
    let theta_nonzero = theta.iter().filter(|&&t| t > 1e-10).count();
    let omega_sum: f64 = omega_sq.iter().sum();
    eprintln!("DEBUG REML spectral: theta.len={}, theta_nonzero={}, theta_sum={:.4}, omega_sq_sum={:.4}",
        theta.len(), theta_nonzero, theta_sum, omega_sum);
    eprintln!("DEBUG REML spectral: theta[0..5]={:?}", &theta[0..5.min(theta.len())]);
    eprintln!("DEBUG REML spectral: omega_sq[0..5]={:?}", &omega_sq[0..5.min(omega_sq.len())]);
    eprintln!("DEBUG REML: n={}, m={}, p={}, branch={}, lambda={:.4}, vu={:.6}, ve={:.6}",
        n, m, p, if use_cholesky { "cholesky" } else { "eigen" },
        lambda_opt, vu_opt, ve_opt);

    Ok(MixedModelCache {
        sample_ids: cache_sample_ids,
        n_obs: n,
        n_ind: m,
        u,
        d_inv_sqrt,
        h_inv,
        y_star,
        x0_star,
        xtx_inv: xtx_star_inv,
        sigma2: ve_opt,
        lambda: lambda_opt,
        vu: vu_opt,
        ve: ve_opt,
    })
}

/// Convert nalgebra DMatrix to faer Mat
fn nalgebra_to_faer(m: &DMatrix<f64>) -> FaerMat<f64> {
    let nrows = m.nrows();
    let ncols = m.ncols();
    FaerMat::from_fn(nrows, ncols, |i, j| m[(i, j)])
}

/// Convert faer MatRef to nalgebra DMatrix
fn faer_to_nalgebra(m: &faer::MatRef<f64>) -> DMatrix<f64> {
    let nrows = m.nrows();
    let ncols = m.ncols();
    DMatrix::from_fn(nrows, ncols, |i, j| m.read(i, j))
}

/// Convert faer Mat to nalgebra DMatrix
fn faer_mat_to_nalgebra(m: &FaerMat<f64>) -> DMatrix<f64> {
    let nrows = m.nrows();
    let ncols = m.ncols();
    DMatrix::from_fn(nrows, ncols, |i, j| m.read(i, j))
}

fn elementwise_scale_vec(v: &mut DVector<f64>, scale: &DVector<f64>) {
    for (val, s) in v.iter_mut().zip(scale.iter()) {
        *val *= *s;
    }
}

fn elementwise_scale_rows(m: &mut DMatrix<f64>, scale: &DVector<f64>) {
    for i in 0..m.nrows() {
        let s = scale[i];
        for j in 0..m.ncols() {
            m[(i, j)] *= s;
        }
    }
}

/// Detect delimiter (comma, tab, space) in a file.
pub fn detect_delimiter<P: AsRef<Path>>(path: P) -> Result<u8> {
    let mut rdr = std::io::BufReader::new(std::fs::File::open(&path)?);
    let mut first_line = String::new();
    std::io::BufRead::read_line(&mut rdr, &mut first_line)?;
    if first_line.contains('\t') {
        Ok(b'\t')
    } else if first_line.contains(',') {
        Ok(b',')
    } else {
        Ok(b' ')
    }
}

/// Load phenotype table from TSV/CSV.
pub fn load_phenotypes_from_tsv<P: AsRef<Path>>(path: P) -> Result<PhenotypeTable> {
    let delim = detect_delimiter(&path)?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(delim)
        .from_path(&path)?;

    let headers = rdr.headers()?.clone();
    if headers.len() < 2 {
        return Err(anyhow!(
            "Phenotype file needs at least 2 columns: sample_id, <trait>"
        ));
    }

    let mut sample_ids = Vec::new();
    let mut raw_columns: Vec<(String, Vec<String>)> = headers
        .iter()
        .skip(1)
        .map(|h| (h.to_string(), Vec::new()))
        .collect();

    for result in rdr.records() {
        let record = result?;
        if record.is_empty() {
            continue;
        }
        sample_ids.push(record.get(0).unwrap_or("").to_string());
        for (i, (_, vals)) in raw_columns.iter_mut().enumerate() {
            vals.push(record.get(i + 1).unwrap_or("").to_string());
        }
    }

    let mut traits = HashMap::new();
    let mut covariates = HashMap::new();
    let mut factor_covariates = HashMap::new();

    for (name, vals) in raw_columns.into_iter() {
        let mut numeric_vals = Vec::with_capacity(vals.len());
        let mut all_numeric = true;
        for v in &vals {
            match v.parse::<f64>() {
                Ok(val) => numeric_vals.push(val),
                Err(_) => {
                    all_numeric = false;
                    break;
                }
            }
        }
        if all_numeric {
            let arr = Array1::from(numeric_vals);
            traits.insert(name.clone(), arr.clone());
            covariates.insert(name, arr);
        } else {
            factor_covariates.insert(name, vals);
        }
    }

    Ok(PhenotypeTable {
        sample_ids,
        traits,
        covariates,
        factor_covariates,
    })
}

/// Load phenotype table with optional row filtering by column value.
pub fn load_phenotypes_filtered<P: AsRef<Path>>(
    path: P,
    filter_column: Option<&str>,
    filter_value: Option<&str>,
) -> Result<PhenotypeTable> {
    load_phenotypes_internal(path, filter_column.zip(filter_value))
}

fn load_phenotypes_internal<P: AsRef<Path>>(
    path: P,
    filter: Option<(&str, &str)>,
) -> Result<PhenotypeTable> {
    let delim = detect_delimiter(&path)?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(delim)
        .from_path(&path)?;

    let headers = rdr.headers()?.clone();
    let col_names: Vec<String> = headers.iter().map(|s| s.to_string()).collect();

    let filter_idx = if let Some((col, _)) = filter {
        headers
            .iter()
            .position(|h| h == col)
            .ok_or_else(|| anyhow!("Filter column '{}' not found in phenotype file", col))?
    } else {
        usize::MAX
    };

    // First column is sample_id/id; the rest are numeric.
    let value_cols: Vec<String> = col_names.iter().skip(1).cloned().collect();

    let mut sample_ids = Vec::new();
    let mut raw_columns: HashMap<String, Vec<String>> = value_cols
        .iter()
        .map(|name| (name.clone(), Vec::new()))
        .collect();

    for result in rdr.records() {
        let record = result?;
        if let Some((_, filter_val)) = filter {
            if record
                .get(filter_idx)
                .map(|v| v != filter_val)
                .unwrap_or(false)
            {
                continue;
            }
        }
        let sample_id = record.get(0).unwrap().to_string();
        sample_ids.push(sample_id);

        for (i, col_name) in value_cols.iter().enumerate() {
            let val_str = record.get(i + 1).unwrap();
            raw_columns
                .get_mut(col_name)
                .unwrap()
                .push(val_str.to_string());
        }
    }

    let mut traits = HashMap::new();
    let mut covariates = HashMap::new();
    let mut factor_covariates = HashMap::new();

    for (name, vals) in raw_columns.into_iter() {
        let mut numeric_vals = Vec::with_capacity(vals.len());
        let mut all_numeric = true;
        for v in &vals {
            match v.parse::<f64>() {
                Ok(val) => numeric_vals.push(val),
                Err(_) => {
                    all_numeric = false;
                    break;
                }
            }
        }
        if all_numeric {
            let arr = Array1::from(numeric_vals);
            traits.insert(name.clone(), arr.clone());
            covariates.insert(name, arr);
        } else {
            factor_covariates.insert(name, vals);
        }
    }

    Ok(PhenotypeTable {
        sample_ids,
        traits,
        covariates,
        factor_covariates,
    })
}

/// Load biallelic dosage matrix from TSV/CSV.
/// Expecting: marker_id, chr, pos, <samples...>
pub fn load_genotypes_biallelic_from_tsv<P: AsRef<Path>>(
    path: P,
    ploidy: u8,
) -> Result<GenotypeMatrixBiallelic> {
    let delim = detect_delimiter(&path)?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(delim)
        .from_path(&path)?;

    let headers = rdr.headers()?.clone();
    let header_strings: Vec<String> = headers.iter().map(|s| s.to_string()).collect();

    if headers.len() < 4 {
        return Err(anyhow!(
            "Genotype file needs at least 4 columns: marker_id, chr, pos, <samples...>"
        ));
    }

    // Assume first 3 columns are marker metadata: id, chr, pos.
    // Sample IDs start at index 3.
    let sample_ids: Vec<String> = header_strings[3..].to_vec();

    let mut marker_ids = Vec::new();
    let mut marker_metadata = Vec::new();
    let mut dosage_rows: Vec<Vec<f64>> = Vec::new();

    for result in rdr.records() {
        let record = result?;
        if record.len() < 4 {
            continue;
        }
        let marker_id = record.get(0).unwrap_or("").to_string();
        let chrom = record.get(1).unwrap_or("").to_string();
        let pos_str = record.get(2).unwrap_or("0");
        let pos: f64 = pos_str.parse().unwrap_or(0.0);

        marker_ids.push(marker_id);
        marker_metadata.push(MarkerMetadata { chrom, pos });

        let row: Vec<f64> = (3..record.len())
            .map(|i| record.get(i).unwrap_or("").parse::<f64>().unwrap_or(f64::NAN))
            .collect();
        dosage_rows.push(row);
    }

    let n_markers = marker_ids.len();
    let n_samples = sample_ids.len();
    if n_markers == 0 {
        return Err(anyhow!("No markers found in genotype file"));
    }

    let mut dosages = Array2::<f64>::zeros((n_markers, n_samples));
    for (i, row) in dosage_rows.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            if j < n_samples {
                dosages[(i, j)] = val;
            }
        }
    }

    Ok(GenotypeMatrixBiallelic {
        ploidy,
        sample_ids,
        marker_ids,
        marker_metadata: Some(marker_metadata),
        dosages,
    })
}

/// Load kinship matrix from TSV.
/// First column is sample_id, remaining columns are the matrix values.
pub fn load_kinship_from_tsv<P: AsRef<Path>>(path: P) -> Result<KinshipMatrix> {
    let delim = detect_delimiter(&path)?;
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(delim)
        .from_path(&path)?;

    let headers = rdr.headers()?.clone();
    if headers.len() < 2 {
        return Err(anyhow!("Kinship file needs at least 2 columns"));
    }

    let col_ids: Vec<String> = headers.iter().skip(1).map(|s| s.to_string()).collect();
    let n = col_ids.len();

    let mut row_ids = Vec::new();
    let mut matrix_rows: Vec<Vec<f64>> = Vec::new();

    for result in rdr.records() {
        let record = result?;
        if record.is_empty() {
            continue;
        }
        row_ids.push(record.get(0).unwrap_or("").to_string());
        let row: Vec<f64> = (1..record.len())
            .map(|i| record.get(i).unwrap_or("0").parse::<f64>().unwrap_or(0.0))
            .collect();
        matrix_rows.push(row);
    }

    if row_ids.len() != n {
        return Err(anyhow!(
            "Kinship matrix not square: {} rows, {} cols",
            row_ids.len(),
            n
        ));
    }

    let mut matrix = Array2::<f64>::zeros((n, n));
    for (i, row) in matrix_rows.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            if j < n {
                matrix[(i, j)] = val;
            }
        }
    }

    Ok(KinshipMatrix {
        sample_ids: row_ids,
        matrix,
    })
}

/// Load principal components from TSV file.
/// First column is sample_id, remaining columns are PC values.
pub fn load_pcs_from_tsv<P: AsRef<Path>>(path: P) -> Result<PcMatrix> {
    let file = std::fs::File::open(&path)?;
    let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').from_reader(file);

    let headers = rdr.headers()?.clone();
    if headers.len() < 2 {
        return Err(anyhow!(
            "PC file must have at least one PC column besides sample_id"
        ));
    }
    let n_pcs = headers.len() - 1;

    let mut sample_ids = Vec::new();
    let mut rows: Vec<Vec<f64>> = Vec::new();

    for result in rdr.records() {
        let record = result?;
        sample_ids.push(
            record
                .get(0)
                .ok_or_else(|| anyhow!("Missing sample_id column in PC file"))?
                .to_string(),
        );

        let mut row = Vec::with_capacity(n_pcs);
        for i in 1..record.len() {
            let val_str = record.get(i).unwrap();
            let val: f64 = val_str.parse().map_err(|e| {
                anyhow!("Failed to parse PC value '{}' at col {}: {}", val_str, i, e)
            })?;
            row.push(val);
        }
        rows.push(row);
    }

    let n_samples = sample_ids.len();
    let mut pcs = Array2::<f64>::zeros((n_samples, n_pcs));
    for (i, row) in rows.into_iter().enumerate() {
        for (j, val) in row.into_iter().enumerate() {
            pcs[[i, j]] = val;
        }
    }

    Ok(PcMatrix { sample_ids, pcs })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_delimiter() {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "a,b,c").unwrap();
        assert_eq!(detect_delimiter(f.path()).unwrap(), b',');
    }
}
