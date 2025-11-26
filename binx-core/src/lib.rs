use anyhow::{anyhow, Result};
use ndarray::{Array1, Array2};
use nalgebra::{DMatrix, DVector, SymmetricEigen};
use serde::Deserialize;
use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

pub type SampleId = String;
pub type MarkerId = String;

/// Simple phenotype table: samples × (traits + covariates).
pub struct PhenotypeTable {
    pub sample_ids: Vec<SampleId>,
    pub traits: HashMap<String, Array1<f64>>, // trait_name -> y
    pub covariates: HashMap<String, Array1<f64>>, // cov_name -> x
}

/// Biallelic dosage matrix: markers × samples, entries 0..ploidy.
pub struct GenotypeMatrixBiallelic {
    pub ploidy: u8,
    pub sample_ids: Vec<SampleId>,
    pub marker_ids: Vec<MarkerId>,
    /// shape: (n_markers, n_samples)
    pub dosages: Array2<f64>,
}

/// Simple dense kinship matrix K (samples × samples).
pub struct KinshipMatrix {
    pub sample_ids: Vec<SampleId>,
    /// shape: (n_samples, n_samples)
    pub matrix: Array2<f64>,
}

/// Placeholder cache for mixed-model computations (Phase 3).
/// Will hold transformed y/X and precomputed matrices for fast per-marker GLS.
pub struct MixedModelCache {
    pub sample_ids: Vec<SampleId>,
    pub u: DMatrix<f64>,
    pub d_inv_sqrt: DVector<f64>,
    pub y_star: DVector<f64>,
    pub x0_star: DMatrix<f64>,
    pub xtx_inv: DMatrix<f64>,
    pub sigma2: f64,
    pub lambda: f64,
}

/// Phase 3 stub: fit null mixed model y ~ X0 with kinship K.
/// Implements a simple grid-search REML over lambda = sigma_g^2 / sigma_e^2.
pub fn fit_null_mixed_model(
    y: &Array1<f64>,
    x0: &Array2<f64>,
    kinship: &KinshipMatrix,
) -> Result<MixedModelCache> {
    let n = y.len();
    if x0.nrows() != n {
        return Err(anyhow!(
            "X0 rows ({}) do not match y length ({})",
            x0.nrows(),
            n
        ));
    }
    if kinship.matrix.nrows() != n || kinship.matrix.ncols() != n {
        return Err(anyhow!(
            "Kinship matrix must be square with size equal to sample count ({})",
            n
        ));
    }
    let p = x0.ncols();
    if p == 0 {
        return Err(anyhow!("X0 must have at least one column (intercept)"));
    }
    if n <= p {
        return Err(anyhow!(
            "Not enough samples (n={}) for fixed effects (p={})",
            n,
            p
        ));
    }

    let y_slice = y
        .as_slice()
        .ok_or_else(|| anyhow!("Trait array is not contiguous"))?;
    let x_slice = x0
        .as_slice()
        .ok_or_else(|| anyhow!("Design matrix is not contiguous"))?;
    let k_slice = kinship
        .matrix
        .as_slice()
        .ok_or_else(|| anyhow!("Kinship matrix is not contiguous"))?;

    let y_vec = DVector::from_row_slice(y_slice);
    let x0_mat = DMatrix::from_row_slice(n, p, x_slice);
    let k_mat = DMatrix::from_row_slice(n, n, k_slice);

    let eig = SymmetricEigen::new(k_mat);
    let vals = eig.eigenvalues;
    let u = eig.eigenvectors;
    let y_t = u.transpose() * &y_vec;
    let x_t = u.transpose() * &x0_mat;

    // Grid search over log10(lambda) in [-4, 4]
    let mut best_ll = f64::NEG_INFINITY;
    let mut best_lambda = 0.0;
    let mut best_sigma2 = None;

    for step in 0..=32 {
        let log10_lambda = -4.0 + (step as f64) * 0.25;
        let lambda = 10f64.powf(log10_lambda);
        if let Some((ll, sigma2)) = reml_loglik(lambda, &vals, &y_t, &x_t) {
            if ll.is_finite() && ll > best_ll {
                best_ll = ll;
                best_lambda = lambda;
                best_sigma2 = Some(sigma2);
            }
        }
    }

    let sigma2 = best_sigma2
        .ok_or_else(|| anyhow!("Failed to optimize REML for mixed model"))?;

    let (_delta, d_inv_sqrt) = build_delta_inv_sqrt(best_lambda, &vals);

    // y* = D^-1/2 U^T y ; X0* = D^-1/2 U^T X0
    let mut y_star = y_t.clone();
    elementwise_scale_vec(&mut y_star, &d_inv_sqrt);

    let mut x0_star = x_t.clone();
    elementwise_scale_rows(&mut x0_star, &d_inv_sqrt);

    let xtx = &x0_star.transpose() * &x0_star;
    let xtx_inv = xtx
        .try_inverse()
        .ok_or_else(|| anyhow!("Failed to invert X0*^T X0*"))?;

    Ok(MixedModelCache {
        sample_ids: kinship.sample_ids.clone(),
        u,
        d_inv_sqrt,
        y_star,
        x0_star,
        xtx_inv,
        sigma2,
        lambda: best_lambda,
    })
}

fn build_delta_inv_sqrt(lambda: f64, vals: &DVector<f64>) -> (Vec<f64>, DVector<f64>) {
    let mut delta = Vec::with_capacity(vals.len());
    let mut inv_sqrt = Vec::with_capacity(vals.len());
    for &d in vals.iter() {
        let val = lambda * d + 1.0;
        delta.push(val);
        inv_sqrt.push(1.0 / val.sqrt());
    }
    let inv_sqrt_vec = DVector::from_row_slice(&inv_sqrt);
    (delta, inv_sqrt_vec)
}

fn reml_loglik(
    lambda: f64,
    vals: &DVector<f64>,
    y_t: &DVector<f64>,
    x_t: &DMatrix<f64>,
) -> Option<(f64, f64)> {
    let n = y_t.len();
    let p = x_t.ncols();
    let mut delta = Vec::with_capacity(n);
    let mut inv = Vec::with_capacity(n);
    for &d in vals.iter() {
        let v = lambda * d + 1.0;
        delta.push(v);
        inv.push(1.0 / v);
    }

    let (xt_dinv_x, xt_dinv_y, y_dinv_y) = accumulate_weighted_moments(x_t, y_t, &inv);
    let lu = xt_dinv_x.lu();
    if !lu.is_invertible() {
        return None;
    }
    let beta = lu.solve(&xt_dinv_y)?;
    let resid = y_dinv_y - xt_dinv_y.dot(&beta);
    if resid <= 0.0 {
        return None;
    }
    let sigma2 = resid / ((n - p) as f64);
    if !sigma2.is_finite() || sigma2 <= 0.0 {
        return None;
    }
    let logdet_xt = lu.determinant().abs().ln();
    let logdet_delta: f64 = delta.iter().map(|v| v.ln()).sum();
    let ll = -0.5 * (logdet_delta + logdet_xt + (n - p) as f64 * (sigma2.ln() + 1.0));
    Some((ll, sigma2))
}

fn accumulate_weighted_moments(
    x_t: &DMatrix<f64>,
    y_t: &DVector<f64>,
    inv: &[f64],
) -> (DMatrix<f64>, DVector<f64>, f64) {
    let n = y_t.len();
    let p = x_t.ncols();
    let mut xt_dinv_x = DMatrix::<f64>::zeros(p, p);
    let mut xt_dinv_y = DVector::<f64>::zeros(p);
    let mut y_dinv_y = 0.0;
    for i in 0..n {
        let w = inv[i];
        let yv = y_t[i];
        y_dinv_y += w * yv * yv;
        for a in 0..p {
            let xa = x_t[(i, a)];
            xt_dinv_y[a] += w * xa * yv;
            for b in 0..p {
                xt_dinv_x[(a, b)] += w * xa * x_t[(i, b)];
            }
        }
    }
    (xt_dinv_x, xt_dinv_y, y_dinv_y)
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

/// Optional PCs: samples × PCs.
pub struct PcMatrix {
    pub sample_ids: Vec<SampleId>,
    /// shape: (n_samples, n_pcs)
    pub pcs: Array2<f64>,
}

/// Row type for phenotype TSV (Phase 1: flexible).
#[allow(dead_code)]
#[derive(Debug, Deserialize)]
struct PhenoRow {
    #[serde(rename = "sample_id")]
    sample_id: String,
    // Remaining columns will be traits/covariates.
    // We can parse them as String->String then cast to f64 as needed.
}

/// Phase 1: load phenotypes from a TSV.
/// Assumes first column is sample_id, others numeric.
pub fn load_phenotypes_from_tsv<P: AsRef<Path>>(path: P) -> Result<PhenotypeTable> {
    let file = File::open(path)?;
    let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').from_reader(file);

    let headers = rdr.headers()?.clone();
    let col_names: Vec<String> = headers.iter().map(|s| s.to_string()).collect();

    // First column is sample_id; the rest are numeric.
    let value_cols: Vec<String> = col_names.iter().skip(1).cloned().collect();

    let mut sample_ids = Vec::new();
    let mut columns: HashMap<String, Vec<f64>> = value_cols
        .iter()
        .map(|name| (name.clone(), Vec::new()))
        .collect();

    for result in rdr.records() {
        let record = result?;
        let sample_id = record.get(0).unwrap().to_string();
        sample_ids.push(sample_id);

        for (i, col_name) in value_cols.iter().enumerate() {
            let val_str = record.get(i + 1).unwrap();
            let val: f64 = val_str.parse().map_err(|e| {
                anyhow!(
                    "Failed to parse value '{}' in column '{}': {}",
                    val_str,
                    col_name,
                    e
                )
            })?;
            columns.get_mut(col_name).unwrap().push(val);
        }
    }

    let mut traits = HashMap::new();
    let covariates = HashMap::new(); // fill later if you want to distinguish

    for (name, vals) in columns {
        let arr = Array1::from(vals);
        traits.insert(name, arr);
    }

    Ok(PhenotypeTable {
        sample_ids,
        traits,
        covariates,
    })
}

/// Phase 1: load biallelic dosage matrix from TSV.
/// Expecting something like:
/// marker_id  chr  pos  S1  S2  S3 ...
pub fn load_genotypes_biallelic_from_tsv<P: AsRef<Path>>(
    path: P,
    ploidy: u8,
) -> Result<GenotypeMatrixBiallelic> {
    let file = File::open(path)?;
    let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').from_reader(file);

    let headers = rdr.headers()?.clone();
    let header_strings: Vec<String> = headers.iter().map(|s| s.to_string()).collect();

    // Assume first 3 columns are marker metadata: id, chr, pos.
    // Sample IDs start at index 3.
    let sample_ids: Vec<SampleId> = header_strings.iter().skip(3).cloned().collect();

    let mut marker_ids = Vec::new();
    let mut dosage_rows: Vec<Vec<f64>> = Vec::new();

    for result in rdr.records() {
        let record = result?;
        let marker_id = record.get(0).unwrap().to_string();
        marker_ids.push(marker_id);

        let mut row = Vec::with_capacity(sample_ids.len());
        for i in 3..record.len() {
            let val_str = record.get(i).unwrap();
            let val: f64 = val_str.parse().map_err(|e| {
                anyhow!("Failed to parse dosage '{}' at col {}: {}", val_str, i, e)
            })?;
            row.push(val);
        }
        dosage_rows.push(row);
    }

    let n_markers = marker_ids.len();
    let n_samples = sample_ids.len();
    let mut dosages = Array2::<f64>::zeros((n_markers, n_samples));

    for (i, row) in dosage_rows.into_iter().enumerate() {
        for (j, val) in row.into_iter().enumerate() {
            dosages[[i, j]] = val;
        }
    }

    Ok(GenotypeMatrixBiallelic {
        ploidy,
        sample_ids,
        marker_ids,
        dosages,
    })
}

/// Phase 1: minimal kinship loader.
/// Expecting:
/// sample_id  S1  S2  S3 ...
pub fn load_kinship_from_tsv<P: AsRef<Path>>(path: P) -> Result<KinshipMatrix> {
    let file = File::open(path)?;
    let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').from_reader(file);

    let headers = rdr.headers()?.clone();
    let header_strings: Vec<String> = headers.iter().map(|s| s.to_string()).collect();
    let sample_ids: Vec<SampleId> = header_strings.iter().skip(1).cloned().collect();

    let mut rows: Vec<Vec<f64>> = Vec::new();

    for result in rdr.records() {
        let record = result?;
        // first column is sample_id; we assume order matches header
        let mut row_vals = Vec::with_capacity(sample_ids.len());
        for i in 1..record.len() {
            let val_str = record.get(i).unwrap();
            let val: f64 = val_str.parse().map_err(|e| {
                anyhow!("Failed to parse kinship value '{}': {}", val_str, e)
            })?;
            row_vals.push(val);
        }
        rows.push(row_vals);
    }

    let n = sample_ids.len();
    let mut mat = Array2::<f64>::zeros((n, n));
    for (i, row) in rows.into_iter().enumerate() {
        for (j, val) in row.into_iter().enumerate() {
            mat[[i, j]] = val;
        }
    }

    Ok(KinshipMatrix {
        sample_ids,
        matrix: mat,
    })
}


/// Phase 2: load PCs from TSV.
/// Expecting rows of: sample_id  PC1  PC2 ...
pub fn load_pcs_from_tsv<P: AsRef<Path>>(path: P) -> Result<PcMatrix> {
    let file = File::open(path)?;
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
    use ndarray::array;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_temp(contents: &str) -> NamedTempFile {
        let mut file = NamedTempFile::new().expect("create temp file");
        write!(file, "{}", contents).expect("write temp file");
        file
    }

    #[test]
    fn load_phenotypes_parses_numeric_columns() {
        let contents = "sample_id\ttrait1\tcov1\nS1\t1.0\t5.0\nS2\t2.0\t6.0\n";
        let file = write_temp(contents);

        let pheno = load_phenotypes_from_tsv(file.path()).expect("load phenotypes");
        assert_eq!(pheno.sample_ids, vec!["S1", "S2"]);
        let trait_vals = pheno.traits.get("trait1").unwrap();
        assert_eq!(trait_vals, &array![1.0, 2.0]);
    }

    #[test]
    fn load_genotypes_parses_dosages() {
        let contents = "marker_id\tchr\tpos\tS1\tS2\nm1\t1\t10\t0\t1\nm2\t1\t20\t2\t1\n";
        let file = write_temp(contents);

        let geno = load_genotypes_biallelic_from_tsv(file.path(), 2).expect("load genotypes");
        assert_eq!(geno.sample_ids, vec!["S1", "S2"]);
        assert_eq!(geno.marker_ids, vec!["m1", "m2"]);
        assert_eq!(geno.dosages.shape(), &[2, 2]);
        assert_eq!(geno.dosages[[0, 1]], 1.0);
        assert_eq!(geno.dosages[[1, 0]], 2.0);
    }

    #[test]
    fn load_kinship_parses_matrix() {
        let contents = "sample_id\tS1\tS2\nS1\t1.0\t0.2\nS2\t0.2\t1.0\n";
        let file = write_temp(contents);

        let kin = load_kinship_from_tsv(file.path()).expect("load kinship");
        assert_eq!(kin.sample_ids, vec!["S1", "S2"]);
        assert_eq!(kin.matrix[[0, 1]], 0.2);
        assert_eq!(kin.matrix[[1, 0]], 0.2);
    }

    #[test]
    fn load_pcs_parses_matrix() {
        let contents = "sample_id\tPC1\tPC2\nS1\t0.1\t0.2\nS2\t0.3\t0.4\n";
        let file = write_temp(contents);

        let pcs = load_pcs_from_tsv(file.path()).expect("load pcs");
        assert_eq!(pcs.sample_ids, vec!["S1", "S2"]);
        assert_eq!(pcs.pcs.shape(), &[2, 2]);
        assert_eq!(pcs.pcs[[0, 1]], 0.2);
    }
}
