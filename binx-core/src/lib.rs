use anyhow::{anyhow, Result};
use ndarray::{Array1, Array2};
use nalgebra::{DMatrix, DVector, SymmetricEigen};
use serde::Deserialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::path::Path;

pub type SampleId = String;
pub type MarkerId = String;

/// Simple phenotype table: samples × (traits + covariates).
pub struct PhenotypeTable {
    pub sample_ids: Vec<SampleId>,
    pub traits: HashMap<String, Array1<f64>>,       // numeric trait_name -> y
    pub covariates: HashMap<String, Array1<f64>>,   // numeric cov_name -> x
    pub factor_covariates: HashMap<String, Vec<String>>, // factor cov_name -> levels per observation
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

/// Cache for mixed-model computations.
pub struct MixedModelCache {
    pub sample_ids: Vec<SampleId>, // matches y length (obs-level if Z provided)
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
}

/// Fit null mixed model y ~ X0 with kinship K, rrBLUP-style spectral REML.
/// If `z` is provided, it is the incidence matrix (n_obs x n_ind). Otherwise, identity is assumed (n_obs must equal n_ind).
pub fn fit_null_mixed_model(
    y: &Array1<f64>,
    x0: &Array2<f64>,
    kinship: &KinshipMatrix,
    z: Option<&Array2<f64>>,
    obs_ids: Option<&[SampleId]>,
) -> Result<MixedModelCache> {
    let n_obs = y.len();
    let p = x0.ncols();
    if x0.nrows() != n_obs {
        return Err(anyhow!(
            "X0 rows ({}) do not match y length ({})",
            x0.nrows(),
            n_obs
        ));
    }
    if p == 0 {
        return Err(anyhow!("X0 must have at least one column (intercept)"));
    }
    if n_obs <= p {
        return Err(anyhow!(
            "Not enough samples (n_obs={}) for fixed effects (p={})",
            n_obs,
            p
        ));
    }

    let y_vec = DVector::from_row_slice(
        y.as_slice()
            .ok_or_else(|| anyhow!("Trait array is not contiguous"))?,
    );
    let x0_mat = DMatrix::from_row_slice(
        n_obs,
        p,
        x0.as_slice()
            .ok_or_else(|| anyhow!("Design matrix is not contiguous"))?,
    );

    // Build Z.
    let n_ind = kinship.matrix.nrows();
    if kinship.matrix.ncols() != n_ind {
        return Err(anyhow!("Kinship matrix must be square"));
    }
    let z_mat = if let Some(z) = z {
        if z.ncols() != n_ind {
            return Err(anyhow!(
                "Z columns ({}) do not match kinship size ({})",
                z.ncols(),
                n_ind
            ));
        }
        if z.nrows() != n_obs {
            return Err(anyhow!(
                "Z rows ({}) do not match observation count ({})",
                z.nrows(),
                n_obs
            ));
        }
        DMatrix::from_row_slice(z.nrows(), z.ncols(), z.as_slice().unwrap())
    } else {
        if n_ind != n_obs {
            return Err(anyhow!(
                "Kinship size ({}) does not match sample count ({})",
                n_ind,
                n_obs
            ));
        }
        DMatrix::identity(n_obs, n_ind)
    };

    // Hb = Z K Z^T + offset * I (offset stabilizes eigen).
    let offset = (n_obs as f64).sqrt();
    let k_mat = DMatrix::from_row_slice(n_ind, n_ind, kinship.matrix.as_slice().unwrap());
    let hb = {
        let zk = &z_mat * &k_mat;
        let mut h = &zk * z_mat.transpose();
        for i in 0..n_obs {
            h[(i, i)] += offset;
        }
        h
    };
    // Projection matrix S = I - X (X'X)^{-1} X'
    let xtx = &x0_mat.transpose() * &x0_mat;
    let xtx_inv = xtx
        .try_inverse()
        .ok_or_else(|| anyhow!("X0 not full rank"))?;
    let proj = {
        let i = DMatrix::<f64>::identity(n_obs, n_obs);
        i - &x0_mat * (&xtx_inv * x0_mat.transpose())
    };

    // SHbS eigen to get theta/Q (projected space).
    let shbs = &proj * &hb * proj.transpose();
    let shbs_eig = SymmetricEigen::new(shbs);

    let hb_eig = SymmetricEigen::new(hb);
    let eigenvectors = hb_eig.eigenvectors.clone();
    let phi: Vec<f64> = hb_eig
        .eigenvalues
        .iter()
        .map(|v| *v - offset)
        .collect();
    if phi.iter().cloned().fold(f64::INFINITY, f64::min) < -1e-6 {
        return Err(anyhow!("K not positive semi-definite (phi < 0)"));
    }
    let u = eigenvectors.clone();

    // SHbS eigen to get theta/Q (projected space).
    let n_theta = n_obs.saturating_sub(p);
    let theta: Vec<f64> = shbs_eig
        .eigenvalues
        .iter()
        .take(n_theta)
        .map(|v| *v - offset)
        .collect();
    let q = shbs_eig.eigenvectors.columns(0, n_theta).into_owned();
    let omega = q.transpose() * &y_vec;
    let omega_sq: Vec<f64> = omega.iter().map(|v| v * v).collect();

    // REML objective and golden-section search on [1e-9, 1e9].
    let reml_obj = |lambda: f64| -> Option<f64> {
        if lambda <= 0.0 {
            return None;
        }
        let denom: f64 = omega_sq
            .iter()
            .zip(theta.iter())
            .map(|(o, t)| o / (t + lambda))
            .sum();
        if denom <= 0.0 {
            return None;
        }
        let term = theta
            .iter()
            .map(|t| (t + lambda).ln())
            .sum::<f64>();
        Some(((n_obs - p) as f64) * denom.ln() + term)
    };
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
    let vu = omega_sq
        .iter()
        .zip(theta.iter())
        .map(|(o, t)| o / (t + lambda_opt))
        .sum::<f64>()
        / ((n_obs - p) as f64);

    // Transforms using phi (length n_obs).
    let d_inv_sqrt = DVector::from_iterator(
        phi.len(),
        phi.iter().map(|v| 1.0 / (v + lambda_opt).sqrt()),
    );
    // Hinv = U diag(1/(phi+lambda)) U^T
    let h_inv = {
        let mut diag_vals = Vec::with_capacity(phi.len());
        for v in phi.iter() {
            diag_vals.push(1.0 / (v + lambda_opt));
        }
        let d_inv = DMatrix::from_diagonal(&DVector::from_row_slice(&diag_vals));
        &u * d_inv * u.transpose()
    };
    let y_t = u.transpose() * &y_vec;
    let mut y_star = y_t.clone();
    elementwise_scale_vec(&mut y_star, &d_inv_sqrt);

    let x_t = u.transpose() * &x0_mat;
    let mut x0_star = x_t.clone();
    elementwise_scale_rows(&mut x0_star, &d_inv_sqrt);

    let xtx_star = &x0_star.transpose() * &x0_star;
    let xtx_inv = xtx_star
        .try_inverse()
        .ok_or_else(|| anyhow!("Failed to invert X0*^T X0*"))?;

    let cache_sample_ids = if let Some(ids) = obs_ids {
        ids.to_vec()
    } else {
        kinship.sample_ids.clone()
    };

    Ok(MixedModelCache {
        sample_ids: cache_sample_ids,
        n_obs,
        n_ind,
        u,
        d_inv_sqrt,
        h_inv,
        y_star,
        x0_star,
        xtx_inv,
        sigma2: lambda_opt * vu,
        lambda: lambda_opt,
    })
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
    load_phenotypes_internal(path, None)
}

/// Phase 1: load phenotypes with optional row filter on a column/value (e.g., env).
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

/// Phase 1: load biallelic dosage matrix from TSV.
/// Expecting something like:
/// marker_id  chr  pos  S1  S2  S3 ...
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
    let delim = detect_delimiter(&path)?;
    let mut rdr = csv::ReaderBuilder::new().delimiter(delim).from_path(&path)?;

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

/// Detect delimiter as tab or comma from the first non-empty line.
fn detect_delimiter<P: AsRef<Path>>(path: P) -> Result<u8> {
    let file = File::open(path.as_ref())?;
    let mut reader = std::io::BufReader::new(file);
    let mut first_line = String::new();
    reader.read_line(&mut first_line)?;
    let tabs = first_line.matches('\t').count();
    let commas = first_line.matches(',').count();
    if commas > tabs {
        Ok(b',')
    } else {
        Ok(b'\t')
    }
}
