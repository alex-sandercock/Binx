use anyhow::{anyhow, Result};
use binx_core::{
    fit_null_mixed_model, load_genotypes_biallelic_from_tsv, load_kinship_from_tsv,
    load_phenotypes_filtered, load_pcs_from_tsv, GenotypeMatrixBiallelic, KinshipMatrix,
    MixedModelCache,
};
use ndarray::{Array1, Array2, Axis};
use nalgebra::{DMatrix, DVector};
use statrs::distribution::{ContinuousCDF, StudentsT};
use std::collections::HashMap;

/// Gene action models (Phase 2: only additive is implemented).
#[derive(Debug, Clone, Copy)]
pub enum GeneActionModel {
    Additive,
    // General, SimplexDomAlt, etc. to be added later.
}

pub struct MarkerResult {
    pub marker_id: String,
    pub beta: f64,
    pub se: f64,
    pub t_stat: f64,
    pub p_value: f64,
}

/// Top-level entry point used by the CLI.
/// Supports LM (default) and will use LMM when --kinship is provided.
pub fn run_gwas(
    geno_path: &str,
    pheno_path: &str,
    trait_name: &str,
    covariate_names: Option<&[String]>,
    pcs_path: Option<&str>,
    kinship_path: Option<&str>,
    allow_missing_samples: bool,
    env_column: Option<&str>,
    env_value: Option<&str>,
    ploidy: u8,
    model_str: &str,
    out_path: &str,
) -> Result<()> {
    let model = parse_model(model_str)?;

    // Load genotype and phenotype data.
    let geno = load_genotypes_biallelic_from_tsv(geno_path, ploidy)?;
    let pheno = load_phenotypes_filtered(pheno_path, env_column, env_value)?;

    // Align samples: choose intersection, optionally allowing missing geno IDs.
    let (geno_keep, pheno_idx) =
        align_samples(&pheno.sample_ids, &geno.sample_ids, allow_missing_samples)?;
    let geno = subset_genotypes(&geno, &geno_keep)?;

    // Extract and align trait vector y.
    let y_full = pheno
        .traits
        .get(trait_name)
        .ok_or_else(|| anyhow!("Trait '{}' not found in phenotype file", trait_name))?;
    let y = align_vector(y_full, &pheno_idx);

    // Build base design: intercept + covariates + PCs (all aligned to genotype order).
    let base_design = build_base_design(&pheno, covariate_names, pcs_path, &geno, &pheno_idx)?;

    if let Some(k_path) = kinship_path {
        let kin = load_kinship(k_path)?;
        let kin_aligned = align_kinship_to_genotypes(kin, &geno.sample_ids)?;
        let x0_array = base_design_to_array2(&base_design)?;
        let cache = fit_null_mixed_model(&y, &x0_array, &kin_aligned)?;
        let results = run_lmm_gwas(&geno, &cache, model)?;
        write_results_tsv(out_path, &results)?;
        return Ok(());
    }

    // Run LM-based GWAS.
    let results = run_lm_gwas(&geno, &y, &base_design, model)?;

    // Write TSV output.
    write_results_tsv(out_path, &results)?;

    Ok(())
}

fn parse_model(model_str: &str) -> Result<GeneActionModel> {
    match model_str.to_lowercase().as_str() {
        "additive" => Ok(GeneActionModel::Additive),
        other => Err(anyhow!("Unsupported model for Phase 2: {}", other)),
    }
}

/// Build alignment index from source IDs (phenotype/PC) to target IDs (genotypes).
fn build_alignment_index(source_ids: &[String], target_ids: &[String]) -> Result<Vec<usize>> {
    let mut index = HashMap::new();
    for (i, sid) in source_ids.iter().enumerate() {
        if index.insert(sid.as_str(), i).is_some() {
            return Err(anyhow!("Duplicate sample_id '{}' in source table", sid));
        }
    }

    let mut aligned = Vec::with_capacity(target_ids.len());
    let mut missing = Vec::new();
    for sid in target_ids {
        match index.get(sid.as_str()) {
            Some(idx) => aligned.push(*idx),
            None => missing.push(sid.clone()),
        }
    }

    if !missing.is_empty() {
        return Err(anyhow!(
            "Missing rows for sample_ids: {}",
            missing.join(", ")
        ));
    }

    Ok(aligned)
}

/// Align samples, returning geno indices to keep (in order) and corresponding pheno indices.
fn align_samples(
    pheno_ids: &[String],
    geno_ids: &[String],
    allow_missing: bool,
) -> Result<(Vec<usize>, Vec<usize>)> {
    let mut pheno_map = HashMap::new();
    for (i, sid) in pheno_ids.iter().enumerate() {
        if pheno_map.insert(sid.as_str(), i).is_some() {
            return Err(anyhow!("Duplicate sample_id '{}' in phenotype file", sid));
        }
    }

    let mut geno_keep = Vec::new();
    let mut pheno_idx = Vec::new();
    let mut missing = Vec::new();
    for (g_idx, sid) in geno_ids.iter().enumerate() {
        if let Some(p_idx) = pheno_map.get(sid.as_str()) {
            geno_keep.push(g_idx);
            pheno_idx.push(*p_idx);
        } else {
            missing.push(sid.clone());
        }
    }

    if geno_keep.is_empty() {
        return Err(anyhow!(
            "No overlapping sample_ids between genotype and phenotype files"
        ));
    }

    if !allow_missing && !missing.is_empty() {
        return Err(anyhow!(
            "Missing rows for sample_ids: {}",
            missing.join(", ")
        ));
    }

    Ok((geno_keep, pheno_idx))
}

fn subset_genotypes(
    geno: &GenotypeMatrixBiallelic,
    keep_indices: &[usize],
) -> Result<GenotypeMatrixBiallelic> {
    let n_markers = geno.marker_ids.len();
    let mut dosages = Array2::<f64>::zeros((n_markers, keep_indices.len()));
    for (row_idx, marker_row) in geno.dosages.outer_iter().enumerate() {
        for (col_out, &col_in) in keep_indices.iter().enumerate() {
            dosages[(row_idx, col_out)] = marker_row[col_in];
        }
    }
    let sample_ids = keep_indices
        .iter()
        .map(|&i| geno.sample_ids[i].clone())
        .collect();
    Ok(GenotypeMatrixBiallelic {
        ploidy: geno.ploidy,
        sample_ids,
        marker_ids: geno.marker_ids.clone(),
        dosages,
    })
}

fn align_vector(values: &Array1<f64>, index: &[usize]) -> Array1<f64> {
    let data: Vec<f64> = index.iter().map(|&i| values[i]).collect();
    Array1::from(data)
}

fn align_vector_to_vec(values: &Array1<f64>, index: &[usize]) -> Vec<f64> {
    index.iter().map(|&i| values[i]).collect()
}

fn align_pcs_to_genotypes(
    pcs: &binx_core::PcMatrix,
    geno_sample_ids: &[String],
) -> Result<Vec<Vec<f64>>> {
    let pcs_to_geno = build_alignment_index(&pcs.sample_ids, geno_sample_ids)?;
    let n = geno_sample_ids.len();
    let n_pcs = pcs.pcs.ncols();
    let mut cols = Vec::with_capacity(n_pcs);
    for pc_idx in 0..n_pcs {
        let mut col = Vec::with_capacity(n);
        for &src_idx in &pcs_to_geno {
            col.push(pcs.pcs[(src_idx, pc_idx)]);
        }
        cols.push(col);
    }
    Ok(cols)
}

fn build_base_design(
    pheno: &binx_core::PhenotypeTable,
    covariate_names: Option<&[String]>,
    pcs_path: Option<&str>,
    geno: &GenotypeMatrixBiallelic,
    pheno_to_geno: &[usize],
) -> Result<Vec<Vec<f64>>> {
    let n_samples = geno.sample_ids.len();
    let mut cols = Vec::new();

    // Intercept
    cols.push(vec![1.0f64; n_samples]);

    // Covariates
    if let Some(names) = covariate_names {
        for name in names {
            let cov = pheno.traits.get(name).ok_or_else(|| {
                anyhow!("Covariate '{}' not found in phenotype file", name)
            })?;
            cols.push(align_vector_to_vec(cov, pheno_to_geno));
        }
    }

    // PCs
    if let Some(pcs_path) = pcs_path {
        let pcs = load_pcs_from_tsv(pcs_path)?;
        let pc_cols = align_pcs_to_genotypes(&pcs, &geno.sample_ids)?;
        cols.extend(pc_cols);
    }

    // Validate column lengths
    for (idx, col) in cols.iter().enumerate() {
        if col.len() != n_samples {
            return Err(anyhow!(
                "Design column {} length {} does not match number of samples {}",
                idx,
                col.len(),
                n_samples
            ));
        }
    }

    Ok(cols)
}

fn base_design_to_array2(cols: &[Vec<f64>]) -> Result<Array2<f64>> {
    if cols.is_empty() {
        return Err(anyhow!("Base design must have at least one column (intercept)"));
    }
    let n_samples = cols[0].len();
    let p = cols.len();
    let mut mat = Array2::<f64>::zeros((n_samples, p));
    for (j, col) in cols.iter().enumerate() {
        if col.len() != n_samples {
            return Err(anyhow!(
                "Design column {} length {} does not match samples {}",
                j,
                col.len(),
                n_samples
            ));
        }
        for (i, val) in col.iter().enumerate() {
            mat[(i, j)] = *val;
        }
    }
    Ok(mat)
}

fn align_kinship_to_genotypes(
    kin: KinshipMatrix,
    geno_ids: &[String],
) -> Result<KinshipMatrix> {
    let map = build_alignment_index(&kin.sample_ids, geno_ids)?;
    let n = geno_ids.len();
    let mut matrix = Array2::<f64>::zeros((n, n));
    for (i_new, &i_old) in map.iter().enumerate() {
        for (j_new, &j_old) in map.iter().enumerate() {
            matrix[(i_new, j_new)] = kin.matrix[(i_old, j_old)];
        }
    }
    Ok(KinshipMatrix {
        sample_ids: geno_ids.to_vec(),
        matrix,
    })
}

/// LM GWAS with arbitrary covariates/PCs: y ~ intercept + covs + marker.
fn run_lm_gwas(
    geno: &GenotypeMatrixBiallelic,
    y: &Array1<f64>,
    base_design_cols: &[Vec<f64>],
    model: GeneActionModel,
) -> Result<Vec<MarkerResult>> {
    let n_samples = geno.sample_ids.len();
    let n_markers = geno.marker_ids.len();

    assert_eq!(y.len(), n_samples);

    let y_slice = y.as_slice().ok_or_else(|| {
        anyhow!("Trait array is not contiguous; cannot obtain slice for LM computation")
    })?;
    let y_vec = DVector::from_row_slice(y_slice);

    // Validate base design columns length
    for (i, col) in base_design_cols.iter().enumerate() {
        if col.len() != n_samples {
            return Err(anyhow!(
                "Base design column {} length {} does not match sample count {}",
                i,
                col.len(),
                n_samples
            ));
        }
    }

    let mut results = Vec::with_capacity(n_markers);

    for (marker_idx, marker_id) in geno.marker_ids.iter().enumerate() {
        // Extract dosage vector for this marker (length n_samples).
        let dosage = geno.dosages.index_axis(Axis(0), marker_idx).to_owned();
        let x_marker = encode_marker(&dosage, geno.ploidy, model)?;

        // Build design matrix X = [base_design_cols | marker].
        let p = base_design_cols.len() + 1;
        let mut data = Vec::with_capacity(n_samples * p);
        for row in 0..n_samples {
            for col in base_design_cols {
                data.push(col[row]);
            }
            data.push(x_marker[row]);
        }
        let x = DMatrix::from_row_slice(n_samples, p, &data);

        // OLS via normal equations solved with LU.
        let xtx = &x.transpose() * &x;
        let xty = &x.transpose() * &y_vec;
        let lu = xtx.lu();
        let beta = lu
            .solve(&xty)
            .ok_or_else(|| anyhow!("Singular design matrix for marker {}", marker_id))?;
        let xtx_inv = lu
            .try_inverse()
            .ok_or_else(|| anyhow!("Failed to invert X^T X for marker {}", marker_id))?;

        // Residuals and variance estimate.
        let y_hat = &x * &beta;
        let resid = &y_vec - y_hat;
        let df = (n_samples as f64) - (p as f64);
        if df <= 0.0 {
            return Err(anyhow!(
                "Degrees of freedom <= 0 (n_samples={}, p={})",
                n_samples,
                p
            ));
        }
        let sigma2 = resid.dot(&resid) / df;

        // Var(beta) = sigma2 * (X^T X)^(-1)
        let se_marker = (sigma2 * xtx_inv[(p - 1, p - 1)]).sqrt();
        let beta_marker = beta[p - 1];
        let (t_stat, p_value) = if se_marker == 0.0 {
            (0.0, 1.0)
        } else {
            let t = beta_marker / se_marker;
            let dist = StudentsT::new(0.0, 1.0, df)?;
            let p = 2.0 * (1.0 - dist.cdf(t.abs()));
            (t, p)
        };

        results.push(MarkerResult {
            marker_id: marker_id.clone(),
            beta: beta_marker,
            se: se_marker,
            t_stat,
            p_value,
        });
    }

    Ok(results)
}

/// LMM GWAS using pre-fit null model cache.
fn run_lmm_gwas(
    geno: &GenotypeMatrixBiallelic,
    cache: &MixedModelCache,
    model: GeneActionModel,
) -> Result<Vec<MarkerResult>> {
    let n_samples = geno.sample_ids.len();
    let n_markers = geno.marker_ids.len();
    let p0 = cache.x0_star.ncols();

    if cache.sample_ids != geno.sample_ids {
        return Err(anyhow!(
            "Kinship sample IDs do not match genotype order after alignment"
        ));
    }

    let mut results = Vec::with_capacity(n_markers);

    for (marker_idx, marker_id) in geno.marker_ids.iter().enumerate() {
        let dosage = geno.dosages.index_axis(Axis(0), marker_idx).to_owned();
        let x_marker = encode_marker(&dosage, geno.ploidy, model)?;

        // Transform marker: x* = D^-1/2 U^T x
        let x_vec = DVector::from_row_slice(&x_marker);
        let x_t = cache.u.transpose() * x_vec;
        let mut x_star_vec = x_t.clone();
        for (xi, s) in x_star_vec.iter_mut().zip(cache.d_inv_sqrt.iter()) {
            *xi *= *s;
        }

        // Build design matrix X* = [X0* | x*]
        let p = p0 + 1;
        let mut data = Vec::with_capacity(n_samples * p);
        for row in 0..n_samples {
            for col in 0..p0 {
                data.push(cache.x0_star[(row, col)]);
            }
            data.push(x_star_vec[row]);
        }
        let x = DMatrix::from_row_slice(n_samples, p, &data);

        let xtx = &x.transpose() * &x;
        let xty = &x.transpose() * &cache.y_star;
        let lu = xtx.lu();
        let beta = match lu.solve(&xty) {
            Some(b) => b,
            None => {
                return Err(anyhow!("Singular design matrix (LMM) for marker {}", marker_id));
            }
        };
        let xtx_inv = match lu.try_inverse() {
            Some(inv) => inv,
            None => {
                return Err(anyhow!("Failed to invert X^T X (LMM) for marker {}", marker_id));
            }
        };

        let y_hat = &x * &beta;
        let resid = &cache.y_star - y_hat;
        let df = (n_samples as f64) - (p as f64);
        if df <= 0.0 {
            return Err(anyhow!(
                "Degrees of freedom <= 0 in LMM (n_samples={}, p={})",
                n_samples,
                p
            ));
        }
        let sigma2 = resid.dot(&resid) / df;

        let se_marker = (sigma2 * xtx_inv[(p - 1, p - 1)]).sqrt();
        let beta_marker = beta[p - 1];
        let (t_stat, p_value) = if se_marker == 0.0 {
            (0.0, 1.0)
        } else {
            let t = beta_marker / se_marker;
            let dist = StudentsT::new(0.0, 1.0, df)?;
            let p = 2.0 * (1.0 - dist.cdf(t.abs()));
            (t, p)
        };

        results.push(MarkerResult {
            marker_id: marker_id.clone(),
            beta: beta_marker,
            se: se_marker,
            t_stat,
            p_value,
        });
    }

    Ok(results)
}

/// Encode marker under a given gene-action model.
/// Phase 2: additive = single column of dosages.
fn encode_marker(
    dosage: &Array1<f64>,
    _ploidy: u8,
    model: GeneActionModel,
) -> Result<Vec<f64>> {
    match model {
        GeneActionModel::Additive => Ok(dosage.to_vec()),
    }
}

fn write_results_tsv(
    path: &str,
    results: &[MarkerResult],
) -> Result<()> {
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(path)?;

    wtr.write_record(&["marker_id", "beta", "se", "t_stat", "p_value"])?;
    for r in results {
        wtr.write_record(&[
            r.marker_id.as_str(),
            &r.beta.to_string(),
            &r.se.to_string(),
            &r.t_stat.to_string(),
            &r.p_value.to_string(),
        ])?;
    }
    wtr.flush()?;
    Ok(())
}

fn load_kinship(path: &str) -> Result<KinshipMatrix> {
    let kin = load_kinship_from_tsv(path)?;
    Ok(kin)
}

#[cfg(test)]
mod tests {
    use super::*;
    use binx_core::{KinshipMatrix, PhenotypeTable};
    use ndarray::array;

    #[test]
    fn aligns_trait_to_genotypes_by_sample_id() {
        let pheno = PhenotypeTable {
            sample_ids: vec!["S1".into(), "S2".into(), "S3".into()],
            traits: vec![("height".into(), array![1.0, 2.0, 3.0])].into_iter().collect(),
            covariates: Default::default(),
        };
        let geno = GenotypeMatrixBiallelic {
            ploidy: 2,
            sample_ids: vec!["S3".into(), "S1".into(), "S2".into()],
            marker_ids: vec!["m1".into()],
            dosages: array![[0.0, 1.0, 2.0]],
        };

        let idx = build_alignment_index(&pheno.sample_ids, &geno.sample_ids).unwrap();
        let aligned = align_vector(pheno.traits.get("height").unwrap(), &idx);
        assert_eq!(aligned.to_vec(), vec![3.0, 1.0, 2.0]);
    }

    #[test]
    fn lm_produces_expected_beta_and_significance() {
        let geno = GenotypeMatrixBiallelic {
            ploidy: 2,
            sample_ids: vec!["S1".into(), "S2".into(), "S3".into()],
            marker_ids: vec!["m1".into()],
            dosages: array![[0.0, 1.0, 2.0]],
        };
        let y = array![1.0, 2.9, 6.1];
        let base = vec![vec![1.0; 3]]; // intercept only

        let results = run_lm_gwas(&geno, &y, &base, GeneActionModel::Additive).unwrap();
        let r = &results[0];
        assert!((r.beta - 2.55).abs() < 1e-6);
        assert!(r.p_value < 0.2);
    }

    #[test]
    fn covariate_absorbs_confounded_signal() {
        // Trait depends only on covariate; marker is correlated but should go to ~0 when covariate is included.
        let geno = GenotypeMatrixBiallelic {
            ploidy: 2,
            sample_ids: vec!["S1".into(), "S2".into(), "S3".into(), "S4".into()],
            marker_ids: vec!["m1".into()],
            dosages: array![[0.0, 1.0, 1.0, 2.0]],
        };
        let y = array![0.0, 0.0, 10.0, 10.0];
        let cov = vec![0.0, 0.0, 1.0, 1.0];
        let base = vec![vec![1.0; 4], cov];

        let results = run_lm_gwas(&geno, &y, &base, GeneActionModel::Additive).unwrap();
        let r = &results[0];
        assert!(r.beta.abs() < 1e-8);
        assert!(r.p_value > 0.9);
    }

    #[test]
    fn lmm_matches_lm_when_kinship_is_identity() {
        let geno = GenotypeMatrixBiallelic {
            ploidy: 2,
            sample_ids: vec!["S1".into(), "S2".into(), "S3".into()],
            marker_ids: vec!["m1".into()],
            dosages: array![[0.0, 1.0, 2.0]],
        };
        let y = array![1.0, 2.9, 6.1];
        let base_cols = vec![vec![1.0; 3]]; // intercept only
        let x0 = base_design_to_array2(&base_cols).unwrap();
        let kin = KinshipMatrix {
            sample_ids: geno.sample_ids.clone(),
            matrix: array![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        };
        let cache = fit_null_mixed_model(&y, &x0, &kin).unwrap();

        let lm = run_lm_gwas(&geno, &y, &base_cols, GeneActionModel::Additive).unwrap();
        let lmm = run_lmm_gwas(&geno, &cache, GeneActionModel::Additive).unwrap();
        let lm_res = &lm[0];
        let lmm_res = &lmm[0];
        assert!((lm_res.beta - lmm_res.beta).abs() < 1e-6);
        assert!((lm_res.se - lmm_res.se).abs() < 1e-6);
    }
}
