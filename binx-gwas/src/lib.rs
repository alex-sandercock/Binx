use anyhow::Result;
use binx_core::{
    load_genotypes_biallelic_from_tsv, load_phenotypes_from_tsv, GenotypeMatrixBiallelic,
};
use ndarray::{s, Array1, Array2, Axis};
use ndarray_linalg::Inverse;
use statrs::distribution::StudentsT;
use std::path::Path;

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
/// Phase 2: simple linear model GWAS, no K yet.
pub fn run_gwas<P: AsRef<Path>>(
    geno_path: P,
    pheno_path: P,
    trait_name: &str,
    ploidy: u8,
    model_str: &str,
    out_path: P,
) -> Result<()> {
    let model = parse_model(model_str)?;

    // Load genotype and phenotype data.
    let geno = load_genotypes_biallelic_from_tsv(&geno_path, ploidy)?;
    let pheno = load_phenotypes_from_tsv(&pheno_path)?;

    // Extract trait vector y, aligned to genotype sample order.
    let (y, sample_ids) = align_trait_to_genotypes(&pheno, trait_name, &geno)?;

    // Phase 2: No covariates/PCs yet: X0 = intercept only.
    let n = sample_ids.len();
    let x0 = Array2::<f64>::ones((n, 1)); // intercept

    // Run LM-based GWAS.
    let results = run_lm_gwas(&geno, &y, &x0, model)?;

    // Write TSV output.
    write_results_tsv(out_path, &results)?;

    Ok(())
}

fn parse_model(model_str: &str) -> Result<GeneActionModel> {
    match model_str.to_lowercase().as_str() {
        "additive" => Ok(GeneActionModel::Additive),
        other => Err(anyhow::anyhow!("Unsupported model for Phase 2: {}", other)),
    }
}

/// Align trait vector to genotype sample order.
/// Phase 2: assumes pheno.sample_ids contain all geno.sample_ids.
fn align_trait_to_genotypes(
    pheno: &binx_core::PhenotypeTable,
    trait_name: &str,
    geno: &GenotypeMatrixBiallelic,
) -> Result<(Array1<f64>, Vec<String>)> {
    let y_full = pheno
        .traits
        .get(trait_name)
        .ok_or_else(|| anyhow::anyhow!("Trait '{}' not found in phenotype file", trait_name))?;

    // For now, assume same order; later, properly align by ID.
    if pheno.sample_ids != geno.sample_ids {
        // TODO: implement proper alignment; this is a placeholder.
        return Err(anyhow::anyhow!(
            "Sample IDs in phenotype and genotype files do not match (Phase 2 placeholder)"
        ));
    }

    Ok((y_full.clone(), geno.sample_ids.clone()))
}

/// Simple LM GWAS: y ~ 1 + marker
fn run_lm_gwas(
    geno: &GenotypeMatrixBiallelic,
    y: &Array1<f64>,
    x0: &Array2<f64>,
    model: GeneActionModel,
) -> Result<Vec<MarkerResult>> {
    let n_samples = geno.sample_ids.len();
    let n_markers = geno.marker_ids.len();

    assert_eq!(y.len(), n_samples);
    assert_eq!(x0.nrows(), n_samples);

    let mut results = Vec::with_capacity(n_markers);

    for (marker_idx, marker_id) in geno.marker_ids.iter().enumerate() {
        // Extract dosage vector for this marker (length n_samples).
        let dosage = geno.dosages.index_axis(Axis(0), marker_idx).to_owned();

        // Encode marker according to gene-action model.
        let marker_design = encode_marker(&dosage, geno.ploidy, model)?;

        // Build full design matrix: [X0 | marker_design]
        // Phase 2: X0 is n×1 (intercept), marker_design is n×1 (additive).
        let x_full = hstack(&[x0.view(), marker_design.view()])?;

        // OLS: beta = (X^T X)^(-1) X^T y
        let xtx = x_full.t().dot(&x_full);
        let xtx_inv = xtx.inv()?;
        let xty = x_full.t().dot(y);
        let beta = xtx_inv.dot(&xty);

        // Residuals and variance estimate.
        let y_hat = x_full.dot(&beta);
        let resid = y - &y_hat;
        let df = (n_samples as f64) - (beta.len() as f64);
        let sigma2 = resid.dot(&resid) / df;

        // Var(beta) = sigma2 * (X^T X)^(-1)
        let var_beta = xtx_inv * sigma2;

        // Marker effect is the last coefficient.
        let beta_marker = beta[beta.len() - 1];
        let var_marker = var_beta[(beta.len() - 1, beta.len() - 1)];
        let se_marker = var_marker.sqrt();
        let t_stat = beta_marker / se_marker;
        let dist = StudentsT::new(0.0, 1.0, df)?;
        let p_value = 2.0 * (1.0 - dist.cdf(t_stat.abs()));

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
) -> Result<Array2<f64>> {
    match model {
        GeneActionModel::Additive => {
            // Return n×1 column.
            let n = dosage.len();
            let mut mat = Array2::<f64>::zeros((n, 1));
            for (i, val) in dosage.iter().enumerate() {
                mat[(i, 0)] = *val;
            }
            Ok(mat)
        }
    }
}

/// Horizontal stack of Array2 views.
/// Phase 2 convenience; can be refactored later.
fn hstack<'a>(arrays: &[ndarray::ArrayView2<'a, f64>]) -> Result<Array2<f64>> {
    let n_rows = arrays
        .first()
        .ok_or_else(|| anyhow::anyhow!("No arrays provided to hstack"))?
        .nrows();

    let total_cols: usize = arrays.iter().map(|a| a.ncols()).sum();

    let mut out = Array2::<f64>::zeros((n_rows, total_cols));

    let mut col_offset = 0;
    for arr in arrays {
        if arr.nrows() != n_rows {
            return Err(anyhow::anyhow!(
                "All arrays must have the same number of rows for hstack"
            ));
        }
        let ncols = arr.ncols();
        out
            .slice_mut(s![.., col_offset..col_offset + ncols])
            .assign(arr);
        col_offset += ncols;
    }

    Ok(out)
}

fn write_results_tsv<P: AsRef<std::path::Path>>(
    path: P,
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
