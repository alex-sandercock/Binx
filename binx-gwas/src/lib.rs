use anyhow::Result;
use binx_core::{
    load_genotypes_biallelic_from_tsv, load_phenotypes_from_tsv, GenotypeMatrixBiallelic,
};
use ndarray::{Array1, Axis};
use statrs::distribution::{ContinuousCDF, StudentsT};
use std::collections::HashMap;
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
    let (y, _) = align_trait_to_genotypes(&pheno, trait_name, &geno)?;

    // Run LM-based GWAS.
    let results = run_lm_gwas(&geno, &y, model)?;

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

    if pheno.sample_ids.len() != y_full.len() {
        return Err(anyhow::anyhow!(
            "Phenotype sample_id count ({}) does not match trait length ({})",
            pheno.sample_ids.len(),
            y_full.len()
        ));
    }

    // Build lookup from sample_id -> phenotype index (fail on duplicates).
    let mut pheno_index: HashMap<&str, usize> = HashMap::new();
    for (idx, sid) in pheno.sample_ids.iter().enumerate() {
        if pheno_index.insert(sid.as_str(), idx).is_some() {
            return Err(anyhow::anyhow!(
                "Duplicate sample_id '{}' in phenotype table",
                sid
            ));
        }
    }

    let mut aligned = Vec::with_capacity(geno.sample_ids.len());
    let mut missing = Vec::new();
    for sid in &geno.sample_ids {
        match pheno_index.get(sid.as_str()) {
            Some(idx) => aligned.push(y_full[*idx]),
            None => missing.push(sid.clone()),
        }
    }

    if !missing.is_empty() {
        return Err(anyhow::anyhow!(
            "Missing phenotype rows for sample_ids: {}",
            missing.join(", ")
        ));
    }

    Ok((Array1::from(aligned), geno.sample_ids.clone()))
}

/// Simple LM GWAS: y ~ 1 + marker
fn run_lm_gwas(
    geno: &GenotypeMatrixBiallelic,
    y: &Array1<f64>,
    model: GeneActionModel,
) -> Result<Vec<MarkerResult>> {
    let n_samples = geno.sample_ids.len();
    let n_markers = geno.marker_ids.len();

    assert_eq!(y.len(), n_samples);

    let y_slice = y.as_slice().ok_or_else(|| {
        anyhow::anyhow!("Trait array is not contiguous; cannot obtain slice for LM computation")
    })?;

    let mut results = Vec::with_capacity(n_markers);

    for (marker_idx, marker_id) in geno.marker_ids.iter().enumerate() {
        // Extract dosage vector for this marker (length n_samples).
        let dosage = geno.dosages.index_axis(Axis(0), marker_idx).to_owned();
        let x = dosage.as_slice().ok_or_else(|| {
            anyhow::anyhow!(
                "Dosage array is not contiguous; cannot obtain slice for marker {}",
                marker_id
            )
        })?;

        // Encode marker according to gene-action model.
        let encoded = encode_marker(x, geno.ploidy, model)?;

        // Simple closed-form OLS for y ~ 1 + x (additive model).
        let n = n_samples as f64;
        let sum_x: f64 = encoded.iter().sum();
        let sum_y: f64 = y_slice.iter().sum();
        let sum_xx: f64 = encoded.iter().map(|v| v * v).sum();
        let sum_xy: f64 = encoded
            .iter()
            .zip(y_slice.iter())
            .map(|(a, b)| a * b)
            .sum();

        let mean_x = sum_x / n;
        let mean_y = sum_y / n;
        let sxx = sum_xx - n * mean_x * mean_x;
        let sxy = sum_xy - n * mean_x * mean_y;

        if sxx == 0.0 {
            return Err(anyhow::anyhow!(
                "Variance of marker {} is zero; cannot fit regression",
                marker_id
            ));
        }

        let beta_marker = sxy / sxx;
        let beta0 = mean_y - beta_marker * mean_x;

        // Residual sum of squares.
        let mut sse = 0.0;
        for (xi, yi) in encoded.iter().zip(y_slice.iter()) {
            let y_hat = beta0 + beta_marker * xi;
            let resid = yi - y_hat;
            sse += resid * resid;
        }

        let df = n - 2.0; // intercept + marker
        let sigma2 = sse / df;

        // Var(beta_marker) = sigma2 / Sxx
        let se_marker = (sigma2 / sxx).sqrt();
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
    dosage: &[f64],
    _ploidy: u8,
    model: GeneActionModel,
) -> Result<Vec<f64>> {
    match model {
        GeneActionModel::Additive => Ok(dosage.to_vec()),
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use binx_core::PhenotypeTable;
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

        let (aligned, order) = align_trait_to_genotypes(&pheno, "height", &geno).unwrap();
        assert_eq!(order, geno.sample_ids);
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
        // Slight noise so variance is non-zero.
        let y = array![1.0, 2.9, 6.1];

        let results = run_lm_gwas(&geno, &y, GeneActionModel::Additive).unwrap();
        let r = &results[0];
        assert!((r.beta - 2.55).abs() < 1e-6);
        assert!(r.p_value < 0.2); // With df=1 heavy tail; just ensure it is reasonably small.
    }
}
