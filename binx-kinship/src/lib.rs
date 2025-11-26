use anyhow::{anyhow, Result};
use binx_core::{load_genotypes_biallelic_from_tsv, GenotypeMatrixBiallelic, KinshipMatrix};
use ndarray::{Array2, Axis};

/// Compute VanRaden kinship matrix from a biallelic genotype TSV and write to TSV.
pub fn run_kinship(geno_path: &str, ploidy: u8, out_path: &str) -> Result<()> {
    let geno = load_genotypes_biallelic_from_tsv(geno_path, ploidy)?;
    let kin = compute_kinship_vanraden(&geno)?;
    write_kinship_tsv(out_path, &kin)
}

/// Compute additive VanRaden kinship matrix for biallelic genotypes.
/// For ploidy m, we center dosages by m * p and scale by m * sum(p(1-p)).
pub fn compute_kinship_vanraden(geno: &GenotypeMatrixBiallelic) -> Result<KinshipMatrix> {
    let n_markers = geno.marker_ids.len();
    let n_samples = geno.sample_ids.len();
    if n_markers == 0 || n_samples == 0 {
        return Err(anyhow!("Genotype matrix is empty"));
    }

    // Compute allele frequencies per marker: p = mean(dosage) / ploidy.
    let mut freqs = Vec::with_capacity(n_markers);
    for i in 0..n_markers {
        let row = geno.dosages.index_axis(Axis(0), i);
        let mean = row.sum() / (n_samples as f64);
        freqs.push(mean / (geno.ploidy as f64));
    }

    // Build centered matrix M (markers x samples): x_ij - m * p_j
    let mut centered = Array2::<f64>::zeros((n_markers, n_samples));
    for i in 0..n_markers {
        let p = freqs[i];
        let offset = (geno.ploidy as f64) * p;
        for j in 0..n_samples {
            centered[(i, j)] = geno.dosages[(i, j)] - offset;
        }
    }

    // Denominator: m * sum(p(1-p)) over markers.
    let denom: f64 = freqs
        .iter()
        .map(|p| (geno.ploidy as f64) * p * (1.0 - p))
        .sum();
    if denom == 0.0 {
        return Err(anyhow!(
            "Kinship denominator is zero; check allele frequencies or ploidy"
        ));
    }

    // K = (M^T M) / denom
    let mtm = centered.t().dot(&centered);
    let kin = mtm / denom;

    Ok(KinshipMatrix {
        sample_ids: geno.sample_ids.clone(),
        matrix: kin,
    })
}

pub fn write_kinship_tsv(path: &str, kin: &KinshipMatrix) -> Result<()> {
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(path)?;

    let mut header = Vec::with_capacity(kin.sample_ids.len() + 1);
    header.push("sample_id".to_string());
    header.extend(kin.sample_ids.iter().cloned());
    wtr.write_record(&header)?;

    for (i, sid) in kin.sample_ids.iter().enumerate() {
        let mut row = Vec::with_capacity(kin.sample_ids.len() + 1);
        row.push(sid.clone());
        for j in 0..kin.sample_ids.len() {
            row.push(kin.matrix[(i, j)].to_string());
        }
        wtr.write_record(&row)?;
    }
    wtr.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use binx_core::GenotypeMatrixBiallelic;
    use ndarray::array;
    use std::fs;
    use tempfile::NamedTempFile;

    #[test]
    fn writes_expected_tsv() {
        let geno = GenotypeMatrixBiallelic {
            ploidy: 2,
            sample_ids: vec!["S1".into(), "S2".into()],
            marker_ids: vec!["m1".into(), "m2".into()],
            dosages: array![[0.0, 2.0], [1.0, 1.0]],
        };
        let kin = compute_kinship_vanraden(&geno).unwrap();
        let file = NamedTempFile::new().unwrap();
        write_kinship_tsv(file.path().to_str().unwrap(), &kin).unwrap();
        let contents = fs::read_to_string(file.path()).unwrap();
        assert!(contents.contains("sample_id\tS1\tS2"));
    }
}
