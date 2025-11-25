use anyhow::Result;
use ndarray::{Array1, Array2};
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
                anyhow::anyhow!(
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
                anyhow::anyhow!("Failed to parse dosage '{}' at col {}: {}", val_str, i, e)
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
                anyhow::anyhow!("Failed to parse kinship value '{}': {}", val_str, e)
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
}
