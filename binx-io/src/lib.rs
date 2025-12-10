//! binx-io: File I/O utilities for the binx genomic toolkit.
//!
//! This crate provides I/O functions for reading and writing genomic data files:
//! - CSV/TSV phenotype files
//! - Genotype dosage matrices
//! - Kinship matrices
//! - Principal components
//!
//! Extended formats (VCF, PLINK) are handled by binx-dosage and binx-convert.

use anyhow::{anyhow, Result};
use ndarray::{Array1, Array2};
use std::collections::HashMap;
use std::path::Path;

pub type SampleId = String;
pub type MarkerId = String;

/// Optional marker metadata for genomic coordinates.
#[derive(Clone, Debug)]
pub struct MarkerMetadata {
    pub chrom: String,
    pub pos: f64,
}

/// Principal components matrix.
#[derive(Clone, Debug)]
pub struct PcMatrix {
    pub sample_ids: Vec<SampleId>,
    /// shape: (n_samples, n_pcs)
    pub pcs: Array2<f64>,
}

/// Simple phenotype table: samples × (traits + covariates).
#[derive(Clone, Debug)]
pub struct PhenotypeTable {
    pub sample_ids: Vec<SampleId>,
    pub traits: HashMap<String, Array1<f64>>,
    pub covariates: HashMap<String, Array1<f64>>,
    pub factor_covariates: HashMap<String, Vec<String>>,
}

/// Biallelic dosage matrix: markers × samples, entries 0..ploidy.
#[derive(Clone, Debug)]
pub struct GenotypeMatrix {
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
pub fn load_phenotypes<P: AsRef<Path>>(path: P) -> Result<PhenotypeTable> {
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
pub fn load_genotypes<P: AsRef<Path>>(path: P, ploidy: u8) -> Result<GenotypeMatrix> {
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

    Ok(GenotypeMatrix {
        ploidy,
        sample_ids,
        marker_ids,
        marker_metadata: Some(marker_metadata),
        dosages,
    })
}

/// Load kinship matrix from TSV.
/// First column is sample_id, remaining columns are the matrix values.
pub fn load_kinship<P: AsRef<Path>>(path: P) -> Result<KinshipMatrix> {
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
pub fn load_pcs<P: AsRef<Path>>(path: P) -> Result<PcMatrix> {
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

/// Write kinship matrix to TSV.
pub fn write_kinship<P: AsRef<Path>>(path: P, kin: &KinshipMatrix) -> Result<()> {
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
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_detect_delimiter_comma() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "a,b,c").unwrap();
        assert_eq!(detect_delimiter(f.path()).unwrap(), b',');
    }

    #[test]
    fn test_detect_delimiter_tab() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "a\tb\tc").unwrap();
        assert_eq!(detect_delimiter(f.path()).unwrap(), b'\t');
    }

    #[test]
    fn test_load_phenotypes() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "sample_id,trait1,trait2").unwrap();
        writeln!(f, "S1,1.0,2.0").unwrap();
        writeln!(f, "S2,3.0,4.0").unwrap();
        f.flush().unwrap();

        let pheno = load_phenotypes(f.path()).unwrap();
        assert_eq!(pheno.sample_ids, vec!["S1", "S2"]);
        assert!(pheno.traits.contains_key("trait1"));
        assert!(pheno.traits.contains_key("trait2"));
    }
}
