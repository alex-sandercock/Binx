//! I/O functions for gwaspoly-rs
//!
//! This module mirrors R/GWASpoly's read.GWASpoly function for loading data.

use crate::types::{GenotypeMatrixBiallelic, KinshipMatrix, MarkerMetadata, PhenotypeTable};
use anyhow::{anyhow, Result};
use ndarray::Array2;
use std::collections::HashMap;
use std::path::Path;

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

/// Read genotype and phenotype data for GWASpoly analysis.
///
/// This function mirrors R/GWASpoly's `read.GWASpoly()` function.
///
/// # Arguments
/// * `geno_path` - Path to genotype file (TSV/CSV with marker_id, chr, pos, <samples...>)
/// * `pheno_path` - Path to phenotype file (TSV/CSV with sample_id, <traits...>)
/// * `ploidy` - Ploidy level (e.g., 2 for diploid, 4 for tetraploid)
///
/// # Returns
/// A tuple of (GenotypeMatrixBiallelic, PhenotypeTable)
pub fn read_gwaspoly(
    geno_path: &str,
    pheno_path: &str,
    ploidy: u8,
) -> Result<(GenotypeMatrixBiallelic, PhenotypeTable)> {
    let geno = load_genotypes(geno_path, ploidy)?;
    let pheno = load_phenotypes(pheno_path)?;
    Ok((geno, pheno))
}

/// Load biallelic dosage matrix from TSV/CSV.
///
/// Expecting format: marker_id, chr, pos, <sample1>, <sample2>, ...
pub fn load_genotypes<P: AsRef<Path>>(path: P, ploidy: u8) -> Result<GenotypeMatrixBiallelic> {
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

/// Load phenotype table from TSV/CSV.
///
/// First column is sample_id, remaining columns are traits/covariates.
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
            traits.insert(name.clone(), numeric_vals.clone());
            covariates.insert(name, numeric_vals);
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

/// Load kinship matrix from TSV/CSV.
///
/// First column is sample_id, remaining columns are the matrix values.
/// The header row should have sample IDs for columns.
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_detect_delimiter_csv() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "a,b,c").unwrap();
        assert_eq!(detect_delimiter(f.path()).unwrap(), b',');
    }

    #[test]
    fn test_detect_delimiter_tsv() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "a\tb\tc").unwrap();
        assert_eq!(detect_delimiter(f.path()).unwrap(), b'\t');
    }
}
