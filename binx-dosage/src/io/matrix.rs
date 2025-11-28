use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};

use ndarray::{Array1, Array2};

use super::{LocusData, MatrixData, TwoLineData};

pub fn parse_two_line_csv(path: &str) -> Result<TwoLineData, Box<dyn Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // First non-empty line must be header with sample names
    let header_line = loop {
        match lines.next() {
            Some(line) => {
                let l = line?;
                if !l.trim().is_empty() {
                    break l;
                }
            }
            None => return Err("CSV file is empty".into()),
        }
    };

    let header_parts: Vec<String> = header_line
        .split(',')
        .map(|s| s.trim().to_string())
        .collect();
    if header_parts.len() < 2 {
        return Err("Header must include locus column plus at least one sample".into());
    }
    let sample_names = header_parts[1..].to_vec();

    let mut loci = Vec::new();
    while let Some(line1_res) = lines.next() {
        let line1 = line1_res?;
        if line1.trim().is_empty() {
            continue;
        }

        let line2 = match lines.next() {
            Some(l) => l?,
            None => return Err("Ref/Total line pairs are incomplete".into()),
        };

        let parts1: Vec<&str> = line1.split(',').map(|s| s.trim()).collect();
        let parts2: Vec<&str> = line2.split(',').map(|s| s.trim()).collect();

        if parts1.len() != header_parts.len() || parts2.len() != header_parts.len() {
            return Err("Row length does not match header sample count".into());
        }

        let id = parts1[0].to_string();

        if parts1[0] != parts2[0] {
            return Err(
                format!(
                    "Locus ID mismatch between Ref and Total lines: {} vs {}",
                    parts1[0], parts2[0]
                )
                .into(),
            );
        }

        let refs: Result<Vec<u32>, _> = parts1[1..].iter().map(|s| s.parse::<u32>()).collect();
        let totals: Result<Vec<u32>, _> = parts2[1..].iter().map(|s| s.parse::<u32>()).collect();

        loci.push(LocusData {
            id,
            ref_counts: Array1::from(refs?),
            total_counts: Array1::from(totals?),
            vcf_chrom: None,
            vcf_pos: None,
            vcf_ref: None,
            vcf_alt: None,
        });
    }

    Ok(TwoLineData { sample_names, loci })
}

fn detect_delimiter(header_line: &str) -> char {
    if header_line.contains('\t') {
        '\t'
    } else {
        ','
    }
}

/// Parse a matrix file with a header row (samples) and first column as marker IDs.
/// Returns marker IDs, rows of counts, and sample names.
fn parse_matrix_file(path: &str) -> Result<(Vec<String>, Vec<Vec<u32>>, Vec<String>), Box<dyn Error>> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut header = String::new();
    let bytes = reader.read_line(&mut header)?;
    if bytes == 0 {
        return Err("Matrix file is empty".into());
    }

    let delim = detect_delimiter(&header);
    let header_fields: Vec<String> = header
        .trim_end_matches(&['\n', '\r'][..])
        .split(delim)
        .map(|s| s.trim().to_string())
        .collect();

    if header_fields.len() < 2 {
        return Err("Matrix header must include marker column plus at least one sample column".into());
    }

    let expected_cols = header_fields.len();
    let sample_names = header_fields[1..].to_vec();
    let mut marker_ids = Vec::new();
    let mut rows: Vec<Vec<u32>> = Vec::new();

    for (idx, line_res) in reader.lines().enumerate() {
        let line = line_res?;
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split(delim).map(|s| s.trim()).collect();
        if parts.len() != expected_cols {
            return Err(format!(
                "Row {} in {} has {} columns, expected {}",
                idx + 2,
                path,
                parts.len(),
                expected_cols
            )
            .into());
        }
        let id = parts[0].to_string();
        let counts: Result<Vec<u32>, _> = parts[1..].iter().map(|s| s.parse::<u32>()).collect();
        let counts = counts.map_err(|_| format!("Failed to parse numeric counts in {}", path))?;
        marker_ids.push(id);
        rows.push(counts);
    }

    if marker_ids.is_empty() {
        return Err("Matrix contains no marker rows".into());
    }

    Ok((marker_ids, rows, sample_names))
}

/// Parse paired ref/total matrices (markers in rows, samples in columns).
pub fn parse_ref_total_matrices(
    ref_path: &str,
    total_path: &str,
) -> Result<MatrixData, Box<dyn Error>> {
    let (ref_ids, ref_rows, ref_names) = parse_matrix_file(ref_path)?;
    let (tot_ids, tot_rows, tot_names) = parse_matrix_file(total_path)?;

    if ref_names != tot_names {
        return Err("Ref and total matrices have different sample column counts".into());
    }
    if ref_rows.len() != tot_rows.len() {
        return Err("Ref and total matrices have different numbers of marker rows".into());
    }
    for (i, (rid, tid)) in ref_ids.iter().zip(tot_ids.iter()).enumerate() {
        if rid != tid {
            return Err(
                format!("Marker ID mismatch at row {}: {} vs {}", i + 1, rid, tid).into(),
            );
        }
    }

    let markers = ref_rows.len();
    let samples = ref_names.len();

    let mut ref_flat = Vec::with_capacity(markers * samples);
    for row in &ref_rows {
        ref_flat.extend_from_slice(row);
    }
    let mut tot_flat = Vec::with_capacity(markers * samples);
    for row in &tot_rows {
        tot_flat.extend_from_slice(row);
    }

    let ref_counts = Array2::from_shape_vec((markers, samples), ref_flat)?;
    let total_counts = Array2::from_shape_vec((markers, samples), tot_flat)?;

    Ok(MatrixData {
        marker_ids: ref_ids,
        sample_names: ref_names,
        ref_counts,
        total_counts,
    })
}
