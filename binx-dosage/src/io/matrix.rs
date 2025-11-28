use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};

use ndarray::{Array1, Array2};

use super::{LocusData, MatrixData};

pub fn parse_two_line_csv(path: &str) -> Result<Vec<LocusData>, Box<dyn Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let mut loci = Vec::new();
    let mut header_skipped = false;

    while let Some(line1_res) = lines.next() {
        let line1 = line1_res?;
        if line1.trim().is_empty() {
            continue;
        }

        if !header_skipped && (line1.contains("Ref") || line1.contains("sample")) {
            let _ = lines.next();
            header_skipped = true;
            continue;
        }

        let line2 = match lines.next() {
            Some(l) => l?,
            None => break,
        };

        let parts1: Vec<&str> = line1.split(',').map(|s| s.trim()).collect();
        let parts2: Vec<&str> = line2.split(',').map(|s| s.trim()).collect();

        if parts1.len() < 2 || parts2.len() < 2 {
            continue;
        }

        let id = parts1[0].to_string();

        if parts1[0] != parts2[0] {
            eprintln!(
                "Warning: Locus ID mismatch between Ref and Total lines: {} vs {}",
                parts1[0], parts2[0]
            );
        }

        let refs: Result<Vec<u32>, _> = parts1[1..].iter().map(|s| s.parse::<u32>()).collect();
        let totals: Result<Vec<u32>, _> = parts2[1..].iter().map(|s| s.parse::<u32>()).collect();

        if let (Ok(r), Ok(t)) = (refs, totals) {
            loci.push(LocusData {
                id,
                ref_counts: Array1::from(r),
                total_counts: Array1::from(t),
            });
        }
    }

    Ok(loci)
}

fn detect_delimiter(header_line: &str) -> char {
    if header_line.contains('\t') {
        '\t'
    } else {
        ','
    }
}

/// Parse a matrix file with a header row (samples) and first column as marker IDs.
/// Returns marker IDs, rows of counts, and the number of samples (columns minus ID).
fn parse_matrix_file(path: &str) -> Result<(Vec<String>, Vec<Vec<u32>>, usize), Box<dyn Error>> {
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
    let sample_count = expected_cols - 1;
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

    Ok((marker_ids, rows, sample_count))
}

/// Parse paired ref/total matrices (markers in rows, samples in columns).
pub fn parse_ref_total_matrices(
    ref_path: &str,
    total_path: &str,
) -> Result<MatrixData, Box<dyn Error>> {
    let (ref_ids, ref_rows, ref_samples) = parse_matrix_file(ref_path)?;
    let (tot_ids, tot_rows, tot_samples) = parse_matrix_file(total_path)?;

    if ref_samples != tot_samples {
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
    let samples = ref_samples;

    let mut ref_flat = Vec::with_capacity(markers * samples);
    for row in &ref_rows {
        if row.len() != samples {
            return Err("Inconsistent ref row length".into());
        }
        ref_flat.extend_from_slice(row);
    }
    let mut tot_flat = Vec::with_capacity(markers * samples);
    for row in &tot_rows {
        if row.len() != samples {
            return Err("Inconsistent total row length".into());
        }
        tot_flat.extend_from_slice(row);
    }

    let ref_counts = Array2::from_shape_vec((markers, samples), ref_flat)?;
    let total_counts = Array2::from_shape_vec((markers, samples), tot_flat)?;

    Ok(MatrixData {
        marker_ids: ref_ids,
        ref_counts,
        total_counts,
    })
}
