use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use ndarray::Array1;

pub struct LocusData {
    pub id: String,
    pub ref_counts: Array1<u32>,
    pub total_counts: Array1<u32>,
}

pub fn parse_two_line_csv(path: &str) -> Result<Vec<LocusData>, Box<dyn Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let mut loci = Vec::new();
    let mut header_skipped = false;

    while let Some(line1_res) = lines.next() {
        let line1 = line1_res?;
        if line1.trim().is_empty() { continue; }
        
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
            eprintln!("Warning: Locus ID mismatch between Ref and Total lines: {} vs {}", parts1[0], parts2[0]);
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
