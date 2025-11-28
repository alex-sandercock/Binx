use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::sync::Arc;

use flate2::read::MultiGzDecoder;
use ndarray::Array1;

use super::VcfRecordCounts;

/// Streaming configuration.
#[derive(Clone, Debug)]
pub struct VcfStreamConfig {
    /// If true, invalid rows/fields raise an error instead of being skipped with a warning.
    pub strict: bool,
    /// Minimum total depth to treat a sample as valid.
    pub min_depth: u32,
}

impl Default for VcfStreamConfig {
    fn default() -> Self {
        Self {
            strict: false,
            min_depth: 0,
        }
    }
}

fn vcf_reader(path: &str) -> Result<Box<dyn BufRead>, Box<dyn Error>> {
    if path == "-" {
        return Ok(Box::new(BufReader::with_capacity(64 * 1024, io::stdin())));
    }

    let file = File::open(path)?;
    if path.to_ascii_lowercase().ends_with(".gz") || path.to_ascii_lowercase().ends_with(".bgz") {
        let decoder = MultiGzDecoder::new(file);
        Ok(Box::new(BufReader::with_capacity(64 * 1024, decoder)))
    } else {
        Ok(Box::new(BufReader::with_capacity(64 * 1024, file)))
    }
}

fn warn_or_err(strict: bool, msg: &str) -> Result<(), Box<dyn Error>> {
    if strict {
        Err(msg.to_string().into())
    } else {
        eprintln!("Warning: {}", msg);
        Ok(())
    }
}

fn parse_uint(field: Option<&str>) -> Result<Option<u32>, ()> {
    match field {
        None => Ok(None),
        Some("") | Some(".") => Ok(None),
        Some(s) => s.parse::<u32>().map(Some).map_err(|_| ()),
    }
}

fn parse_ad(field: Option<&str>) -> Result<Option<(u32, u32)>, ()> {
    let Some(ad_str) = field else { return Ok(None); };
    if ad_str.is_empty() || ad_str == "." {
        return Ok(None);
    }
    let mut parts = ad_str.split(',');
    let ref_depth = parts.next().ok_or(())?.parse::<u32>().map_err(|_| ())?;
    let total = parts.try_fold(ref_depth, |acc, part| {
        part.parse::<u32>()
            .map(|v| acc.saturating_add(v))
            .map_err(|_| ())
    })?;
    Ok(Some((ref_depth, total)))
}

fn format_indices(format_str: &str) -> (Option<usize>, Option<usize>, Option<usize>) {
    let mut ad_idx = None;
    let mut ra_idx = None;
    let mut dp_idx = None;
    for (i, key) in format_str.split(':').enumerate() {
        match key {
            "AD" => ad_idx = Some(i),
            "RA" => ra_idx = Some(i),
            "DP" => dp_idx = Some(i),
            _ => {}
        }
    }
    (ad_idx, ra_idx, dp_idx)
}

fn parse_sample(
    sample_str: &str,
    ad_idx: Option<usize>,
    ra_idx: Option<usize>,
    dp_idx: Option<usize>,
    min_depth: u32,
    strict: bool,
) -> Result<(u32, u32, bool), Box<dyn Error>> {
    let tokens: Vec<&str> = sample_str.split(':').collect();
    let field_at = |idx: usize| -> Option<&str> { tokens.get(idx).copied() };

    // AD preferred
    if let Some(idx) = ad_idx {
        match parse_ad(field_at(idx)) {
            Ok(Some((r, t))) => return Ok((r, t, t >= min_depth && t > 0)),
            Ok(None) => {}
            Err(_) => warn_or_err(strict, "Malformed AD value")?,
        }
    }

    // Fallback: RA + DP
    let ref_depth = match ra_idx {
        Some(idx) => match parse_uint(field_at(idx)) {
            Ok(v) => v.unwrap_or(0),
            Err(_) => {
                warn_or_err(strict, "Malformed RA value")?;
                0
            }
        },
        None => 0,
    };
    let total = match dp_idx {
        Some(idx) => match parse_uint(field_at(idx)) {
            Ok(v) => v.unwrap_or(0),
            Err(_) => {
                warn_or_err(strict, "Malformed DP value")?;
                0
            }
        },
        None => 0,
    };

    let is_valid = total >= min_depth && total > 0;
    Ok((ref_depth, total, is_valid))
}

fn parse_record_line(
    line: &str,
    expected_samples: &mut Option<usize>,
    config: &VcfStreamConfig,
) -> Result<Option<VcfRecordCounts>, Box<dyn Error>> {
    // Parse fixed fields without allocating Vec
    let mut fields = line.split('\t');

    let chrom = fields.next().ok_or("Missing CHROM field")?;
    let pos_str = fields.next().ok_or("Missing POS field")?;
    let id = fields.next().ok_or("Missing ID field")?;
    let ref_allele = fields.next().ok_or("Missing REF field")?;
    let alt_allele = fields.next().ok_or("Missing ALT field")?;
    let _qual = fields.next().ok_or("Missing QUAL field")?;
    let _filter = fields.next().ok_or("Missing FILTER field")?;
    let _info = fields.next().ok_or("Missing INFO field")?;
    let format_str = fields.next().ok_or("Missing FORMAT field")?;

    // Remaining fields are samples - process without collecting

    let (ad_idx, ra_idx, dp_idx) = format_indices(format_str);
    if ad_idx.is_none() && ra_idx.is_none() && dp_idx.is_none() {
        warn_or_err(config.strict, "No AD/RA/DP in FORMAT; skipping locus")?;
        return Ok(None);
    }

    // Pre-allocate vectors if we know expected sample count
    let initial_capacity = expected_samples.unwrap_or(100);
    let mut ref_counts = Vec::with_capacity(initial_capacity);
    let mut total_counts = Vec::with_capacity(initial_capacity);
    let mut valid_mask = Vec::with_capacity(initial_capacity);
    let mut any_valid = false;
    let mut sample_count = 0;

    // Process samples from iterator without collecting
    for sample_str in fields {
        sample_count += 1;
        match parse_sample(
            sample_str,
            ad_idx,
            ra_idx,
            dp_idx,
            config.min_depth,
            config.strict,
        ) {
            Ok((r, t, is_valid)) => {
                ref_counts.push(r);
                total_counts.push(t);
                valid_mask.push(is_valid);
                any_valid |= is_valid;
            }
            Err(e) => {
                if config.strict {
                    return Err(e);
                }
                ref_counts.push(0);
                total_counts.push(0);
                valid_mask.push(false);
            }
        }
    }

    // Validate sample count matches expected
    if let Some(expected) = expected_samples {
        if sample_count != *expected {
            warn_or_err(
                config.strict,
                &format!(
                    "Sample count mismatch (expected {}, got {})",
                    expected,
                    sample_count
                ),
            )?;
            return Ok(None);
        }
    } else {
        *expected_samples = Some(sample_count);
    }

    if !any_valid {
        return Ok(None);
    }

    // Parse position as integer
    let pos: u64 = pos_str.parse().unwrap_or(0);

    // Use VCF ID field if present and not ".", otherwise use chrom:pos
    let locus_id = if id == "." {
        format!("{}:{}", chrom, pos)
    } else {
        id.to_string()
    };

    Ok(Some(VcfRecordCounts {
        id: locus_id,
        ref_counts: Array1::from(ref_counts),
        total_counts: Array1::from(total_counts),
        chrom: Arc::new(chrom.to_string()),
        pos,
        ref_allele: Arc::new(ref_allele.to_string()),
        alt_allele: Arc::new(alt_allele.to_string()),
    }))
}

/// Quickly counts the number of data records in a VCF file (excluding headers).
/// Used to enable progress bars with ETA.
pub fn count_vcf_records(path: &str) -> Result<usize, Box<dyn Error>> {
    let mut reader = vcf_reader(path)?;
    let mut line = String::new();
    let mut count = 0;
    let mut in_header = true;

    loop {
        line.clear();
        let bytes = reader.read_line(&mut line)?;
        if bytes == 0 {
            break;
        }
        let trimmed = line.trim_end_matches(&['\n', '\r'][..]);
        if trimmed.is_empty() {
            continue;
        }
        if in_header {
            if trimmed.starts_with("#") {
                continue;
            } else {
                in_header = false;
            }
        }
        count += 1;
    }
    Ok(count)
}

/// Stream VCF records using a fast text parser (supports stdin, gzip/bgzip).
/// AD is preferred; if missing, RA+DP is used. Invalid fields warn by default.
pub fn stream_vcf_records_with_config<F, H>(
    path: &str,
    mut on_header: H,
    mut on_record: F,
    config: VcfStreamConfig,
)-> Result<(), Box<dyn Error>>
where
    F: FnMut(VcfRecordCounts),
    H: FnMut(&[String]),
{
    let mut reader = vcf_reader(path)?;
    let mut line = String::with_capacity(8192); // Start with default, will resize after header
    let mut sample_count: Option<usize> = None;
    let mut line_sized = false;

    loop {
        line.clear();
        let bytes = reader.read_line(&mut line)?;
        if bytes == 0 {
            break;
        }
        let trimmed = line.trim_end_matches(&['\n', '\r'][..]);
        if trimmed.is_empty() || trimmed.starts_with("##") {
            continue;
        }
        if trimmed.starts_with("#CHROM") {
            // Parse header without collecting all fields into Vec
            let mut header_fields = trimmed.split('\t');

            // Skip first 9 fixed column headers
            let mut valid_header = true;
            for _ in 0..9 {
                if header_fields.next().is_none() {
                    warn_or_err(config.strict, "Header has fewer than 9 columns")?;
                    valid_header = false;
                    break;
                }
            }

            if !valid_header {
                continue;
            }

            // Collect only sample names (not all fields)
            let names: Vec<String> = header_fields.map(|s| s.to_string()).collect();
            sample_count = Some(names.len());

            // Resize line buffer based on expected line size
            // Formula: ~100 bytes for fixed fields + (samples * 20 bytes average per sample field)
            if !line_sized && names.len() > 0 {
                let expected_line_size = 100 + (names.len() * 20);
                if expected_line_size > line.capacity() {
                    line.reserve(expected_line_size - line.capacity());
                }
                line_sized = true;
            }

            on_header(&names);
            continue;
        }

        match parse_record_line(trimmed, &mut sample_count, &config)? {
            Some(rec) => on_record(rec),
            None => continue,
        }
    }

    Ok(())
}

/// Default configuration: non-strict, no depth filtering.
pub fn stream_vcf_records<F, H>(path: &str, on_header: H, on_record: F) -> Result<(), Box<dyn Error>>
where
    F: FnMut(VcfRecordCounts),
    H: FnMut(&[String]),
{
    stream_vcf_records_with_config(path, on_header, on_record, VcfStreamConfig::default())
}
