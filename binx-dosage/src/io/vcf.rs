use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

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
    let parts: Vec<&str> = line.split('\t').collect();
    if parts.len() < 9 {
        warn_or_err(
            config.strict,
            &format!("Line has {} columns (expected >=9)", parts.len()),
        )?;
        return Ok(None);
    }

    let chrom = parts[0];
    let pos = parts[1];
    let format_str = parts[8];
    let sample_fields = &parts[9..];

    let (ad_idx, ra_idx, dp_idx) = format_indices(format_str);
    if ad_idx.is_none() && ra_idx.is_none() && dp_idx.is_none() {
        warn_or_err(config.strict, "No AD/RA/DP in FORMAT; skipping locus")?;
        return Ok(None);
    }

    if let Some(expected) = expected_samples {
        if sample_fields.len() != *expected {
            warn_or_err(
                config.strict,
                &format!(
                    "Sample count mismatch (expected {}, got {})",
                    expected,
                    sample_fields.len()
                ),
            )?;
            return Ok(None);
        }
    } else {
        *expected_samples = Some(sample_fields.len());
    }

    let mut ref_counts = Vec::with_capacity(sample_fields.len());
    let mut total_counts = Vec::with_capacity(sample_fields.len());
    let mut valid_mask = Vec::with_capacity(sample_fields.len());
    let mut any_valid = false;

    for sample_str in sample_fields {
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

    if !any_valid {
        return Ok(None);
    }

    let locus_id = format!("{}:{}", chrom, pos);
    Ok(Some(VcfRecordCounts {
        id: locus_id,
        ref_counts: Array1::from(ref_counts),
        total_counts: Array1::from(total_counts),
    }))
}

/// Stream VCF records using a fast text parser (supports stdin, gzip/bgzip).
/// AD is preferred; if missing, RA+DP is used. Invalid fields warn by default.
pub fn stream_vcf_records_with_config<F>(
    path: &str,
    mut on_record: F,
    config: VcfStreamConfig,
) -> Result<(), Box<dyn Error>>
where
    F: FnMut(VcfRecordCounts),
{
    let mut reader = vcf_reader(path)?;
    let mut line = String::with_capacity(8192);
    let mut sample_count: Option<usize> = None;

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
            let cols = trimmed.split('\t').count();
            if cols < 9 {
                warn_or_err(config.strict, "Header has fewer than 9 columns")?;
            } else {
                sample_count = cols.checked_sub(9);
            }
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
pub fn stream_vcf_records<F>(path: &str, on_record: F) -> Result<(), Box<dyn Error>>
where
    F: FnMut(VcfRecordCounts),
{
    stream_vcf_records_with_config(path, on_record, VcfStreamConfig::default())
}
