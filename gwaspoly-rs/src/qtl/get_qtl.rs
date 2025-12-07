//! QTL detection from GWAS results
//!
//! Implements R/GWASpoly's get.QTL function for identifying significant QTL
//! and optionally pruning nearby redundant signals.

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use crate::MarkerResult;

/// Result from QTL detection
#[derive(Debug, Clone)]
pub struct QtlResult {
    /// Marker identifier
    pub marker_id: String,
    /// Chromosome
    pub chrom: String,
    /// Position in base pairs
    pub pos: f64,
    /// Genetic model (e.g., "additive", "general")
    pub model: String,
    /// -log10(p-value) score
    pub score: f64,
    /// Effect size (if available)
    pub effect: Option<f64>,
    /// Significance threshold used
    pub threshold: f64,
}

/// Detect QTLs from GWAS marker results
///
/// This function:
/// 1. Filters markers to those exceeding their model's threshold
/// 2. Optionally prunes nearby markers using bp_window (keeping most significant)
///
/// # Arguments
/// * `results` - Slice of MarkerResult from GWAS
/// * `bp_window` - Optional window size for pruning nearby signals (in bp)
///
/// # Returns
/// Vector of QtlResult for significant markers
pub fn get_qtl(results: &[MarkerResult], bp_window: Option<u64>) -> Result<Vec<QtlResult>> {
    // Group results by model
    let mut by_model: HashMap<String, Vec<&MarkerResult>> = HashMap::new();
    for r in results {
        by_model.entry(r.model.clone()).or_default().push(r);
    }

    let mut qtls = Vec::new();

    for (_model, markers) in by_model {
        // Filter to significant markers (those with threshold and score >= threshold)
        let mut significant: Vec<&MarkerResult> = markers
            .into_iter()
            .filter(|m| {
                if let Some(thresh) = m.threshold {
                    m.score >= thresh
                } else {
                    false
                }
            })
            .collect();

        if significant.is_empty() {
            continue;
        }

        // Sort by score descending (most significant first)
        significant.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());

        // Apply bp_window pruning if specified
        let pruned = if let Some(window) = bp_window {
            prune_by_window(&significant, window)
        } else {
            significant
        };

        // Convert to QtlResult
        for m in pruned {
            qtls.push(QtlResult {
                marker_id: m.marker_id.clone(),
                chrom: m.chrom.clone().unwrap_or_else(|| "NA".to_string()),
                pos: m.pos.unwrap_or(0.0),
                model: m.model.clone(),
                score: m.score,
                effect: m.effect,
                threshold: m.threshold.unwrap_or(0.0),
            });
        }
    }

    // Sort output by model, then chromosome, then position
    qtls.sort_by(|a, b| {
        a.model
            .cmp(&b.model)
            .then_with(|| natural_chrom_cmp(&a.chrom, &b.chrom))
            .then_with(|| a.pos.partial_cmp(&b.pos).unwrap())
    });

    Ok(qtls)
}

/// Prune markers within a bp window, keeping the most significant
///
/// R/GWASpoly algorithm:
/// 1. Group by chromosome
/// 2. Sort by score (descending)
/// 3. For each marker, check distance to all previously retained markers
/// 4. Keep only if distance > window to all retained markers
fn prune_by_window<'a>(markers: &[&'a MarkerResult], window: u64) -> Vec<&'a MarkerResult> {
    // Group by chromosome
    let mut by_chrom: HashMap<String, Vec<&'a MarkerResult>> = HashMap::new();
    for m in markers {
        let chrom = m.chrom.clone().unwrap_or_else(|| "NA".to_string());
        by_chrom.entry(chrom).or_default().push(m);
    }

    let mut retained = Vec::new();

    for (_chrom, mut chrom_markers) in by_chrom {
        // Sort by score descending (should already be sorted, but ensure)
        chrom_markers.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());

        let mut chrom_retained: Vec<&MarkerResult> = Vec::new();

        for m in chrom_markers {
            let pos = m.pos.unwrap_or(0.0) as i64;

            // Check if this marker is far enough from all retained markers
            let is_far_enough = chrom_retained.iter().all(|r| {
                let r_pos = r.pos.unwrap_or(0.0) as i64;
                (pos - r_pos).unsigned_abs() > window
            });

            if is_far_enough {
                chrom_retained.push(m);
            }
        }

        retained.extend(chrom_retained);
    }

    retained
}

/// Natural chromosome comparison (chr1 < chr2 < chr10 < chrX)
fn natural_chrom_cmp(a: &str, b: &str) -> std::cmp::Ordering {
    // Strip common prefixes
    let a_num = a.trim_start_matches("chr").trim_start_matches("Chr").trim_start_matches("CHR");
    let b_num = b.trim_start_matches("chr").trim_start_matches("Chr").trim_start_matches("CHR");

    // Try to parse as numbers
    match (a_num.parse::<u32>(), b_num.parse::<u32>()) {
        (Ok(an), Ok(bn)) => an.cmp(&bn),
        (Ok(_), Err(_)) => std::cmp::Ordering::Less, // Numbers before letters
        (Err(_), Ok(_)) => std::cmp::Ordering::Greater,
        (Err(_), Err(_)) => a.cmp(b), // Fall back to string comparison
    }
}

/// Load GWAS results from a file and detect QTLs
///
/// The input file should be a CSV with columns:
/// marker_id, chrom, pos, model, score, p_value, effect, n_obs, threshold
pub fn get_qtl_from_file(path: &str, bp_window: Option<u64>) -> Result<Vec<QtlResult>> {
    let results = load_marker_results(path)?;
    get_qtl(&results, bp_window)
}

/// Load MarkerResult from a GWAS results CSV file
fn load_marker_results(path: &str) -> Result<Vec<MarkerResult>> {
    let file = File::open(path).with_context(|| format!("Failed to open {}", path))?;
    let reader = BufReader::new(file);
    let mut results = Vec::new();

    let mut lines = reader.lines();

    // Parse header
    let header = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("Empty file"))??;
    let columns: Vec<&str> = header.split(',').map(|s| s.trim()).collect();

    // Find column indices
    let find_col = |name: &str| -> Option<usize> {
        columns.iter().position(|&c| c == name)
    };

    let idx_marker = find_col("marker_id").ok_or_else(|| anyhow::anyhow!("Missing marker_id column"))?;
    let idx_chrom = find_col("chrom");
    let idx_pos = find_col("pos");
    let idx_model = find_col("model").ok_or_else(|| anyhow::anyhow!("Missing model column"))?;
    let idx_score = find_col("score").ok_or_else(|| anyhow::anyhow!("Missing score column"))?;
    let idx_pvalue = find_col("p_value");
    let idx_effect = find_col("effect");
    let idx_nobs = find_col("n_obs");
    let idx_threshold = find_col("threshold");

    for line in lines {
        let line = line?;
        let fields: Vec<&str> = line.split(',').map(|s| s.trim()).collect();

        if fields.len() <= idx_marker {
            continue;
        }

        let parse_optional_f64 = |idx: Option<usize>| -> Option<f64> {
            idx.and_then(|i| fields.get(i))
                .and_then(|s| {
                    if *s == "NA" || s.is_empty() {
                        None
                    } else {
                        s.parse().ok()
                    }
                })
        };

        let parse_optional_str = |idx: Option<usize>| -> Option<String> {
            idx.and_then(|i| fields.get(i))
                .and_then(|s| {
                    if *s == "NA" || s.is_empty() {
                        None
                    } else {
                        Some(s.to_string())
                    }
                })
        };

        let score = fields.get(idx_score)
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.0);

        let p_value = parse_optional_f64(idx_pvalue).unwrap_or_else(|| {
            // Calculate p_value from score if not present
            10_f64.powf(-score)
        });

        let n_obs = idx_nobs
            .and_then(|i| fields.get(i))
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        results.push(MarkerResult {
            marker_id: fields[idx_marker].to_string(),
            chrom: parse_optional_str(idx_chrom),
            pos: parse_optional_f64(idx_pos),
            model: fields[idx_model].to_string(),
            score,
            p_value,
            effect: parse_optional_f64(idx_effect),
            n_obs,
            threshold: parse_optional_f64(idx_threshold),
        });
    }

    Ok(results)
}

/// Write QTL results to a file
pub fn write_qtl_results(qtls: &[QtlResult], path: &str) -> Result<()> {
    let file = File::create(path).with_context(|| format!("Failed to create {}", path))?;
    let mut writer = BufWriter::new(file);

    // Write header
    writeln!(writer, "marker_id,chrom,pos,model,score,effect,threshold")?;

    // Write results
    for qtl in qtls {
        let effect_str = qtl.effect.map(|e| format!("{:.6}", e)).unwrap_or_else(|| "NA".to_string());
        writeln!(
            writer,
            "{},{},{:.0},{},{:.4},{},{}",
            qtl.marker_id,
            qtl.chrom,
            qtl.pos,
            qtl.model,
            qtl.score,
            effect_str,
            qtl.threshold
        )?;
    }

    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_marker(id: &str, chrom: &str, pos: f64, model: &str, score: f64, threshold: f64) -> MarkerResult {
        MarkerResult {
            marker_id: id.to_string(),
            chrom: Some(chrom.to_string()),
            pos: Some(pos),
            model: model.to_string(),
            score,
            p_value: 10_f64.powf(-score),
            effect: Some(0.5),
            n_obs: 100,
            threshold: Some(threshold),
        }
    }

    #[test]
    fn test_get_qtl_filters_significant() {
        let results = vec![
            make_marker("m1", "1", 1000.0, "additive", 6.0, 5.0), // Significant
            make_marker("m2", "1", 2000.0, "additive", 4.0, 5.0), // Not significant
            make_marker("m3", "1", 3000.0, "additive", 5.5, 5.0), // Significant
        ];

        let qtls = get_qtl(&results, None).unwrap();
        assert_eq!(qtls.len(), 2);
        assert!(qtls.iter().any(|q| q.marker_id == "m1"));
        assert!(qtls.iter().any(|q| q.marker_id == "m3"));
    }

    #[test]
    fn test_get_qtl_bp_window_pruning() {
        let results = vec![
            make_marker("m1", "1", 1000.0, "additive", 7.0, 5.0),     // Keep (highest)
            make_marker("m2", "1", 1500.0, "additive", 6.0, 5.0),     // Prune (within 1kb of m1)
            make_marker("m3", "1", 3000.0, "additive", 5.5, 5.0),     // Keep (>1kb from m1)
            make_marker("m4", "1", 3200.0, "additive", 5.2, 5.0),     // Prune (within 1kb of m3)
        ];

        let qtls = get_qtl(&results, Some(1000)).unwrap();
        assert_eq!(qtls.len(), 2);
        assert!(qtls.iter().any(|q| q.marker_id == "m1"));
        assert!(qtls.iter().any(|q| q.marker_id == "m3"));
    }

    #[test]
    fn test_get_qtl_multiple_models() {
        let results = vec![
            make_marker("m1", "1", 1000.0, "additive", 6.0, 5.0),
            make_marker("m1", "1", 1000.0, "general", 5.5, 5.0),
            make_marker("m2", "1", 2000.0, "additive", 4.0, 5.0), // Not significant
        ];

        let qtls = get_qtl(&results, None).unwrap();
        assert_eq!(qtls.len(), 2);
        assert!(qtls.iter().any(|q| q.marker_id == "m1" && q.model == "additive"));
        assert!(qtls.iter().any(|q| q.marker_id == "m1" && q.model == "general"));
    }

    #[test]
    fn test_get_qtl_different_chromosomes() {
        // Markers on different chromosomes should not be pruned by each other
        let results = vec![
            make_marker("m1", "1", 1000.0, "additive", 7.0, 5.0),
            make_marker("m2", "2", 1000.0, "additive", 6.0, 5.0), // Different chrom, same pos
        ];

        let qtls = get_qtl(&results, Some(1000)).unwrap();
        assert_eq!(qtls.len(), 2);
    }

    #[test]
    fn test_natural_chrom_cmp() {
        assert_eq!(natural_chrom_cmp("chr1", "chr2"), std::cmp::Ordering::Less);
        assert_eq!(natural_chrom_cmp("chr2", "chr10"), std::cmp::Ordering::Less);
        assert_eq!(natural_chrom_cmp("chr10", "chrX"), std::cmp::Ordering::Less);
        assert_eq!(natural_chrom_cmp("1", "2"), std::cmp::Ordering::Less);
    }
}
