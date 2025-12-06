//! binx-plotting: GWAS visualization tools for Binx
//!
//! This crate provides plotting functionality for genome-wide association study results,
//! including Manhattan plots and QQ plots.
//!
//! ## Features
//! - Manhattan plots with customizable themes and significance thresholds
//! - QQ plots with confidence bands
//! - SVG output (default)
//! - PNG output (optional, requires `png` feature)
//!
//! ## Example
//! ```ignore
//! use binx_plotting::{manhattan_plot, qq_plot, PlotConfig};
//!
//! // Load GWAS results and create Manhattan plot
//! let results = binx_plotting::load_gwas_results("gwas_results.csv")?;
//! manhattan_plot(&results, "manhattan.svg", PlotConfig::default())?;
//!
//! // Create QQ plot
//! qq_plot(&results, "qq.svg", PlotConfig::default())?;
//! ```

pub mod manhattan;
pub mod qq;
pub mod themes;
pub mod output;

use anyhow::Result;
use serde::Deserialize;
use std::path::Path;

/// A GWAS result point for plotting
#[derive(Debug, Clone, Deserialize)]
pub struct GwasPoint {
    pub marker_id: String,
    #[serde(default)]
    pub chrom: Option<String>,
    #[serde(default)]
    pub pos: Option<f64>,
    pub model: String,
    pub score: f64,      // -log10(p_value)
    pub p_value: f64,
    #[serde(default)]
    pub effect: Option<f64>,
    pub n_obs: usize,
}

/// Configuration for plot appearance
#[derive(Debug, Clone)]
pub struct PlotConfig {
    /// Plot width in pixels
    pub width: u32,
    /// Plot height in pixels
    pub height: u32,
    /// Significance threshold as -log10(p), default 5.0 (p = 1e-5)
    pub significance_threshold: f64,
    /// Suggestive threshold as -log10(p), default 3.0 (p = 1e-3)
    pub suggestive_threshold: Option<f64>,
    /// Plot title
    pub title: Option<String>,
    /// Color theme
    pub theme: themes::Theme,
    /// Point size
    pub point_size: u32,
    /// Show chromosome labels on x-axis
    pub show_chrom_labels: bool,
    /// Filter to specific chromosomes (None = show all)
    pub chromosomes: Option<Vec<String>>,
}

impl Default for PlotConfig {
    fn default() -> Self {
        Self {
            width: 1200,
            height: 600,
            significance_threshold: 5.0,
            suggestive_threshold: Some(3.0),
            title: None,
            theme: themes::Theme::default(),
            point_size: 3,
            show_chrom_labels: true,
            chromosomes: None,
        }
    }
}

/// Load GWAS results from a CSV file
pub fn load_gwas_results<P: AsRef<Path>>(path: P) -> Result<Vec<GwasPoint>> {
    let mut reader = csv::Reader::from_path(path)?;
    let mut results = Vec::new();

    for record in reader.deserialize() {
        let point: GwasPoint = record?;
        results.push(point);
    }

    Ok(results)
}

/// Filter results to a single model (useful when multiple models are present)
pub fn filter_by_model(results: &[GwasPoint], model: &str) -> Vec<GwasPoint> {
    results
        .iter()
        .filter(|p| p.model == model)
        .cloned()
        .collect()
}

/// Get unique models present in the results
pub fn get_models(results: &[GwasPoint]) -> Vec<String> {
    let mut models: Vec<String> = results
        .iter()
        .map(|p| p.model.clone())
        .collect();
    models.sort();
    models.dedup();
    models
}

/// Get unique chromosomes present in the results
pub fn get_chromosomes(results: &[GwasPoint]) -> Vec<String> {
    let mut chroms: Vec<String> = results
        .iter()
        .filter_map(|p| p.chrom.clone())
        .collect();
    chroms.sort_by(|a, b| {
        // Natural sort: numeric chromosomes first, then alphabetic
        let a_num: Option<u32> = a.trim_start_matches(|c: char| !c.is_numeric()).parse().ok();
        let b_num: Option<u32> = b.trim_start_matches(|c: char| !c.is_numeric()).parse().ok();
        match (a_num, b_num) {
            (Some(an), Some(bn)) => an.cmp(&bn),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => a.cmp(b),
        }
    });
    chroms.dedup();
    chroms
}

/// Filter results to specific chromosomes
pub fn filter_by_chromosomes(results: &[GwasPoint], chromosomes: &[String]) -> Vec<GwasPoint> {
    results
        .iter()
        .filter(|p| {
            p.chrom.as_ref().map_or(false, |c| chromosomes.contains(c))
        })
        .cloned()
        .collect()
}

// Re-export main functions
pub use manhattan::manhattan_plot;
pub use qq::qq_plot;
