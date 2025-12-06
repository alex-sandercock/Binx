//! QQ (Quantile-Quantile) plot generation for GWAS results

use crate::{GwasPoint, PlotConfig};
use anyhow::{Result, Context};
use plotters::prelude::*;
use std::path::Path;

/// Generate a QQ plot from GWAS results
///
/// The QQ plot compares the observed -log10(p-values) against the expected
/// values under the null hypothesis (uniform distribution).
///
/// # Arguments
/// * `results` - GWAS results to plot
/// * `output_path` - Path for output file (SVG or PNG based on extension)
/// * `config` - Plot configuration
///
/// # Example
/// ```ignore
/// use binx_plotting::{qq_plot, PlotConfig, load_gwas_results};
///
/// let results = load_gwas_results("gwas_results.csv")?;
/// qq_plot(&results, "qq.svg", PlotConfig::default())?;
/// ```
pub fn qq_plot<P: AsRef<Path>>(
    results: &[GwasPoint],
    output_path: P,
    config: PlotConfig,
) -> Result<()> {
    let output_path = output_path.as_ref();

    if results.is_empty() {
        anyhow::bail!("No GWAS results to plot");
    }

    // Calculate expected and observed values
    let (expected, observed) = calculate_qq_values(results);

    let max_val = expected.iter()
        .chain(observed.iter())
        .cloned()
        .fold(0.0_f64, f64::max);
    let axis_max = (max_val * 1.1).max(1.0);

    // Determine output format from extension
    let ext = output_path.extension()
        .and_then(|e| e.to_str())
        .unwrap_or("svg")
        .to_lowercase();

    match ext.as_str() {
        "svg" => draw_qq_svg(output_path, &expected, &observed, &config, axis_max),
        #[cfg(feature = "png")]
        "png" => draw_qq_png(output_path, &expected, &observed, &config, axis_max),
        _ => anyhow::bail!("Unsupported output format: {}", ext),
    }
}

/// Calculate expected and observed -log10(p) values for QQ plot
fn calculate_qq_values(results: &[GwasPoint]) -> (Vec<f64>, Vec<f64>) {
    let n = results.len();

    // Sort p-values in ascending order
    let mut p_values: Vec<f64> = results.iter().map(|r| r.p_value).collect();
    p_values.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Calculate expected p-values under null (uniform distribution)
    // Expected p-value for rank i out of n is: i / (n + 1)
    let expected: Vec<f64> = (1..=n)
        .map(|i| -((i as f64) / ((n + 1) as f64)).log10())
        .collect();

    // Convert observed p-values to -log10
    let observed: Vec<f64> = p_values
        .iter()
        .map(|p| {
            if *p > 0.0 {
                -p.log10()
            } else {
                // Handle p = 0 by capping at a reasonable value
                16.0
            }
        })
        .collect();

    (expected, observed)
}

/// Calculate confidence band for QQ plot (95% CI)
fn calculate_confidence_band(n: usize, expected: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut lower = Vec::with_capacity(expected.len());
    let mut upper = Vec::with_capacity(expected.len());

    for (i, _exp) in expected.iter().enumerate() {
        let rank = i + 1;
        let p_expected = rank as f64 / (n + 1) as f64;

        // Beta distribution confidence intervals
        // Using normal approximation for simplicity
        let se = (p_expected * (1.0 - p_expected) / n as f64).sqrt();
        let z = 1.96; // 95% CI

        let p_lower = (p_expected - z * se).max(1e-10);
        let p_upper = (p_expected + z * se).min(1.0 - 1e-10);

        lower.push(-p_upper.log10()); // Note: inverted because -log10
        upper.push(-p_lower.log10());
    }

    (lower, upper)
}

fn draw_qq_svg(
    output_path: &Path,
    expected: &[f64],
    observed: &[f64],
    config: &PlotConfig,
    axis_max: f64,
) -> Result<()> {
    let root = SVGBackend::new(output_path, (config.height, config.height)) // Square plot
        .into_drawing_area();

    draw_qq_impl(&root, expected, observed, config, axis_max)
        .context("Failed to draw QQ plot")?;

    root.present().context("Failed to write SVG")?;
    Ok(())
}

#[cfg(feature = "png")]
fn draw_qq_png(
    output_path: &Path,
    expected: &[f64],
    observed: &[f64],
    config: &PlotConfig,
    axis_max: f64,
) -> Result<()> {
    let root = BitMapBackend::new(output_path, (config.height, config.height))
        .into_drawing_area();

    draw_qq_impl(&root, expected, observed, config, axis_max)
        .context("Failed to draw QQ plot")?;

    root.present().context("Failed to write PNG")?;
    Ok(())
}

fn draw_qq_impl<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    expected: &[f64],
    observed: &[f64],
    config: &PlotConfig,
    axis_max: f64,
) -> Result<(), DrawingAreaErrorKind<DB::ErrorType>> {
    root.fill(&config.theme.background)?;

    let title = config.title.as_deref().unwrap_or("QQ Plot");

    let mut chart = ChartBuilder::on(root)
        .caption(title, ("sans-serif", 24).into_font().color(&config.theme.text))
        .margin(10)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(0.0..axis_max, 0.0..axis_max)?;

    chart
        .configure_mesh()
        .x_desc("Expected -log₁₀(p)")
        .y_desc("Observed -log₁₀(p)")
        .x_label_style(("sans-serif", 14).into_font().color(&config.theme.text))
        .y_label_style(("sans-serif", 14).into_font().color(&config.theme.text))
        .axis_style(&config.theme.axis)
        .draw()?;

    // Draw confidence band
    let (ci_lower, ci_upper) = calculate_confidence_band(observed.len(), expected);

    // Draw as filled area
    let ci_points: Vec<_> = expected.iter()
        .zip(ci_lower.iter())
        .zip(ci_upper.iter())
        .map(|((e, l), u)| (*e, *l, *u))
        .collect();

    for (exp, lower, upper) in &ci_points {
        chart.draw_series(std::iter::once(Rectangle::new(
            [(*exp - 0.02, *lower), (*exp + 0.02, *upper)],
            config.theme.confidence_band.filled(),
        )))?;
    }

    // Draw diagonal reference line (y = x)
    chart.draw_series(LineSeries::new(
        vec![(0.0, 0.0), (axis_max, axis_max)],
        config.theme.reference_line.stroke_width(2),
    ))?;

    // Draw points
    let point_color = &config.theme.chromosome_colors[0];
    let point_size = config.point_size;

    chart.draw_series(
        expected.iter()
            .zip(observed.iter())
            .map(|(e, o)| Circle::new((*e, *o), point_size, point_color.filled()))
    )?;

    Ok(())
}
