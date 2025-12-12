//! QQ (Quantile-Quantile) plot generation for GWAS results

use crate::{GwasPoint, PlotConfig};
use anyhow::{Result, Context};
use plotters::prelude::*;
use std::path::Path;

/// Data for a single model's QQ plot
struct ModelQQData {
    name: String,
    expected: Vec<f64>,
    observed: Vec<f64>,
}

/// Prepare QQ data grouped by model
fn prepare_qq_data(results: &[GwasPoint]) -> Vec<ModelQQData> {
    // Collect unique models in order of first appearance
    let mut models: Vec<String> = Vec::new();
    for point in results {
        if !models.contains(&point.model) {
            models.push(point.model.clone());
        }
    }

    // Calculate QQ values for each model
    models.iter().map(|model_name| {
        let model_results: Vec<&GwasPoint> = results
            .iter()
            .filter(|p| &p.model == model_name)
            .collect();

        let (expected, observed) = calculate_qq_values(&model_results);

        ModelQQData {
            name: model_name.clone(),
            expected,
            observed,
        }
    }).collect()
}

/// Generate a QQ plot from GWAS results
///
/// The QQ plot compares the observed -log10(p-values) against the expected
/// values under the null hypothesis (uniform distribution).
/// When multiple models are present, each model is shown in a different color.
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

    // Prepare data grouped by model
    let model_data = prepare_qq_data(results);

    // Find max value across all models
    let max_val = model_data.iter()
        .flat_map(|m| m.expected.iter().chain(m.observed.iter()))
        .cloned()
        .fold(0.0_f64, f64::max);
    let axis_max = (max_val * 1.1).max(1.0);

    // Determine output format from extension
    let ext = output_path.extension()
        .and_then(|e| e.to_str())
        .unwrap_or("svg")
        .to_lowercase();

    match ext.as_str() {
        "svg" => draw_qq_svg(output_path, &model_data, &config, axis_max),
        #[cfg(feature = "png")]
        "png" => draw_qq_png(output_path, &model_data, &config, axis_max),
        _ => anyhow::bail!("Unsupported output format: {}", ext),
    }
}

/// Calculate expected and observed -log10(p) values for QQ plot
fn calculate_qq_values(results: &[&GwasPoint]) -> (Vec<f64>, Vec<f64>) {
    let n = results.len();

    if n == 0 {
        return (vec![], vec![]);
    }

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
    model_data: &[ModelQQData],
    config: &PlotConfig,
    axis_max: f64,
) -> Result<()> {
    let root = SVGBackend::new(output_path, (config.height, config.height)) // Square plot
        .into_drawing_area();

    draw_qq_impl(&root, model_data, config, axis_max)
        .context("Failed to draw QQ plot")?;

    root.present().context("Failed to write SVG")?;
    Ok(())
}

#[cfg(feature = "png")]
fn draw_qq_png(
    output_path: &Path,
    model_data: &[ModelQQData],
    config: &PlotConfig,
    axis_max: f64,
) -> Result<()> {
    let root = BitMapBackend::new(output_path, (config.height, config.height))
        .into_drawing_area();

    draw_qq_impl(&root, model_data, config, axis_max)
        .context("Failed to draw QQ plot")?;

    root.present().context("Failed to write PNG")?;
    Ok(())
}

fn draw_qq_impl<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    model_data: &[ModelQQData],
    config: &PlotConfig,
    axis_max: f64,
) -> Result<(), DrawingAreaErrorKind<DB::ErrorType>> {
    root.fill(&config.theme.background)?;

    let title = config.title.as_deref().unwrap_or("QQ Plot");
    let has_multiple_models = model_data.len() > 1;

    // Reserve space for legend on the right if multiple models
    let right_margin = if has_multiple_models { 120 } else { 10 };

    let mut chart = ChartBuilder::on(root)
        .caption(title, ("sans-serif", 24).into_font().color(&config.theme.text))
        .margin(10)
        .margin_right(right_margin)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(0.0..axis_max, 0.0..axis_max)?;

    chart
        .configure_mesh()
        .disable_mesh()
        .x_desc("Expected -log₁₀(p)")
        .y_desc("Observed -log₁₀(p)")
        .x_label_style(("sans-serif", 14).into_font().color(&config.theme.text))
        .y_label_style(("sans-serif", 14).into_font().color(&config.theme.text))
        .axis_desc_style(("sans-serif", 18).into_font().color(&config.theme.text))
        .axis_style(&config.theme.axis)
        .draw()?;

    // Draw confidence band (use first model's data for the band, or combined if preferred)
    // For simplicity, draw band based on total number of points
    let total_points: usize = model_data.iter().map(|m| m.expected.len()).sum();
    if total_points > 0 {
        // Use evenly spaced expected values for the confidence band
        let band_n = 100.min(total_points);
        let band_expected: Vec<f64> = (1..=band_n)
            .map(|i| -((i as f64) / ((band_n + 1) as f64)).log10())
            .collect();
        let (ci_lower, ci_upper) = calculate_confidence_band(band_n, &band_expected);

        for (exp, (lower, upper)) in band_expected.iter().zip(ci_lower.iter().zip(ci_upper.iter())) {
            let width = axis_max * 0.008;
            chart.draw_series(std::iter::once(Rectangle::new(
                [(*exp - width, *lower), (*exp + width, *upper)],
                config.theme.confidence_band.filled(),
            )))?;
        }
    }

    // Draw diagonal reference line (y = x)
    chart.draw_series(LineSeries::new(
        vec![(0.0, 0.0), (axis_max, axis_max)],
        config.theme.reference_line.stroke_width(2),
    ))?;

    // Draw points for each model with different colors
    let colors = &config.theme.chromosome_colors;
    let point_size = config.point_size;

    // Define a palette for models (extend if needed)
    let model_colors: Vec<RGBColor> = vec![
        RGBColor(31, 119, 180),   // Blue
        RGBColor(255, 127, 14),   // Orange
        RGBColor(44, 160, 44),    // Green
        RGBColor(214, 39, 40),    // Red
        RGBColor(148, 103, 189),  // Purple
        RGBColor(140, 86, 75),    // Brown
        RGBColor(227, 119, 194),  // Pink
        RGBColor(127, 127, 127),  // Gray
    ];

    for (model_idx, model) in model_data.iter().enumerate() {
        let color = if has_multiple_models {
            &model_colors[model_idx % model_colors.len()]
        } else {
            &colors[0]
        };

        chart.draw_series(
            model.expected.iter()
                .zip(model.observed.iter())
                .map(|(e, o)| Circle::new((*e, *o), point_size, color.filled()))
        )?;
    }

    // Draw legend if multiple models
    if has_multiple_models {
        let legend_x = axis_max * 1.05;
        let legend_y_start = axis_max * 0.95;
        let legend_spacing = axis_max * 0.08;

        for (idx, model) in model_data.iter().enumerate() {
            let y_pos = legend_y_start - (idx as f64 * legend_spacing);
            let color = &model_colors[idx % model_colors.len()];

            // Draw color sample
            chart.draw_series(std::iter::once(Circle::new(
                (legend_x, y_pos),
                point_size,
                color.filled(),
            )))?;

            // Draw model name
            chart.draw_series(std::iter::once(Text::new(
                model.name.clone(),
                (legend_x + axis_max * 0.03, y_pos),
                ("sans-serif", 10).into_font().color(&config.theme.text),
            )))?;
        }
    }

    Ok(())
}
