//! Manhattan plot generation for GWAS results

use crate::{GwasPoint, PlotConfig};
use anyhow::{Result, Context};
use plotters::prelude::*;
use std::collections::BTreeMap;
use std::path::Path;

/// Processed data ready for Manhattan plot rendering
struct ManhattanData {
    /// Points with cumulative x positions: (cumulative_pos, score, chrom_index)
    points: Vec<(f64, f64, usize)>,
    /// Chromosome boundaries: (name, start_pos, end_pos, mid_pos)
    chrom_info: Vec<(String, f64, f64, f64)>,
    /// Maximum cumulative position
    max_x: f64,
    /// Maximum score
    max_y: f64,
}

/// Prepare GWAS results for Manhattan plot
fn prepare_manhattan_data(results: &[GwasPoint]) -> ManhattanData {
    // Group points by chromosome
    let mut by_chrom: BTreeMap<String, Vec<&GwasPoint>> = BTreeMap::new();

    for point in results {
        let chrom = point.chrom.clone().unwrap_or_else(|| "Unknown".to_string());
        by_chrom.entry(chrom).or_default().push(point);
    }

    // Sort chromosomes naturally (1, 2, ..., 10, 11, ... X, Y)
    let mut chrom_order: Vec<String> = by_chrom.keys().cloned().collect();
    chrom_order.sort_by(|a, b| {
        let a_num: Option<u32> = a.trim_start_matches(|c: char| !c.is_numeric())
            .parse().ok();
        let b_num: Option<u32> = b.trim_start_matches(|c: char| !c.is_numeric())
            .parse().ok();
        match (a_num, b_num) {
            (Some(an), Some(bn)) => an.cmp(&bn),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => a.cmp(b),
        }
    });

    // Calculate cumulative positions
    let mut cumulative_offset = 0.0;
    let mut points = Vec::new();
    let mut chrom_info = Vec::new();
    let mut max_y = 0.0_f64;

    for (chrom_idx, chrom_name) in chrom_order.iter().enumerate() {
        let chrom_points = &by_chrom[chrom_name];

        // Sort by position within chromosome
        let mut sorted_points: Vec<_> = chrom_points.iter().collect();
        sorted_points.sort_by(|a, b| {
            a.pos.unwrap_or(0.0).partial_cmp(&b.pos.unwrap_or(0.0)).unwrap()
        });

        let chrom_start = cumulative_offset;
        let mut chrom_max_pos = 0.0_f64;

        for point in sorted_points {
            let pos = point.pos.unwrap_or(0.0);
            let cum_pos = cumulative_offset + pos;
            chrom_max_pos = chrom_max_pos.max(pos);
            max_y = max_y.max(point.score);
            points.push((cum_pos, point.score, chrom_idx));
        }

        let chrom_end = cumulative_offset + chrom_max_pos;
        let chrom_mid = (chrom_start + chrom_end) / 2.0;
        chrom_info.push((chrom_name.clone(), chrom_start, chrom_end, chrom_mid));

        // Add gap between chromosomes
        cumulative_offset = chrom_end + chrom_max_pos * 0.02;
    }

    let max_x = cumulative_offset;

    ManhattanData {
        points,
        chrom_info,
        max_x,
        max_y,
    }
}

/// Generate a Manhattan plot from GWAS results
///
/// # Arguments
/// * `results` - GWAS results to plot
/// * `output_path` - Path for output file (SVG or PNG based on extension)
/// * `config` - Plot configuration
///
/// # Example
/// ```ignore
/// use binx_plotting::{manhattan_plot, PlotConfig, load_gwas_results};
///
/// let results = load_gwas_results("gwas_results.csv")?;
/// manhattan_plot(&results, "manhattan.svg", PlotConfig::default())?;
/// ```
pub fn manhattan_plot<P: AsRef<Path>>(
    results: &[GwasPoint],
    output_path: P,
    config: PlotConfig,
) -> Result<()> {
    let output_path = output_path.as_ref();

    if results.is_empty() {
        anyhow::bail!("No GWAS results to plot");
    }

    let data = prepare_manhattan_data(results);

    // Add padding to y-axis
    let y_max = (data.max_y * 1.1).max(config.significance_threshold + 1.0);

    // Determine output format from extension
    let ext = output_path.extension()
        .and_then(|e| e.to_str())
        .unwrap_or("svg")
        .to_lowercase();

    match ext.as_str() {
        "svg" => draw_manhattan_svg(output_path, &data, &config, y_max),
        #[cfg(feature = "png")]
        "png" => draw_manhattan_png(output_path, &data, &config, y_max),
        _ => anyhow::bail!("Unsupported output format: {}", ext),
    }
}

fn draw_manhattan_svg(
    output_path: &Path,
    data: &ManhattanData,
    config: &PlotConfig,
    y_max: f64,
) -> Result<()> {
    let root = SVGBackend::new(output_path, (config.width, config.height))
        .into_drawing_area();

    draw_manhattan_impl(&root, data, config, y_max)
        .context("Failed to draw Manhattan plot")?;

    root.present().context("Failed to write SVG")?;
    Ok(())
}

#[cfg(feature = "png")]
fn draw_manhattan_png(
    output_path: &Path,
    data: &ManhattanData,
    config: &PlotConfig,
    y_max: f64,
) -> Result<()> {
    let root = BitMapBackend::new(output_path, (config.width, config.height))
        .into_drawing_area();

    draw_manhattan_impl(&root, data, config, y_max)
        .context("Failed to draw Manhattan plot")?;

    root.present().context("Failed to write PNG")?;
    Ok(())
}

fn draw_manhattan_impl<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    data: &ManhattanData,
    config: &PlotConfig,
    y_max: f64,
) -> Result<(), DrawingAreaErrorKind<DB::ErrorType>> {
    root.fill(&config.theme.background)?;

    let title = config.title.as_deref().unwrap_or("Manhattan Plot");

    let mut chart = ChartBuilder::on(root)
        .caption(title, ("sans-serif", 24).into_font().color(&config.theme.text))
        .margin(10)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(0.0..data.max_x, 0.0..y_max)?;

    // Configure mesh
    chart
        .configure_mesh()
        .disable_x_mesh()
        .y_desc("-log₁₀(p)")
        .y_label_style(("sans-serif", 14).into_font().color(&config.theme.text))
        .x_label_style(("sans-serif", 12).into_font().color(&config.theme.text))
        .axis_style(&config.theme.axis)
        .draw()?;

    // Draw significance threshold line
    chart.draw_series(LineSeries::new(
        vec![(0.0, config.significance_threshold), (data.max_x, config.significance_threshold)],
        config.theme.significance_line.stroke_width(2),
    ))?;

    // Draw suggestive threshold line if configured
    if let Some(suggestive) = config.suggestive_threshold {
        chart.draw_series(LineSeries::new(
            vec![(0.0, suggestive), (data.max_x, suggestive)],
            config.theme.suggestive_line.stroke_width(1),
        ))?;
    }

    // Draw points colored by chromosome
    let colors = &config.theme.chromosome_colors;
    let point_size = config.point_size;

    for (cum_pos, score, chrom_idx) in &data.points {
        let color = &colors[*chrom_idx % colors.len()];
        chart.draw_series(std::iter::once(Circle::new(
            (*cum_pos, *score),
            point_size,
            color.filled(),
        )))?;
    }

    // Draw chromosome labels on x-axis
    if config.show_chrom_labels {
        for (chrom_name, _start, _end, mid) in &data.chrom_info {
            chart.draw_series(std::iter::once(Text::new(
                chrom_name.clone(),
                (*mid, -y_max * 0.02),
                ("sans-serif", 11).into_font().color(&config.theme.text),
            )))?;
        }
    }

    Ok(())
}
