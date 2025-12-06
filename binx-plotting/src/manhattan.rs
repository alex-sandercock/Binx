//! Manhattan plot generation for GWAS results

use crate::{GwasPoint, PlotConfig};
use anyhow::{Result, Context};
use plotters::prelude::*;
use std::collections::BTreeMap;
use std::path::Path;

/// Processed data ready for Manhattan plot rendering
struct ManhattanData {
    /// Points with cumulative x positions: (cumulative_pos, score, chrom_index, model_index)
    points: Vec<(f64, f64, usize, usize)>,
    /// Chromosome boundaries: (name, start_pos, end_pos, mid_pos)
    chrom_info: Vec<ChromInfo>,
    /// Unique model names in order
    models: Vec<String>,
    /// Maximum cumulative position
    max_x: f64,
    /// Maximum score
    max_y: f64,
}

/// Information about a chromosome for plotting
struct ChromInfo {
    name: String,
    /// Midpoint for label placement (cumulative position)
    mid: f64,
}

/// Prepare GWAS results for Manhattan plot
fn prepare_manhattan_data(results: &[GwasPoint], chrom_filter: Option<&[String]>) -> ManhattanData {
    // Collect unique models in order of first appearance
    let mut models: Vec<String> = Vec::new();
    for point in results {
        if !models.contains(&point.model) {
            models.push(point.model.clone());
        }
    }

    // Group points by chromosome
    let mut by_chrom: BTreeMap<String, Vec<&GwasPoint>> = BTreeMap::new();

    for point in results {
        let chrom = point.chrom.clone().unwrap_or_else(|| "Unknown".to_string());

        // Skip if chromosome filter is set and this chromosome is not in the list
        if let Some(filter) = chrom_filter {
            if !filter.contains(&chrom) {
                continue;
            }
        }

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
            let model_idx = models.iter().position(|m| m == &point.model).unwrap_or(0);
            points.push((cum_pos, point.score, chrom_idx, model_idx));
        }

        let chrom_end = cumulative_offset + chrom_max_pos;
        let chrom_mid = (chrom_start + chrom_end) / 2.0;
        chrom_info.push(ChromInfo {
            name: chrom_name.clone(),
            mid: chrom_mid,
        });

        // Add gap between chromosomes
        cumulative_offset = chrom_end + chrom_max_pos * 0.02;
    }

    let max_x = cumulative_offset;

    ManhattanData {
        points,
        chrom_info,
        models,
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

    let data = prepare_manhattan_data(results, config.chromosomes.as_deref());

    if data.points.is_empty() {
        anyhow::bail!("No GWAS results to plot after chromosome filtering");
    }

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

/// Shape types for different models
#[derive(Clone, Copy)]
enum PointShape {
    CircleFilled,
    TriangleFilled,
    CircleOutline,
    TriangleOutline,
    Cross,
    Plus,
}

impl PointShape {
    fn from_index(idx: usize) -> Self {
        match idx % 6 {
            0 => PointShape::CircleFilled,
            1 => PointShape::TriangleFilled,
            2 => PointShape::CircleOutline,
            3 => PointShape::Cross,
            4 => PointShape::TriangleOutline,
            _ => PointShape::Plus,
        }
    }
}

fn draw_manhattan_impl<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    data: &ManhattanData,
    config: &PlotConfig,
    y_max: f64,
) -> Result<(), DrawingAreaErrorKind<DB::ErrorType>> {
    root.fill(&config.theme.background)?;

    let title = config.title.as_deref().unwrap_or("Manhattan Plot");
    let has_multiple_models = data.models.len() > 1;

    // Reserve space for legend on the right if multiple models
    let right_margin = if has_multiple_models { 120 } else { 10 };

    let mut chart = ChartBuilder::on(root)
        .caption(title, ("sans-serif", 24).into_font().color(&config.theme.text))
        .margin(10)
        .margin_right(right_margin)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(0.0..data.max_x, 0.0..y_max)?;

    // Configure mesh - disable grid lines and default x-axis
    chart
        .configure_mesh()
        .disable_mesh()
        .disable_x_axis()
        .y_desc("-log₁₀(p)")
        .y_label_style(("sans-serif", 14).into_font().color(&config.theme.text))
        .axis_desc_style(("sans-serif", 18).into_font().color(&config.theme.text))
        .axis_style(&config.theme.axis)
        .draw()?;

    // Draw x-axis line
    chart.draw_series(LineSeries::new(
        vec![(0.0, 0.0), (data.max_x, 0.0)],
        config.theme.axis.stroke_width(1),
    ))?;

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

    // Draw points colored by chromosome, shaped by model
    let colors = &config.theme.chromosome_colors;
    let point_size = config.point_size as i32;

    for (cum_pos, score, chrom_idx, model_idx) in &data.points {
        let color = &colors[*chrom_idx % colors.len()];
        let shape = PointShape::from_index(*model_idx);

        draw_point(&mut chart, *cum_pos, *score, point_size, color, shape)?;
    }

    // Draw chromosome labels centered below x-axis
    if config.show_chrom_labels {
        for chrom in &data.chrom_info {
            chart.draw_series(std::iter::once(Text::new(
                chrom.name.clone(),
                (chrom.mid, -y_max * 0.20),
                ("sans-serif", 14)
                    .into_font()
                    .color(&config.theme.text),
            )))?;
        }
    }

    // Draw legend if multiple models
    if has_multiple_models {
        let legend_x = data.max_x * 1.02;
        let legend_y_start = y_max * 0.95;
        let legend_spacing = y_max * 0.08;

        for (idx, model_name) in data.models.iter().enumerate() {
            let y_pos = legend_y_start - (idx as f64 * legend_spacing);
            let shape = PointShape::from_index(idx);

            // Draw shape sample (use first chromosome color for legend)
            let sample_color = &colors[0];
            draw_point(&mut chart, legend_x, y_pos, point_size, sample_color, shape)?;

            // Draw model name
            chart.draw_series(std::iter::once(Text::new(
                model_name.clone(),
                (legend_x + data.max_x * 0.02, y_pos),
                ("sans-serif", 10).into_font().color(&config.theme.text),
            )))?;
        }
    }

    Ok(())
}

/// Draw a point with the specified shape
fn draw_point<DB: DrawingBackend>(
    chart: &mut ChartContext<DB, Cartesian2d<plotters::coord::types::RangedCoordf64, plotters::coord::types::RangedCoordf64>>,
    x: f64,
    y: f64,
    size: i32,
    color: &RGBColor,
    shape: PointShape,
) -> Result<(), DrawingAreaErrorKind<DB::ErrorType>> {
    match shape {
        PointShape::CircleFilled => {
            chart.draw_series(std::iter::once(Circle::new((x, y), size as u32, color.filled())))?;
        }
        PointShape::TriangleFilled => {
            chart.draw_series(std::iter::once(TriangleMarker::new((x, y), size, color.filled())))?;
        }
        PointShape::CircleOutline => {
            chart.draw_series(std::iter::once(Circle::new((x, y), size as u32, color.stroke_width(2))))?;
        }
        PointShape::TriangleOutline => {
            chart.draw_series(std::iter::once(TriangleMarker::new((x, y), size, color.stroke_width(2))))?;
        }
        PointShape::Cross => {
            chart.draw_series(std::iter::once(Cross::new((x, y), size, color.stroke_width(2))))?;
        }
        PointShape::Plus => {
            // Draw a plus sign using two lines would require different approach
            // Use a smaller cross rotated - actually Cross is already an X, use it differently
            chart.draw_series(std::iter::once(Cross::new((x, y), size + 1, color.stroke_width(1))))?;
        }
    }
    Ok(())
}
