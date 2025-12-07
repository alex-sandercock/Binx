//! LD (Linkage Disequilibrium) decay plot
//!
//! Plots r² vs physical distance, showing the decay of LD with distance.
//! Implements functionality similar to R/GWASpoly's LD.plot function.

use anyhow::{Context, Result};
use ndarray::{Array1, Array2};
use plotters::prelude::*;
use rand::seq::SliceRandom;
use std::collections::HashMap;

use crate::PlotConfig;

/// Configuration specific to LD plots
#[derive(Clone, Debug)]
pub struct LDPlotConfig {
    /// Maximum number of r² pairs to use for plotting (default: 10000, matches R/GWASpoly)
    pub max_pairs: usize,
    /// Maximum number of markers to use per chromosome (None = all)
    pub max_loci_per_chrom: Option<usize>,
    /// Number of distance bins for smoothing (default: 50)
    pub n_bins: usize,
    /// Whether to show individual points (can be slow for many pairs)
    pub show_points: bool,
    /// R² threshold to mark on plot (e.g., 0.2) - draws vertical line at distance where r² drops to this value
    pub r2_threshold: Option<f64>,
    /// Filter to specific chromosomes (None = all)
    pub chromosomes: Option<Vec<String>>,
    /// Base plot config
    pub plot_config: PlotConfig,
}

impl Default for LDPlotConfig {
    fn default() -> Self {
        Self {
            max_pairs: 10_000, // Match R/GWASpoly default
            max_loci_per_chrom: None,
            n_bins: 100, // More bins for smoother curve
            show_points: false,
            r2_threshold: None,
            chromosomes: None,
            plot_config: PlotConfig::default(),
        }
    }
}

/// A single LD data point (distance, r²)
#[derive(Clone, Debug)]
pub struct LDPoint {
    pub distance_mb: f64,
    pub r2: f64,
}

/// Processed LD data for plotting
#[derive(Clone, Debug)]
pub struct LDData {
    /// Raw LD points (distance in Mb, r²)
    pub points: Vec<LDPoint>,
    /// Binned data: (bin_center_mb, mean_r2) after isotonic regression
    pub binned: Vec<(f64, f64)>,
    /// Smooth curve points for plotting (generated via cubic spline)
    pub smooth_curve: Vec<(f64, f64)>,
    /// Maximum distance in Mb
    pub max_distance: f64,
    /// Distance at which r² drops to threshold (if threshold specified)
    pub threshold_distance: Option<f64>,
}

/// Calculate LD decay data from genotype matrix
///
/// # Arguments
/// * `geno` - Genotype matrix with marker metadata
/// * `config` - LD plot configuration
///
/// # Returns
/// LDData containing raw points and binned averages
pub fn calculate_ld_data(
    geno: &binx_core::GenotypeMatrixBiallelic,
    config: &LDPlotConfig,
) -> Result<LDData> {
    let metadata = geno
        .marker_metadata
        .as_ref()
        .context("Genotype data must include marker metadata (chrom, pos)")?;

    // Group markers by chromosome
    let mut by_chrom: HashMap<String, Vec<usize>> = HashMap::new();
    for (i, meta) in metadata.iter().enumerate() {
        // Apply chromosome filter if specified
        if let Some(ref chroms) = config.chromosomes {
            if !chroms.contains(&meta.chrom) {
                continue;
            }
        }
        by_chrom.entry(meta.chrom.clone()).or_default().push(i);
    }

    if by_chrom.is_empty() {
        anyhow::bail!("No markers found for specified chromosomes");
    }

    let mut all_points: Vec<LDPoint> = Vec::new();
    let mut rng = rand::thread_rng();

    for (chrom, indices) in by_chrom.iter() {
        // Optionally subsample markers per chromosome
        let working_indices: Vec<usize> = if let Some(max_loci) = config.max_loci_per_chrom {
            if indices.len() > max_loci {
                let mut sampled = indices.clone();
                sampled.shuffle(&mut rng);
                sampled.truncate(max_loci);
                sampled.sort(); // Keep in position order
                sampled
            } else {
                indices.clone()
            }
        } else {
            indices.clone()
        };

        let m = working_indices.len();
        if m < 2 {
            continue;
        }

        eprintln!("  Processing {} with {} markers...", chrom, m);

        // Extract genotype submatrix for this chromosome
        let mut geno_sub = Array2::<f64>::zeros((m, geno.dosages.ncols()));
        for (new_i, &orig_i) in working_indices.iter().enumerate() {
            geno_sub.row_mut(new_i).assign(&geno.dosages.row(orig_i));
        }

        // Calculate pairwise r² and distances
        for i in 0..m {
            for j in (i + 1)..m {
                let orig_i = working_indices[i];
                let orig_j = working_indices[j];

                // Calculate distance in Mb
                let pos_i = metadata[orig_i].pos;
                let pos_j = metadata[orig_j].pos;
                let distance_mb = (pos_j - pos_i).abs() / 1_000_000.0;

                // Calculate r² (squared correlation)
                let r2 = calculate_r2(&geno_sub.row(i).to_owned(), &geno_sub.row(j).to_owned());

                if r2.is_finite() {
                    all_points.push(LDPoint { distance_mb, r2 });
                }
            }
        }
    }

    // Subsample if too many pairs
    if all_points.len() > config.max_pairs {
        all_points.shuffle(&mut rng);
        all_points.truncate(config.max_pairs);
    }

    // Calculate max distance
    let max_distance = all_points
        .iter()
        .map(|p| p.distance_mb)
        .fold(0.0_f64, |a, b| a.max(b));

    // Bin the data
    let binned = bin_ld_data(&all_points, config.n_bins, max_distance);

    // Generate smooth curve using monotone cubic spline (500 points for smooth rendering)
    let smooth_curve = generate_smooth_curve(&binned, 500);

    // Find distance where r² drops to threshold using smooth curve for accuracy
    let threshold_distance = config.r2_threshold.and_then(|thresh| {
        find_threshold_distance_smooth(&smooth_curve, thresh)
    });

    Ok(LDData {
        points: all_points,
        binned,
        smooth_curve,
        max_distance,
        threshold_distance,
    })
}

/// Calculate r² (squared Pearson correlation) between two genotype vectors
fn calculate_r2(x: &Array1<f64>, y: &Array1<f64>) -> f64 {
    let n = x.len();
    if n == 0 {
        return f64::NAN;
    }

    // Filter out pairs where either is NaN
    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut sum_xy = 0.0;
    let mut sum_x2 = 0.0;
    let mut sum_y2 = 0.0;
    let mut count = 0;

    for i in 0..n {
        let xi = x[i];
        let yi = y[i];
        if xi.is_finite() && yi.is_finite() {
            sum_x += xi;
            sum_y += yi;
            sum_xy += xi * yi;
            sum_x2 += xi * xi;
            sum_y2 += yi * yi;
            count += 1;
        }
    }

    if count < 2 {
        return f64::NAN;
    }

    let n_f = count as f64;
    let mean_x = sum_x / n_f;
    let mean_y = sum_y / n_f;

    let var_x = sum_x2 / n_f - mean_x * mean_x;
    let var_y = sum_y2 / n_f - mean_y * mean_y;
    let cov_xy = sum_xy / n_f - mean_x * mean_y;

    if var_x <= 0.0 || var_y <= 0.0 {
        return 0.0; // No variance = no correlation
    }

    let r = cov_xy / (var_x.sqrt() * var_y.sqrt());
    r * r // Return r²
}

/// Bin LD data by distance and calculate mean r² per bin
fn bin_ld_data(points: &[LDPoint], n_bins: usize, max_distance: f64) -> Vec<(f64, f64)> {
    if points.is_empty() || max_distance <= 0.0 {
        return Vec::new();
    }

    let bin_width = max_distance / n_bins as f64;
    let mut bin_sums: Vec<f64> = vec![0.0; n_bins];
    let mut bin_counts: Vec<usize> = vec![0; n_bins];

    for p in points {
        let bin_idx = ((p.distance_mb / bin_width) as usize).min(n_bins - 1);
        bin_sums[bin_idx] += p.r2;
        bin_counts[bin_idx] += 1;
    }

    let mut result = Vec::new();
    for i in 0..n_bins {
        if bin_counts[i] > 0 {
            let bin_center = (i as f64 + 0.5) * bin_width;
            let mean_r2 = bin_sums[i] / bin_counts[i] as f64;
            result.push((bin_center, mean_r2));
        }
    }

    // Apply isotonic regression for monotone decreasing (optional but recommended)
    isotonic_decreasing(&mut result);

    result
}

/// Apply isotonic regression to ensure monotone decreasing values
fn isotonic_decreasing(data: &mut [(f64, f64)]) {
    if data.len() <= 1 {
        return;
    }

    // Pool Adjacent Violators Algorithm (PAVA) for decreasing
    let n = data.len();
    let mut i = 0;

    while i < n - 1 {
        // If current value is less than next (violates decreasing), pool them
        if data[i].1 < data[i + 1].1 {
            // Find the block that needs to be pooled
            let mut j = i + 1;
            let mut sum = data[i].1;
            let mut count = 1;

            while j < n && (count == 1 || data[j].1 > sum / count as f64) {
                sum += data[j].1;
                count += 1;
                j += 1;
            }

            // Set all values in the block to the average
            let avg = sum / count as f64;
            for k in i..j {
                data[k].1 = avg;
            }

            // Go back to check if this created new violations
            if i > 0 {
                i -= 1;
            }
        } else {
            i += 1;
        }
    }
}

/// Generate smooth curve using monotone cubic spline interpolation (Fritsch-Carlson method)
/// This preserves monotonicity while creating a smooth curve through the data points
fn generate_smooth_curve(data: &[(f64, f64)], n_points: usize) -> Vec<(f64, f64)> {
    if data.len() < 2 {
        return data.to_vec();
    }

    let n = data.len();
    let x: Vec<f64> = data.iter().map(|(x, _)| *x).collect();
    let y: Vec<f64> = data.iter().map(|(_, y)| *y).collect();

    // Calculate slopes between consecutive points
    let mut delta: Vec<f64> = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        let dx = x[i + 1] - x[i];
        if dx.abs() < 1e-10 {
            delta.push(0.0);
        } else {
            delta.push((y[i + 1] - y[i]) / dx);
        }
    }

    // Initialize tangents using Fritsch-Carlson method
    let mut m: Vec<f64> = vec![0.0; n];

    // Endpoints
    m[0] = delta[0];
    m[n - 1] = delta[n - 2];

    // Interior points: average of adjacent secants, but only if same sign
    for i in 1..n - 1 {
        if delta[i - 1] * delta[i] <= 0.0 {
            // Different signs or zero - set tangent to zero for monotonicity
            m[i] = 0.0;
        } else {
            // Same sign - use harmonic mean for better monotonicity preservation
            m[i] = (delta[i - 1] + delta[i]) / 2.0;
        }
    }

    // Fritsch-Carlson modification to ensure monotonicity
    for i in 0..n - 1 {
        if delta[i].abs() < 1e-10 {
            // Flat segment
            m[i] = 0.0;
            m[i + 1] = 0.0;
        } else {
            let alpha = m[i] / delta[i];
            let beta = m[i + 1] / delta[i];

            // Check if we're in the monotonicity region
            let tau = alpha * alpha + beta * beta;
            if tau > 9.0 {
                // Restrict to monotonicity region
                let tau_sqrt = tau.sqrt();
                m[i] = 3.0 * delta[i] * alpha / tau_sqrt;
                m[i + 1] = 3.0 * delta[i] * beta / tau_sqrt;
            }
        }
    }

    // Generate smooth curve points
    let x_min = x[0];
    let x_max = x[n - 1];
    let mut result = Vec::with_capacity(n_points);

    for i in 0..n_points {
        let t = i as f64 / (n_points - 1) as f64;
        let xi = x_min + t * (x_max - x_min);

        // Find the interval containing xi
        let mut seg = 0;
        for j in 0..n - 1 {
            if xi >= x[j] && xi <= x[j + 1] {
                seg = j;
                break;
            }
            if j == n - 2 {
                seg = j; // Last segment
            }
        }

        // Evaluate cubic Hermite spline
        let h = x[seg + 1] - x[seg];
        if h.abs() < 1e-10 {
            result.push((xi, y[seg]));
            continue;
        }

        let t_local = (xi - x[seg]) / h;
        let t2 = t_local * t_local;
        let t3 = t2 * t_local;

        // Hermite basis functions
        let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
        let h10 = t3 - 2.0 * t2 + t_local;
        let h01 = -2.0 * t3 + 3.0 * t2;
        let h11 = t3 - t2;

        let yi = h00 * y[seg] + h10 * h * m[seg] + h01 * y[seg + 1] + h11 * h * m[seg + 1];

        // Clamp to valid range [0, max_y] and ensure monotonicity
        let yi_clamped = yi.max(0.0).min(y[0]);
        result.push((xi, yi_clamped));
    }

    // Final pass: enforce strict monotone decreasing
    for i in 1..result.len() {
        if result[i].1 > result[i - 1].1 {
            result[i].1 = result[i - 1].1;
        }
    }

    result
}

/// Find threshold distance using smooth curve (more accurate than binned data)
fn find_threshold_distance_smooth(curve: &[(f64, f64)], threshold: f64) -> Option<f64> {
    if curve.is_empty() {
        return None;
    }

    // Find first point where r² <= threshold
    for i in 0..curve.len() {
        if curve[i].1 <= threshold {
            if i == 0 {
                return Some(curve[0].0);
            }
            // Linear interpolate for precise distance
            let (x0, y0) = curve[i - 1];
            let (x1, y1) = curve[i];
            if (y0 - y1).abs() < 1e-10 {
                return Some(x0);
            }
            let x = x0 + (threshold - y0) * (x1 - x0) / (y1 - y0);
            return Some(x);
        }
    }

    // r² never drops to threshold
    None
}

/// Generate LD decay plot
pub fn ld_plot(
    geno: &binx_core::GenotypeMatrixBiallelic,
    output_path: &str,
    config: LDPlotConfig,
) -> Result<LDData> {
    eprintln!("Calculating LD data...");
    let ld_data = calculate_ld_data(geno, &config)?;
    eprintln!(
        "Calculated {} LD pairs, max distance: {:.2} Mb",
        ld_data.points.len(),
        ld_data.max_distance
    );

    // Report threshold distance if calculated
    if let (Some(thresh), Some(dist)) = (config.r2_threshold, ld_data.threshold_distance) {
        eprintln!("Distance at r² = {:.2}: {:.3} Mb", thresh, dist);
    } else if let Some(thresh) = config.r2_threshold {
        eprintln!("r² never drops to {:.2} within the data range", thresh);
    }

    if output_path.ends_with(".png") {
        #[cfg(feature = "png")]
        draw_ld_png(&ld_data, output_path, &config)?;
        #[cfg(not(feature = "png"))]
        anyhow::bail!("PNG output requires the 'png' feature");
    } else {
        draw_ld_svg(&ld_data, output_path, &config)?;
    }

    Ok(ld_data)
}

fn draw_ld_svg(ld_data: &LDData, path: &str, config: &LDPlotConfig) -> Result<()> {
    let root = SVGBackend::new(path, (config.plot_config.width, config.plot_config.height))
        .into_drawing_area();
    draw_ld_impl(&root, ld_data, config)?;
    root.present()?;
    Ok(())
}

#[cfg(feature = "png")]
fn draw_ld_png(ld_data: &LDData, path: &str, config: &LDPlotConfig) -> Result<()> {
    let root = BitMapBackend::new(path, (config.plot_config.width, config.plot_config.height))
        .into_drawing_area();
    draw_ld_impl(&root, ld_data, config)?;
    root.present()?;
    Ok(())
}

fn draw_ld_impl<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    ld_data: &LDData,
    config: &LDPlotConfig,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    let theme = &config.plot_config.theme;

    root.fill(&theme.background)?;

    let max_x = ld_data.max_distance * 1.05; // 5% padding

    // Auto-scale y-axis based on actual data (like GWASpoly)
    let data_max_r2 = ld_data
        .smooth_curve
        .iter()
        .map(|(_, y)| *y)
        .fold(0.0_f64, |a, b| a.max(b));
    let max_y = (data_max_r2 * 1.15).max(0.1); // 15% padding, minimum 0.1

    // Build chart - title is optional (default to None for cleaner look like GWASpoly)
    let mut chart_builder = ChartBuilder::on(root);
    chart_builder
        .margin(15)
        .x_label_area_size(50)
        .y_label_area_size(60);

    // Only add caption if title is explicitly provided
    if let Some(ref title) = config.plot_config.title {
        chart_builder.caption(title, ("sans-serif", 20).into_font().color(&theme.text));
    }

    let mut chart = chart_builder.build_cartesian_2d(0.0..max_x, 0.0..max_y)?;

    // Light gray gridlines for professional look (like GWASpoly/ggplot2)
    let grid_color = RGBColor(220, 220, 220);

    chart
        .configure_mesh()
        .x_desc("Distance (Mb)")
        .y_desc("r²")
        .axis_style(&theme.axis)
        .light_line_style(grid_color)
        .bold_line_style(grid_color.mix(0.8))
        .x_label_style(("sans-serif", 14).into_font().color(&theme.text))
        .y_label_style(("sans-serif", 14).into_font().color(&theme.text))
        .axis_desc_style(("sans-serif", 18).into_font().color(&theme.text))
        .draw()?;

    // Draw individual points if requested (semi-transparent)
    if config.show_points && !ld_data.points.is_empty() {
        let point_color = theme.chromosome_colors[0].mix(0.1);
        chart.draw_series(ld_data.points.iter().map(|p| {
            Circle::new((p.distance_mb, p.r2), 2, point_color.filled())
        }))?;
    }

    // Draw r² threshold line if specified
    if let (Some(thresh), Some(dist)) = (config.r2_threshold, ld_data.threshold_distance) {
        // Vertical line at threshold distance
        chart.draw_series(LineSeries::new(
            vec![(dist, 0.0), (dist, thresh)],
            theme.significance_line.stroke_width(2),
        ))?;

        // Horizontal line at threshold r²
        chart.draw_series(LineSeries::new(
            vec![(0.0, thresh), (dist, thresh)],
            theme.significance_line.stroke_width(2),
        ))?;

        // Annotation
        let label = format!("r²={:.2} @ {:.3} Mb", thresh, dist);
        chart.draw_series(std::iter::once(Text::new(
            label,
            (dist + max_x * 0.02, thresh + 0.05),
            ("sans-serif", 12).into_font().color(&theme.text),
        )))?;
    }

    // Draw the smooth spline curve - thicker line, no markers (like GWASpoly)
    if !ld_data.smooth_curve.is_empty() {
        let line_color = BLACK; // Clean black line like GWASpoly

        chart.draw_series(LineSeries::new(
            ld_data.smooth_curve.iter().map(|(x, y)| (*x, *y)),
            line_color.stroke_width(3), // Thicker line
        ))?;
    }

    Ok(())
}

/// Load genotype data and create LD plot
pub fn ld_plot_from_file(
    geno_path: &str,
    output_path: &str,
    ploidy: u8,
    config: LDPlotConfig,
) -> Result<LDData> {
    eprintln!("Loading genotype data from {}...", geno_path);
    let geno = binx_core::load_genotypes_biallelic_from_tsv(geno_path, ploidy)?;
    eprintln!(
        "Loaded {} markers x {} samples",
        geno.marker_ids.len(),
        geno.sample_ids.len()
    );

    if let Some(ref chroms) = config.chromosomes {
        eprintln!("Filtering to chromosomes: {:?}", chroms);
    }

    ld_plot(&geno, output_path, config)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_r2() {
        let x = Array1::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let y = Array1::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let r2 = calculate_r2(&x, &y);
        assert!((r2 - 1.0).abs() < 1e-10); // Perfect correlation

        let z = Array1::from_vec(vec![5.0, 4.0, 3.0, 2.0, 1.0]);
        let r2_neg = calculate_r2(&x, &z);
        assert!((r2_neg - 1.0).abs() < 1e-10); // Perfect negative correlation, r² = 1
    }

    #[test]
    fn test_isotonic_decreasing() {
        let mut data = vec![(0.5, 0.8), (1.5, 0.9), (2.5, 0.7), (3.5, 0.5)];
        isotonic_decreasing(&mut data);
        // After isotonic regression, should be monotone decreasing
        for i in 1..data.len() {
            assert!(data[i].1 <= data[i - 1].1);
        }
    }

    #[test]
    fn test_bin_ld_data() {
        let points = vec![
            LDPoint { distance_mb: 0.1, r2: 0.9 },
            LDPoint { distance_mb: 0.2, r2: 0.85 },
            LDPoint { distance_mb: 0.5, r2: 0.6 },
            LDPoint { distance_mb: 1.0, r2: 0.3 },
        ];
        let binned = bin_ld_data(&points, 5, 1.0);
        assert!(!binned.is_empty());
        // Should be monotone decreasing after isotonic regression
        for i in 1..binned.len() {
            assert!(binned[i].1 <= binned[i - 1].1);
        }
    }

    #[test]
    fn test_find_threshold_distance() {
        let curve = vec![
            (0.1, 0.9),
            (0.2, 0.7),
            (0.3, 0.5),
            (0.4, 0.3),
            (0.5, 0.1),
        ];

        // Find distance where r² = 0.4 (interpolate between 0.5 and 0.3)
        let dist = find_threshold_distance_smooth(&curve, 0.4);
        assert!(dist.is_some());
        let d = dist.unwrap();
        assert!(d > 0.3 && d < 0.4);
    }

    #[test]
    fn test_generate_smooth_curve() {
        let data = vec![
            (0.0, 1.0),
            (1.0, 0.5),
            (2.0, 0.2),
            (3.0, 0.1),
        ];
        let smooth = generate_smooth_curve(&data, 100);
        assert_eq!(smooth.len(), 100);
        // Should be monotone decreasing
        for i in 1..smooth.len() {
            assert!(smooth[i].1 <= smooth[i - 1].1);
        }
        // Should start and end near original values
        assert!((smooth[0].1 - 1.0).abs() < 0.01);
        assert!((smooth[99].1 - 0.1).abs() < 0.05);
    }
}
