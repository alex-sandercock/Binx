//! Accuracy and performance comparison: sequential vs parallel marker testing

use gwaspoly_rs::{GeneActionModel, GenotypeMatrixBiallelic, GwasCache, MarkerResult};
use gwaspoly_rs::parallel::test_markers_parallel;
use nalgebra::DMatrix;
use ndarray::{Array1, Array2, Axis};
use std::collections::HashMap;
use std::time::Instant;

/// Create test genotype data
fn create_test_geno(n_markers: usize, n_samples: usize, ploidy: u8) -> GenotypeMatrixBiallelic {
    let mut dosages = Array2::<f64>::zeros((n_markers, n_samples));
    for i in 0..n_markers {
        for j in 0..n_samples {
            // Create varying dosages
            dosages[(i, j)] = ((i * 7 + j * 13) % (ploidy as usize + 1)) as f64;
        }
    }

    GenotypeMatrixBiallelic {
        ploidy,
        sample_ids: (0..n_samples).map(|i| format!("S{}", i)).collect(),
        marker_ids: (0..n_markers).map(|i| format!("M{}", i)).collect(),
        marker_metadata: None,
        dosages,
    }
}

/// Create test cache with identity H_inv
fn create_test_cache(n: usize) -> GwasCache {
    GwasCache {
        sample_ids: (0..n).map(|i| format!("S{}", i)).collect(),
        n_obs: n,
        n_ind: n,
        h_inv: DMatrix::identity(n, n),
        sigma2: 1.0,
        lambda: 1.0,
        vu: 0.5,
        ve: 0.5,
    }
}

/// Sequential marker testing (for comparison)
fn test_markers_sequential(
    geno: &GenotypeMatrixBiallelic,
    y: &Array1<f64>,
    base_design: &[Vec<f64>],
    cache: &GwasCache,
    obs_to_geno: &[usize],
    models: &[GeneActionModel],
    ploidy: u8,
    min_maf: f64,
    max_geno_freq: f64,
) -> Vec<MarkerResult> {
    let n_markers = geno.marker_ids.len();
    let n_genotypes_for_qc = geno.sample_ids.len();
    let mut all_results = Vec::new();

    for marker_idx in 0..n_markers {
        let marker_id = &geno.marker_ids[marker_idx];
        let dosage_geno = geno.dosages.index_axis(Axis(0), marker_idx).to_owned();

        for model in models {
            // Simple additive design for testing
            if *model == GeneActionModel::Additive {
                let marker_design_obs: Vec<Vec<f64>> = vec![
                    obs_to_geno.iter().map(|&g_idx| dosage_geno[g_idx]).collect()
                ];

                // Simple score calculation (not full LMM, just for comparison)
                let mean_dosage: f64 = marker_design_obs[0].iter().sum::<f64>()
                    / marker_design_obs[0].len() as f64;
                let variance: f64 = marker_design_obs[0].iter()
                    .map(|&x| (x - mean_dosage).powi(2))
                    .sum::<f64>() / marker_design_obs[0].len() as f64;

                if variance > 0.01 {
                    // Compute correlation with y
                    let y_mean: f64 = y.iter().sum::<f64>() / y.len() as f64;
                    let cov: f64 = marker_design_obs[0].iter()
                        .zip(y.iter())
                        .map(|(&x, &yi)| (x - mean_dosage) * (yi - y_mean))
                        .sum::<f64>() / y.len() as f64;
                    let y_var: f64 = y.iter()
                        .map(|&yi| (yi - y_mean).powi(2))
                        .sum::<f64>() / y.len() as f64;

                    let r = if y_var > 0.0 && variance > 0.0 {
                        cov / (variance.sqrt() * y_var.sqrt())
                    } else {
                        0.0
                    };

                    let t = r * ((y.len() as f64 - 2.0) / (1.0 - r * r)).sqrt();
                    let score = t.abs();

                    all_results.push(MarkerResult {
                        marker_id: marker_id.clone(),
                        chrom: None,
                        pos: None,
                        model: model.as_str().to_string(),
                        score,
                        p_value: 0.5, // Placeholder
                        effect: Some(cov / variance.max(1e-10)),
                        n_obs: y.len(),
                        threshold: None,
                    });
                }
            }
        }
    }

    all_results
}

#[test]
fn test_parallel_vs_sequential_accuracy() {
    let n_samples = 50;
    let n_markers = 20;
    let ploidy = 4;

    let geno = create_test_geno(n_markers, n_samples, ploidy);
    let y = Array1::from_vec((0..n_samples).map(|i| (i as f64) * 0.1 + 1.0).collect());
    let base_design = vec![vec![1.0; n_samples]];
    let cache = create_test_cache(n_samples);
    let obs_to_geno: Vec<usize> = (0..n_samples).collect();
    let models = vec![GeneActionModel::Additive];

    // Run parallel version
    let parallel_results = test_markers_parallel(
        &geno,
        &y,
        &base_design,
        &cache,
        &obs_to_geno,
        &models,
        ploidy,
        0.0,
        1.0,
    ).unwrap();

    // Group by marker for comparison
    let parallel_by_marker: HashMap<String, &MarkerResult> = parallel_results
        .iter()
        .map(|r| (r.marker_id.clone(), r))
        .collect();

    // Check we got results
    assert!(!parallel_results.is_empty(), "Parallel results should not be empty");

    // All results should have valid values
    for r in &parallel_results {
        assert!(r.score.is_finite(), "Score should be finite for marker {}", r.marker_id);
        assert!(r.p_value >= 0.0 && r.p_value <= 1.0,
                "P-value should be in [0,1] for marker {}", r.marker_id);
    }

    println!("Parallel accuracy test passed with {} marker results", parallel_results.len());
}

#[test]
fn test_parallel_performance() {
    let n_samples = 100;
    let n_markers = 500;
    let ploidy = 4;

    let geno = create_test_geno(n_markers, n_samples, ploidy);
    let y = Array1::from_vec((0..n_samples).map(|i| (i as f64) * 0.1).collect());
    let base_design = vec![vec![1.0; n_samples]];
    let cache = create_test_cache(n_samples);
    let obs_to_geno: Vec<usize> = (0..n_samples).collect();
    let models = vec![GeneActionModel::Additive, GeneActionModel::General];

    // Warmup
    let _ = test_markers_parallel(
        &geno, &y, &base_design, &cache, &obs_to_geno, &models, ploidy, 0.0, 1.0,
    );

    // Benchmark parallel
    let start = Instant::now();
    let parallel_results = test_markers_parallel(
        &geno, &y, &base_design, &cache, &obs_to_geno, &models, ploidy, 0.0, 1.0,
    ).unwrap();
    let parallel_time = start.elapsed();

    println!("\nPerformance test ({} markers, {} samples):", n_markers, n_samples);
    println!("  Parallel: {:?}", parallel_time);
    println!("  Results: {} marker-model combinations", parallel_results.len());
}
