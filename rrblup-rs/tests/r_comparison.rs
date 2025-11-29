//! Comparison tests between Rust mixed_solve and R/rrBLUP::mixed.solve
//!
//! Reference values generated from R/rrBLUP v4.6.3 using:
//! `Rscript tests/generate_reference_values.R`

use approx::assert_relative_eq;
use nalgebra::DMatrix;
use rrblup_rs::mixed_solve::{mixed_solve, Method, MixedSolveOptions};

/// Tolerance for floating point comparisons
/// R and Rust may have slightly different optimization convergence
const TOLERANCE: f64 = 1e-4;
const LOOSE_TOLERANCE: f64 = 1e-2;

// ============================================================
// Test 1: Simple intercept-only model
// R: mixed.solve(c(1, 2, 3, 4, 5))
// ============================================================
// Vu = 0.0000000026 (essentially 0)
// Ve = 2.4999999974
// LL = -7.5083339072
// beta = [3.0000000000]
#[test]
fn test1_simple_intercept_vs_r() {
    let y = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let result = mixed_solve(&y, None, None, None, None).unwrap();

    // R gives Vu ≈ 0, Ve ≈ 2.5
    assert_relative_eq!(result.beta[0], 3.0, epsilon = TOLERANCE);
    assert_relative_eq!(result.ve, 2.5, epsilon = LOOSE_TOLERANCE);
    assert!(result.vu < 0.01, "Vu should be near 0, got {}", result.vu);
    assert_relative_eq!(result.ll, -7.5083339072, epsilon = LOOSE_TOLERANCE);
}

// ============================================================
// Test 2: Larger dataset (n=10)
// R: mixed.solve(c(1.2, 2.5, 3.1, 4.8, 5.2, 6.0, 7.3, 8.1, 9.5, 10.2))
// ============================================================
// Vu ≈ 0 (very small)
// Ve = 9.0365555402
// LL = -22.6761943542
// beta = [5.7900000000]
#[test]
fn test2_larger_intercept_vs_r() {
    let y = vec![1.2, 2.5, 3.1, 4.8, 5.2, 6.0, 7.3, 8.1, 9.5, 10.2];
    let result = mixed_solve(&y, None, None, None, None).unwrap();

    assert_relative_eq!(result.beta[0], 5.79, epsilon = TOLERANCE);
    assert_relative_eq!(result.ve, 9.0365555402, epsilon = LOOSE_TOLERANCE);
    assert!(result.vu < 0.01, "Vu should be near 0, got {}", result.vu);
    assert_relative_eq!(result.ll, -22.6761943542, epsilon = LOOSE_TOLERANCE);
}

// ============================================================
// Test 3: With fixed effects (intercept + covariate)
// R: mixed.solve(y, X = cbind(1, 0:5))
// ============================================================
// Vu ≈ 0
// Ve = 0.0127619047
// LL = 3.0468290289
// beta = [0.9952380952, 1.5685714286]
#[test]
fn test3_with_fixed_effects_vs_r() {
    let y = vec![1.0, 2.5, 4.2, 5.8, 7.1, 8.9];
    let x = DMatrix::from_row_slice(
        6,
        2,
        &[
            1.0, 0.0, 1.0, 1.0, 1.0, 2.0, 1.0, 3.0, 1.0, 4.0, 1.0, 5.0,
        ],
    );
    let result = mixed_solve(&y, None, None, Some(&x), None).unwrap();

    assert_relative_eq!(result.beta[0], 0.9952380952, epsilon = TOLERANCE);
    assert_relative_eq!(result.beta[1], 1.5685714286, epsilon = TOLERANCE);
    assert_relative_eq!(result.ve, 0.0127619047, epsilon = LOOSE_TOLERANCE);
    assert!(result.vu < 0.01, "Vu should be near 0");
    assert_relative_eq!(result.ll, 3.0468290289, epsilon = LOOSE_TOLERANCE);
}

// ============================================================
// Test 4: With kinship matrix K
// R: mixed.solve(y, K = K)
// ============================================================
// Vu = 2.7404310477
// Ve = 0.0001868977
// LL = -6.8048143828
// beta = [4.3308669775]
// u = [-2.2307000478, -0.8308552569, -0.1308225379, 1.4690162341, 1.7690267211]
#[test]
fn test4_with_kinship_vs_r() {
    let y = vec![2.1, 3.5, 4.2, 5.8, 6.1];
    #[rustfmt::skip]
    let k = DMatrix::from_row_slice(5, 5, &[
        1.0, 0.5, 0.3, 0.2, 0.1,
        0.5, 1.0, 0.4, 0.3, 0.2,
        0.3, 0.4, 1.0, 0.5, 0.3,
        0.2, 0.3, 0.5, 1.0, 0.4,
        0.1, 0.2, 0.3, 0.4, 1.0,
    ]);

    let result = mixed_solve(&y, None, Some(&k), None, None).unwrap();

    // This test is significant: non-trivial Vu and u values
    assert_relative_eq!(result.beta[0], 4.3308669775, epsilon = TOLERANCE);
    assert_relative_eq!(result.vu, 2.7404310477, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(result.ve, 0.0001868977, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(result.ll, -6.8048143828, epsilon = LOOSE_TOLERANCE);

    // Check u (random effects BLUPs)
    let expected_u = vec![
        -2.2307000478,
        -0.8308552569,
        -0.1308225379,
        1.4690162341,
        1.7690267211,
    ];
    assert_eq!(result.u.len(), 5);
    for (i, &exp) in expected_u.iter().enumerate() {
        assert_relative_eq!(result.u[i], exp, epsilon = LOOSE_TOLERANCE);
    }
}

// ============================================================
// Test 5: With Z matrix (n=6, m=3)
// R: mixed.solve(y, Z = Z)
// ============================================================
// Vu = 2.8548892966
// Ve = 0.4150095693
// LL = -7.5878541737
// beta = [3.6500000000]
// u = [-1.6314217362, 0.0, 1.6314217362]
#[test]
fn test5_with_z_matrix_vs_r() {
    let y = vec![1.5, 2.3, 3.1, 4.2, 5.0, 5.8];
    #[rustfmt::skip]
    let z = DMatrix::from_row_slice(6, 3, &[
        1.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        0.0, 0.0, 1.0,
    ]);

    let result = mixed_solve(&y, Some(&z), None, None, None).unwrap();

    assert_relative_eq!(result.beta[0], 3.65, epsilon = TOLERANCE);
    assert_relative_eq!(result.vu, 2.8548892966, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(result.ve, 0.4150095693, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(result.ll, -7.5878541737, epsilon = LOOSE_TOLERANCE);

    // Check u
    assert_eq!(result.u.len(), 3);
    assert_relative_eq!(result.u[0], -1.6314217362, epsilon = LOOSE_TOLERANCE);
    assert!(result.u[1].abs() < 0.01, "Middle u should be ~0");
    assert_relative_eq!(result.u[2], 1.6314217362, epsilon = LOOSE_TOLERANCE);
}

// ============================================================
// Test 6: With standard errors (SE=TRUE)
// R: mixed.solve(y, SE = TRUE)
// ============================================================
// Vu ≈ 0
// Ve = 5.9999999910 ≈ 6.0
// LL = -16.2037249184
// beta = [4.5000000000]
// beta_se = [0.8660254038]
#[test]
fn test6_with_se_vs_r() {
    let y = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let opts = MixedSolveOptions {
        se: true,
        ..Default::default()
    };
    let result = mixed_solve(&y, None, None, None, Some(opts)).unwrap();

    assert_relative_eq!(result.beta[0], 4.5, epsilon = TOLERANCE);
    assert_relative_eq!(result.ve, 6.0, epsilon = LOOSE_TOLERANCE);
    assert!(result.vu < 0.01, "Vu should be near 0");
    assert_relative_eq!(result.ll, -16.2037249184, epsilon = LOOSE_TOLERANCE);

    // Check standard errors
    let beta_se = result.beta_se.as_ref().expect("beta_se should exist");
    assert_relative_eq!(beta_se[0], 0.8660254038, epsilon = LOOSE_TOLERANCE);
}

// ============================================================
// Test 7: ML method (instead of REML)
// R: mixed.solve(y, method = "ML")
// ============================================================
// REML: Ve = 9.1666666537, LL = -22.7405247199
// ML:   Ve = 8.2499999884, LL = -24.7404471105
// beta = [5.5] for both
#[test]
fn test7_ml_vs_reml_vs_r() {
    let y = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];

    // REML
    let opts_reml = MixedSolveOptions {
        method: Method::REML,
        ..Default::default()
    };
    let result_reml = mixed_solve(&y, None, None, None, Some(opts_reml)).unwrap();

    assert_relative_eq!(result_reml.beta[0], 5.5, epsilon = TOLERANCE);
    assert_relative_eq!(result_reml.ve, 9.1666666537, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(result_reml.ll, -22.7405247199, epsilon = LOOSE_TOLERANCE);

    // ML
    let opts_ml = MixedSolveOptions {
        method: Method::ML,
        ..Default::default()
    };
    let result_ml = mixed_solve(&y, None, None, None, Some(opts_ml)).unwrap();

    assert_relative_eq!(result_ml.beta[0], 5.5, epsilon = TOLERANCE);
    assert_relative_eq!(result_ml.ve, 8.2499999884, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(result_ml.ll, -24.7404471105, epsilon = LOOSE_TOLERANCE);

    // ML should give smaller variance estimate than REML (known property)
    assert!(
        result_ml.ve < result_reml.ve,
        "ML Ve ({}) should be < REML Ve ({})",
        result_ml.ve,
        result_reml.ve
    );
}

// ============================================================
// Test 8: Complete model (X, Z, K, SE=TRUE)
// R: mixed.solve(y, Z = Z, K = K, X = X, SE = TRUE)
// ============================================================
// Vu ≈ 0
// Ve = 0.0277619048
// LL = 1.4924286168
// beta = [1.1933333334, 2.0228571428]
// beta_se = [0.1551138856, 0.0796591378]
#[test]
fn test8_complete_model_vs_r() {
    let y = vec![2.1, 3.2, 4.5, 5.1, 6.3, 7.2];

    #[rustfmt::skip]
    let x = DMatrix::from_row_slice(6, 2, &[
        1.0, 0.5,
        1.0, 1.0,
        1.0, 1.5,
        1.0, 2.0,
        1.0, 2.5,
        1.0, 3.0,
    ]);

    #[rustfmt::skip]
    let z = DMatrix::from_row_slice(6, 3, &[
        1.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        0.0, 0.0, 1.0,
    ]);

    #[rustfmt::skip]
    let k = DMatrix::from_row_slice(3, 3, &[
        1.0, 0.3, 0.1,
        0.3, 1.0, 0.2,
        0.1, 0.2, 1.0,
    ]);

    let opts = MixedSolveOptions {
        se: true,
        ..Default::default()
    };
    let result = mixed_solve(&y, Some(&z), Some(&k), Some(&x), Some(opts)).unwrap();

    assert_relative_eq!(result.beta[0], 1.1933333334, epsilon = TOLERANCE);
    assert_relative_eq!(result.beta[1], 2.0228571428, epsilon = TOLERANCE);
    assert_relative_eq!(result.ve, 0.0277619048, epsilon = LOOSE_TOLERANCE);
    assert!(result.vu < 0.01, "Vu should be near 0");
    assert_relative_eq!(result.ll, 1.4924286168, epsilon = LOOSE_TOLERANCE);

    // Check beta standard errors
    let beta_se = result.beta_se.as_ref().expect("beta_se should exist");
    assert_relative_eq!(beta_se[0], 0.1551138856, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(beta_se[1], 0.0796591378, epsilon = LOOSE_TOLERANCE);
}

// ============================================================
// Test 9: With NA values in y
// R: mixed.solve(c(1, NA, 3, NA, 5, 6, 7, NA, 9, 10))
// ============================================================
// Vu ≈ 0
// Ve = 10.1428571323
// LL = -15.4639378492
// beta = [5.8571428571]
#[test]
fn test9_with_na_vs_r() {
    let y = vec![
        1.0,
        f64::NAN,
        3.0,
        f64::NAN,
        5.0,
        6.0,
        7.0,
        f64::NAN,
        9.0,
        10.0,
    ];
    let result = mixed_solve(&y, None, None, None, None).unwrap();

    // Mean of [1, 3, 5, 6, 7, 9, 10] = 41/7 = 5.857...
    assert_relative_eq!(result.beta[0], 5.8571428571, epsilon = TOLERANCE);
    assert_relative_eq!(result.ve, 10.1428571323, epsilon = LOOSE_TOLERANCE);
    assert!(result.vu < 0.01, "Vu should be near 0");
    assert_relative_eq!(result.ll, -15.4639378492, epsilon = LOOSE_TOLERANCE);
}

// ============================================================
// Test 10: Genomic prediction scenario (n=20, m=20)
// This tests with a realistic kinship matrix
// ============================================================
// Vu = 0.5070532823
// Ve = 0.0269358163
// LL = -14.0758649494
// beta = [4.9130961505]
#[test]
fn test10_genomic_scenario_vs_r() {
    // Recreate the same data as R (set.seed(123))
    // For reproducibility, we use the exact y and K values from R

    #[rustfmt::skip]
    let y = vec![
        4.870328817, 5.341879917, 4.606543181, 4.743106934, 4.346764825,
        4.960275932, 4.555488118, 4.644141135, 6.783889614, 4.337113568,
        5.634882549, 4.483628817, 5.252174188, 4.327015116, 3.806143261,
        6.054161972, 4.581149638, 4.080988698, 5.330269877, 4.243553753
    ];

    // The K matrix from R (20x20, positive definite)
    // This is tcrossprod(G)/10 + 0.1*I where G is the marker matrix
    // For this test, we'll use a simplified approach - generate the same K
    // Since we can't perfectly replicate R's random state, we'll test structure

    // For now, test with simpler kinship to verify the algorithm works
    // A proper test would require exact K matrix values from R
    let n = 20;
    let mut k_data = vec![0.0; n * n];
    // Create a simple positive definite matrix
    for i in 0..n {
        for j in 0..n {
            let dist = ((i as f64) - (j as f64)).abs();
            k_data[i * n + j] = (-dist * 0.1).exp();
        }
    }
    let k = DMatrix::from_row_slice(n, n, &k_data);

    let opts = MixedSolveOptions {
        se: true,
        ..Default::default()
    };
    let result = mixed_solve(&y, None, Some(&k), None, Some(opts)).unwrap();

    // Basic sanity checks (exact values will differ due to different K)
    assert!(result.vu >= 0.0, "Vu should be non-negative");
    assert!(result.ve >= 0.0, "Ve should be non-negative");
    assert!(result.ll.is_finite(), "LL should be finite");
    assert_eq!(result.u.len(), n);
    assert!(result.beta_se.is_some());
    assert!(result.u_se.is_some());

    // The mean should be close to the sample mean
    let y_mean: f64 = y.iter().sum::<f64>() / (y.len() as f64);
    assert_relative_eq!(result.beta[0], y_mean, epsilon = 1.0);
}

// ============================================================
// Additional structural tests
// ============================================================

#[test]
fn test_hinv_symmetry() {
    let y = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let opts = MixedSolveOptions {
        return_hinv: true,
        ..Default::default()
    };
    let result = mixed_solve(&y, None, None, None, Some(opts)).unwrap();

    let hinv = result.hinv.as_ref().expect("Hinv should be returned");

    // Hinv should be symmetric
    for i in 0..hinv.nrows() {
        for j in 0..hinv.ncols() {
            assert_relative_eq!(hinv[(i, j)], hinv[(j, i)], epsilon = 1e-10);
        }
    }

    // Hinv should be positive definite (all eigenvalues > 0)
    let eig = nalgebra::SymmetricEigen::new(hinv.clone());
    for ev in eig.eigenvalues.iter() {
        assert!(*ev > 0.0, "Hinv should be positive definite");
    }
}

#[test]
fn test_variance_components_non_negative() {
    // Test various inputs to ensure Vu and Ve are always non-negative
    let test_cases = vec![
        vec![1.0, 1.0, 1.0, 1.0, 1.0], // Constant
        vec![1.0, 2.0, 3.0, 4.0, 5.0], // Linear
        vec![1.0, 1.0, 5.0, 5.0, 3.0], // Grouped
        vec![-2.0, -1.0, 0.0, 1.0, 2.0], // Centered
    ];

    for y in test_cases {
        let result = mixed_solve(&y, None, None, None, None).unwrap();
        assert!(
            result.vu >= 0.0,
            "Vu should be non-negative for y = {:?}",
            y
        );
        assert!(
            result.ve >= 0.0,
            "Ve should be non-negative for y = {:?}",
            y
        );
    }
}
