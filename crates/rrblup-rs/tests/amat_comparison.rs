//! Comparison tests between Rust a_mat and R/rrBLUP::A.mat
//!
//! Reference values generated from R/rrBLUP v4.6.3 using:
//! `Rscript tests/generate_amat_reference.R`

use approx::assert_relative_eq;
use nalgebra::DMatrix;
use rrblup_rs::a_mat::{a_mat, AMatOptions, ImputeMethod, ShrinkConfig, ShrinkMethod};

/// Tolerance for floating point comparisons
const TOLERANCE: f64 = 1e-6;

// ============================================================
// Test 1: Simple case (3 individuals x 5 markers)
// ============================================================
// Input X:
//      [,1] [,2] [,3] [,4] [,5]
// [1,]   -1    0    1    0   -1
// [2,]    0    1   -1    1    0
// [3,]    1   -1    0   -1    1
//
// Output A:
// [1.2, -0.4, -0.8]
// [-0.4, 1.2, -0.8]
// [-0.8, -0.8, 1.6]
#[test]
fn test1_simple_case_vs_r() {
    #[rustfmt::skip]
    let x = DMatrix::from_row_slice(3, 5, &[
        -1.0,  0.0,  1.0,  0.0, -1.0,
         0.0,  1.0, -1.0,  1.0,  0.0,
         1.0, -1.0,  0.0, -1.0,  1.0,
    ]);

    let result = a_mat(&x, None).unwrap();

    #[rustfmt::skip]
    let expected = [
        1.2, -0.4, -0.8,
       -0.4,  1.2, -0.8,
       -0.8, -0.8,  1.6,
    ];

    assert_eq!(result.a.nrows(), 3);
    assert_eq!(result.a.ncols(), 3);

    for i in 0..3 {
        for j in 0..3 {
            assert_relative_eq!(
                result.a[(i, j)],
                expected[i * 3 + j],
                epsilon = TOLERANCE
            );
        }
    }
}

// ============================================================
// Test 2: Larger case (5 individuals x 10 markers)
// ============================================================
// Input X (from R with set.seed(123)):
//      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
// [1,]    1    0    0   -1   -1    0   -1    0    0    -1
// [2,]    1    0    0    1   -1   -1    1   -1    1     1
// [3,]    1    0   -1    1    1    0    1    1    1     0
// [4,]    0    1    0   -1    0    1   -1   -1   -1    -1
// [5,]    1   -1    1   -1    1    0    1   -1    1     0
#[test]
fn test2_larger_case_vs_r() {
    #[rustfmt::skip]
    let x = DMatrix::from_row_slice(5, 10, &[
        1.0,  0.0,  0.0, -1.0, -1.0,  0.0, -1.0,  0.0,  0.0, -1.0,
        1.0,  0.0,  0.0,  1.0, -1.0, -1.0,  1.0, -1.0,  1.0,  1.0,
        1.0,  0.0, -1.0,  1.0,  1.0,  0.0,  1.0,  1.0,  1.0,  0.0,
        0.0,  1.0,  0.0, -1.0,  0.0,  1.0, -1.0, -1.0, -1.0, -1.0,
        1.0, -1.0,  1.0, -1.0,  1.0,  0.0,  1.0, -1.0,  1.0,  0.0,
    ]);

    let result = a_mat(&x, None).unwrap();

    // Expected from R
    #[rustfmt::skip]
    let expected = [
        0.9439252336, -0.5514018692, -0.6448598131,  0.7102803738, -0.4579439252,
       -0.5514018692,  1.4579439252,  0.1962616822, -1.0186915888, -0.0841121495,
       -0.6448598131,  0.1962616822,  1.5046728972, -0.8785046729, -0.1775700935,
        0.7102803738, -1.0186915888, -0.8785046729,  1.6448598131, -0.4579439252,
       -0.4579439252, -0.0841121495, -0.1775700935, -0.4579439252,  1.1775700935,
    ];

    assert_eq!(result.a.nrows(), 5);
    assert_eq!(result.a.ncols(), 5);

    for i in 0..5 {
        for j in 0..5 {
            assert_relative_eq!(
                result.a[(i, j)],
                expected[i * 5 + j],
                epsilon = TOLERANCE
            );
        }
    }
}

// ============================================================
// Test 3: With MAF filtering (min.MAF = 0.2)
// ============================================================
// Some markers should be filtered due to low MAF
#[test]
fn test3_with_maf_filtering_vs_r() {
    #[rustfmt::skip]
    let x = DMatrix::from_row_slice(4, 5, &[
        -1.0, -1.0, -1.0,  0.0,  1.0,  // marker 1 has low MAF
         0.0,  1.0, -1.0,  1.0,  0.0,
         1.0, -1.0,  0.0, -1.0,  1.0,
         0.0,  0.0,  1.0,  0.0, -1.0,
    ]);

    let opts = AMatOptions {
        min_maf: Some(0.2),
        ..Default::default()
    };

    let result = a_mat(&x, Some(opts)).unwrap();

    // Expected from R
    #[rustfmt::skip]
    let expected = [
        1.1168831169, -0.2337662338, -0.0259740260, -0.8571428571,
       -0.2337662338,  1.3246753247, -0.9610389610, -0.1298701299,
       -0.0259740260, -0.9610389610,  1.3246753247, -0.3376623377,
       -0.8571428571, -0.1298701299, -0.3376623377,  1.3246753247,
    ];

    assert_eq!(result.a.nrows(), 4);
    assert_eq!(result.a.ncols(), 4);

    for i in 0..4 {
        for j in 0..4 {
            assert_relative_eq!(
                result.a[(i, j)],
                expected[i * 4 + j],
                epsilon = TOLERANCE
            );
        }
    }
}

// ============================================================
// Test 4: With missing values (mean imputation)
// ============================================================
// Input X (with NA):
//      [,1] [,2] [,3] [,4]
// [1,]   -1    0    1    0
// [2,]    0   NA   -1    1
// [3,]    1   -1   NA   -1
// [4,]    0    1    0   NA
//
// Output A:
// [1.0, -0.5, -0.5, 0.0]
// [-0.5, 1.0, -0.5, 0.0]
// [-0.5, -0.5, 1.5, -0.5]
// [0.0, 0.0, -0.5, 0.5]
#[test]
fn test4_with_missing_mean_imputation_vs_r() {
    #[rustfmt::skip]
    let x = DMatrix::from_row_slice(4, 4, &[
        -1.0,      0.0,      1.0,  0.0,
         0.0,  f64::NAN,    -1.0,  1.0,
         1.0,     -1.0,  f64::NAN, -1.0,
         0.0,      1.0,      0.0, f64::NAN,
    ]);

    let opts = AMatOptions {
        impute_method: ImputeMethod::Mean,
        ..Default::default()
    };

    let result = a_mat(&x, Some(opts)).unwrap();

    #[rustfmt::skip]
    let expected = [
        1.0, -0.5, -0.5,  0.0,
       -0.5,  1.0, -0.5,  0.0,
       -0.5, -0.5,  1.5, -0.5,
        0.0,  0.0, -0.5,  0.5,
    ];

    assert_eq!(result.a.nrows(), 4);
    assert_eq!(result.a.ncols(), 4);

    for i in 0..4 {
        for j in 0..4 {
            assert_relative_eq!(
                result.a[(i, j)],
                expected[i * 4 + j],
                epsilon = TOLERANCE
            );
        }
    }
}

// ============================================================
// Test 5: With shrinkage (EJ method)
// ============================================================
// Note: Shrinkage involves randomness, so we check properties rather than exact values
#[test]
fn test5_with_shrinkage_ej_vs_r() {
    #[rustfmt::skip]
    let x = DMatrix::from_row_slice(8, 10, &[
        -1.0,  1.0,  0.0,  0.0, -1.0,  1.0,  0.0,  0.0, -1.0,  1.0,
        -1.0,  0.0,  0.0, -1.0,  0.0,  0.0,  1.0,  1.0,  0.0, -1.0,
         1.0, -1.0, -1.0, -1.0, -1.0,  1.0,  1.0,  1.0,  1.0,  0.0,
         0.0,  1.0, -1.0,  0.0,  1.0,  1.0,  0.0, -1.0,  1.0,  1.0,
        -1.0,  0.0,  0.0,  0.0,  0.0, -1.0,  0.0,  1.0, -1.0, -1.0,
         1.0,  1.0,  1.0,  1.0,  1.0, -1.0, -1.0, -1.0, -1.0,  1.0,
        -1.0,  1.0, -1.0,  0.0,  0.0, -1.0,  1.0,  0.0,  0.0,  0.0,
         0.0,  1.0, -1.0,  1.0,  0.0,  0.0,  1.0, -1.0,  0.0, -1.0,
    ]);

    let opts = AMatOptions {
        shrink: Some(ShrinkConfig {
            method: ShrinkMethod::EJ,
            ..Default::default()
        }),
        ..Default::default()
    };

    let result = a_mat(&x, Some(opts)).unwrap();

    // Check basic properties
    assert_eq!(result.a.nrows(), 8);
    assert_eq!(result.a.ncols(), 8);

    // Should have shrinkage intensity
    assert!(result.shrink_intensity.is_some());
    let delta = result.shrink_intensity.unwrap();
    assert!(delta >= 0.0 && delta <= 1.0);

    // Matrix should be symmetric
    for i in 0..8 {
        for j in 0..8 {
            assert_relative_eq!(result.a[(i, j)], result.a[(j, i)], epsilon = 1e-10);
        }
    }

    // R gives shrinkage intensity ≈ 0.62
    // Allow some tolerance due to numerical differences
    assert!(delta > 0.4 && delta < 0.8, "Shrinkage intensity {} not in expected range", delta);
}

// ============================================================
// Test 6: Return imputed matrix
// ============================================================
#[test]
fn test6_return_imputed_vs_r() {
    #[rustfmt::skip]
    let x = DMatrix::from_row_slice(3, 3, &[
        -1.0,      0.0,  1.0,
         0.0,  f64::NAN, -1.0,
         1.0,     -1.0,  0.0,
    ]);

    let opts = AMatOptions {
        impute_method: ImputeMethod::Mean,
        return_imputed: true,
        ..Default::default()
    };

    let result = a_mat(&x, Some(opts)).unwrap();

    // Expected A from R
    #[rustfmt::skip]
    let expected_a = [
        1.6363636364, -0.7272727273, -0.9090909091,
       -0.7272727273,  0.7272727273,  0.0000000000,
       -0.9090909091,  0.0000000000,  0.9090909091,
    ];

    for i in 0..3 {
        for j in 0..3 {
            assert_relative_eq!(
                result.a[(i, j)],
                expected_a[i * 3 + j],
                epsilon = TOLERANCE
            );
        }
    }

    // Check imputed matrix exists
    assert!(result.imputed.is_some());
    let imputed = result.imputed.unwrap();

    // Expected imputed X from R:
    // [-1.0,  0.0,  1.0]
    // [ 0.0, -0.5, -1.0]
    // [ 1.0, -1.0,  0.0]
    #[rustfmt::skip]
    let expected_imputed = [
        -1.0,  0.0,  1.0,
         0.0, -0.5, -1.0,
         1.0, -1.0,  0.0,
    ];

    for i in 0..3 {
        for j in 0..3 {
            assert_relative_eq!(
                imputed[(i, j)],
                expected_imputed[i * 3 + j],
                epsilon = TOLERANCE
            );
        }
    }
}

// ============================================================
// Test 7: All heterozygotes (edge case)
// ============================================================
// When all genotypes are 0 (heterozygous), A should be all zeros
// because there's no genetic variance
#[test]
fn test7_all_heterozygotes_vs_r() {
    let x = DMatrix::from_element(3, 4, 0.0); // All zeros

    let result = a_mat(&x, None).unwrap();

    // R gives all zeros
    for i in 0..3 {
        for j in 0..3 {
            assert_relative_eq!(result.a[(i, j)], 0.0, epsilon = TOLERANCE);
        }
    }
}

// ============================================================
// Test 8: Symmetry check
// ============================================================
#[test]
fn test8_a_matrix_symmetry() {
    #[rustfmt::skip]
    let x = DMatrix::from_row_slice(5, 8, &[
        -1.0,  0.0,  1.0,  0.0, -1.0,  1.0,  0.0, -1.0,
         0.0,  1.0, -1.0,  1.0,  0.0, -1.0,  1.0,  0.0,
         1.0, -1.0,  0.0, -1.0,  1.0,  0.0, -1.0,  1.0,
        -1.0,  1.0,  1.0,  0.0,  0.0,  1.0, -1.0,  0.0,
         0.0,  0.0, -1.0,  1.0, -1.0, -1.0,  0.0,  1.0,
    ]);

    let result = a_mat(&x, None).unwrap();

    // A should be symmetric
    for i in 0..5 {
        for j in 0..5 {
            assert_relative_eq!(result.a[(i, j)], result.a[(j, i)], epsilon = 1e-10);
        }
    }
}

// ============================================================
// Test 9: Diagonal values check
// ============================================================
// Diagonal of A represents individual inbreeding + 1
// Should be positive for realistic data
#[test]
fn test9_diagonal_positive() {
    #[rustfmt::skip]
    let x = DMatrix::from_row_slice(4, 6, &[
        -1.0,  0.0,  1.0,  0.0, -1.0,  1.0,
         0.0,  1.0, -1.0,  1.0,  0.0, -1.0,
         1.0, -1.0,  0.0, -1.0,  1.0,  0.0,
        -1.0,  1.0,  1.0,  0.0,  0.0,  1.0,
    ]);

    let result = a_mat(&x, None).unwrap();

    // All diagonal values should be positive
    for i in 0..4 {
        assert!(
            result.a[(i, i)] > 0.0,
            "Diagonal element A[{},{}] = {} should be positive",
            i,
            i,
            result.a[(i, i)]
        );
    }
}

// ============================================================
// Test 10: VanRaden formula verification
// ============================================================
// Manually compute A using VanRaden formula and compare
#[test]
fn test10_vanraden_formula_verification() {
    #[rustfmt::skip]
    let x = DMatrix::from_row_slice(3, 5, &[
        -1.0,  0.0,  1.0,  0.0, -1.0,
         0.0,  1.0, -1.0,  1.0,  0.0,
         1.0, -1.0,  0.0, -1.0,  1.0,
    ]);

    let n = 3;
    let m = 5;

    // Compute allele frequencies: freq = mean(X + 1) / 2
    let mut freq = vec![0.0; m];
    for j in 0..m {
        let sum: f64 = (0..n).map(|i| x[(i, j)] + 1.0).sum();
        freq[j] = sum / (2.0 * n as f64);
    }

    // For this data, all frequencies should be 0.5
    for &f in &freq {
        assert_relative_eq!(f, 0.5, epsilon = 1e-10);
    }

    // Center: W = X + 1 - 2*freq
    let mut w = DMatrix::zeros(n, m);
    for i in 0..n {
        for j in 0..m {
            w[(i, j)] = x[(i, j)] + 1.0 - 2.0 * freq[j];
        }
    }

    // var.A = 2 * mean(p*(1-p))
    let var_a: f64 = 2.0 * freq.iter().map(|&p| p * (1.0 - p)).sum::<f64>() / m as f64;
    assert_relative_eq!(var_a, 0.5, epsilon = 1e-10);

    // A = WW' / (var_a * m)
    let a_manual = (&w * w.transpose()) / (var_a * m as f64);

    // Compare with a_mat result
    let result = a_mat(&x, None).unwrap();

    for i in 0..n {
        for j in 0..n {
            assert_relative_eq!(result.a[(i, j)], a_manual[(i, j)], epsilon = 1e-10);
        }
    }
}
