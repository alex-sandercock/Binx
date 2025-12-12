//! Comparison tests between Rust kin_blup and R/rrBLUP::kin.blup
//!
//! Reference values generated from R/rrBLUP v4.6.3 using:
//! `Rscript tests/generate_kinblup_reference.R`

use approx::assert_relative_eq;
use nalgebra::DMatrix;
use rrblup_rs::kin_blup::{kin_blup, KinBlupData, KinBlupOptions};

/// Tolerance for floating point comparisons
const TOLERANCE: f64 = 1e-4;
const LOOSE_TOLERANCE: f64 = 1e-2;

// ============================================================
// Test 1: No kinship matrix (genotype as random factor)
// ============================================================
// Vg = 0.9375653009
// Ve = 0.1249947765
// g: A = -0.9375065, B ≈ 0, C = 0.9375065
// pred: A = 1.312493, B = 2.25, C = 3.187507
#[test]
fn test1_no_kinship_vs_r() {
    let data = KinBlupData {
        geno_ids: vec![
            "A".into(),
            "B".into(),
            "C".into(),
            "A".into(),
            "B".into(),
            "C".into(),
        ],
        pheno: vec![1.0, 2.0, 3.0, 1.5, 2.5, 3.5],
        fixed: None,
        covariates: None,
    };

    let result = kin_blup(&data, None).unwrap();

    assert_relative_eq!(result.vg, 0.9375653009, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(result.ve, 0.1249947765, epsilon = LOOSE_TOLERANCE);

    // Find indices for A, B, C
    let idx_a = result.g_ids.iter().position(|x| x == "A").unwrap();
    let idx_b = result.g_ids.iter().position(|x| x == "B").unwrap();
    let idx_c = result.g_ids.iter().position(|x| x == "C").unwrap();

    assert_relative_eq!(result.g[idx_a], -0.9375065, epsilon = LOOSE_TOLERANCE);
    assert!(result.g[idx_b].abs() < 0.01, "B should be near 0");
    assert_relative_eq!(result.g[idx_c], 0.9375065, epsilon = LOOSE_TOLERANCE);

    // Check predictions
    assert_relative_eq!(result.pred[idx_a], 1.312493, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(result.pred[idx_b], 2.25, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(result.pred[idx_c], 3.187507, epsilon = LOOSE_TOLERANCE);
}

// ============================================================
// Test 2: With kinship matrix
// ============================================================
// Vg = 2.0636058069
// Ve = 0.0001407381
// g: A = -1.6698503, B = -0.1700292, C = 0.5299801, D = 1.4298831
#[test]
fn test2_with_kinship_vs_r() {
    let data = KinBlupData {
        geno_ids: vec!["A".into(), "B".into(), "C".into(), "D".into()],
        pheno: vec![1.0, 2.5, 3.2, 4.1],
        fixed: None,
        covariates: None,
    };

    #[rustfmt::skip]
    let k = DMatrix::from_row_slice(4, 4, &[
        1.0, 0.5, 0.3, 0.2,
        0.5, 1.0, 0.4, 0.3,
        0.3, 0.4, 1.0, 0.5,
        0.2, 0.3, 0.5, 1.0,
    ]);

    let opts = KinBlupOptions {
        k: Some(k),
        k_ids: Some(vec!["A".into(), "B".into(), "C".into(), "D".into()]),
        ..Default::default()
    };

    let result = kin_blup(&data, Some(opts)).unwrap();

    assert_relative_eq!(result.vg, 2.0636058069, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(result.ve, 0.0001407381, epsilon = LOOSE_TOLERANCE);

    // Find indices
    let idx_a = result.g_ids.iter().position(|x| x == "A").unwrap();
    let idx_b = result.g_ids.iter().position(|x| x == "B").unwrap();
    let idx_c = result.g_ids.iter().position(|x| x == "C").unwrap();
    let idx_d = result.g_ids.iter().position(|x| x == "D").unwrap();

    assert_relative_eq!(result.g[idx_a], -1.6698503, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(result.g[idx_b], -0.1700292, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(result.g[idx_c], 0.5299801, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(result.g[idx_d], 1.4298831, epsilon = LOOSE_TOLERANCE);
}

// ============================================================
// Test 3: With missing phenotypes (genomic prediction)
// ============================================================
// Individuals C and E have missing phenotypes
// Vg = 2.9029583725
// Ve = 0.0001979821
// g: A=-1.5966, B=-0.0968, C=0.3871, D=1.4031, E=0.5564
#[test]
fn test3_missing_phenotypes_vs_r() {
    let data = KinBlupData {
        geno_ids: vec![
            "A".into(),
            "B".into(),
            "C".into(),
            "D".into(),
            "E".into(),
        ],
        pheno: vec![1.0, 2.5, f64::NAN, 4.0, f64::NAN],
        fixed: None,
        covariates: None,
    };

    #[rustfmt::skip]
    let k = DMatrix::from_row_slice(5, 5, &[
        1.0, 0.5, 0.3, 0.2, 0.1,
        0.5, 1.0, 0.4, 0.3, 0.2,
        0.3, 0.4, 1.0, 0.5, 0.3,
        0.2, 0.3, 0.5, 1.0, 0.4,
        0.1, 0.2, 0.3, 0.4, 1.0,
    ]);

    let opts = KinBlupOptions {
        k: Some(k),
        k_ids: Some(vec![
            "A".into(),
            "B".into(),
            "C".into(),
            "D".into(),
            "E".into(),
        ]),
        ..Default::default()
    };

    let result = kin_blup(&data, Some(opts)).unwrap();

    assert_relative_eq!(result.vg, 2.9029583725, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(result.ve, 0.0001979821, epsilon = LOOSE_TOLERANCE);

    // Should have predictions for all 5 genotypes including C and E
    assert_eq!(result.g.len(), 5);
    assert_eq!(result.g_ids.len(), 5);

    // Find indices
    let idx_a = result.g_ids.iter().position(|x| x == "A").unwrap();
    let idx_b = result.g_ids.iter().position(|x| x == "B").unwrap();
    let idx_c = result.g_ids.iter().position(|x| x == "C").unwrap();
    let idx_d = result.g_ids.iter().position(|x| x == "D").unwrap();
    let idx_e = result.g_ids.iter().position(|x| x == "E").unwrap();

    assert_relative_eq!(result.g[idx_a], -1.5966209, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(result.g[idx_b], -0.0968024, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(result.g[idx_c], 0.3870611, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(result.g[idx_d], 1.4031151, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(result.g[idx_e], 0.5564034, epsilon = LOOSE_TOLERANCE);

    // Residuals for C and E should be NaN
    assert!(result.resid[2].is_nan()); // C
    assert!(result.resid[4].is_nan()); // E
}

// ============================================================
// Test 4: With PEV=TRUE
// ============================================================
// PEV: A=1.073168, B=1.073196, C=1.073196, D=1.073168
#[test]
fn test4_with_pev_vs_r() {
    let data = KinBlupData {
        geno_ids: vec!["A".into(), "B".into(), "C".into(), "D".into()],
        pheno: vec![1.0, 2.5, 3.2, 4.1],
        fixed: None,
        covariates: None,
    };

    #[rustfmt::skip]
    let k = DMatrix::from_row_slice(4, 4, &[
        1.0, 0.5, 0.3, 0.2,
        0.5, 1.0, 0.4, 0.3,
        0.3, 0.4, 1.0, 0.5,
        0.2, 0.3, 0.5, 1.0,
    ]);

    let opts = KinBlupOptions {
        k: Some(k),
        k_ids: Some(vec!["A".into(), "B".into(), "C".into(), "D".into()]),
        pev: true,
        ..Default::default()
    };

    let result = kin_blup(&data, Some(opts)).unwrap();

    assert!(result.pev.is_some());
    let pev = result.pev.unwrap();
    assert_eq!(pev.len(), 4);

    // Find indices
    let idx_a = result.g_ids.iter().position(|x| x == "A").unwrap();
    let idx_b = result.g_ids.iter().position(|x| x == "B").unwrap();
    let idx_c = result.g_ids.iter().position(|x| x == "C").unwrap();
    let idx_d = result.g_ids.iter().position(|x| x == "D").unwrap();

    // R gives PEV ≈ 1.073
    assert_relative_eq!(pev[idx_a], 1.073168, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(pev[idx_b], 1.073196, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(pev[idx_c], 1.073196, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(pev[idx_d], 1.073168, epsilon = LOOSE_TOLERANCE);
}

// ============================================================
// Test 5: With fixed effects
// ============================================================
// Vg = 0.4999829506
// Ve = 0.0000340989
// g: A ≈ -1, B ≈ 0, C ≈ 1
#[test]
fn test5_with_fixed_effects_vs_r() {
    let data = KinBlupData {
        geno_ids: vec![
            "A".into(),
            "B".into(),
            "C".into(),
            "A".into(),
            "B".into(),
            "C".into(),
        ],
        pheno: vec![1.0, 2.0, 3.0, 2.0, 3.0, 4.0],
        fixed: Some(vec![vec![
            "E1".into(),
            "E1".into(),
            "E1".into(),
            "E2".into(),
            "E2".into(),
            "E2".into(),
        ]]),
        covariates: None,
    };

    let result = kin_blup(&data, None).unwrap();

    assert_relative_eq!(result.vg, 0.4999829506, epsilon = LOOSE_TOLERANCE);
    assert_relative_eq!(result.ve, 0.0000340989, epsilon = LOOSE_TOLERANCE);

    // Find indices
    let idx_a = result.g_ids.iter().position(|x| x == "A").unwrap();
    let idx_b = result.g_ids.iter().position(|x| x == "B").unwrap();
    let idx_c = result.g_ids.iter().position(|x| x == "C").unwrap();

    assert_relative_eq!(result.g[idx_a], -1.0, epsilon = LOOSE_TOLERANCE);
    assert!(result.g[idx_b].abs() < 0.01, "B should be near 0");
    assert_relative_eq!(result.g[idx_c], 1.0, epsilon = LOOSE_TOLERANCE);
}

// ============================================================
// Test 6: With covariates
// ============================================================
// When covariate explains all variance, Vg ≈ 0, g ≈ 0
#[test]
fn test6_with_covariates_vs_r() {
    let data = KinBlupData {
        geno_ids: vec!["A".into(), "B".into(), "C".into(), "D".into()],
        pheno: vec![1.0, 2.2, 3.5, 4.8],
        fixed: None,
        covariates: Some(vec![vec![0.0, 1.0, 2.0, 3.0]]),
    };

    #[rustfmt::skip]
    let k = DMatrix::from_row_slice(4, 4, &[
        1.0, 0.5, 0.3, 0.2,
        0.5, 1.0, 0.4, 0.3,
        0.3, 0.4, 1.0, 0.5,
        0.2, 0.3, 0.5, 1.0,
    ]);

    let opts = KinBlupOptions {
        k: Some(k),
        k_ids: Some(vec!["A".into(), "B".into(), "C".into(), "D".into()]),
        ..Default::default()
    };

    let result = kin_blup(&data, Some(opts)).unwrap();

    // R gives Vg ≈ 0 because covariate explains all variation
    assert!(result.vg < 0.01, "Vg should be near 0, got {}", result.vg);

    // Genetic values should be near 0
    for &g in &result.g {
        assert!(g.abs() < 0.01, "g should be near 0, got {}", g);
    }
}

// ============================================================
// Test 7: Check residuals sum to near zero
// ============================================================
#[test]
fn test7_residuals_sum_to_zero() {
    let data = KinBlupData {
        geno_ids: vec![
            "A".into(),
            "B".into(),
            "C".into(),
            "A".into(),
            "B".into(),
            "C".into(),
        ],
        pheno: vec![1.0, 2.0, 3.0, 1.5, 2.5, 3.5],
        fixed: None,
        covariates: None,
    };

    let result = kin_blup(&data, None).unwrap();

    // Residuals should sum to approximately zero
    let resid_sum: f64 = result.resid.iter().filter(|x| x.is_finite()).sum();
    assert!(
        resid_sum.abs() < 0.01,
        "Residuals should sum to ~0, got {}",
        resid_sum
    );
}

// ============================================================
// Test 8: Gaussian kernel (basic structure check)
// ============================================================
#[test]
fn test8_gaussian_kernel_structure() {
    let data = KinBlupData {
        geno_ids: vec!["A".into(), "B".into(), "C".into(), "D".into()],
        pheno: vec![1.0, 1.8, 2.5, 3.5],
        fixed: None,
        covariates: None,
    };

    // Use a distance-like matrix
    #[rustfmt::skip]
    let k = DMatrix::from_row_slice(4, 4, &[
        0.0, 1.0, 2.0, 3.0,
        1.0, 0.0, 1.0, 2.0,
        2.0, 1.0, 0.0, 1.0,
        3.0, 2.0, 1.0, 0.0,
    ]);

    let opts = KinBlupOptions {
        k: Some(k),
        k_ids: Some(vec!["A".into(), "B".into(), "C".into(), "D".into()]),
        gauss: true,
        theta_seq: Some(vec![0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0]),
        ..Default::default()
    };

    let result = kin_blup(&data, Some(opts)).unwrap();

    // Should have profile likelihood
    assert!(result.profile.is_some());
    let profile = result.profile.unwrap();
    assert_eq!(profile.len(), 10);

    // Variance components should be positive
    assert!(result.vg >= 0.0);
    assert!(result.ve >= 0.0);

    // Should have genetic values for all 4 genotypes
    assert_eq!(result.g.len(), 4);
}

// ============================================================
// Test 9: Symmetric K matrix handling
// ============================================================
#[test]
fn test9_k_matrix_symmetry() {
    let data = KinBlupData {
        geno_ids: vec!["A".into(), "B".into(), "C".into()],
        pheno: vec![1.0, 2.0, 3.0],
        fixed: None,
        covariates: None,
    };

    // Symmetric K
    #[rustfmt::skip]
    let k = DMatrix::from_row_slice(3, 3, &[
        1.0, 0.5, 0.3,
        0.5, 1.0, 0.4,
        0.3, 0.4, 1.0,
    ]);

    let opts = KinBlupOptions {
        k: Some(k),
        k_ids: Some(vec!["A".into(), "B".into(), "C".into()]),
        ..Default::default()
    };

    let result = kin_blup(&data, Some(opts)).unwrap();

    assert!(result.vg >= 0.0);
    assert!(result.ve >= 0.0);
    assert_eq!(result.g.len(), 3);
}

// ============================================================
// Test 10: Different K ordering than data
// ============================================================
#[test]
fn test10_k_different_ordering() {
    let data = KinBlupData {
        geno_ids: vec!["C".into(), "A".into(), "B".into()], // Different order
        pheno: vec![3.0, 1.0, 2.0],
        fixed: None,
        covariates: None,
    };

    // K has standard ordering A, B, C
    #[rustfmt::skip]
    let k = DMatrix::from_row_slice(3, 3, &[
        1.0, 0.5, 0.3,
        0.5, 1.0, 0.4,
        0.3, 0.4, 1.0,
    ]);

    let opts = KinBlupOptions {
        k: Some(k),
        k_ids: Some(vec!["A".into(), "B".into(), "C".into()]),
        ..Default::default()
    };

    let result = kin_blup(&data, Some(opts)).unwrap();

    // Should correctly map data genotypes to K
    assert_eq!(result.g.len(), 3);
    assert!(result.vg >= 0.0);
}
