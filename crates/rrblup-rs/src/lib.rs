//! # rrblup-rs
//!
//! Rust implementation of the [R/rrBLUP](https://cran.r-project.org/package=rrBLUP) package
//! for mixed model analysis and genomic prediction.
//!
//! ## Features
//!
//! - [`mixed_solve`](mixed_solve::mixed_solve) - REML-based mixed model solver
//! - [`a_mat()`] - Additive relationship matrix from marker data
//! - [`kin_blup()`] - Genomic BLUP with kinship matrix
//!
//! ## Example
//!
//! ```
//! use rrblup_rs::mixed_solve::{mixed_solve, MixedSolveOptions};
//!
//! let y = vec![1.0, 2.0, 3.0, 4.0, 5.0];
//! let result = mixed_solve(&y, None, None, None, None).unwrap();
//! println!("Vu = {}, Ve = {}", result.vu, result.ve);
//! ```
//!
//! ## References
//!
//! Endelman, J.B. 2011. Ridge regression and other kernels for genomic selection
//! with R package rrBLUP. Plant Genome 4:250-255.

/// REML-based mixed model solver.
pub mod mixed_solve;

/// Additive relationship matrix computation from marker data.
pub mod a_mat;

/// Genomic BLUP with kinship matrix.
pub mod kin_blup;

// Re-export main types and functions from mixed_solve
pub use mixed_solve::{
    mixed_solve as mixed_solve_reml, Method, MixedSolveOptions, MixedSolveResult,
};

// Re-export main types and functions from a_mat
pub use a_mat::{a_mat, AMatOptions, AMatResult, ImputeMethod, ShrinkConfig, ShrinkMethod};

// Re-export main types and functions from kin_blup
pub use kin_blup::{kin_blup, KinBlupData, KinBlupOptions, KinBlupResult};
