//! gwaspoly-rs: Rust implementation of R/GWASpoly package
//!
//! This crate provides a faithful implementation of the GWASpoly package for
//! genome-wide association studies in polyploid organisms.
//!
//! ## Features
//! - Eight genetic models (additive, general, 1-dom-ref, 1-dom-alt, 2-dom-ref, 2-dom-alt, diplo-general, diplo-additive)
//! - LOCO (Leave-One-Chromosome-Out) support
//! - Multiple model testing per marker
//! - Q+K mixed model framework
//! - Support for repeated observations per genotype (multi-environment trials)
//! - Parallel marker testing via rayon
//!
//! ## Module Organization
//! - `gwaspoly`: Main GWASpoly function and types (GwasCache, GeneActionModel, MarkerResult)
//! - `set_k`: Kinship matrix computation (faithful implementation of R/GWASpoly set.K)
//! - `parallel`: Parallel marker testing using rayon
//!
//! ## Example
//! ```ignore
//! use gwaspoly_rs::{run_gwaspoly, GeneActionModel};
//!
//! run_gwaspoly(
//!     "genotypes.csv",
//!     "phenotypes.csv",
//!     "yield",
//!     None,           // covariates
//!     None,           // kinship path
//!     true,           // allow missing
//!     None,           // env column (deprecated)
//!     None,           // env value (deprecated)
//!     4,              // ploidy
//!     &[GeneActionModel::Additive],
//!     false,          // LOCO
//!     0.05,           // min MAF
//!     0.95,           // max geno freq
//!     "results.csv",
//!     true,           // parallel
//! )?;
//! ```

// Main GWASpoly module
pub mod gwaspoly;

// Kinship matrix computation (set.K implementation)
pub mod set_k;

// Parallel marker testing module
pub mod parallel;

// Re-export main types and functions from gwaspoly
pub use gwaspoly::{
    GwasCache,
    GeneActionModel,
    MarkerResult,
    run_gwaspoly,
    write_gwaspoly_scores_wide,
};
