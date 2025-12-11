//! gwaspoly-rs: Rust implementation of R/GWASpoly package
//!
//! This crate provides a Rust implementation of the GWASpoly package for
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
//! - `types`: Core data types (GenotypeMatrixBiallelic, KinshipMatrix, etc.)
//! - `io`: I/O functions including read_gwaspoly (mirrors R's read.GWASpoly)
//! - `gwaspoly`: Main GWASpoly function and types (GwasCache, GeneActionModel, MarkerResult)
//! - `set_k`: Kinship matrix computation (Rust implementation of R/GWASpoly set.K)
//! - `parallel`: Parallel marker testing using rayon
//! - `threshold`: Significance threshold calculation (set.threshold)
//! - `qtl`: QTL detection (get.QTL)
//!
//! ## Example
//! ```ignore
//! use gwaspoly_rs::{gwaspoly, read_gwaspoly, GeneActionModel};
//!
//! // Load data (mirrors R's read.GWASpoly)
//! let (geno, pheno) = read_gwaspoly("genotypes.csv", "phenotypes.csv", 4)?;
//!
//! // Run GWAS (mirrors R's GWASpoly())
//! gwaspoly(
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

// Core types (standalone, no binx dependencies)
pub mod types;

// I/O functions (mirrors R's read.GWASpoly)
pub mod io;

// Main GWASpoly module
pub mod gwaspoly;

// Kinship matrix computation (set.K implementation)
pub mod set_k;

// Parallel marker testing module
pub mod parallel;

// Threshold calculation module (set.threshold implementation)
pub mod threshold;

// QTL detection module (get.QTL implementation)
pub mod qtl;

// Re-export main types and functions from gwaspoly (mirrors R's GWASpoly())
pub use gwaspoly::{
    GwasCache,
    GeneActionModel,
    MarkerResult,
    gwaspoly,
    write_gwaspoly_scores_wide,
};

// Re-export threshold types and functions (mirrors R's set.threshold())
pub use threshold::{
    ThresholdMethod,
    ThresholdResult,
    set_threshold,
    set_threshold_simple,
    load_gwas_results_for_threshold,
    print_thresholds,
};

// Re-export QTL types and functions
pub use qtl::{
    QtlResult,
    get_qtl,
    get_qtl_from_file,
    write_qtl_results,
};

// Re-export types (standalone, mirrors R/GWASpoly data structures)
pub use types::{
    GenotypeMatrixBiallelic,
    KinshipMatrix,
    MarkerMetadata,
    PhenotypeTable,
    SampleId,
    MarkerId,
};

// Re-export I/O functions (mirrors R's read.GWASpoly)
pub use io::{
    read_gwaspoly,
    load_genotypes,
    load_phenotypes,
    load_kinship,
};

// Re-export kinship functions (mirrors R's set.K())
pub use set_k::{
    set_k,
    compute_loco_kinship_gwaspoly,
};
