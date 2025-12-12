//! QTL detection and analysis module
//!
//! This module provides functionality for identifying QTL (Quantitative Trait Loci)
//! from GWAS results. It implements the get.QTL function from R/GWASpoly.
//!
//! ## Features
//! - Filter significant markers based on threshold
//! - Prune redundant signals using bp-window (keeps most significant per region)
//! - Support for multiple genetic models
//!
//! ## Example
//! ```ignore
//! use gwaspoly_rs::qtl::{get_qtl, QtlResult};
//!
//! let qtls = get_qtl(&marker_results, Some(1_000_000))?;
//! for qtl in &qtls {
//!     println!("{}: {} @ {}:{}", qtl.model, qtl.marker_id, qtl.chrom, qtl.pos);
//! }
//! ```

mod get_qtl;

pub use get_qtl::{
    QtlResult,
    get_qtl,
    get_qtl_from_file,
    write_qtl_results,
};
