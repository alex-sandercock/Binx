//! Core data types for gwaspoly-rs
//!
//! These types mirror the data structures used in R/GWASpoly.

use ndarray::Array2;
use std::collections::HashMap;

pub type SampleId = String;
pub type MarkerId = String;

/// Optional marker metadata for genomic coordinates.
#[derive(Clone, Debug)]
pub struct MarkerMetadata {
    pub chrom: String,
    pub pos: f64,
}

/// Biallelic dosage matrix: markers × samples, entries 0..ploidy.
#[derive(Clone, Debug)]
pub struct GenotypeMatrixBiallelic {
    pub ploidy: u8,
    pub sample_ids: Vec<SampleId>,
    pub marker_ids: Vec<MarkerId>,
    pub marker_metadata: Option<Vec<MarkerMetadata>>,
    /// Shape: (n_markers, n_samples)
    pub dosages: Array2<f64>,
}

/// Kinship matrix with sample labels.
#[derive(Clone, Debug)]
pub struct KinshipMatrix {
    pub sample_ids: Vec<SampleId>,
    pub matrix: Array2<f64>,
}

/// Simple phenotype table: samples × (traits + covariates).
#[derive(Clone, Debug)]
pub struct PhenotypeTable {
    pub sample_ids: Vec<SampleId>,
    pub traits: HashMap<String, Vec<f64>>,
    pub covariates: HashMap<String, Vec<f64>>,
    pub factor_covariates: HashMap<String, Vec<String>>,
}
