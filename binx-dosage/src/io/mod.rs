use ndarray::{Array1, Array2};

pub mod matrix;
pub mod vcf;

#[cfg(feature = "bcf")]
pub mod bcf;

/// CSV/two-line inputs produce per-locus data in this shape.
#[derive(Debug, Clone)]
pub struct LocusData {
    pub id: String,
    pub ref_counts: Array1<u32>,
    pub total_counts: Array1<u32>,
}

/// Matrix inputs (markers x samples).
pub struct MatrixData {
    pub marker_ids: Vec<String>,
    pub ref_counts: Array2<u32>,
    pub total_counts: Array2<u32>,
}

/// VCF/BCF streaming output.
pub struct VcfRecordCounts {
    pub id: String,
    pub ref_counts: Array1<u32>,
    pub total_counts: Array1<u32>,
}

pub use matrix::{parse_ref_total_matrices, parse_two_line_csv};
pub use vcf::stream_vcf_records;

#[cfg(feature = "bcf")]
pub use bcf::stream_bcf_records;

// Placeholder when BCF is not built in.
#[cfg(not(feature = "bcf"))]
pub fn stream_bcf_records<F>(
    _path: &str,
    _on_record: F,
) -> Result<(), Box<dyn std::error::Error>>
where
    F: FnMut(VcfRecordCounts),
{
    Err("BCF support is not enabled (compile with feature \"bcf\")".into())
}
