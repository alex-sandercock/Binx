mod matrix;
mod stats;
mod beagle;
mod vcf;
mod plink;
mod gwaspoly;

use crate::{CompressMode, OutputFormat};
use flate2::write::GzEncoder;
use flate2::Compression;
use ndarray::Array1;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::Arc;

/// Aggregated result for a single locus across all samples.
#[derive(Clone)]
pub struct LocusOutput {
    pub id: String,
    pub best: Vec<usize>,
    /// Flattened probabilities in row-major order: [sample0_gt0, sample0_gt1, ..., sample1_gt0, sample1_gt1, ...]
    /// Access pattern: probs[sample * num_genotypes + genotype]
    pub probs: Vec<f64>,
    pub num_genotypes: usize,
    pub bias: f64,
    pub rho: f64,
    pub mu: f64,
    pub sigma: f64,
    pub loglik: f64,
    pub ref_counts: Arc<Array1<u32>>,
    pub total_counts: Arc<Array1<u32>>,
    /// VCF metadata (only present when input is VCF)
    pub vcf_chrom: Option<Arc<String>>,
    pub vcf_pos: Option<u64>,
    pub vcf_ref: Option<Arc<String>>,
    pub vcf_alt: Option<Arc<String>>,
}

impl LocusOutput {
    /// Helper to get probabilities for a specific sample
    #[inline]
    pub fn sample_probs(&self, sample_idx: usize) -> &[f64] {
        let start = sample_idx * self.num_genotypes;
        let end = start + self.num_genotypes;
        &self.probs[start..end]
    }
}

/// Streaming writer for formats that can be written incrementally (locus by locus).
/// Supports matrix, stats, beagle, and vcf formats.
pub struct StreamingWriter {
    writer: BufWriter<Box<dyn Write>>,
    format: OutputFormat,
    ploidy: usize,
    sample_names: Vec<String>,
    header_written: bool,
}

impl StreamingWriter {
    /// Creates a new streaming writer for the specified format.
    pub fn new(
        format: OutputFormat,
        compress: CompressMode,
        output_path: Option<&str>,
        ploidy: usize,
    ) -> anyhow::Result<Self> {
        // Open output writer
        let writer: Box<dyn Write> = match output_path {
            Some(path) => Box::new(File::create(path)?),
            None => Box::new(std::io::stdout()),
        };

        // Wrap in gzip encoder if requested
        let buf_writer: BufWriter<Box<dyn Write>> = match compress {
            CompressMode::Gzip => {
                let encoder = GzEncoder::new(writer, Compression::default());
                BufWriter::with_capacity(64 * 1024, Box::new(encoder))
            }
            CompressMode::None => BufWriter::with_capacity(64 * 1024, writer),
        };

        Ok(Self {
            writer: buf_writer,
            format,
            ploidy,
            sample_names: Vec::new(),
            header_written: false,
        })
    }

    /// Writes the header line. Must be called once before writing any chunks.
    pub fn write_header(&mut self, sample_names: &[String]) -> anyhow::Result<()> {
        if self.header_written {
            return Ok(());
        }

        self.sample_names = sample_names.to_vec();

        match self.format {
            OutputFormat::Matrix => matrix::write_header(&mut self.writer, sample_names)?,
            OutputFormat::Stats => stats::write_header(&mut self.writer)?,
            OutputFormat::Beagle => {
                beagle::write_header(&mut self.writer, sample_names, self.ploidy)?
            }
            OutputFormat::Vcf => vcf::write_header(&mut self.writer, sample_names)?,
            OutputFormat::GwasPoly => gwaspoly::write_header(&mut self.writer, sample_names)?,
            OutputFormat::PlinkRaw => {
                return Err(anyhow::anyhow!("PLINK format requires collection, not streaming"))
            }
        }

        self.header_written = true;
        Ok(())
    }

    /// Writes a chunk of locus results.
    pub fn write_chunk(&mut self, chunk: &[LocusOutput]) -> anyhow::Result<()> {
        if !self.header_written {
            return Err(anyhow::anyhow!("Header must be written before chunks"));
        }

        // Validate sample counts
        let expected_samples = self.sample_names.len();
        for res in chunk {
            if res.best.len() != expected_samples {
                return Err(anyhow::anyhow!(
                    "Locus {} has {} samples, expected {}",
                    res.id,
                    res.best.len(),
                    expected_samples
                ));
            }
        }

        match self.format {
            OutputFormat::Matrix => matrix::write_chunk(&mut self.writer, chunk)?,
            OutputFormat::Stats => stats::write_chunk(&mut self.writer, chunk)?,
            OutputFormat::Beagle => beagle::write_chunk(&mut self.writer, chunk)?,
            OutputFormat::Vcf => vcf::write_chunk(&mut self.writer, chunk, self.ploidy)?,
            OutputFormat::GwasPoly => gwaspoly::write_chunk(&mut self.writer, chunk, self.ploidy)?,
            OutputFormat::PlinkRaw => unreachable!("PLINK uses collection path"),
        }

        Ok(())
    }

    /// Flushes and finishes writing.
    pub fn finish(mut self) -> anyhow::Result<()> {
        self.writer.flush()?;
        Ok(())
    }
}

/// Writes collected locus results to the specified output in the requested format.
/// Used only for PLINK .raw format which requires transposition.
/// Handles optional gzip compression.
pub fn write_output(
    results: Vec<LocusOutput>,
    sample_names: &[String],
    format: OutputFormat,
    compress: CompressMode,
    output_path: Option<&str>,
    _ploidy: usize,
) -> anyhow::Result<()> {
    // Validate sample count consistency
    let expected_samples = sample_names.len();
    for res in &results {
        if res.best.len() != expected_samples {
            return Err(anyhow::anyhow!(
                "Locus {} has {} samples, expected {}",
                res.id,
                res.best.len(),
                expected_samples
            ));
        }
    }

    // Open output writer
    let writer: Box<dyn Write> = match output_path {
        Some(path) => Box::new(File::create(path)?),
        None => Box::new(std::io::stdout()),
    };

    // Wrap in gzip encoder if requested
    let mut buf_writer: BufWriter<Box<dyn Write>> = match compress {
        CompressMode::Gzip => {
            let encoder = GzEncoder::new(writer, Compression::default());
            BufWriter::with_capacity(64 * 1024, Box::new(encoder))
        }
        CompressMode::None => BufWriter::with_capacity(64 * 1024, writer),
    };

    // Only PLINK uses this path now
    match format {
        OutputFormat::PlinkRaw => plink::write(&mut buf_writer, &results, sample_names)?,
        _ => {
            return Err(anyhow::anyhow!(
                "Format {:?} should use streaming writer, not collection",
                format
            ))
        }
    }

    buf_writer.flush()?;
    Ok(())
}
