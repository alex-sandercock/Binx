use super::LocusOutput;
use std::io::Write;

/// Writes the GWASpoly header.
/// Format: Marker, Chrom, Position, followed by sample names.
pub fn write_header(writer: &mut dyn Write, sample_names: &[String]) -> anyhow::Result<()> {
    write!(writer, "Marker\tChrom\tPosition")?;
    for name in sample_names {
        write!(writer, "\t{}", name)?;
    }
    writeln!(writer)?;
    Ok(())
}

/// Writes a chunk of results in GWASpoly format.
/// Each row: marker_id, chromosome, position, dosage1, dosage2, ...
/// Uses VCF metadata (chrom, pos) when available, otherwise uses placeholders.
pub fn write_chunk(
    writer: &mut dyn Write,
    chunk: &[LocusOutput],
    _ploidy: usize,
) -> anyhow::Result<()> {
    for res in chunk {
        // Use VCF metadata if available, otherwise use placeholder values
        let chrom = res.vcf_chrom.as_deref().unwrap_or("chr1");
        let pos = res.vcf_pos.unwrap_or(0);

        // Write marker ID, chromosome, position
        write!(writer, "{}\t{}\t{}", res.id, chrom, pos)?;

        // Write dosages for each sample (best genotype: 0, 1, 2, ..., ploidy)
        for &dosage in &res.best {
            write!(writer, "\t{}", dosage)?;
        }
        writeln!(writer)?;
    }
    Ok(())
}
