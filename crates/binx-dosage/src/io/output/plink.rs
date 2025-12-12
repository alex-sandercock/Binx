use super::LocusOutput;
use std::io::Write;

/// Writes results in PLINK .raw format:
/// Header: FID IID PAT MAT SEX PHENOTYPE followed by SNP names
/// Rows: one per sample with dummy IDs and best-call dosages for each SNP
pub fn write(
    writer: &mut dyn Write,
    results: &[LocusOutput],
    sample_names: &[String],
) -> anyhow::Result<()> {
    if sample_names.is_empty() {
        return Ok(());
    }

    // Write header
    write!(writer, "FID\tIID\tPAT\tMAT\tSEX\tPHENOTYPE")?;
    for res in results {
        write!(writer, "\t{}", res.id)?;
    }
    writeln!(writer)?;

    // Write one row per sample
    for (sample_idx, sample_name) in sample_names.iter().enumerate() {
        // Use sample name as both FID and IID; 0 for missing values
        write!(writer, "{}\t{}\t0\t0\t0\t-9", sample_name, sample_name)?;

        // Write best genotype for each locus
        for res in results {
            write!(writer, "\t{}", res.best[sample_idx])?;
        }
        writeln!(writer)?;
    }

    Ok(())
}
