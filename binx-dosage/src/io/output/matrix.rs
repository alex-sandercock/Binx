use super::LocusOutput;
use std::io::Write;

/// Writes the matrix format header.
pub fn write_header(writer: &mut dyn Write, sample_names: &[String]) -> anyhow::Result<()> {
    write!(writer, "SNP_ID")?;
    for name in sample_names {
        write!(writer, "\t{}", name)?;
    }
    writeln!(writer)?;
    Ok(())
}

/// Writes a chunk of results in matrix format.
pub fn write_chunk(writer: &mut dyn Write, chunk: &[LocusOutput]) -> anyhow::Result<()> {
    for res in chunk {
        write!(writer, "{}", res.id)?;
        for &genotype in &res.best {
            write!(writer, "\t{}", genotype)?;
        }
        writeln!(writer)?;
    }
    Ok(())
}
