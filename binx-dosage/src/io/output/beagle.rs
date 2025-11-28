use super::LocusOutput;
use std::io::Write;

/// Writes the Beagle format header.
/// For ploidy n, there are n+1 possible genotypes.
pub fn write_header(
    writer: &mut dyn Write,
    sample_names: &[String],
    ploidy: usize,
) -> anyhow::Result<()> {
    let num_genotypes = ploidy + 1;

    write!(writer, "marker\tallele1\tallele2")?;
    for name in sample_names {
        for _ in 0..num_genotypes {
            write!(writer, "\t{}", name)?;
        }
    }
    writeln!(writer)?;
    Ok(())
}

/// Writes a chunk of results in Beagle format.
pub fn write_chunk(writer: &mut dyn Write, chunk: &[LocusOutput]) -> anyhow::Result<()> {
    for res in chunk {
        write!(writer, "{}\tA\tT", res.id)?;
        for sample_probs in &res.probs {
            for prob in sample_probs {
                write!(writer, " {:.6}", prob)?;
            }
        }
        writeln!(writer)?;
    }
    Ok(())
}
