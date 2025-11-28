use super::LocusOutput;
use std::io::Write;

/// Writes the stats format header.
pub fn write_header(writer: &mut dyn Write) -> anyhow::Result<()> {
    writeln!(writer, "SNP_ID\tBias\tRho\tMu\tSigma\tLogLik\tBestGenotypes")?;
    Ok(())
}

/// Writes a chunk of results in stats format.
pub fn write_chunk(writer: &mut dyn Write, chunk: &[LocusOutput]) -> anyhow::Result<()> {
    for res in chunk {
        // Write fixed fields
        write!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t",
            res.id, res.bias, res.rho, res.mu, res.sigma, res.loglik
        )?;

        // Write genotypes directly without intermediate allocation
        if let Some((first, rest)) = res.best.split_first() {
            write!(writer, "{}", first)?;
            for gt in rest {
                write!(writer, ",{}", gt)?;
            }
        }
        writeln!(writer)?;
    }
    Ok(())
}
