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
        let genotypes: Vec<String> = res.best.iter().map(|g| g.to_string()).collect();
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            res.id,
            res.bias,
            res.rho,
            res.mu,
            res.sigma,
            res.loglik,
            genotypes.join(",")
        )?;
    }
    Ok(())
}
