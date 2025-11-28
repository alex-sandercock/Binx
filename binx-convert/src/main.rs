use anyhow::Result;
use binx_dosage::io::{stream_vcf_records, VcfRecordCounts};
use clap::Parser;
use std::cell::RefCell;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;

#[derive(Parser)]
#[command(name = "binx-convert", version, about = "Convert VCF (AD) to Binx two-line CSV")]
struct Cli {
    /// Input VCF (plain or gzipped)
    #[arg(long)]
    vcf: String,

    /// Output CSV path (two-line format)
    #[arg(long)]
    output: String,

    /// Emit progress to stderr
    #[arg(long, default_value_t = false)]
    verbose: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    convert_vcf_to_csv(&cli.vcf, &cli.output, cli.verbose)
}

fn convert_vcf_to_csv(vcf_path: &str, output_path: &str, verbose: bool) -> Result<()> {
    let writer = RefCell::new(BufWriter::with_capacity(64 * 1024, File::create(output_path)?));
    let header_written = RefCell::new(false);
    let write_error = RefCell::new(None::<anyhow::Error>);
    let processed = RefCell::new(0usize);
    let last_report = RefCell::new(Instant::now());

    stream_vcf_records(
        vcf_path,
        |names| {
            let mut w = writer.borrow_mut();
            let header = std::iter::once("locus".to_string())
                .chain(names.iter().cloned())
                .collect::<Vec<_>>()
                .join(",");
            if let Err(e) = writeln!(w, "{}", header) {
                *write_error.borrow_mut() = Some(anyhow::anyhow!(e));
            } else {
                *header_written.borrow_mut() = true;
            }
        },
        |rec: VcfRecordCounts| {
            if write_error.borrow().is_some() {
                return;
            }
            if !*header_written.borrow() {
                *write_error.borrow_mut() =
                    Some(anyhow::anyhow!("VCF header with sample names was not found"));
                return;
            }

            if let Err(e) = write_record(&mut writer.borrow_mut(), &rec) {
                *write_error.borrow_mut() = Some(e);
                return;
            }

            if verbose && last_report.borrow().elapsed().as_secs_f64() >= 2.0 {
                let mut lr = last_report.borrow_mut();
                *lr = Instant::now();
                let count = {
                    let mut p = processed.borrow_mut();
                    *p += 1;
                    *p
                };
                eprint!("\rConverted {} loci...", count);
                let _ = std::io::stderr().flush();
            } else {
                let mut p = processed.borrow_mut();
                *p += 1;
            }
        },
    )
    .map_err(|e| anyhow::anyhow!("Failed to read VCF: {}", e))?;

    if let Some(e) = write_error.into_inner() {
        return Err(e);
    }

    writer.into_inner().flush()?;
    if verbose {
        let count = processed.into_inner();
        eprintln!("\n✓ Wrote {} loci to {}", count, output_path);
    }
    Ok(())
}

fn write_record(writer: &mut BufWriter<File>, rec: &VcfRecordCounts) -> Result<()> {
    // Ref line
    write!(writer, "{}", rec.id)?;
    for &c in rec.ref_counts.iter() {
        write!(writer, ",{}", c)?;
    }
    writeln!(writer)?;

    // Total line
    write!(writer, "{}", rec.id)?;
    for &c in rec.total_counts.iter() {
        write!(writer, ",{}", c)?;
    }
    writeln!(writer)?;

    Ok(())
}
