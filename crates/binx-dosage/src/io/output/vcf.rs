use super::LocusOutput;
use std::io::Write;

/// Writes the VCF header.
pub fn write_header(writer: &mut dyn Write, sample_names: &[String]) -> anyhow::Result<()> {
    writeln!(writer, "##fileformat=VCFv4.2")?;
    writeln!(writer, "##source=binx-dosage")?;
    writeln!(writer, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">")?;
    writeln!(writer, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")?;
    writeln!(writer, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">")?;
    writeln!(writer, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">")?;
    writeln!(writer, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for REF and ALT\">")?;

    write!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    for name in sample_names {
        write!(writer, "\t{}", name)?;
    }
    writeln!(writer)?;

    Ok(())
}

/// Writes a chunk of results in VCF format.
pub fn write_chunk(
    writer: &mut dyn Write,
    chunk: &[LocusOutput],
    ploidy: usize,
) -> anyhow::Result<()> {
    for res in chunk {
        let num_samples = res.best.len();

        // Use real VCF metadata if available, otherwise use dummy values
        let chrom = res.vcf_chrom.as_ref().map(|s| s.as_str()).unwrap_or("chr1");
        let pos = res.vcf_pos.unwrap_or(1);
        let ref_allele = res.vcf_ref.as_ref().map(|s| s.as_str()).unwrap_or("A");
        let alt_allele = res.vcf_alt.as_ref().map(|s| s.as_str()).unwrap_or("T");

        // CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
        write!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t.\tPASS\tNS={}\tGT:PL:DP:AD",
            chrom, pos, res.id, ref_allele, alt_allele, num_samples
        )?;

        // Sample fields
        for (sample_idx, &best_gt) in res.best.iter().enumerate() {
            // GT field - convert genotype to VCF format
            let gt_str = genotype_to_vcf(best_gt, ploidy);

            // PL field - convert probabilities to Phred-scaled likelihoods
            let sample_probs = res.sample_probs(sample_idx);
            let pl_values: Vec<String> = sample_probs
                .iter()
                .map(|&prob| {
                    let phred = if prob > 0.0 {
                        (-10.0 * prob.log10()).round() as i32
                    } else {
                        255
                    };
                    phred.min(255).to_string()
                })
                .collect();

            // DP field - total depth
            let dp = res.total_counts[[sample_idx]];

            // AD field - ref and alt depths
            let ref_depth = res.ref_counts[[sample_idx]];
            let alt_depth = dp - ref_depth;

            write!(
                writer,
                "\t{}:{}:{}:{},{}",
                gt_str,
                pl_values.join(","),
                dp,
                ref_depth,
                alt_depth
            )?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

/// Converts a genotype (reference allele count: 0, 1, 2, ...) to VCF GT format.
/// For diploid (ploidy=2): 0 -> 1/1, 1 -> 0/1, 2 -> 0/0
/// For tetraploid (ploidy=4): 0 -> 1/1/1/1, 1 -> 0/1/1/1, etc.
fn genotype_to_vcf(genotype: usize, ploidy: usize) -> String {
    let num_ref = genotype;  // genotype is the reference allele count
    let num_alt = ploidy - num_ref;

    let mut alleles = vec!["0"; num_ref];
    alleles.extend(vec!["1"; num_alt]);

    alleles.join("/")
}
