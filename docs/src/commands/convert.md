# binx convert

Convert VCF files to Binx-compatible formats.

## Synopsis

```bash
binx convert --vcf <FILE> --output <FILE> [OPTIONS]
```

## Description

The `convert` command transforms VCF (Variant Call Format) files into tabular formats used by Binx. It provides two output formats:

- **csv**: Extracts allele depths (AD field) as two-line ref/total counts for use with `binx dosage`
- **gwaspoly**: Extracts genotype dosages (from GT field) for direct use with `binx gwas`

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--vcf <FILE>` | Input VCF file (plain or gzipped) |
| `--output <FILE>` | Output file path |

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--format <FMT>` | csv | Output format: `csv` or `gwaspoly` |
| `--verbose` | false | Enable verbose progress output |

## Output Formats

### csv (default)

Outputs allele depths in a two-line format suitable for `binx dosage`:

- Reads the `AD` (Allelic Depths) field from the VCF
- First line for each locus: reference allele counts
- Second line for each locus: total read counts

```
locus,Sample1,Sample2,Sample3
SNP001,10,2,15
SNP001,15,20,18
SNP002,8,12,5
SNP002,16,25,10
```

Use this format when you want to estimate dosages with `binx dosage`.

### gwaspoly

Outputs genotype dosages in GWASpoly format suitable for `binx gwas`:

- Reads the `GT` (Genotype) field from the VCF
- Converts genotypes to dosage values (count of alternate alleles)
- Handles missing genotypes as `NA`

```
Marker	Chrom	Position	Sample1	Sample2	Sample3
SNP001	chr1	1000	0	2	4
SNP002	chr1	2000	1	1	3
```

Use this format when your VCF has reliable genotype calls.

## Examples

### Convert VCF to Allele Depths (for dosage estimation)

```bash
binx convert \
  --vcf variants.vcf.gz \
  --format csv \
  --output allele_depths.csv
```

Then estimate dosages:

```bash
binx dosage \
  --csv allele_depths.csv \
  --ploidy 4 \
  --output genotypes.tsv
```

### Convert VCF to GWASpoly Format (direct use)

```bash
binx convert \
  --vcf variants.vcf.gz \
  --format gwaspoly \
  --output genotypes.tsv
```

Then run GWAS directly:

```bash
binx gwas \
  --geno genotypes.tsv \
  --pheno phenotypes.csv \
  --trait yield \
  --ploidy 4 \
  --out results.csv
```

### With Verbose Output

```bash
binx convert \
  --vcf variants.vcf.gz \
  --format gwaspoly \
  --output genotypes.tsv \
  --verbose
```

## VCF Requirements

### For csv format (allele depths)

The VCF must contain the `AD` (Allelic Depths) field in the FORMAT column:

```
#CHROM  POS     ID      REF  ALT  QUAL  FILTER  INFO  FORMAT        Sample1       Sample2
chr1    1000    SNP001  A    T    .     .       .     GT:AD:DP      0/1:10,5:15   1/1:2,18:20
```

### For gwaspoly format (genotype dosages)

The VCF must contain the `GT` (Genotype) field:

```
#CHROM  POS     ID      REF  ALT  QUAL  FILTER  INFO  FORMAT  Sample1     Sample2
chr1    1000    SNP001  A    T    .     .       .     GT      0/0/0/1     0/0/1/1
```

Polyploid genotypes are supported (e.g., `0/0/0/1` for tetraploid).

## Choosing Between Formats

| Use Case | Recommended Format |
|----------|-------------------|
| Raw sequencing data with uncertain genotypes | `csv` → `binx dosage` |
| Imputed or high-confidence genotypes | `gwaspoly` → `binx gwas` |
| Low-depth sequencing | `csv` → `binx dosage` |
| Array genotyping data | `gwaspoly` → `binx gwas` |

## See Also

- [binx dosage](dosage.md) - Estimate dosages from allele depths
- [binx gwas](gwas.md) - Run GWAS with genotype data
- [Input Formats](../getting-started/input-formats.md)
