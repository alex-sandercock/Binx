# binx convert

Convert VCF files to Binx-compatible formats.

## Synopsis

```bash
binx convert [OPTIONS] --vcf <FILE>
```

## Description

The `convert` command transforms VCF (Variant Call Format) files into the tabular formats used by Binx. It extracts genotype dosages and marker information while optionally filtering variants.

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--vcf <FILE>` | Path to VCF file (.vcf or .vcf.gz) |

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--output <FILE>` | stdout | Output file path |
| `--format <FMT>` | gwaspoly | Output format |
| `--min-maf <FLOAT>` | 0.0 | Minimum minor allele frequency |
| `--max-missing <FLOAT>` | 1.0 | Maximum missing rate |
| `--biallelic-only` | false | Keep only biallelic variants |

## Output Formats

| Format | Description |
|--------|-------------|
| `gwaspoly` | Standard Binx/GWASpoly format (default) |
| `numeric` | Simple numeric matrix |
| `dosage` | Dosage values with uncertainty |

## Examples

### Basic Conversion

```bash
binx convert \
  --vcf variants.vcf.gz \
  --format gwaspoly \
  --output genotypes.tsv
```

### With Quality Filters

```bash
binx convert \
  --vcf variants.vcf.gz \
  --min-maf 0.05 \
  --max-missing 0.20 \
  --biallelic-only \
  --output genotypes.tsv
```

## VCF Requirements

The VCF file should contain:
- **GT field**: Genotype calls (required)
- **Standard headers**: #CHROM, POS, ID, REF, ALT, etc.

### Supported Genotype Formats

| Ploidy | GT Examples |
|--------|-------------|
| Diploid | 0/0, 0/1, 1/1 |
| Tetraploid | 0/0/0/0, 0/0/0/1, etc. |

## Output Format

```
marker_id	chrom	pos	Sample1	Sample2	Sample3
chr1_1000_A_T	chr1	1000	0	1	2
chr1_2500_G_C	chr1	2500	1	1	0
```

## See Also

- [Input Formats](../getting-started/input-formats.md)
- [binx gwas](gwas.md)
