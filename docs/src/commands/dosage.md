# binx dosage

Estimate genotype dosages from sequencing read count data.

## Synopsis

```bash
binx dosage [OPTIONS] --input <FILE>
```

## Description

The `dosage` command estimates genotype dosages from read count data using algorithms based on the R/Updog package (Gerard et al., 2018). This is useful when working with genotyping-by-sequencing (GBS) or similar data where discrete genotype calls may be uncertain.

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--input <FILE>` | Path to read count file |

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--output <FILE>` | stdout | Output dosage file |
| `--ploidy <INT>` | 2 | Ploidy level |
| `--model <MODEL>` | norm | Model for dosage estimation |
| `--min-depth <INT>` | 5 | Minimum read depth |
| `--error-rate <FLOAT>` | 0.01 | Sequencing error rate |

## Input Format

Read count file with reference and alternate allele counts:

```
marker_id	chrom	pos	Sample1_ref	Sample1_alt	Sample2_ref	Sample2_alt
SNP001	1	1000	25	5	0	30
SNP002	1	2000	15	15	20	10
```

## Examples

### Basic Dosage Estimation

```bash
binx dosage \
  --input read_counts.tsv \
  --ploidy 4 \
  --output dosages.tsv
```

### With Quality Filters

```bash
binx dosage \
  --input read_counts.tsv \
  --ploidy 4 \
  --min-depth 10 \
  --error-rate 0.005 \
  --output dosages.tsv
```

## Output

Standard Binx genotype format with estimated dosages:

```
marker_id	chrom	pos	Sample1	Sample2
SNP001	1	1000	1	4
SNP002	1	2000	2	3
```

## See Also

- [binx convert](convert.md) - Alternative genotype input methods
- [Input Formats](../getting-started/input-formats.md)
