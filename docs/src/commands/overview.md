# Command Overview

Binx provides a suite of commands for genomic analysis. Each command is designed to handle a specific task in the analysis pipeline.

## Command Structure

All Binx commands follow the pattern:

```bash
binx <command> [options]
```

Get help for any command with `--help`:

```bash
binx --help           # List all commands
binx gwas --help      # Help for specific command
```

## Available Commands

### Analysis Commands

| Command | Description | Primary Use |
|---------|-------------|-------------|
| [`gwas`](gwas.md) | Genome-wide association study | Identify trait-associated markers |
| [`kinship`](kinship.md) | Compute kinship matrix | Account for population structure |
| [`dosage`](dosage.md) | Estimate genotype dosages | Process read count data |

### Utility Commands

| Command | Description | Primary Use |
|---------|-------------|-------------|
| [`convert`](convert.md) | Convert file formats | Prepare VCF data for analysis |
| [`plot`](plot.md) | Generate visualizations | Create Manhattan/QQ plots |
| [`qtl`](qtl.md) | Extract QTLs | Identify significant loci |
| [`threshold`](threshold.md) | Calculate thresholds | Determine significance cutoffs |

## Typical Workflows

### Basic GWAS Pipeline

```bash
# 1. Convert VCF to Binx format
binx convert --vcf data.vcf.gz --output geno.tsv --format gwaspoly

# 2. Compute kinship matrix
binx kinship --geno geno.tsv --ploidy 4 --out kinship.tsv

# 3. Run GWAS
binx gwas --geno geno.tsv --pheno pheno.csv --trait yield \
          --kinship kinship.tsv --ploidy 4 --out results.csv

# 4. Visualize results
binx plot --input results.csv --output manhattan.svg --plot-type manhattan
```

### With Dosage Estimation from VCF

```bash
# 1. Estimate dosages from VCF with allele depths
binx dosage --vcf data.vcf.gz --ploidy 4 --output geno.tsv --format gwaspoly

# 2. Continue with GWAS...
binx gwas --geno geno.tsv --pheno pheno.csv --out results.csv --trait yield --ploidy 4
```

## Common Options

These options are available across multiple commands:

| Option | Description |
|--------|-------------|
| `--help`, `-h` | Display help information |
| `--version`, `-V` | Display version information (top-level only) |
| `--verbose` | Enable verbose output (where applicable) |
| `--threads` | Number of threads (where applicable) |
| `--output` or `--out` | Output file path (varies by command) |

## Exit Codes

| Code | Meaning |
|------|---------|
| `0` | Success |
| `1` | Error (invalid arguments, file not found, processing error, etc.) |

## Next Steps

Explore individual command documentation:

- [binx gwas](gwas.md) - The main GWAS command
- [binx kinship](kinship.md) - Kinship matrix computation
- [binx convert](convert.md) - File format conversion
