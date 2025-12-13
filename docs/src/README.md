<p align="center">
  <img src="assets/binx-logo.png" alt="Binx Logo" width="200">
</p>

# Welcome to Binx!

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CI](https://github.com/alex-sandercock/Binx/actions/workflows/release.yml/badge.svg)](https://github.com/alex-sandercock/Binx/actions/workflows/release.yml)

**Binx** is a Rust command-line genomics workbench for diploid and polyploid species. It targets GWAS and related analyses with a familiar UX: fast defaults, explicit inputs, and clear outputs.

## What can Binx do?

Binx provides a suite of tools for genomic analysis:

| Command | Description |
|---------|-------------|
| [`binx gwas`](commands/gwas.md) | GWASpoly-style GWAS with multiple genetic models |
| [`binx kinship`](commands/kinship.md) | Compute kinship matrix (VanRaden method) |
| [`binx dosage`](commands/dosage.md) | Estimate genotype dosages from read counts |
| [`binx convert`](commands/convert.md) | Convert VCF to other formats |
| [`binx plot`](commands/plot.md) | Generate Manhattan, QQ, or LD decay plots |
| [`binx qtl`](commands/qtl.md) | Identify significant QTLs from GWAS results |
| [`binx threshold`](commands/threshold.md) | Calculate significance thresholds |

## Key Features

- **GWASpoly-style GWAS** with eight genetic models for polyploids, validated against R/GWASpoly
- **Accurate mixed model fitting** via rrblup-rs, a Rust implementation of R/rrBLUP's `mixed.solve`
- **Genotype dosage estimation** from VCF or read count data using R/Updog-based algorithms
- **Polyploid-aware**: supports ploidy levels 2, 4, 6, etc.
- **LOCO support**: Leave-One-Chromosome-Out analysis
- **Multi-environment trials**: handles repeated phenotype IDs

## Quick Example

```bash
# Convert VCF to GWASpoly format
binx convert --vcf samples.vcf.gz --format gwaspoly --output genotypes.tsv

# Run GWAS with multiple genetic models
binx gwas \
  --geno genotypes.tsv \
  --pheno phenotypes.csv \
  --trait yield \
  --ploidy 4 \
  --models additive,general \
  --out gwas_results.csv

# Create a Manhattan plot
binx plot \
  --input gwas_results.csv \
  --plot-type manhattan \
  --model additive \
  --output gwas_manhattan.svg
```

## Getting Started

New to Binx? Start here:

1. **[Installation](getting-started/installation.md)** - Download and install Binx
2. **[Quick Start](getting-started/quickstart.md)** - Run your first analysis in minutes
3. **[Input Formats](getting-started/input-formats.md)** - Understand the data formats

## Tutorials

Learn Binx through practical examples:

- [Your First GWAS Analysis](tutorials/first-gwas.md) - Step-by-step GWAS walkthrough
- [Working with Polyploids](tutorials/polyploid-analysis.md) - Tetraploid and hexaploid analysis
- [Multi-Environment Trials](tutorials/multi-environment.md) - Handling complex experimental designs
- [From VCF to Results](tutorials/vcf-to-results.md) - Complete pipeline example

## Why Binx?

Binx was created to bring the power of R/GWASpoly and R/rrBLUP to the command line with:

- **Speed**: Written in Rust for fast execution
- **Reproducibility**: Explicit parameters and deterministic outputs
- **Validation**: Results match R implementations to 4-6 decimal places
- **Simplicity**: No R environment or dependencies required

## Getting Help

- **Questions?** Open an issue on [GitHub](https://github.com/alex-sandercock/Binx/issues)
- **Found a bug?** Please [report it](https://github.com/alex-sandercock/Binx/issues/new)

## License

Binx is released under the [GPL-3.0 license](https://www.gnu.org/licenses/gpl-3.0).
