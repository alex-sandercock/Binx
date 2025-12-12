# Binx <img src="docs/assets/binx-logo.png" align="right" width="150"/>

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CI](https://github.com/alex-sandercock/Binx/actions/workflows/release.yml/badge.svg)](https://github.com/alex-sandercock/Binx/actions/workflows/release.yml)

Rust command-line genomics workbench for diploid and polyploid species. `binx` targets GWAS and related analyses with a familiar UX: fast defaults, explicit inputs, and clear TSV/CSV outputs.

This repo contains:
- rrblup-rs: a faithful Rust port of R/rrBLUP’s mixed.solve and related routines
- gwaspoly-rs: a faithful Rust port of R/GWASpoly
- binx-*: crates that compose these into a multi-method GWAS CLI (binx gwas)

## Highlights

- **GWASpoly-style GWAS** (`binx gwas`) with eight genetic models for polyploids, validated against R/GWASpoly
- **Accurate mixed model fitting** via rrblup-rs, a Rust implementation of R/rrBLUP's `mixed.solve`
- **Genotype dosage estimation** (`binx dosage`) from VCF or read count data using R/Updog-based algorithms
- **Kinship matrix computation** (`binx kinship`) via VanRaden methods
- **VCF conversion** (`binx convert`) to Binx two-line CSV format
- Polyploid-aware: supports ploidy levels 2, 4, 6, etc.
- Handles repeated phenotype IDs (multi-environment trials) and LOCO (Leave-One-Chromosome-Out)

## Installation

### Pre-built Binaries (Recommended)
Download the latest release for your platform from [Releases](https://github.com/alex-sandercock/Binx/releases):
- `binx-linux-x86_64.tar.gz` — Linux
- `binx-macos-x86_64.tar.gz` — macOS (Intel)
- `binx-macos-aarch64.tar.gz` — macOS (Apple Silicon)

```bash
# Download and extract (example: Linux)
curl -LO https://github.com/alex-sandercock/Binx/releases/latest/download/binx-linux-x86_64.tar.gz
tar -xzf binx-linux-x86_64.tar.gz
./binx --help

# Optional: install to PATH
mkdir -p ~/bin && mv binx ~/bin/
# Add to PATH if not already (add to ~/.bashrc or ~/.zshrc):
# export PATH="$HOME/bin:$PATH"
```

### Build from Source
Requires Rust toolchain (`cargo` + `rustc`).

```bash
# Build locally (binary at target/release/binx)
cargo build --release

# Or install into $CARGO_HOME/bin
cargo install --path binx-cli
```

## Quick Start

```bash
# 1) Compute a kinship matrix from biallelic dosages
binx kinship \
  --geno data/genotypes.tsv \
  --ploidy 4 \
  --out kinship.tsv

# 2) Run GWASpoly-style GWAS with multiple genetic models
binx gwas \
  --geno data/genotypes.tsv \
  --pheno data/phenotypes.csv \
  --trait yield \
  --ploidy 4 \
  --models additive,general \
  --out gwas_results.csv

# 3) Estimate genotype dosages from a VCF file
binx dosage \
  --vcf data/samples.vcf.gz \
  --ploidy 4 \
  --output dosages.tsv

# 4) Convert VCF to Binx two-line CSV format
binx convert \
  --vcf data/samples.vcf.gz \
  --output counts.csv
```

## Commands

| Command | Description |
|---------|-------------|
| `binx gwas` | GWASpoly-style GWAS with multiple genetic models |
| `binx kinship` | Compute kinship matrix (VanRaden method) |
| `binx dosage` | Estimate genotype dosages from read counts |
| `binx convert` | Convert VCF to other formats |
| `binx plot` | Generate Manhattan, QQ, or LD decay plots |
| `binx qtl` | Identify significant QTLs from GWAS results |
| `binx threshold` | Calculate significance thresholds (Bonferroni, M.eff, FDR) |

Run `binx <command> --help` for full options.

## Input Formats

### Genotypes (TSV/CSV)
First three columns: `marker_id`, `chrom`, `pos`, followed by sample dosage columns (values 0 to ploidy).
```
marker_id   chrom   pos     Sample1 Sample2 Sample3
SNP001      1       1000    0       2       4
SNP002      1       2000    1       1       3
```

### Phenotypes (TSV/CSV)
First column `sample_id`, remaining columns are traits/covariates. Factor covariates kept as strings.
```
sample_id   yield   env
Sample1     45.2    field_A
Sample2     52.1    field_A
Sample3     48.7    field_B
```

### Kinship (TSV/CSV)
Square symmetric matrix with sample IDs as header and row names.
```
sample_id   Sample1 Sample2 Sample3
Sample1     1.0     0.25    0.15
Sample2     0.25    1.0     0.20
Sample3     0.15    0.20    1.0
```

## Architecture

Binx is organized as a Cargo workspace with specialized crates:

| Crate | Description |
|-------|-------------|
| `binx-cli` | Main CLI binary (`binx`) |
| `binx-gwas` | GWAS orchestration layer |
| `binx-kinship` | Kinship matrix computation |
| `binx-dosage` | Genotype dosage estimation (updog-style) |
| `binx-convert` | VCF conversion utilities |
| `binx-types` | Core data structures and I/O utilities |
| `binx-io` | File I/O helpers |
| `binx-plotting` | Plotting utilities |
| `binx-multigwas` | Multiallelic GWAS (planned) |
| `gwaspoly-rs` | GWASpoly GWAS implementation |
| `rrblup-rs` | R/rrBLUP mixed model solver (validated against R) |

## Validation

### rrblup-rs
Validated against R/rrBLUP with 52 test cases covering:
- Variance component estimation (REML)
- Fixed and random effect predictions
- Kinship matrix handling
- Missing data handling

Results match R to 5-6 decimal places.

### gwaspoly-rs
Validated against R/GWASpoly across all configurations:
- LOCO and non-LOCO kinship
- With and without principal components (P+K model)
- Multiple genetic models (additive, general)
- With and without covariates

Results match R/GWASpoly to 4-5 decimal places.

## Citation

If Binx is useful in your work, please cite the original methods:

> Endelman, J.B. (2011). Ridge regression and other kernels for genomic selection with R package rrBLUP. *The Plant Genome* 4:250-255.

> Rosyara, U.R., De Jong, W.S., Douches, D.S., & Endelman, J.B. (2016). Software for genome-wide association studies in autopolyploids and its application to potato. *The Plant Genome* 9(2).

> Gerard, D., Ferrão, L.F.V., Garcia, A.A.F., & Stephens, M. (2018). Genotyping polyploids from messy sequencing data. *Genetics* 210(3):789-807.

You may also cite the Binx implementation:

> Sandercock, A.M. (2025). Binx: A Rust-based CLI tool for polyploid and diploid genomic analysis. GitHub repository: https://github.com/alex-sandercock/Binx

*(A formal article/DOI is planned; this section will be updated when available.)*

## License

GPL-3.0-or-later
