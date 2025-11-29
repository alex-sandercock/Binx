
# Binx <img src="docs/assets/binx-logo.png" align="right" width="150"/>

Rust command-line genomics workbench for diploid and polyploid species. `binx` targets GWAS-style analyses with a familiar UX: fast defaults, explicit inputs, and clear TSV/CSV outputs.

## Highlights

- **GWASpoly-style GWAS** (`binx gwaspoly`) with eight genetic models for polyploids, validated against R/GWASpoly
- **Accurate mixed model fitting** via rrblup-rs, a faithful Rust implementation of R/rrBLUP's `mixed.solve`
- **Genotype dosage estimation** (`binx dosage`) from VCF or read count data using updog-style algorithms
- **Kinship matrix computation** (`binx kinship`) via VanRaden or GWASpoly methods
- **VCF conversion** (`binx convert`) to Binx two-line CSV format
- Polyploid-aware: supports ploidy levels 2, 4, 6, etc.
- Handles repeated phenotype IDs (multi-environment trials) and LOCO (Leave-One-Chromosome-Out)

## Installation

Requires a Rust toolchain (`cargo` + `rustc`).

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
binx gwaspoly \
  --geno data/genotypes.tsv \
  --pheno data/phenotypes.csv \
  --trait-name yield \
  --ploidy 4 \
  --kinship kinship.tsv \
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

## Command Reference

### binx gwaspoly
GWASpoly-style GWAS for polyploids with multiple genetic models. Uses validated rrblup-rs mixed model solver.

**Required:**
- `--geno <file>` — Genotype dosage file (TSV/CSV)
- `--pheno <file>` — Phenotype file (TSV/CSV)
- `--trait-name <name>` — Trait column to analyze
- `--ploidy <int>` — Ploidy level (2, 4, 6, etc.)
- `--out <file>` — Output results file

**Optional:**
- `--kinship <file>` — Pre-computed kinship matrix (auto-computed if omitted)
- `--models <list>` — Comma-separated models: `additive`, `general`, `1-dom-ref`, `1-dom-alt`, `2-dom-ref`, `2-dom-alt`, `diplo-general`, `diplo-additive` (default: `additive,general`)
- `--covariates <list>` — Comma-separated covariate column names
- `--loco` — Enable Leave-One-Chromosome-Out kinship
- `--min-maf <float>` — Minimum minor allele frequency filter
- `--max-geno-freq <float>` — Maximum genotype frequency for QC
- `--env-column <col>` + `--env-value <val>` — Filter phenotype rows by environment
- `--allow-missing-samples` — Allow samples in genotypes missing from phenotypes

**Output:** CSV with `marker_id`, `chrom`, `pos`, `model`, `score` (-log10 p), `p_value`, `effect`, `n_obs`

### binx kinship
Compute kinship matrix from biallelic dosages.

**Required:**
- `--geno <file>` — Genotype dosage file
- `--ploidy <int>` — Ploidy level
- `--out <file>` — Output kinship matrix

**Optional:**
- `--method <method>` — `vanraden` (default) or `gwaspoly`

**Output:** Symmetric TSV with header `sample_id <S1> <S2> ...`

### binx dosage
Estimate genotype dosages from sequencing read counts.

**Required:**
- `--ploidy <int>` — Ploidy level
- One of:
  - `--vcf <file>` — VCF file with FORMAT/AD field
  - `--csv <file>` — Two-line CSV (ref counts, total counts per locus)
  - `--counts --ref-path <file> --total-path <file>` — Separate count matrices

**Optional:**
- `--output <file>` — Output path (defaults to stdout)
- `--format <fmt>` — Output format: `matrix`, `stats`, `beagle`, `vcf`, `plink`, `gwaspoly`
- `--mode <mode>` — Optimization mode: `auto`, `updog`, `updog-fast`, `updog-exact`, `fast`, `turbo`, `turboauto`, `turboauto-safe`
- `--threads <int>` — Number of parallel threads
- `--compress <mode>` — `none` or `gzip`
- `--verbose` — Show progress

### binx convert
Convert VCF (with AD field) to Binx two-line CSV format.

**Required:**
- `--vcf <file>` — Input VCF file (plain or gzipped)
- `--output <file>` — Output CSV path

**Optional:**
- `--verbose` — Show progress

### binx multigwas
Multiallelic GWAS (stub for future development).

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
| `binx-core` | Core data structures and I/O utilities |
| `binx-kinship` | Kinship matrix computation |
| `binx-dosage` | Genotype dosage estimation (updog-style) |
| `gwaspoly-rs` | GWASpoly GWAS implementation |
| `rrblup-rs` | R/rrBLUP mixed model solver (validated against R) |

## Validation

The mixed model implementation in `rrblup-rs` has been validated against R/rrBLUP with 52 test cases covering:
- Variance component estimation (REML)
- Fixed and random effect predictions
- Kinship matrix handling
- Missing data handling

Results match R to 5-6 decimal places.

## License

GPL-3.0-or-later
