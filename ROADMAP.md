# binx Development Roadmap

This document outlines the main phases for developing **binx**, a Rust-based genomic analysis toolkit with a Sourmash/bcftools-style CLI.

The initial target (`v0.1`) is:

- `binx gwas`: GWASpoly-style **biallelic** GWAS (LM first, then Q+K LMM).
- `binx multigwas`: **multiallelic** GWAS (built on shared mixed-model core).
- Shared core crate (`binx-core`) for data structures, IO, and numerical routines.

---

## Phase 0 – Workspace Skeleton & Tooling

**Goals**

- Create a Rust workspace with separate crates:
  - `binx-core` (core data types + IO)
  - `binx-gwas` (biallelic GWAS)
  - `binx-multigwas` (multiallelic GWAS)
  - `binx-cli` (user-facing binary)
- Establish a bcftools/sourmash-like CLI with subcommands.

**Tasks**

- Set up `Cargo.toml` workspace and crate manifests.
- Add baseline dependencies:
  - `anyhow`, `ndarray`, `ndarray-linalg`, `csv`, `serde`
  - `clap` for CLI
- Implement a minimal `binx` CLI:
  - `binx gwas --help`
  - `binx multigwas --help`
- Configure basic formatting and linting:
  - `cargo fmt`, `cargo clippy`
- Add GPL-3 license and initial README.

---

## Phase 1 – Core Data Model & IO

**Goals**

- Implement shared data structures and IO utilities that both `gwas` and `multigwas` can use.

**Tasks**

- In `binx-core`:
  - Define:
    - `PhenotypeTable` (sample IDs, traits, covariates)
    - `GenotypeMatrixBiallelic` (markers × samples, dosages 0..ploidy)
    - `KinshipMatrix` (samples × samples)
    - `PcMatrix` (optional PCs)
- Implement TSV readers:
  - `load_phenotypes_from_tsv()`
  - `load_genotypes_biallelic_from_tsv()`
  - `load_kinship_from_tsv()`
- Implement basic sample alignment logic between phenotypes and genotypes.
- Add unit tests with tiny toy datasets to verify shapes, IDs, and parsing.

---

## Phase 2 – Simple Linear Model GWAS (No Kinship)

**Goals**

- Implement a first working GWAS engine for **biallelic** markers using a simple linear model (`y ~ intercept + marker`).
- Use this to validate:
  - IO
  - marker-wise loops
  - gene-action encoding
  - integration between `binx-cli`, `binx-gwas`, and `binx-core`.

**Tasks**

- In `binx-gwas`:
  - Define `GeneActionModel` enum (start with `Additive`).
  - Implement `encode_marker()` for additive model.
  - Implement `run_lm_gwas()`:
    - Build design matrix X = [intercept | marker].
    - Fit OLS: `beta = (X^T X)^(-1) X^T y`.
    - Compute residuals, residual variance, standard errors, t-statistics, and p-values.
  - Add TSV writer for GWAS results (`marker_id`, `beta`, `se`, `t_stat`, `p_value`).
  - Implement `run_gwas()` as the main entry point used by the CLI:
    - Load phenotypes and genotypes.
    - Extract the desired trait.
    - For now, ignore K and PCs (LM only).
    - Run GWAS across all markers.
- Wire `binx gwas` subcommand to `binx-gwas::run_gwas()`.
- Validate against R:
  - Build a small synthetic dataset in R.
  - Run per-marker `lm(y ~ dosage)`.
  - Compare `beta` and `p_value` with `binx gwas` output.

---

## Phase 3 – Mixed Model Engine (Q+K) for Biallelic GWAS

**Goals**

- Implement a Q+K linear mixed model engine (LMM) and integrate it with `binx-gwas`.
- Reuse this core later for `multigwas`.

**Tasks**

- In `binx-core`:
  - Implement `fit_null_mixed_model(y, X0, K) -> MixedModelCache`:
    - Estimate variance components (e.g., EMMA/AI-REML style).
    - Build covariance matrix `V = σ_g^2 K + σ_e^2 I`.
    - Compute Cholesky `V = L L^T`.
    - Transform `y* = L^-1 y`, `X0* = L^-1 X0`.
    - Precompute `X0*^T X0*` and `X0*^T y*`.
  - Implement per-marker GLS routines using the cache.
- In `binx-gwas`:
  - Extend `run_gwas()` to support kinship:
    - If `--kinship` is provided, use LMM (Q+K).
    - Otherwise, fall back to LM.
- Validate vs GWASpoly:
  - On a small and a medium dataset, compare:
    - −log10(p) correlation
    - Manhattan and QQ plots.

---

## Phase 4 – Multiallelic GWAS (`binx multigwas`)

**Goals**

- Enable mixed-model GWAS for multiallelic loci (3+ alleles per marker).
- Build this on top of the same mixed-model engine used in `binx-gwas`.

**Tasks**

- Extend `binx-core` genotype layer to support multiallelic sites:
  - Per-site allele list.
  - Per-sample allele counts.
- Define a `MarkerEncoder` trait to separate:
  - Biallelic encoders (for `gwas`).
  - Multiallelic encoders (for `multigwas`).
- Implement multiallelic encoding strategies:
  - Reference-based dosage (K–1 dosage columns for K alleles).
  - Optional collapsed/allelic-series modes.
- In `binx-multigwas`:
  - Implement `run_multigwas()` using the shared mixed-model engine and multiallelic encoders.
  - Expose options like `--encoding`, `--min-allele-freq`, `--collapse-rare-alleles`.
- Validate with simulated multiallelic data.

---

## Phase 5 – Filtering, Kinship, PCA, and Basic Plots

**Goals**

- Provide key supporting functionality to make `binx` usable in end-to-end workflows.

**Tasks**

- `binx filter`:
  - Filter markers/samples based on missingness, MAF, and optional quality metrics.
- `binx kinship`:
  - Compute kinship matrices (e.g., VanRaden) from biallelic dosages.
- `binx pca`:
  - Compute PCA on genotype or kinship matrices.
- `binx plot`:
  - Generate Manhattan and QQ plots from GWAS results (e.g., via SVG/PNG output or JSON for downstream plotting).

---

## Phase 6 – Documentation, Examples, and Release

**Goals**

- Make `binx` easy to install, understand, and cite.

**Tasks**

- Write documentation:
  - Installation instructions.
  - Input format specifications.
  - End-to-end tutorials for `gwas` and `multigwas`.
- Include example datasets and scripts for:
  - Simple GWAS.
  - Multiallelic GWAS.
  - Comparison vs GWASpoly.
- Add a citation file and a short “how to cite” section.
- Tag a stable `v0.1.0` release.
