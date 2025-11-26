# GWAS & FastGWAS Requirements

This document defines the requirements for the `gwas` and `fastgwas` crates in the binx toolkit.

- `gwas` = GWASpoly-faithful engine (R GWASpoly -> Rust, crate: `binx-gwas`)
- `fastgwas` = enhanced engine with performance and feature improvements (crate: `binx-fastgwas`)
- `binx-kinship` remains a separate crate used by both engines.
- `binx-multigwas` is a distinct crate for multiallelic data GWAS and is out of scope here.

The goal is to guarantee that `gwas` is scientifically equivalent to GWASpoly under the same inputs and settings, while `fastgwas` is free to innovate.

---

## 0. Terminology & Scope

- **GWASpoly**: The R package for polyploid GWAS, treated as the reference implementation.
- **binx-core**: Shared code (genotype/phenotype IO, math utils, etc.).
- **gwas crate**: Rust implementation that mirrors GWASpoly behavior (`binx-gwas`).
- **fastgwas crate**: Rust implementation that builds on binx-core but may diverge from GWASpoly internals and outputs (`binx-fastgwas`).
- **multigwas crate**: Separate crate for multiallelic data GWAS (`binx-multigwas`), not covered by this spec.

Scope of this spec:

- Requirements for v1.0 of `gwas` and `fastgwas`.
- Focused on mixed-model GWAS for quantitative traits in (tetra)ploids, with 5k tetraploid alfalfa as the main benchmark scale.

---

## 1. High-Level Goals

### 1.1 `gwas` (GWASpoly-faithful)

- Primary goal: For a given dataset and set of options, `binx gwas` reproduces GWASpoly results (marker set plus statistics) within numerical tolerance.
- Secondary goal: Be at least several times faster and more memory efficient than the R implementation, but never at the expense of correctness.

### 1.2 `fastgwas`

- Primary goal: Provide a high-performance GWAS engine for large polyploid panels using more advanced algorithms and multithreading.
- Secondary goal: Maintain similar scientific conclusions but allow small, documented deviations from GWASpoly (for example slightly different p-values, different default filters, extra models).

---

## 2. `gwas` Crate Requirements (GWASpoly-Faithful)

### 2.1 Functional Requirements

#### 2.1.1 Input Data & Semantics

- [ ] Genotype input
  - [ ] Accept a dense dosage matrix or GWASpoly-style dosage file.
  - [ ] Support integer allele dosages from `0` to `ploidy` (for example `0..=4` for tetraploid).
  - [ ] Honor sample/marker IDs and maintain consistent order across geno/pheno/kinship.
- [ ] Phenotype input
  - [ ] Accept a table (CSV/TSV) with:
    - [ ] One row per individual.
    - [ ] One or more trait columns (numeric).
    - [ ] Optional covariates (for example PCs, factors).
  - [ ] Provide a way to choose:
    - [ ] Trait column (`--trait`).
    - [ ] Covariates (`--covariates`, comma-separated).
- [ ] Kinship input
  - [ ] Accept an optional precomputed kinship matrix (dense).
  - [ ] If not provided, compute kinship internally using a method equivalent to GWASpoly default (for example VanRaden-style).
- [ ] Ploidy
  - [ ] Require an explicit `--ploidy` argument.
  - [ ] Use ploidy to interpret dosages and construct design matrices exactly as GWASpoly does.

#### 2.1.1a Repeated Measures / Multi-Environment Phenotypes (GWASpoly-Style)

Goal: `gwas` must handle phenotype files in the same way as GWASpoly, including cases where the same genotype is measured in multiple environments or trials.

- [ ] Allow repeated genotype IDs
  - Phenotype tables may contain multiple rows with the same `id`; each row is a distinct observation.
- [ ] Unique individuals vs observations
  - Distinguish `n_ind` = number of unique genotype IDs from `n_obs` = total rows/observations.
  - Genotypes (dosage matrix, kinship) are at the individual level; response vector and fixed-effect design are at the observation level.
- [ ] Model structure
  - Use an incidence matrix `Z` of size `n_obs x n_ind` mapping observations to genotypes.
  - Build the mixed model so it matches GWASpoly for repeated measures, for example:
    ```
    y_ij = mu + env_j + g_i + beta_m * x_im + e_ij
    ```
    where `env_j` is a fixed effect, `g_i ~ N(0, sigma_g^2 K)` is the random polygenic effect, and `e_ij ~ N(0, sigma_e^2)` is the residual.
- [ ] Fixed-effect detection
  - After the first column (`id`) and the configured number of trait columns, treat remaining columns as candidate fixed effects (mirroring GWASpoly).
  - Respect factor vs numeric typing in params/CLI (for example `fixed=env`, `fixed.type=factor` in GWASpoly should be reproducible via CLI).
- [ ] Reference parity
  - Use the potato example phenotypes/genotypes (new_potato_pheno / new_potato_geno) as a parity test:
    - Confirm the same number of unique individuals and fixed effects.
    - Confirm numerically equivalent GWAS results within tolerance when using the same kinship and model settings.

`fastgwas` must also accept these repeated-measure phenotypes. It may either mirror the observation-level model or use a two-stage approach (for example precomputed means/BLUEs) so long as behavior differences are documented and the CLI does not force users to pre-collapse repeated measures.

#### 2.1.2 Models & Mixed-Model Structure

- [ ] Implement at least the following models:
  - [ ] `additive`
  - [ ] `general`
- [ ] Match GWASpoly treatment of:
  - [ ] Fixed effects: intercept, covariates, marker effects.
  - [ ] Random effects: polygenic term with kinship matrix.
  - [ ] Residual variance.
- [ ] Implement variance component estimation that is numerically compatible with GWASpoly approach (for example restricted maximum likelihood; behavior matched via tests even if internals differ).

#### 2.1.3 Marker Filtering / Curation

- [ ] Implement QC steps equivalent to GWASpoly:
  - [ ] Minor allele frequency filter (`min_maf`).
  - [ ] Maximum missing rate per marker.
  - [ ] Any model-specific filtering (for example monomorphic markers removal).
- [ ] Defaults:
  - [ ] Mirror GWASpoly default thresholds where known.
  - [ ] Expose CLI flags to override (for example `--min-maf`, `--max-missing`).
- [ ] Requirement: Given the same raw data and thresholds, `gwas` must include/exclude the same markers as GWASpoly.

#### 2.1.4 Kinship & LOCO

- [ ] Implement a kinship estimator equivalent to GWASpoly default (for example from centered genotype matrix).
- [ ] Support LOCO (leave-one-chromosome-out) if GWASpoly does:
  - [ ] Ability to specify chromosome per marker.
  - [ ] Use chromosome-specific kinship matrices when requested.
- [ ] Allow passing a precomputed kinship; in that case, do not recompute.

#### 2.1.5 Per-Marker Tests & Outputs

- [ ] For each marker passing QC:
  - [ ] Fit the model or compute a test equivalent to GWASpoly per-marker logic:
    - [ ] Effect size(s) for marker.
    - [ ] Standard error(s).
    - [ ] Test statistic (for example LRT, score).
    - [ ] p-value.
- [ ] Output:
  - [ ] One row per marker plus model.
  - [ ] Columns must include, at minimum:
    - [ ] `marker_id`
    - [ ] `chrom` (if available)
    - [ ] `pos` (if available)
    - [ ] `model` (for example `additive`, `general`)
    - [ ] `n_obs`
    - [ ] `effect` (or `effect_1`, `effect_2`, etc., for multi-parameter models)
    - [ ] `se`
    - [ ] `stat` (for example score/LRT)
    - [ ] `p_value`
  - [ ] Column naming should be clearly mappable to GWASpoly output schema (document mapping in README).

#### 2.1.6 CLI & Config

- [ ] Provide a `binx gwas` subcommand that:
  - [ ] Accepts:
    - [ ] `--geno`
    - [ ] `--pheno`
    - [ ] `--kinship` (optional)
    - [ ] `--ploidy`
    - [ ] `--trait`
    - [ ] `--model` (or multiple)
    - [ ] `--covariates`
    - [ ] `--min-maf`, `--max-missing`
    - [ ] `--threads` (conservative use)
    - [ ] `--out`
  - [ ] Runs a full GWAS scan using the gwas crate.
  - [ ] Writes result tables to `${out}.gwas.tsv` (or similar).

---

### 2.2 Non-Functional Requirements

#### 2.2.1 GWASpoly Parity

- [ ] Define at least two reference datasets:
  - [ ] A small toy dataset (for example ~100 individuals, ~1k SNPs).
  - [ ] A medium dataset representative of 5k tetraploid alfalfa (or a subset, for example 1k individuals, 10k SNPs).
- [ ] Provide R scripts that:
  - [ ] Run GWASpoly with specified options on these datasets.
  - [ ] Write:
    - [ ] Filtered marker list.
    - [ ] GWAS results (one table per trait/model).
- [ ] Provide Rust tests or integration harness that:
  - [ ] Run `binx gwas` on the same datasets and options.
  - [ ] Compare:
    - [ ] Marker IDs in output vs GWASpoly (no missing or extra markers).
    - [ ] Per-marker stats:
      - [ ] Maximum absolute difference in -log10(p) < `epsilon_p` (for example `1e-4`).
      - [ ] Maximum absolute difference in effect size < `epsilon_beta`.
  - [ ] Fail tests if tolerances are exceeded.
- [ ] Add a CI job that runs these parity tests.

#### 2.2.2 Determinism & Reproducibility

- [ ] GWAS runs must be deterministic given:
  - [ ] Same input files.
  - [ ] Same CLI options.
  - [ ] Same number of threads.
- [ ] Any stochastic components (if any) must accept a fixed `--seed` and use it.

#### 2.2.3 Performance Baseline

- [ ] Document a baseline performance target (non-binding but aspirational):
  - [ ] On the medium reference dataset, `binx gwas` should be at least 2x faster than GWASpoly in R on the same single core.
  - [ ] On multi-core (for example 8-16 threads), it should show near-linear scaling in the per-SNP scan stage.

---

## 3. `fastgwas` Crate Requirements (Enhanced Engine)

### 3.1 Functional Requirements

#### 3.1.1 Input / Output

- [ ] Use the same input semantics as `gwas` (reusing `binx-core` IO):
  - [ ] `--geno`, `--pheno`, `--kinship`, `--ploidy`, `--trait`, `--covariates`.
- [ ] Output:
  - [ ] Similar result table schema (marker-level stats) to `gwas` so downstream plotting/scripts can be shared.
  - [ ] May include additional columns (for example convergence diagnostics, approximation flags).

#### 3.1.2 Algorithms & Modes

- [ ] Implement at least two modes:
  - [ ] `exact`: Results should be numerically close to `gwas` (may use more efficient solvers, but not major approximations).
  - [ ] `approx`: Use faster algorithms (for example score tests, approximations) that may yield small deviations but improve runtime.
- [ ] Extend the model support to at least:
  - [ ] `additive` (v1).
  - [ ] Optionally `general` when feasible.
- [ ] Implement improved variance component estimation:
  - [ ] AI-REML or similar (details can evolve).
- [ ] Use aggressive multithreading for:
  - [ ] Per-marker scans.
  - [ ] Possibly per-trait or per-chromosome loops.

#### 3.1.3 Performance / Scalability

- [ ] Design target scale: 5k tetraploid individuals, 100k-1M SNPs.
- [ ] Requirements:
  - [ ] Must be able to process 5k x 100k SNPs on a modern 16-core node with <= 32-64 GB RAM.
  - [ ] Provide CLI controls:
    - [ ] `--threads N` to set number of threads.
    - [ ] `--chunk-size` or similar to cap memory usage.
- [ ] Performance target (non-binding but guiding):
  - [ ] In `approx` mode, aim for >= 5x faster than `gwas` on the same data and thread count, with similar conclusions (for example top hits, Manhattan shape).

---

### 3.2 Non-Functional Requirements

#### 3.2.1 Robustness & Diagnostics

- [ ] Provide clear logging for:
  - [ ] Convergence status for variance component estimation.
  - [ ] Any markers skipped or dropped due to numerical issues.
- [ ] If using approximations:
  - [ ] Document the trade-offs in the README.
  - [ ] Optionally provide a way to cross-check a subset of markers with the `gwas` engine for validation.

#### 3.2.2 Separation from `gwas`

- [ ] `fastgwas` must not change the behavior of `gwas`:
  - [ ] Shared code should live in `binx-core`.
  - [ ] `gwas` parity tests must pass independently of `fastgwas` changes.

---

## 4. CLI Requirements Summary

### 4.1 `binx gwas`

- [ ] Provides GWASpoly-like behavior.
- [ ] Minimal interface, focused on faithfulness and simplicity.

Example:

```bash
binx gwas \
  --geno alfalfa.geno.tsv \
  --pheno alfalfa.pheno.tsv \
  --ploidy 4 \
  --trait Height2024 \
  --model additive \
  --covariates PC1,PC2 \
  --min-maf 0.05 \
  --max-missing 0.1 \
  --threads 4 \
  --out alfalfa_additive
```

### 4.2 `binx fastgwas`

- [ ] Exposes performance and advanced options.

Example:

```bash
binx fastgwas \
  --geno alfalfa.geno.tsv \
  --pheno alfalfa.pheno.tsv \
  --ploidy 4 \
  --trait Height2024 \
  --model additive \
  --covariates PC1,PC2 \
  --mode approx \
  --threads 16 \
  --chunk-size 5000 \
  --out alfalfa_fastgwas
```

---

## 5. Milestones / Checklists

### 5.1 `gwas` MVP (GWASpoly-faithful)

- [ ] Read geno/pheno; basic ploidy-aware dosage matrix.
- [ ] Compute kinship equivalent to GWASpoly default.
- [ ] Fit null mixed model for one trait.
- [ ] Implement additive model per-marker tests.
- [ ] Implement QC filters (MAF, missingness).
- [ ] Implement CLI `binx gwas`.
- [ ] Create one toy dataset plus parity test against GWASpoly.
- [ ] CI job for parity test.

### 5.2 `gwas` v1.1

- [ ] Add `general` model support.
- [ ] Add LOCO support (if desired).
- [ ] Add second, medium-sized reference dataset for parity testing.

### 5.3 `fastgwas` MVP

- [ ] Reuse same IO, data models from `binx-core`.
- [ ] Implement additive model in `exact` mode using improved solvers.
- [ ] Implement multithreaded per-marker scanning.
- [ ] Provide `binx fastgwas` CLI with `--mode exact` and `--threads`.

### 5.4 `fastgwas` v1.1

- [ ] Add `approx` mode with faster algorithms.
- [ ] Add basic benchmarking script comparing `gwas`, `fastgwas`, and GWASpoly on 5k tetra alfalfa subset.
- [ ] Document performance and trade-offs.

---
