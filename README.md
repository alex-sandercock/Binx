
# Binx <img src="docs/assets/binx-logo.png" align="right" width="250"/>

Rust command-line genomics workbench for diploid and polyploid species. `binx` targets GWAS-style analyses with a familiar UX: fast defaults, explicit inputs, and clear TSV outputs.

-
-
-
-
## Highlights
- Biallelic GWAS (`binx gwas`) with linear models by default; switch to LMM when a kinship matrix is provided.
- Polyploid-aware: pass ploidy for genotype scaling; supports repeated phenotype IDs and optional environment filtering.
- Covariates and PCs handled as additional fixed effects.
- Kinship generation (`binx kinship`) via additive VanRaden on dosage matrices.
- Multiallelic GWAS is planned (`binx multigwas`) and currently stubbed with a clear error.

## Installation
- Requires a Rust toolchain (`cargo` + `rustc`).
- Build locally: `cargo build --release -p binx-cli` (binary at `target/release/binx`).
- Or install into `$CARGO_HOME/bin`: `cargo install --path binx-cli`.

## Quick start
```bash
# 1) Compute a VanRaden kinship from a biallelic dosage matrix.
cargo run -p binx-cli -- \
  kinship \
  --geno tests/parity/data/toy/toy.geno.tsv \
  --ploidy 4 \
  --out /tmp/toy_kinship.tsv

# 2) Run GWAS (LMM because we pass --kinship). Without --kinship it runs an LM.
cargo run -p binx-cli -- \
  gwas \
  --geno tests/parity/data/toy/toy.geno.tsv \
  --pheno tests/parity/data/toy/toy.pheno.tsv \
  --trait-name Trait1 \
  --ploidy 4 \
  --kinship /tmp/toy_kinship.tsv \
  --out /tmp/toy_gwas.tsv

# 3) Inspect results (TSV): marker_id, beta, se, t_stat, p_value.
head /tmp/toy_gwas.tsv
```

## Command reference
**binx gwas** — Biallelic GWAS (LM or LMM)
- Required: `--geno <dosage.tsv>`; `--pheno <phenotypes.tsv>`; `--trait-name <trait>`; `--ploidy <int>`; `--out <results.tsv>`.
- Optional: `--covariates <name1,name2>`; `--pcs <pcs.tsv>`; `--kinship <kinship.tsv>`; `--allow-missing-samples`; `--env-column <col>` + `--env-value <val>`; `--model <additive>` (only additive supported today).
- Output: tab-delimited results with `marker_id`, `beta`, `se`, `t_stat`, `p_value`.

**binx kinship** — Additive VanRaden kinship from biallelic dosages
- Required: `--geno <dosage.tsv>`; `--ploidy <int>`; `--out <kinship.tsv>`.
- Output: symmetric TSV with header `sample_id <S1> <S2> ...` and matching rows.

**binx multigwas** — Multiallelic GWAS (stub)
- Exists to reserve the interface; currently returns a “not implemented” error.

## Input formats
- **Genotypes (TSV/CSV auto-detected)**: first three columns `marker_id`, `chr`, `pos`, followed by sample dosage columns (`0..ploidy` for biallelic). Rows = markers, columns = samples.
- **Phenotypes (TSV/CSV auto-detected)**: first column `sample_id`, remaining columns are numeric traits/covariates. Factor covariates are kept as strings. Use `--env-column`/`--env-value` to subset rows (e.g., environment).
- **PCs (TSV)**: `sample_id`, `PC1`, `PC2`, … (header required).
- **Kinship (TSV/CSV auto-detected)**: header `sample_id <S1> <S2> ...`; each row starts with the row `sample_id` and contains the square matrix values.

## Notes and tips
- Pass `--kinship` to `binx gwas` to enable the mixed-model path; otherwise it runs an LM.
- Use `--allow-missing-samples` when phenotypes include IDs missing from the genotype table.
- The test harness under `tests/parity/` contains tiny toy datasets that can be used for local smoke tests.
