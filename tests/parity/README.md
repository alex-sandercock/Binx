# Parity Harness: GWASpoly vs binx gwas

This harness is the starter pack for proving `binx gwas` matches GWASpoly. It keeps R and Rust steps side-by-side and produces fixtures the Rust tests consume.

## Layout

- `tests/parity/data/toy/`
  - `toy.geno.tsv` — small dosage matrix (dense, GWASpoly-style).
  - `toy.pheno.tsv` — phenotype with a single trait and optional covariates.
  - `toy.kinship.tsv` — optional precomputed kinship (dense).
- `tests/parity/data/potato/`
  - `new_potato_geno.csv`, `new_potato_pheno.csv` (from GWASpoly repo examples; includes repeated IDs and `env` factor).
  - `new_potato_kinship.tsv` (optional precomputed).
- `tests/parity/fixtures/`
  - Output tables from GWASpoly runs, e.g. `toy_additive.tsv`, `potato_additive_env.tsv`. (Empty by default; populate via the R step.)
- `scripts/gwaspoly/run_parity.R`
  - R script that loads GWASpoly and writes the fixture tables.
- `binx-cli/tests/parity.rs`
  - Rust integration test that runs `binx gwas` and compares to the fixtures.

## R step: produce reference outputs

Script sketch (`scripts/gwaspoly/run_parity.R`):

```r
# Inputs: args or env vars
geno      <- Sys.getenv("GENO")
pheno     <- Sys.getenv("PHENO")
trait     <- Sys.getenv("TRAIT")
ploidy    <- as.integer(Sys.getenv("PLOIDY"))
out       <- Sys.getenv("OUT")
fixed     <- strsplit(Sys.getenv("FIXED"), ",")[[1]] # optional, e.g. "env"
fixedType <- strsplit(Sys.getenv("FIXED_TYPE"), ",")[[1]] # e.g. "factor"

library(GWASpoly)
data <- read.GWASpoly(pheno, geno, ploidy)
if (!is.na(fixed[1])) {
  params <- set.params(geno.freq = 1 - 5 / nrow(data@geno), fixed = fixed, fixed.type = fixedType)
} else {
  params <- set.params(geno.freq = 1 - 5 / nrow(data@geno))
}
data <- set.K(data) # or set.K.from.file if provided
data <- set.pheno(data, trait)
data <- GWASpoly(data, models = c("additive"), traits = trait, params = params)
res  <- get.results(data, model = "additive", trait = trait)
write.table(res, file = out, sep = "\t", quote = FALSE, row.names = FALSE)
```

Example invocations (run from repo root):

```bash
# toy, additive only
Rscript scripts/gwaspoly/run_parity.R \
  --geno=tests/parity/data/toy/toy.geno.tsv \
  --pheno=tests/parity/data/toy/toy.pheno.tsv \
  --trait=Trait1 \
  --ploidy=4 \
  --models=additive \
  --out=tests/parity/fixtures/toy_additive.tsv

# potato example with env fixed effect and repeated IDs
Rscript scripts/gwaspoly/run_parity.R \
  --geno=tests/parity/data/potato/new_potato_geno.csv \
  --pheno=tests/parity/data/potato/new_potato_pheno.csv \
  --trait=vine.maturity \
  --ploidy=4 \
  --fixed=env \
  --fixed_type=factor \
  --models=additive \
  --out=tests/parity/fixtures/potato_additive_env.tsv
```

Notes:
- Use GWASpoly defaults for QC thresholds unless explicitly set.
- Keep chromosome and position columns if available; if absent, marker ID suffices for matching.
- Capture `n.ind` and `n.obs` in the output if possible; otherwise record them in a sidecar YAML/JSON for the Rust test.

## Rust integration test

Implemented in `binx-cli/tests/parity.rs`:

- Runs `binx gwas` with the same inputs as the R call.
- Reads the GWASpoly fixture and the binx output.
- Asserts:
  - Marker IDs match exactly (after QC).
  - For matched markers, `max | -log10(p_binx) - -log10(p_ref) | < epsilon_p` (e.g. 1e-4).
  - `max | beta_binx - beta_ref | < epsilon_beta`.
  - `n_obs` (if present) matches (important for repeated IDs).
- Column names are normalized with common synonyms (GWASpoly `Marker`, `Effect`, `pValue`, etc.).

The tests auto-skip unless `BINX_PARITY=1` and fixtures + data files are present. Run locally with:

```bash
BINX_PARITY=1 cargo test -p binx-cli parity
```

The potato parity test is currently gated because `binx gwas` does not yet support repeated genotype IDs with an `env` fixed effect. Enable it (once supported) via:

```bash
BINX_PARITY=1 BINX_PARITY_POTATO=1 cargo test -p binx-cli parity -- --nocapture
```

## Expected columns (minimum)

- `marker_id`
- `chrom` (optional)
- `pos` (optional)
- `model`
- `n_obs`
- `effect` (or `effect_1`, `effect_2`, … for multi-parameter models)
- `se`
- `stat`
- `p_value`

Mapping from GWASpoly typical output:
- `Marker` -> `marker_id`
- `Chr` -> `chrom`
- `Pos` -> `pos`
- `nObs` -> `n_obs`
- `Effect`/`additive` -> `effect`
- `Score`/`Stat` -> `stat`
- `pValue` -> `p_value`

## Refreshing fixtures

1) Ensure GWASpoly is installed (R).  
2) Run the R script for each dataset to regenerate `tests/parity/fixtures/*.tsv`.  
3) Commit both the fixture updates and any harness code changes.  
4) CI job should run the Rust tests under `cargo test --tests parity` (or similar) and fail if tolerances are exceeded.

## TODOs to finish the harness

- Populate `tests/parity/fixtures/` by running `scripts/gwaspoly/run_parity.R` (toy data and GWASpoly example potato data are already placed under `tests/parity/data/`).
- Optional: enable the GitHub Actions workflow `.github/workflows/parity.yml` (manual trigger, installs R+GWASpoly) to run the parity tests in CI.
