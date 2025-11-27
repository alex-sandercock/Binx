# GWASpoly math + behavior notes

Quick notes from inspecting GWASpoly R code (jendelman/GWASpoly) to guide Rust parity.

## Core pipeline

- `read.GWASpoly(ploidy, pheno.file, geno.file, format, n.traits, delim)`: loads phenotypes (can include repeated genotype IDs and fixed-effect columns) and genotype dosages. Phenotypes keep all rows; genotype IDs are per individual.
- `set.K(data, LOCO=FALSE)`: computes kinship. LOCO optional (not used in fixtures). Requires an explicit `LOCO` flag; defaults to error if missing.
- `set.params(...)`: returns a **list**, not an S4 slot object. Key fields used in GWASpoly body:
  - `fixed`, `fixed.type` (factor/numeric)
  - `min.MAF`, `max.geno.freq`, `P3D` (defaults TRUE), `n.PC`
  - When neither `min.MAF` nor `max.geno.freq` is set, GWASpoly sets `max.geno.freq = 1 - 5/N` (N = genotypes) and `min.MAF = 0`.
- `GWASpoly(data, models, traits, params, n.core=1, quiet=F)`:
  - Validates inputs, fills defaults as above.
  - Builds per-trait model:
    - `not.miss`: non-missing observations for the trait (observation-level).
    - `y`: trait vector (length `n_obs`).
    - `pheno.gid`: genotype IDs per observation.
    - `Z`: incidence matrix (`n_obs x n_genotypes`), with 1 linking each observation row to its genotype.
    - `X`: intercept; augments with fixed effects from `params$fixed` (factor -> dummy, numeric -> numeric). Adds PCs if `n.PC > 0` (from eigenvectors of K or LOCO).
    - `X2 <- .make.full(X)` enforces full rank.
  - P3D path (default): fits variance components once via `mixed.solve` to obtain `Hinv` (the inverse of `V = Z K Z^T + I * sigma_e^2`); if LOCO is requested, computes chromosome-specific `Hinv`.
  - For each model and chromosome, calls `.score.calc(...)`.
  - Returns `GWASpoly.fitted` with slots: `map`, `pheno`, `fixed`, `geno`, `ploidy`, `K`, `scores`, `effects`, `params`.

## Marker encoding `.design.score`

Source: `getAnywhere(.design.score)` in GWASpoly.

```
.design.score <- function(Mi, model, ploidy, min.MAF, max.geno.freq) {
  n <- length(Mi)
  freq <- mean(Mi, na.rm = TRUE) / ploidy
  if (min(freq, 1 - freq) < min.MAF) return(NULL)
  if (model == "additive") {
    geno.freq <- table(round(Mi)) / n
    if (max(geno.freq) <= max.geno.freq) return(matrix(Mi)) else return(NULL)
  }
  if (model == "general") {
    Mi <- round(Mi)
    geno.freq <- table(Mi) / n
    if (max(geno.freq) <= max.geno.freq) {
      tmp <- model.matrix(~x, data.frame(x=factor(Mi)))[, -1]
      return(if (is.null(dim(tmp))) matrix(tmp) else tmp)
    } else return(NULL)
  }
  # other models (diplo*, dominance) follow similar pattern
}
```

Key points for parity:
- MAF uses mean dosage / ploidy vs `min.MAF`.
- `max.geno.freq` is the maximum genotype-class frequency; defaults to `1 - 5/N`.
- Additive model uses **raw dosages (not centered)** as the single predictor column.
- General model dummy-codes rounded dosages (reference level dropped).

## Score test `.score.calc`

Source: `getAnywhere(.score.calc)` in GWASpoly (trimmed summary):

```
.score.calc(geno, y, Z, X, K, Hinv, ploidy, model, min.MAF, max.geno.freq)
  P3D if Hinv provided, otherwise fits each marker with mixed.solve
  For each marker:
    S <- .design.score(Mi, ...)
    X2 <- cbind(X, Z %*% S)
    p <- ncol(X2); v1 <- ncol(S); v2 <- nrow(Z) - p
    if (!P3D) Hinv <- mixed.solve(..., return.Hinv=TRUE)$Hinv
    W <- t(X2) %*% Hinv %*% X2
    Winv <- solve(W)
    beta <- Winv %*% t(X2) %*% Hinv %*% y
    resid <- y - X2 %*% beta
    s2 <- (t(resid) %*% Hinv %*% resid) / v2
    Q <- s2 * Winv[last v1 cols/rows]
    Tt <- solve(Q); V <- beta[last v1]
    Fstat <- t(V) %*% Tt %*% V / v1
    score <- -log10(pbeta(v2/(v2+v1*Fstat), v2/2, v1/2))
    effect (beta) reported only for non-general models.
```

Key parity requirements:
- Uses a score/F-test style with `Hinv` from P3D null fit.
- Degrees of freedom: `v2 = n_obs - p` (observation-level), `v1 = ncol(S)` (marker parameters).
- P-values are **-log10(pbeta(...))** in the score test; GWASpoly reports this in `scores`. Our fixture converter converts back to `p_value = 10^(-score)`.
- For additive, effect is the last coefficient in `beta`; for `general`, effects are not populated (NA).

## Fixed effects and repeated measures

- Phenotypes may have repeated genotype IDs (multi-environment). Each row is an observation.
- `Z` incidence matrix maps observations to genotype-level random effects.
- Fixed effects: `params$fixed` names columns in `data@fixed`; `fixed.type` controls factor vs numeric; factors are dummy-coded with reference level dropped.
- PCs: if `n.PC > 0`, derived from eigenvectors of K (or LOCO K), appended to X.

## Defaults and thresholds

- When neither `min.MAF` nor `max.geno.freq` is provided, GWASpoly uses `max.geno.freq = 1 - 5/N` and `min.MAF = 0`.
- Missingness: enforced via `max.geno.freq` (monomorphic/high-frequency filter). No explicit per-marker missingness in `.design.score` (but `min.MAF` and `max.geno.freq` gate markers).
- P3D is default TRUE; if set to FALSE, requires LOCO=TRUE in `set.K`.

## Implications for Rust parity (binx-gwas)

- Add observation-level support: allow repeated phenotype IDs, build `Z` (n_obs x n_ind), and fit null mixed model with `V = Z K Z^T + I * sigma_e^2`.
- Implement P3D path: fit variance components once per trait (and per chromosome if LOCO) using REML with `Hinv`.
- Use marker design per `.design.score`: additive = raw dosages (rounded only in MAF/genotype checks), general = dummy-coded rounded dosages.
- Filtering: implement `min.MAF` and `max.geno.freq = 1 - 5/N` default; enforce marker exclusion if either fails.
- Score test: match GWASpoly equations above; report -log10(p) or convert to p-value consistently.
- Fixed effects: support factor vs numeric with dummy-coding; accept `--covariates` and factor typing equivalent to `fixed`/`fixed.type`.
- LOCO: optional; when enabled, compute chromosome-specific K (not yet in Rust).

These notes are intended as a reference while bringing `binx-gwas` to GWASpoly parity (including repeated-ID phenotypes and env fixed effects).***
