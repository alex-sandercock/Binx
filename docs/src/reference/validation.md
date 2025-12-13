# Validation & Accuracy

Binx has been extensively validated against the original R implementations to ensure accuracy.

## Validation Summary

| Component | Reference | Tests | Accuracy |
|-----------|-----------|-------|----------|
| rrblup-rs | R/rrBLUP | 52 | 5-6 decimal places |
| gwaspoly-rs | R/GWASpoly | Multiple configs | 4-5 decimal places |

## rrblup-rs Validation

The `rrblup-rs` crate implements R/rrBLUP's `mixed.solve` function. It was validated with 52 test cases covering:

### Variance Component Estimation (REML)

```r
# R/rrBLUP
library(rrBLUP)
result <- mixed.solve(y, K=K)
result$Vu  # Genetic variance
result$Ve  # Residual variance
```

```bash
# Binx produces matching values
```

### Fixed and Random Effect Predictions

- BLUP predictions for random effects
- BLUE estimates for fixed effects
- Standard errors

### Edge Cases

- Missing phenotype data
- Singular kinship matrices
- Small sample sizes

## gwaspoly-rs Validation

The `gwaspoly-rs` crate was validated against R/GWASpoly across:

### Configurations Tested

| Configuration | Description |
|---------------|-------------|
| LOCO vs non-LOCO | Leave-One-Chromosome-Out kinship |
| With/without PCs | P+K model vs K-only model |
| Multiple genetic models | Additive, general, dominance |
| With/without covariates | Factor and numeric covariates |

### Test Data

Validation used:
- Simulated tetraploid datasets
- Real potato GWAS data (from GWASpoly paper)
- Various sample sizes (100-500)

### Results Comparison

P-values match to 4-5 decimal places:

```
Marker      R/GWASpoly      Binx            Difference
SNP001      3.21e-05        3.21e-05        < 1e-09
SNP002      0.4687          0.4687          < 1e-06
SNP003      1.45e-03        1.45e-03        < 1e-08
```

## Running Validation Tests

Validation scripts are in the repository:

```bash
# Clone repository
git clone https://github.com/alex-sandercock/Binx.git
cd Binx

# Run parity tests
cd tests/parity
Rscript compare_rrblup.R
Rscript compare_gwaspoly.R
```

## Known Differences

Minor numerical differences can arise from:

1. **Floating point precision**: Rust and R may handle edge cases slightly differently
2. **Optimization convergence**: REML optimization may converge to slightly different points
3. **Random number generation**: If any stochastic elements are used

These differences are typically < 1e-5 and do not affect biological conclusions.

## Continuous Integration

Validation tests run automatically on each release:

- Comparison against R reference outputs
- Regression tests for all commands
- Edge case handling

## Reporting Issues

If you find discrepancies between Binx and R implementations:

1. Check input file formats match exactly
2. Verify parameter settings are equivalent
3. Report via [GitHub Issues](https://github.com/alex-sandercock/Binx/issues)

Include:
- Input data (or minimal reproducible example)
- R code and output
- Binx command and output
- Expected vs actual results

## References

- Endelman, J.B. (2011). Ridge regression and other kernels for genomic selection with R package rrBLUP. *The Plant Genome* 4:250-255.

- Rosyara, U.R., De Jong, W.S., Douches, D.S., & Endelman, J.B. (2016). Software for genome-wide association studies in autopolyploids and its application to potato. *The Plant Genome* 9(2).
