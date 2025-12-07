# gwaspoly-rs TODO

## Score discrepancy with R/GWASpoly (2024-12-06)

When comparing binx gwaspoly scores to R/GWASpoly on the potato tutorial dataset, some markers show score differences:

| Marker | R Score | Binx Score | Difference |
|--------|---------|------------|------------|
| solcap_snp_c2_23032 | 5.38 | 4.45 | -0.93 |
| solcap_snp_c2_48597 | 5.42 | 4.93 | -0.49 |

These differences cause 2 borderline-significant QTLs to be missed (8 in R vs 6 in binx).

### Possible causes to investigate:
- [ ] `env` fixed effect handling (params had `fixed = "env", fixed.type = "factor"`)
- [ ] Kinship matrix computation differences
- [ ] Z matrix handling for repeated observations
- [ ] Numerical precision in mixed model solving

### Test case:
```r
# R/GWASpoly
params <- set.params(geno.freq = 1 - 5/N, fixed = "env", fixed.type = "factor")
data.original.scan <- GWASpoly(data.original, models=c("additive"),
                               traits=c("vine.maturity"), params=params, n.core=10)
```

```bash
# binx
binx gwaspoly --geno potato_geno.tsv --pheno potato_pheno.csv \
    --trait_name vine.maturity --ploidy 4 --models additive \
    --threshold m.eff --out results.csv
```

## LD Plot curve smoothness (2025-12-06)

**Problem:** LD decay curve is not as smooth as GWASpoly's output. Current implementation uses binning + isotonic regression + cubic spline interpolation, but results don't match GWASpoly quality.

### To investigate:
- [ ] Benchmark r² threshold distance against GWASpoly output
- [ ] Compare curve shape visually with GWASpoly LD.plot
- [ ] Evaluate if scam (shape-constrained additive model) approach is needed
- [ ] Check if `n_bins` parameter still affects results inappropriately
