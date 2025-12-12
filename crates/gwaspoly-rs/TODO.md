# gwaspoly-rs TODO

## LD Plot curve smoothness (2025-12-06)

**Problem:** LD decay curve is not as smooth as GWASpoly's output. Current implementation uses binning + isotonic regression + cubic spline interpolation, but results don't match GWASpoly quality.

### To investigate:
- [ ] Benchmark r² threshold distance against GWASpoly output
- [ ] Compare curve shape visually with GWASpoly LD.plot
- [ ] Evaluate if scam (shape-constrained additive model) approach is needed
- [ ] Check if `n_bins` parameter still affects results inappropriately

## Crate extraction

These crates are "engines/methods" - CLI stays in Binx.

### rrblup-rs
- [ ] Extract rrblup-rs as a standalone crate/GitHub repo
- [ ] Should remain faithful to R/rrBLUP originals
- [ ] Binx will depend on this crate; new features/optimizations stay on Binx side

### gwaspoly-rs
- [ ] Extract gwaspoly-rs as a standalone crate/GitHub repo
- [ ] Should remain faithful to R/GWASpoly originals
- [ ] Binx will depend on this crate; new features/optimizations stay on Binx side
