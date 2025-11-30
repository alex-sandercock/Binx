Preliminary benchmarks of Binx gwaspoly-rs vs R/GWASpoly show a massive speed improvement.

### Potato dataset
Dataset:
 - 958 Samples
 - 9.8k SNPs
 - Ploidy: 4

Analysis:
 - Two-environment (fixed effects: env)
 - No LOCO

Table 1.

  | method | R (s) | Rust (s) | Rust Speedup |
  |--------|--------|-----------|--------------|
  | non-LOCO    | 63.90 | 10.72     | 6.0x         |
  | xx      | xx    | xx        | xx           |

  Summary

  For the GWASpoly runs:
  - Rust is 6x faster than R when not using LOCO
