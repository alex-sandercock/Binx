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

  | method   | Threads  | R (s)   | Rust (s)  | Rust Speedup |
  |----------|----------|---------|-----------|--------------|
  | non-LOCO | 8-10     | 63.90   | 10.72     | 6.0x         |
  | LOCO*     | 8-10     | 125.21  | 46.1      | 2.7x         |
  
  *LOCO has not been optimized in Rust for parallelization
  
  **rrblup-rs also will be optimized with faer in Rust, which
  will improve the LOCO performance

  --------------------------
  For the GWASpoly runs:
  - Rust is 6x faster than R when not using LOCO
  - Rust is ~2.7x faster than R when using LOCO, however, additional optimizations to rrblup-rs and improvements to the LOCO implementation should improve performance.
