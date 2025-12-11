# gwaspoly-rs

A Rust implementation of [GWASpoly](https://github.com/jendelman/GWASpoly) for genome-wide association studies in polyploid species.

> **For CLI workflows, use [Binx](https://github.com/Breeding-Insight/Binx)**
>
> This crate provides a reference Rust implementation of GWASpoly as a library.
> For end users and CLI-based workflows, we recommend the **Binx CLI** (`binx gwas`),
> which wraps this crate and adds I/O handling, plotting, multi-method support, and more.

## Features

- Faithful port of R/GWASpoly algorithms
- Multiple genetic models: additive, general, dominance variants, diploidized
- LOCO (Leave-One-Chromosome-Out) support
- Parallel marker testing
- Threshold calculation (Bonferroni, M.eff, FDR)
- QTL detection with window-based pruning

## Usage as a Library

```rust
use gwaspoly_rs::{gwaspoly, GeneActionModel, ThresholdMethod};

// Run GWAS
gwaspoly(
    "genotypes.tsv",
    "phenotypes.csv",
    "yield",                           // trait name
    None,                              // covariates
    Some("kinship.tsv"),               // kinship matrix
    false,                             // allow_missing_samples
    None, None,                        // env_column, env_value (deprecated)
    4,                                 // ploidy
    &[GeneActionModel::Additive],      // models
    false,                             // LOCO
    0.05,                              // min_maf
    0.0,                               // max_geno_freq
    "results.csv",                     // output
    true,                              // parallel
    Some(ThresholdMethod::Meff),       // threshold method
    0.05,                              // alpha
)?;
```

## Benchmarks

Preliminary benchmarks vs R/GWASpoly on a potato dataset (958 samples, 9.8k SNPs, ploidy 4):

| Method   | Threads | R (s)  | Rust (s) | Speedup |
|----------|---------|--------|----------|---------|
| non-LOCO | 8-10    | 63.90  | 10.72    | **6.0x**    |
| LOCO     | 8-10    | 125.21 | 46.1     | **2.7x**    |

*LOCO parallelization is not yet fully optimized. Future improvements to rrblup-rs with faer will improve LOCO performance.*

## Related Crates

- **[rrblup-rs](../rrblup-rs)** - Mixed model solving (dependency)
- **[binx-gwas](../binx-gwas)** - GWAS orchestration layer
- **[binx-cli](../binx-cli)** - User-facing CLI (recommended for most users)

## License

GPL-3.0-or-later
