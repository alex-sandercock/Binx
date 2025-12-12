# binx-dosage

A Rust implementation of genotype dosage estimation using [updog](https://cran.r-project.org/package=updog)'s norm model for polyploid species. Based on R/updog version **2.1.6**.

> **This is a library crate.**
>
> For dosage estimation workflows, see the **Binx CLI** (`binx dosage`) for a complete
> command-line interface with VCF/CSV input and multiple output formats.

## Features

- `fit_norm()` - EM-based genotype calling with beta-binomial likelihood
- Multiple fitting modes: `Auto`, `Updog`, `UpdogExact`, `UpdogFast`, `Fast`, `Turbo`, `TurboAuto`
- Parallel processing with rayon
- Streaming I/O for memory-efficient processing of large VCF files
- Output formats: matrix, stats, BEAGLE, VCF, PLINK .raw, GWASpoly

## Usage

```rust
use binx_dosage::{run_norm_model, FitMode};
use ndarray::Array1;

let ref_counts = Array1::from(vec![5, 10, 0, 15]);
let total_counts = Array1::from(vec![10, 20, 5, 30]);
let ploidy = 4;

let result = run_norm_model(&ref_counts, &total_counts, ploidy, FitMode::Auto)?;
println!("Best genotypes: {:?}", result.best_genotypes);
println!("Bias: {}, Seq error: {}", result.bias, result.seq_error);
```

## Fitting Modes

| Mode | Description |
|------|-------------|
| `Auto` | Hybrid sprint: 3 bias starts, 8 iters each, pick winner (default) |
| `Updog` | Full multi-start: 5 starts to convergence (matches R/updog) |
| `UpdogExact` | Updog with relaxed bounds for exact R parity |
| `Fast` | Single start at bias=1.0 |
| `Turbo` | Single start with parallel E-step (fastest) |
| `TurboAuto` | 3-start sprint with parallel E-step |

## Citation

If you use `binx-dosage`, please cite the original updog paper:

> Gerard, D., Ferrão, L.F.V., Garcia, A.A.F., & Stephens, M. (2018). Genotyping polyploids from messy sequencing data. *Genetics* 210(3):789-807.

This crate is a Rust reimplementation of updog's `flexdog` norm model.
For a command-line interface, see [Binx](https://github.com/alex-sandercock/Binx).

## Related Crates

- **[gwaspoly-rs](../gwaspoly-rs)** - GWASpoly GWAS implementation
- **[rrblup-rs](../rrblup-rs)** - Mixed model solving
- **[binx-cli](../binx-cli)** - User-facing CLI (recommended for most users)

## License

GPL-3.0-or-later
