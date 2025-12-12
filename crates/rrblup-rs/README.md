# rrblup-rs

A Rust implementation of [rrBLUP](https://cran.r-project.org/package=rrBLUP) for mixed model solving and genomic prediction. Based on R/rrBLUP version **4.6.3**.

> **This is a library crate.**
>
> For GWAS workflows, see [gwaspoly-rs](../gwaspoly-rs) (which uses this crate) or the
> **Binx CLI** (`binx gwas`) for a complete command-line interface.

## Features

- `mixed_solve()` - REML-based mixed model solver
- `a_mat()` - Additive relationship matrix computation
- `kin_blup()` - Genomic BLUP with kinship matrix
- Faithful port of R/rrBLUP algorithms

## Usage

```rust
use rrblup_rs::{mixed_solve_reml, MixedSolveOptions};
use nalgebra::DMatrix;

let y = vec![1.0, 2.0, 3.0, 4.0, 5.0];
let x = DMatrix::from_row_slice(5, 1, &[1.0, 1.0, 1.0, 1.0, 1.0]); // intercept

let result = mixed_solve_reml(&y, Some(&x), None, MixedSolveOptions::default())?;
println!("Vu: {}, Ve: {}", result.vu, result.ve);
```

## Benchmarks

### mixed.solve with marker matrix Z

**Varying markers (n=500 genotypes):**

| Markers (m) | R (ms) | Rust (ms) | Speedup |
|-------------|--------|-----------|---------|
| 1,000       | 1,014  | 532       | **1.9x**    |
| 5,000       | 2,447  | 706       | **3.5x**    |
| 10,000      | 4,222  | 1,000     | **4.2x**    |
| 20,000      | 7,841  | 1,604     | **4.9x**    |

**Varying genotypes (m=10,000 markers):**

| Genotypes (n) | R (ms) | Rust (ms) | Speedup |
|---------------|--------|-----------|---------|
| 100           | 159    | 31        | **5.1x**    |
| 200           | 638    | 117       | **5.5x**    |
| 500           | 4,237  | 936       | **4.5x**    |
| 1,000         | 19,937 | 5,960     | **3.3x**    |

### kin.blup

| n_geno | R (ms) | Rust (ms) | Speedup |
|--------|--------|-----------|---------|
| 100    | 93     | 54        | **1.7x**    |
| 200    | 704    | 474       | **1.5x**    |

**Summary:** Rust is 3-5x faster for marker-based workflows, with advantages growing as problem size increases.

## Citation

If you use `rrblup-rs`, please cite the original R/rrBLUP paper:

> Endelman, J.B. (2011). Ridge regression and other kernels for genomic selection with R package rrBLUP. *The Plant Genome* 4:250-255.

This crate is a Rust reimplementation of `mixed.solve` and related functionality.
For a command-line interface, see [Binx](https://github.com/alex-sandercock/Binx).

## Related Crates

- **[gwaspoly-rs](../gwaspoly-rs)** - GWASpoly implementation (uses this crate)
- **[binx-cli](../binx-cli)** - User-facing CLI for GWAS workflows

## License

GPL-3.0-or-later
