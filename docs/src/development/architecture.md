# Architecture

Binx is organized as a Cargo workspace with specialized crates.

## Crate Structure

```
Binx/
├── binx-cli/          # Main CLI binary
├── binx-gwas/         # GWAS orchestration
├── binx-kinship/      # Kinship computation
├── binx-dosage/       # Dosage estimation
├── binx-convert/      # File conversion
├── binx-plotting/     # Visualization
├── binx-types/        # Core data structures
├── binx-io/           # I/O utilities
├── gwaspoly-rs/       # GWASpoly implementation
└── rrblup-rs/         # rrBLUP implementation
```

## Crate Descriptions

| Crate | Description |
|-------|-------------|
| `binx-cli` | Command-line interface and argument parsing |
| `binx-gwas` | GWAS workflow orchestration |
| `binx-kinship` | VanRaden kinship matrix computation |
| `binx-dosage` | Updog-style dosage estimation |
| `binx-convert` | VCF parsing and format conversion |
| `binx-types` | Shared types: `GenotypeMatrix`, `Phenotype`, etc. |
| `binx-io` | File I/O helpers |
| `binx-plotting` | SVG/PNG plot generation |
| `gwaspoly-rs` | Core GWASpoly statistical methods |
| `rrblup-rs` | REML mixed model solver |

## Data Flow

```
Input Files
    │
    ▼
binx-io (parsing)
    │
    ▼
binx-types (data structures)
    │
    ├──► binx-kinship ──► Kinship Matrix
    │
    ▼
binx-gwas
    │
    ├──► gwaspoly-rs (genetic models)
    │
    └──► rrblup-rs (mixed model)
            │
            ▼
        Results CSV
            │
            ▼
        binx-plotting
            │
            ▼
        Plots (SVG/PNG)
```

## Key Dependencies

- `nalgebra`: Linear algebra
- `ndarray`: N-dimensional arrays
- `csv`: CSV parsing
- `clap`: CLI argument parsing
- `plotters`: Visualization

## See Also

- [Contributing](contributing.md) for development setup
- [GitHub Repository](https://github.com/alex-sandercock/Binx)
