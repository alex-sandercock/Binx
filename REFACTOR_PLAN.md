# Binx Crate Refactoring Plan

## Goal

Make `gwaspoly-rs` and `rrblup-rs` standalone libraries (publishable to crates.io) while `binx` serves as the opinionated, breeder-facing orchestration layer.

## Final Crate Structure

### Standalone Crates (no binx dependencies)

```
rrblup-rs/
├── a_mat()           # R/rrBLUP::A.mat
├── mixed_solve()     # R/rrBLUP::mixed.solve
├── kin_blup()        # R/rrBLUP::kin.blup
└── (own I/O if needed)

gwaspoly-rs/
├── read_gwaspoly()   # R/GWASpoly::read.GWASpoly
├── set_k()           # R/GWASpoly::set.K (GWASpoly-style kinship)
├── gwaspoly()        # R/GWASpoly::GWASpoly
├── set_threshold()   # R/GWASpoly::set.threshold
├── get_qtl()         # R/GWASpoly::get.QTL
└── depends on: rrblup-rs (for mixed model)
```

### Binx Ecosystem Crates

```
binx-types/           # Core data structures
├── GenotypeMatrixBiallelic
├── KinshipMatrix
├── PhenotypeTable
├── MarkerMetadata
└── PcMatrix

binx-io/              # File I/O beyond R packages
├── VCF parsing
├── PLINK format support
├── Extended CSV/TSV handling
└── Output format conversions

binx-gwas/            # GWAS orchestration
├── Sample/marker alignment
├── Multi-method comparison
├── LOCO coordination
└── Workflow management

binx-kinship/         # Kinship computation
├── compute_kinship_vanraden()  # True VanRaden (2008): ZZ'/(ploidy*Σpq)
├── compute_kinship_rrblup()    # Wrapper for rrblup-rs::a_mat
└── compute_kinship_gwaspoly()  # Wrapper for gwaspoly-rs::set_k

binx-dosage/          # Dosage estimation (keep as-is)
binx-plotting/        # Visualization (keep as-is)
binx-cli/             # User-facing CLI
```

## API Naming Convention (mirrors R)

Rust APIs mirror their R counterparts with underscores instead of dots:

| R Function | Rust Function | Crate |
|------------|---------------|-------|
| `read.GWASpoly()` | `read_gwaspoly()` | gwaspoly-rs |
| `set.K()` | `set_k()` | gwaspoly-rs |
| `GWASpoly()` | `gwaspoly()` | gwaspoly-rs |
| `set.threshold()` | `set_threshold()` | gwaspoly-rs |
| `get.QTL()` | `get_qtl()` | gwaspoly-rs |
| `A.mat()` | `a_mat()` | rrblup-rs |
| `mixed.solve()` | `mixed_solve_reml()` | rrblup-rs |
| `kin.blup()` | `kin_blup()` | rrblup-rs |

## Kinship Methods Summary

| Method | Formula | Location |
|--------|---------|----------|
| VanRaden (2008) | `ZZ' / (ploidy × Σpq)` | `binx-kinship` |
| rrBLUP A.mat | rrBLUP centering/scaling | `rrblup-rs` |
| GWASpoly set.K | `MM' / mean(diag)` | `gwaspoly-rs` |

## Dependency Flow

```
                    ┌─────────────────────────────────────┐
                    │         STANDALONE CRATES           │
                    │         (no binx deps)              │
                    │                                     │
                    │  rrblup-rs                          │
                    │    ├── a_mat                        │
                    │    ├── mixed_solve                  │
                    │    └── kin_blup                     │
                    │                                     │
                    │  gwaspoly-rs                        │
                    │    ├── read_gwaspoly                │
                    │    ├── set_k                        │
                    │    ├── gwaspoly                     │
                    │    └── depends on: rrblup-rs        │
                    └─────────────────────────────────────┘
                                     │
                                     ▼
                    ┌─────────────────────────────────────┐
                    │         BINX ECOSYSTEM              │
                    │                                     │
                    │  binx-types (data structures)       │
                    │         │                           │
                    │         ▼                           │
                    │  binx-io (VCF, PLINK, CSV)          │
                    │         │                           │
                    │         ▼                           │
                    │  binx-kinship (wraps methods)       │
                    │  binx-gwas (orchestration)          │
                    │  binx-dosage, binx-plotting         │
                    │         │                           │
                    │         ▼                           │
                    │  binx-cli (user interface)          │
                    └─────────────────────────────────────┘
```

## Refactoring Phases

### Phase 1: Create New Crates (No Breaking Changes) ✅ COMPLETED

1. ✅ Create `binx-io/` with I/O functions and types
2. ✅ Create `binx-gwas/` as placeholder with sample alignment
3. ✅ Run `cargo test --workspace` - all 103 tests pass

### Phase 2: Make Standalone Crates Independent

#### 2a. rrblup-rs ✅ COMPLETED
- ✅ Removed `binx-core` and `binx-kinship` dependencies from Cargo.toml
- ✅ Removed wrapper functions that used binx types from lib.rs
- ✅ Kept standalone modules: `a_mat`, `mixed_solve`, `mixed_solve_fast`, `kin_blup`
- ✅ Renamed `mixed_solve_new` → `mixed_solve_reml` for clarity
- ✅ All 57 rrblup-rs tests pass (including R comparison tests)

#### 2b. gwaspoly-rs ✅ COMPLETED
- ✅ Created `types.rs` module with standalone types (GenotypeMatrixBiallelic, KinshipMatrix, PhenotypeTable, MarkerMetadata)
- ✅ Created `io.rs` module with `read_gwaspoly()` function (mirrors R's read.GWASpoly)
- ✅ Renamed functions to mirror R API:
  - `compute_kinship_all_markers()` → `set_k()` (mirrors R's set.K)
  - `run_gwaspoly()` → `gwaspoly()` (mirrors R's GWASpoly)
  - `calculate_thresholds()` → `set_threshold()` (mirrors R's set.threshold)
- ✅ Removed `binx-core` and `binx-kinship` dependencies from Cargo.toml
- ✅ Kept only `rrblup-rs` dependency (mirrors R architecture where GWASpoly wraps rrBLUP)
- ✅ All 26 gwaspoly-rs tests pass (24 unit + 2 integration)
- ✅ Updated binx-cli to use new API names
- ✅ All 101+ workspace tests pass

### Phase 3: Update Binx Ecosystem ✅ COMPLETED

**Key Achievement:** gwaspoly-rs and rrblup-rs are now fully standalone and can be published to crates.io independently.

**Type Architecture Decision:**
- Standalone crates (gwaspoly-rs, rrblup-rs) have their own types for zero dependencies
- Binx ecosystem crates (binx-core, binx-kinship, etc.) have their own compatible types
- This allows standalone crates to be used outside of binx while binx can orchestrate multiple tools

**Changes Made:**
1. ✅ binx-cli uses gwaspoly_rs directly for GWAS operations (gwaspoly, set_threshold, load_genotypes)
2. ✅ binx-kinship keeps its own kinship implementations (VanRaden, GWASpoly-style) for non-GWAS workflows
3. ✅ All 101+ workspace tests pass

**Note:** binx-kinship's `compute_kinship_gwaspoly` and gwaspoly-rs's `set_k` have identical algorithms. The duplication is intentional - binx-kinship uses binx-types while gwaspoly-rs is standalone.

### Phase 4: Cleanup ✅ COMPLETED

1. ✅ Renamed `binx-core` to `binx-types`
   - Directory renamed
   - Package name updated in Cargo.toml
   - All dependent crates updated (binx-multigwas, binx-plotting, binx-cli, binx-kinship)
   - All source imports updated (binx_core → binx_types)
2. ✅ All 101+ workspace tests pass
3. ✅ Final validation complete

## Test Coverage (Safety Net)

| Test File | Validates |
|-----------|-----------|
| `rrblup-rs/tests/r_comparison.rs` | `mixed_solve` vs R (Vu, Ve, beta, u, LL) |
| `rrblup-rs/tests/amat_comparison.rs` | `a_mat` vs R (matrix values, MAF, shrinkage) |
| `rrblup-rs/tests/kinblup_comparison.rs` | `kin_blup` vs R (Vg, Ve, BLUPs, PEV) |
| `gwaspoly-rs/tests/parallel_accuracy.rs` | Parallel vs sequential marker testing |
| `binx-cli/tests/parity.rs` | End-to-end GWAS vs R/GWASpoly |

## Key Principles

1. **Tests must pass at every step** - if any R comparison or parity test fails, stop and investigate
2. **No loss of accuracy** - numerical results must match R packages within tolerance
3. **No loss of functionality** - all current features remain available
4. **Incremental migration** - duplicate first, then migrate, then clean up

## What Makes Binx Special

With standalone method crates, Binx provides:

- **Opinionated workflows**: `binx init`, `binx check`, `binx gwas`, `binx predict`
- **Unified data model**: `binx.toml` project files, standardized output folders
- **Multi-tool orchestration**: Run multiple methods, compare results, consensus tables
- **Reproducibility**: Run metadata, version tracking, "Methods" text generation
- **HPC integration**: Chunking, checkpointing, SLURM/LSF support
- **Ecosystem integration**: Backend for Weaver, Mandala, BIGapp GUIs

The crates are "replaceable engines." Binx is the platform.
