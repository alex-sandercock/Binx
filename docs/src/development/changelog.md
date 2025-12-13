# Changelog

All notable changes to Binx are documented here.

## [0.1.0] - 2025-12-12

### Added

- Initial release
- `binx gwas` - GWASpoly-style GWAS with multiple genetic models
- `binx kinship` - VanRaden kinship matrix computation
- `binx dosage` - Genotype dosage estimation
- `binx convert` - VCF to Binx format conversion
- `binx plot` - Manhattan and QQ plot generation
- `binx qtl` - QTL extraction from GWAS results
- `binx threshold` - Significance threshold calculation

### Validated

- rrblup-rs validated against R/rrBLUP (52 test cases, 5-6 decimal accuracy)
- gwaspoly-rs validated against R/GWASpoly (4-5 decimal accuracy)

---

## Version Format

Binx follows [Semantic Versioning](https://semver.org/):

- MAJOR: Incompatible API/format changes
- MINOR: New features (backward compatible)
- PATCH: Bug fixes (backward compatible)
