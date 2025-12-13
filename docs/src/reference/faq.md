# Frequently Asked Questions

## General

### What species can I analyze with Binx?

Binx works with any diploid or polyploid species. It has been validated with potato (tetraploid), but works equally well for wheat (hexaploid), strawberry (octoploid), or any diploid crop.

### How does Binx compare to PLINK/GCTA/TASSEL?

Binx focuses specifically on polyploid GWAS and implements the GWASpoly framework. For diploid-only analysis, tools like PLINK may be more efficient. Binx's strength is its native polyploid support.

### Why Rust?

Rust provides memory safety and performance comparable to C/C++, while being easier to maintain. This makes Binx fast and reliable.

## Data & Input

### What's the maximum dataset size?

Binx has been tested with:
- 10,000+ samples
- 1,000,000+ markers
- Memory usage scales with sample count

### Can I use imputed genotypes?

Yes, dosages can be fractional (e.g., 1.7 instead of 2) for imputed data.

### Why do my sample IDs not match?

Sample IDs are case-sensitive and must match exactly. Check for:
- Leading/trailing spaces
- Different naming conventions
- Underscore vs hyphen

## Analysis

### Which genetic model should I use?

Start with `additive` and `general`. The general model catches complex patterns, while additive has more power for linear effects. See [Genetic Models](genetic-models.md).

### What causes QQ plot inflation?

Common causes:
- Population structure (add PCs with `--n-pc`)
- Cryptic relatedness (use kinship matrix)
- Technical artifacts (check data quality)

### How do I handle related individuals?

Include a kinship matrix computed with `binx kinship`. The mixed model accounts for relatedness.

## Troubleshooting

### "Out of memory" error

Try:
- Process chromosomes separately
- Use a machine with more RAM
- Filter markers more stringently

### Results don't match R/GWASpoly

Check:
- Same input data format
- Same model parameters
- Same kinship matrix

Minor differences (<1e-5) are expected due to numerical precision.

## See Also

- [GitHub Issues](https://github.com/alex-sandercock/Binx/issues) for bug reports
- [Tutorials](../tutorials/first-gwas.md) for worked examples
