# binx kinship

Compute kinship (genomic relationship) matrices from genotype data.

## Synopsis

```bash
binx kinship --geno <FILE> --ploidy <INT> --out <FILE> [OPTIONS]
```

## Description

The `kinship` command computes a genomic relationship matrix (GRM) from marker dosage data. The kinship matrix captures genetic similarity between individuals and is used in GWAS to account for population structure and relatedness.

> **When to use:** While `binx gwas` auto-generates a kinship matrix if not provided (using gwaspoly-rs's `set_k()`), pre-computing with `binx kinship` is recommended when:
> - Running GWAS on multiple traits (avoids recomputation)
> - You need a specific kinship method (VanRaden vs GWASpoly)
> - You want to inspect or reuse the kinship matrix

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--geno <FILE>` | Path to genotype file (TSV/CSV with dosages) |
| `--ploidy <INT>` | Ploidy level (e.g., 2, 4, 6) |
| `--out <FILE>` | Output file path |

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--method <METHOD>` | vanraden | Kinship method: `vanraden` or `gwaspoly` |

## Methods

### VanRaden (default)

The standard VanRaden (2008) Method 1 additive relationship matrix, extended for polyploids:

```
K = M'M / (ploidy × Σ pq)
```

Where:
- M is the centered genotype matrix (markers × samples)
- Centering: dosage - (ploidy × p)
- p = allele frequency, q = 1-p

### GWASpoly

GWASpoly-style kinship matching R/GWASpoly's `set.K()` function:

```
K = MM' / mean(diag(K))
```

Where M is centered by column means and normalized to have unit diagonal mean.

## Examples

### Basic Usage (Tetraploid)

```bash
binx kinship \
  --geno genotypes.tsv \
  --ploidy 4 \
  --out kinship.tsv
```

### Using GWASpoly Method

```bash
binx kinship \
  --geno genotypes.tsv \
  --ploidy 4 \
  --method gwaspoly \
  --out kinship.tsv
```

### For Diploids

```bash
binx kinship \
  --geno genotypes.tsv \
  --ploidy 2 \
  --out kinship.tsv
```

## Output Format

A symmetric matrix with sample IDs as row and column headers:

```
sample_id	Sample1	Sample2	Sample3
Sample1	1.0000	0.2534	0.1256
Sample2	0.2534	1.0000	0.1892
Sample3	0.1256	0.1892	1.0000
```

## Tips

1. **Pre-compute for multiple traits**: Compute kinship once and reuse across multiple GWAS runs

2. **Check diagonal values**: Diagonal values should be close to 1.0; much higher values may indicate inbreeding or data issues

3. **Method selection**: Use `vanraden` (default) for standard GWAS, or `gwaspoly` for compatibility with R/GWASpoly workflows

## See Also

- [binx gwas](gwas.md) - Use kinship in GWAS
- [Input Formats](../getting-started/input-formats.md) - Kinship file format
