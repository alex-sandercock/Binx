# binx kinship

Compute kinship (genomic relationship) matrices from genotype data.

## Synopsis

```bash
binx kinship [OPTIONS] --geno <FILE>
```

## Description

The `kinship` command computes a genomic relationship matrix (GRM) from marker dosage data. The kinship matrix captures genetic similarity between individuals and is used in GWAS to account for population structure and relatedness.

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--geno <FILE>` | Path to genotype file (TSV/CSV with dosages) |

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--method <METHOD>` | vanraden | Kinship method (vanraden, gwaspoly) |
| `--out <FILE>` | - | Output file path (required) |
| `--ploidy <INT>` | 2 | Ploidy level |

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

### Basic Usage

```bash
binx kinship --geno genotypes.tsv --out kinship.tsv
```

### Using GWASpoly Method

```bash
binx kinship \
  --geno genotypes.tsv \
  --method gwaspoly \
  --out kinship.tsv
```

### For Tetraploids

```bash
binx kinship \
  --geno genotypes.tsv \
  --ploidy 4 \
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
