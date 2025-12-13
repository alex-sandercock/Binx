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
| `--method <METHOD>` | vanraden | Kinship method (vanraden, vanraden2) |
| `--output <FILE>` | stdout | Output file path |
| `--ploidy <INT>` | 2 | Ploidy level |
| `--min-maf <FLOAT>` | 0.0 | Minimum MAF filter |

## Methods

### VanRaden Method 1 (default)

The standard VanRaden (2008) method:

```
K = ZZ' / (2 * Σ p(1-p))
```

Where Z is the centered and scaled genotype matrix.

### VanRaden Method 2

Alternative scaling:

```
K = ZZ' / n
```

Where n is the number of markers.

## Examples

### Basic Usage

```bash
binx kinship --geno genotypes.tsv --output kinship.tsv
```

### With MAF Filter

```bash
binx kinship \
  --geno genotypes.tsv \
  --min-maf 0.05 \
  --output kinship.tsv
```

### For Tetraploids

```bash
binx kinship \
  --geno genotypes.tsv \
  --ploidy 4 \
  --output kinship.tsv
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

2. **MAF filtering**: Consider filtering low-MAF markers which can distort the kinship matrix

3. **Check diagonal values**: Diagonal values should be close to 1.0; much higher values may indicate inbreeding or data issues

## See Also

- [binx gwas](gwas.md) - Use kinship in GWAS
- [Input Formats](../getting-started/input-formats.md) - Kinship file format
