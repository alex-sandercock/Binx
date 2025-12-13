# binx qtl

Identify and filter significant QTLs from GWAS results.

## Synopsis

```bash
binx qtl [OPTIONS] --input <FILE>
```

## Description

The `qtl` command processes GWAS results to identify significant quantitative trait loci (QTLs). It applies significance thresholds and groups nearby markers into single QTL peaks.

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--input <FILE>` | Path to GWAS results file |

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--output <FILE>` | stdout | Output file path |
| `--threshold <FLOAT>` | - | -log10(p) significance threshold |
| `--bp-window <INT>` | 1000000 | Window size (bp) for grouping markers |
| `--model <MODEL>` | - | Filter to specific model |
| `--best-per-window` | false | Keep only best marker per window |

## Examples

### Basic QTL Extraction

```bash
binx qtl \
  --input gwas_results.csv \
  --threshold 5 \
  --output qtls.csv
```

### With Window-based Grouping

```bash
binx qtl \
  --input gwas_results.csv \
  --threshold 5 \
  --bp-window 10000000 \
  --best-per-window \
  --output qtls.csv
```

### Model-specific QTLs

```bash
binx qtl \
  --input gwas_results.csv \
  --model additive \
  --threshold 5 \
  --output additive_qtls.csv
```

## Output Format

```csv
marker_id,chrom,pos,model,effect,pvalue,log10p,qtl_id
SNP_1_1500,1,1500,additive,0.82,1.2e-08,7.92,QTL_1
SNP_3_8200,3,8200,additive,0.45,3.1e-06,5.51,QTL_2
```

The `qtl_id` column groups markers belonging to the same QTL region.

## QTL Grouping Algorithm

1. Sort markers by chromosome and position
2. Apply significance threshold
3. Group markers within `--bp-window` of each other
4. Optionally select best marker per group

## See Also

- [binx gwas](gwas.md) - Generate GWAS results
- [binx threshold](threshold.md) - Calculate thresholds
