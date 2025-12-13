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

| Column | Description |
|--------|-------------|
| `marker_id` | Peak marker identifier |
| `chrom` | Chromosome |
| `pos` | Base pair position |
| `model` | Genetic model used |
| `score` | -log10(p-value) |
| `effect` | Effect size estimate |
| `threshold` | Significance threshold used |

### Example Output

```csv
marker_id,chrom,pos,model,score,effect,threshold
SNP_1_1500,1,1500,additive,7.92,0.82,5.0
SNP_3_8200,3,8200,additive,5.51,0.45,5.0
```

## QTL Grouping Algorithm

1. Sort markers by chromosome and position
2. Apply significance threshold
3. Group markers within `--bp-window` of each other
4. Optionally select best marker per group

## See Also

- [binx gwas](gwas.md) - Generate GWAS results
- [binx threshold](threshold.md) - Calculate thresholds
