# binx threshold

Calculate significance thresholds for GWAS results.

## Synopsis

```bash
binx threshold [OPTIONS] --results <FILE>
```

## Description

The `threshold` command calculates significance thresholds for GWAS using various multiple testing correction methods. It accounts for the number of tests performed and correlation structure among markers.

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--results <FILE>` | Path to GWAS results file |

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--method <METHOD>` | bonferroni | Correction method |
| `--alpha <FLOAT>` | 0.05 | Desired significance level |
| `--output <FILE>` | stdout | Output file |

## Methods

### Bonferroni

The most conservative approach:

```
threshold = α / n_tests
```

```bash
binx threshold \
  --results gwas.csv \
  --method bonferroni \
  --alpha 0.05
```

### Effective Number of Tests (M.eff)

Accounts for LD between markers:

```bash
binx threshold \
  --results gwas.csv \
  --method meff \
  --alpha 0.05
```

### FDR (Benjamini-Hochberg)

False Discovery Rate control:

```bash
binx threshold \
  --results gwas.csv \
  --method fdr \
  --alpha 0.05
```

### Permutation

Empirical threshold via permutation:

```bash
binx threshold \
  --results gwas.csv \
  --method permutation \
  --n-perm 1000 \
  --alpha 0.05
```

## Output

```
Method: Bonferroni
Alpha: 0.05
Number of tests: 50000
P-value threshold: 1.00e-06
-log10(p) threshold: 6.00
```

## Examples

### Compare Methods

```bash
for method in bonferroni meff fdr; do
  echo "=== $method ==="
  binx threshold --results gwas.csv --method $method --alpha 0.05
done
```

### Use Threshold in QTL Extraction

```bash
# Calculate threshold
THRESHOLD=$(binx threshold --results gwas.csv --method bonferroni | grep log10 | awk '{print $NF}')

# Apply to QTL extraction
binx qtl --input gwas.csv --threshold $THRESHOLD --output qtls.csv
```

## See Also

- [binx qtl](qtl.md) - Apply thresholds to extract QTLs
- [binx gwas](gwas.md) - Generate GWAS results
