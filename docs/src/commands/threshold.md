# binx threshold

Calculate significance thresholds for GWAS results.

## Synopsis

```bash
binx threshold --results <FILE> --method <METHOD> [OPTIONS]
```

## Description

The `threshold` command calculates significance thresholds for GWAS using various multiple testing correction methods. It accounts for the number of tests performed and optionally the correlation structure among markers.

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--results <FILE>` | GWAS results CSV file |
| `--method <METHOD>` | Threshold method: bonferroni, m.eff, or fdr |

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--alpha <FLOAT>` | 0.05 | Significance level |
| `--geno <FILE>` | - | Genotype file (required for m.eff) |
| `--ploidy <INT>` | - | Ploidy level (required for m.eff) |

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

Accounts for LD between markers using the method of Moskvina & Schmidt (2008). Requires genotype data to calculate marker correlations.

```bash
binx threshold \
  --results gwas.csv \
  --method m.eff \
  --geno geno.tsv \
  --ploidy 4 \
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
# Bonferroni and FDR (don't require genotype data)
for method in bonferroni fdr; do
  echo "=== $method ==="
  binx threshold --results gwas.csv --method $method --alpha 0.05
done

# M.eff requires genotype data
binx threshold --results gwas.csv --method m.eff --geno geno.tsv --ploidy 4
```

### Integrated GWAS Workflow (Recommended)

The recommended approach is to use `binx gwas --threshold` to calculate thresholds during GWAS, which adds the threshold to each result row for use with `binx qtl`:

```bash
# Run GWAS with threshold calculation
binx gwas \
  --geno geno.tsv \
  --pheno pheno.csv \
  --trait yield \
  --ploidy 4 \
  --threshold m.eff \
  --out results.csv

# Extract QTLs (uses threshold column from results)
binx qtl --input results.csv --bp-window 1000000 --output qtls.csv
```

### Adding Thresholds to Existing Results

If you have GWAS results without thresholds, you can calculate them separately and add them to the results file.

**Step 1: Calculate thresholds**

```bash
binx threshold --results gwas_results.csv --method bonferroni
```

Output:
```
Thresholds (Bonferroni):
Model                   Threshold       M.eff   n_markers
------------------------------------------------------------
additive                     5.30           -        1000
general                      5.60           -        1000
```

**Step 2: Add threshold column matching by model**

```bash
# Define thresholds per model (from Step 1 output)
awk -F',' -v OFS=',' '
BEGIN {
    thresh["additive"] = 5.30
    thresh["general"] = 5.60
}
NR==1 { print $0",threshold" }
NR>1 {
    model = $4
    t = (model in thresh) ? thresh[model] : "NA"
    print $0","t
}' gwas_results.csv > gwas_with_threshold.csv
```

**Step 3: Extract QTLs**

```bash
binx qtl --input gwas_with_threshold.csv --bp-window 1000000 --output qtls.csv
```

## See Also

- [binx qtl](qtl.md) - Apply thresholds to extract QTLs
- [binx gwas](gwas.md) - Generate GWAS results
