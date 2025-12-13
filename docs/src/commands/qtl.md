# binx qtl

Identify and filter significant QTLs from GWAS results.

## Synopsis

```bash
binx qtl [OPTIONS]
```

## Description

The `qtl` command processes GWAS results to identify significant quantitative trait loci (QTLs). It filters markers where `score >= threshold` and optionally prunes nearby signals within a specified window.

> **Important:** The input file must contain a `threshold` column. Use `binx gwas --threshold` to generate results with thresholds, or `binx threshold` to calculate thresholds separately.

> **Tip:** If your results file doesn't have a threshold column, see [Adding Thresholds to Existing Results](threshold.md#adding-thresholds-to-existing-results) for instructions on how to add one.

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--input <FILE>` | stdin | Input GWAS results file |
| `--output <FILE>` | stdout | Output file path |
| `--bp-window <INT>` | - | Prune signals within this window (bp) |

## Examples

### Basic QTL Extraction

```bash
binx qtl \
  --input gwas_results.csv \
  --output qtls.csv
```

### With Window-based Pruning

Prune nearby signals within a 1 Mb window:

```bash
binx qtl \
  --input gwas_results.csv \
  --bp-window 1000000 \
  --output qtls.csv
```

### Pipeline from GWAS

Pipe directly from GWAS with threshold calculation:

```bash
binx gwas \
  --geno geno.tsv \
  --pheno pheno.csv \
  --trait yield \
  --ploidy 4 \
  --threshold m.eff \
  --out /dev/stdout 2>/dev/null | \
binx qtl --bp-window 1000000 --output qtls.csv
```

### Reading from stdin

```bash
cat gwas_results.csv | binx qtl --bp-window 1000000
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

## Algorithm

The QTL detection algorithm (matching R/GWASpoly's `get.QTL`):

1. **Filter significant markers**: Keep only markers where `score >= threshold`
2. **Group by model**: Process each genetic model separately
3. **Sort by significance**: Order markers by score (descending)
4. **Window-based pruning** (if `--bp-window` specified):
   - For each chromosome, iterate through markers from most to least significant
   - Keep a marker only if it's more than `bp-window` away from all previously retained markers
   - This ensures the most significant marker in each region is retained

## See Also

- [binx gwas](gwas.md) - Generate GWAS results
- [binx threshold](threshold.md) - Calculate thresholds
