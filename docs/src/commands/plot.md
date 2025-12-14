# binx plot

Generate publication-quality visualizations from GWAS results.

## Synopsis

```bash
binx plot --input <FILE> --output <FILE> [OPTIONS]
```

## Description

The `plot` command creates visualizations from GWAS output files, including Manhattan plots, QQ plots, and LD decay plots.

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--input <FILE>` | Input file: GWAS results CSV (manhattan/qq) or genotype TSV (ld) |
| `--output <FILE>` | Output file path (.svg or .png) |

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--plot-type <TYPE>` | manhattan | Type of plot: manhattan, qq, or ld |
| `--model <MODEL>` | - | Filter to specific genetic model |
| `--threshold <FLOAT>` | 5.0 | Significance threshold as -log10(p) |
| `--suggestive <FLOAT>` | 3.0 | Suggestive threshold as -log10(p) (0 to disable) |
| `--theme <THEME>` | classic | Visual theme |
| `--width <INT>` | 1200 | Plot width in pixels |
| `--height <INT>` | 600 | Plot height in pixels |
| `--title <TEXT>` | - | Plot title |
| `--chromosomes <LIST>` | - | Filter to specific chromosomes (comma-separated) |

> **Threshold Recommendation:** For accurate significance thresholds, use the value calculated by `binx gwas --threshold` or `binx threshold`. These commands compute thresholds using Bonferroni correction, M.eff (effective number of tests), or FDR methods appropriate for your dataset.

## Plot Types

### Manhattan Plot

Classic GWAS visualization showing -log10(p) across chromosomes:

```bash
binx plot \
  --input results.csv \
  --plot-type manhattan \
  --threshold 5 \
  --output manhattan.svg
```

> **Beta Feature:** When `--model` is not specified, all models from the results file are plotted together with different colors. This multi-model visualization is currently in beta.

### QQ Plot

Quantile-quantile plot for assessing genomic inflation:

```bash
binx plot \
  --input results.csv \
  --plot-type qq \
  --model additive \
  --output qq.svg
```

The plot includes:
- Expected vs observed -log10(p)
- 95% confidence band
- Diagonal reference line for visual inflation assessment

> **Beta Feature:** When `--model` is not specified, all models from the results file are plotted together with different colors. This multi-model visualization is currently in beta.

### LD Plot

Linkage disequilibrium decay over distance:

```bash
binx plot \
  --input geno.tsv \
  --plot-type ld \
  --ploidy 4 \
  --output ld_decay.svg
```

> **Important:** LD plots require **genotype data** (dosage matrix), not GWAS results. The input file must have the format: `marker_id, chr, pos, sample1, sample2, ...` where sample columns contain dosage values. If you accidentally use a GWAS results file, all r² values will appear as 1.0 because the statistical columns are misinterpreted as samples.

LD plot specific options:

| Option | Default | Description |
|--------|---------|-------------|
| `--ploidy <INT>` | - | Ploidy level (required for LD plot) |
| `--r2-threshold <FLOAT>` | - | R² threshold to mark on plot |
| `--max-pairs <INT>` | 10000 | Maximum marker pairs to sample |
| `--max-loci <INT>` | - | Maximum markers per chromosome |
| `--n-bins <INT>` | 50 | Number of distance bins for smoothing |

## Themes

| Theme | Description |
|-------|-------------|
| `classic` | Blue/orange alternating chromosomes (default) |
| `nature` | Muted gray tones for publication |
| `colorful` | Multi-color distinct chromosomes |
| `dark` | Dark background for presentations |
| `high_contrast` | High contrast for accessibility |

## Examples

### Basic Manhattan Plot

```bash
binx plot \
  --input gwas_results.csv \
  --plot-type manhattan \
  --output manhattan.svg
```

### Styled Manhattan Plot

```bash
binx plot \
  --input gwas_results.csv \
  --output manhattan.svg \
  --plot-type manhattan \
  --model additive \
  --threshold 7.3 \
  --theme nature \
  --title "Yield GWAS - Additive Model"
```

### QQ Plot for Model Comparison

```bash
# Generate QQ plots for each model
for model in additive general 1-dom-alt; do
  binx plot \
    --input gwas_results.csv \
    --output qq_${model}.svg \
    --plot-type qq \
    --model $model
done
```

### LD Plot with Threshold

```bash
binx plot \
  --input geno.tsv \
  --output ld.svg \
  --plot-type ld \
  --ploidy 4 \
  --r2-threshold 0.2
```

### LD Plot for Specific Chromosomes

```bash
binx plot \
  --input geno.tsv \
  --output ld.svg \
  --plot-type ld \
  --ploidy 4 \
  --chromosomes chr05,chr09
```

### Multi-panel Figure

```bash
# Create individual plots, then combine externally
binx plot --input results.csv --output manhattan.svg --plot-type manhattan
binx plot --input results.csv --output qq.svg --plot-type qq
```

## Output Formats

The output format is determined by file extension:

| Extension | Format |
|-----------|--------|
| `.svg` | Scalable Vector Graphics (recommended) |
| `.png` | PNG raster image |

## See Also

- [binx gwas](gwas.md) - Generate GWAS results
- [binx qtl](qtl.md) - Extract significant peaks
