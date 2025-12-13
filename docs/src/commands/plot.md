# binx plot

Generate publication-quality visualizations from GWAS results.

## Synopsis

```bash
binx plot [OPTIONS] --input <FILE> --plot-type <TYPE>
```

## Description

The `plot` command creates visualizations from GWAS output files, including Manhattan plots, QQ plots, and LD decay plots.

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--input <FILE>` | Path to GWAS results file |
| `--plot-type <TYPE>` | Type of plot to generate |

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--output <FILE>` | plot.svg | Output file path |
| `--model <MODEL>` | - | Filter to specific genetic model |
| `--threshold <FLOAT>` | - | -log10(p) significance threshold line |
| `--theme <THEME>` | default | Visual theme |
| `--width <INT>` | 1200 | Plot width in pixels |
| `--height <INT>` | 600 | Plot height in pixels |
| `--title <TEXT>` | - | Plot title |

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

Options specific to Manhattan plots:
- `--highlight <FILE>`: File with markers to highlight
- `--annotate-top <INT>`: Label top N peaks
- `--colors <LIST>`: Alternating chromosome colors

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
- Genomic inflation factor (λ)
- 95% confidence band

### LD Decay Plot

Linkage disequilibrium decay over distance:

```bash
binx plot \
  --input ld_results.csv \
  --plot-type ld-decay \
  --output ld_decay.svg
```

## Themes

| Theme | Description |
|-------|-------------|
| `default` | Clean, minimal style |
| `classic` | Traditional publication style |
| `dark` | Dark background |
| `colorblind` | Colorblind-friendly palette |

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
  --plot-type manhattan \
  --model additive \
  --threshold 7.3 \
  --theme classic \
  --title "Yield GWAS - Additive Model" \
  --annotate-top 5 \
  --output manhattan.svg
```

### QQ Plot for Model Comparison

```bash
# Generate QQ plots for each model
for model in additive general simplex-dom; do
  binx plot \
    --input gwas_results.csv \
    --plot-type qq \
    --model $model \
    --output qq_${model}.svg
done
```

### Multi-panel Figure

```bash
# Create individual plots, then combine externally
binx plot --input results.csv --plot-type manhattan --output manhattan.svg
binx plot --input results.csv --plot-type qq --output qq.svg
```

## Output Formats

The output format is determined by file extension:

| Extension | Format |
|-----------|--------|
| `.svg` | Scalable Vector Graphics (default) |
| `.png` | PNG raster image |
| `.pdf` | PDF document |

## See Also

- [binx gwas](gwas.md) - Generate GWAS results
- [binx qtl](qtl.md) - Extract significant peaks
