# Output Formats

Reference for all output file formats produced by Binx.

## GWAS Results

Output from `binx gwas`:

| Column | Type | Description |
|--------|------|-------------|
| `marker_id` | string | Marker identifier |
| `chrom` | string | Chromosome |
| `pos` | integer | Base pair position |
| `model` | string | Genetic model used |
| `effect` | float | Effect size estimate |
| `stderr` | float | Standard error |
| `pvalue` | float | P-value |
| `log10p` | float | -log10(p-value) |
| `maf` | float | Minor allele frequency |
| `n` | integer | Sample size |

## Kinship Matrix

Output from `binx kinship`:

- Tab-separated values
- Square symmetric matrix
- Sample IDs as header and first column
- Values represent genetic relatedness (typically 0-2)

## QTL Results

Output from `binx qtl`:

| Column | Type | Description |
|--------|------|-------------|
| `marker_id` | string | Peak marker |
| `chrom` | string | Chromosome |
| `pos` | integer | Position |
| `model` | string | Best model |
| `effect` | float | Effect size |
| `pvalue` | float | P-value |
| `log10p` | float | -log10(p) |
| `qtl_id` | string | QTL group identifier |

## Plot Outputs

`binx plot` supports:
- `.svg` - Scalable Vector Graphics
- `.png` - PNG raster image
- `.pdf` - PDF document
