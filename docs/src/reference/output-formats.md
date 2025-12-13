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
| `score` | float | -log10(p-value) |
| `p_value` | float | Association p-value |
| `effect` | float | Effect size estimate |
| `n_obs` | integer | Sample size (non-missing) |
| `threshold` | float | Significance threshold used |

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
| `score` | float | -log10(p-value) |
| `effect` | float | Effect size |
| `threshold` | float | Significance threshold used |

## Plot Outputs

`binx plot` supports:
- `.svg` - Scalable Vector Graphics
- `.png` - PNG raster image
- `.pdf` - PDF document
