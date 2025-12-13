# binx gwas

Perform genome-wide association studies (GWAS) using GWASpoly-style methods with support for multiple genetic models.

## Synopsis

```bash
binx gwas [OPTIONS] --geno <FILE> --pheno <FILE> --trait <NAME>
```

## Description

The `gwas` command performs association analysis between genetic markers and phenotypic traits. It implements the statistical methods from GWASpoly (Rosyara et al., 2016) and rrBLUP (Endelman, 2011), supporting both diploid and polyploid species.

Key features:
- Multiple genetic models for polyploid analysis
- Mixed model framework (K model, P+K model)
- Leave-One-Chromosome-Out (LOCO) kinship
- Support for covariates and multi-environment trials

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--geno <FILE>` | Path to genotype file (TSV/CSV with dosages) |
| `--pheno <FILE>` | Path to phenotype file (TSV/CSV) |
| `--trait <NAME>` | Name of the trait column to analyze |

## Options

### Analysis Options

| Option | Default | Description |
|--------|---------|-------------|
| `--ploidy <INT>` | 2 | Ploidy level (2, 4, 6, etc.) |
| `--models <LIST>` | additive | Genetic models to test (comma-separated) |
| `--kinship <FILE>` | - | Pre-computed kinship matrix |
| `--loco` | false | Use Leave-One-Chromosome-Out kinship |
| `--n-pc <INT>` | 0 | Number of principal components to include |
| `--covariates <LIST>` | - | Covariate column names (comma-separated) |

### Output Options

| Option | Default | Description |
|--------|---------|-------------|
| `--out <FILE>` | stdout | Output file path |
| `--format <FMT>` | csv | Output format (csv, tsv) |

### Filtering Options

| Option | Default | Description |
|--------|---------|-------------|
| `--min-maf <FLOAT>` | 0.0 | Minimum minor allele frequency |
| `--max-missing <FLOAT>` | 1.0 | Maximum missing rate per marker |

## Genetic Models

The `--models` option accepts the following values. These match R/GWASpoly's gene action models (Rosyara et al., 2016).

### For Diploids (ploidy=2)

| Model | Description | Encoding | df |
|-------|-------------|----------|-----|
| `additive` | Linear dosage effect | 0, 1, 2 | 1 |
| `general` | Separate effect per dosage | dummy coded | 2 |
| `1-dom` | Dominant (tests both ref and alt) | — | 1 each |
| `1-dom-ref` | Dominant (ref group distinct) | 0, 1, 1 | 1 |
| `1-dom-alt` | Dominant (alt group distinct) | 0, 0, 1 | 1 |

### For Tetraploids (ploidy=4)

| Model | Description | Encoding | df |
|-------|-------------|----------|-----|
| `additive` | Linear dosage effect | 0, 1, 2, 3, 4 | 1 |
| `general` | Separate effect per dosage | dummy coded | 4 |
| `1-dom` | Simplex dominant (tests both ref and alt) | — | 1 each |
| `1-dom-ref` | Simplex dominant (ref group distinct) | 0, 1, 1, 1, 1 | 1 |
| `1-dom-alt` | Simplex dominant (alt group distinct) | 0, 0, 0, 0, 1 | 1 |
| `2-dom` | Duplex dominant (tests both ref and alt) | — | 1 each |
| `2-dom-ref` | Duplex dominant (ref side distinct) | 0, 0, 1, 1, 1 | 1 |
| `2-dom-alt` | Duplex dominant (alt side distinct) | 0, 0, 0, 1, 1 | 1 |
| `diplo-general` | Diploidized general (hets collapsed) | dummy coded | 2 |
| `diplo-additive` | Diploidized additive (hets = 0.5) | 0, 0.5, 0.5, 0.5, 1 | 1 |

### Model Expansion

Like R/GWASpoly, specifying `1-dom` automatically tests both `1-dom-ref` and `1-dom-alt`. Similarly, `2-dom` expands to both `2-dom-ref` and `2-dom-alt`. Use the specific `-ref` or `-alt` variants if you only want one direction.

### Using Multiple Models

Specify multiple models separated by commas:

```bash
binx gwas --models additive,general,1-dom,2-dom ...
```

See [Genetic Models Reference](../reference/genetic-models.md) for detailed explanations.

## Examples

### Basic GWAS

```bash
binx gwas \
  --geno genotypes.tsv \
  --pheno phenotypes.csv \
  --trait yield \
  --ploidy 4 \
  --out results.csv
```

### With Pre-computed Kinship

```bash
# First compute kinship
binx kinship --geno genotypes.tsv --out kinship.tsv

# Then run GWAS
binx gwas \
  --geno genotypes.tsv \
  --pheno phenotypes.csv \
  --trait yield \
  --kinship kinship.tsv \
  --ploidy 4 \
  --out results.csv
```

### LOCO Analysis

Leave-One-Chromosome-Out reduces proximal contamination:

```bash
binx gwas \
  --geno genotypes.tsv \
  --pheno phenotypes.csv \
  --trait yield \
  --ploidy 4 \
  --loco \
  --out results.csv
```

### With Covariates and PCs

```bash
binx gwas \
  --geno genotypes.tsv \
  --pheno phenotypes.csv \
  --trait yield \
  --covariates environment,block \
  --n-pc 3 \
  --ploidy 4 \
  --out results.csv
```

### Multiple Genetic Models

```bash
binx gwas \
  --geno genotypes.tsv \
  --pheno phenotypes.csv \
  --trait yield \
  --ploidy 4 \
  --models additive,general,1-dom,2-dom \
  --out results.csv
```

## Output Format

The output file contains the following columns:

| Column | Description |
|--------|-------------|
| `marker_id` | Marker identifier |
| `chrom` | Chromosome |
| `pos` | Base pair position |
| `model` | Genetic model used |
| `score` | -log10(p-value) |
| `p_value` | Association p-value |
| `effect` | Effect size estimate |
| `n_obs` | Sample size (non-missing) |
| `threshold` | Significance threshold used |

### Example Output

```csv
marker_id,chrom,pos,model,score,p_value,effect,n_obs,threshold
SNP_1_1000,1,1000,additive,4.49,3.21e-05,0.523,198,5.0
SNP_1_2000,1,2000,additive,0.33,0.469,0.081,200,5.0
SNP_1_3500,1,3500,additive,2.84,1.45e-03,-0.312,195,5.0
```

## Statistical Details

### Mixed Model

The GWAS uses a linear mixed model:

```
y = Xβ + Zu + e
```

Where:
- `y` = phenotype vector
- `X` = fixed effects design matrix (intercept, covariates, marker)
- `β` = fixed effects
- `Z` = random effects design matrix
- `u ~ N(0, Kσ²ᵤ)` = random polygenic effects
- `K` = kinship matrix
- `e ~ N(0, Iσ²ₑ)` = residual errors

### P+K Model

When `--n-pc` is specified, principal components are included as fixed effects to account for population structure (P+K model).

### LOCO

With `--loco`, the kinship matrix is recalculated for each chromosome, excluding markers on the chromosome being tested. This prevents the tested marker from influencing its own significance through the kinship matrix.

## Tips and Best Practices

1. **Choose appropriate models**: For autopolyploids, start with `additive` and `general`. The general model captures complex dominance patterns but uses more degrees of freedom.

2. **Use LOCO for accurate p-values**: LOCO prevents proximal contamination and generally provides better-calibrated p-values.

3. **Pre-compute kinship for efficiency**: If running multiple traits, compute the kinship matrix once and reuse it.

4. **Filter markers**: Use `--min-maf` to remove rare variants that have low power.

5. **Check QQ plots**: After analysis, generate QQ plots to assess genomic inflation.

## See Also

- [binx kinship](kinship.md) - Compute kinship matrices
- [binx plot](plot.md) - Visualize GWAS results
- [binx qtl](qtl.md) - Extract significant QTLs
- [Genetic Models Reference](../reference/genetic-models.md)
