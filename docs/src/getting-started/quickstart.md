# Quick Start

This guide will walk you through your first Binx analysis in under 5 minutes.

## Overview

A typical Binx workflow looks like this:

```
VCF file → binx convert → Genotype file ─┐
                                         ├─→ binx gwas → Results → binx plot
                     Phenotype file ─────┘
```

## Step 1: Prepare Your Data

Binx requires two main input files:

### Genotype File

A tab-separated file with marker information and sample dosages:

```
marker_id   chrom   pos     Sample1 Sample2 Sample3
SNP001      1       1000    0       2       4
SNP002      1       2000    1       1       3
SNP003      2       1500    2       2       2
```

- First three columns: `marker_id`, `chrom`, `pos`
- Remaining columns: sample dosage values (0 to ploidy)

### Phenotype File

A CSV/TSV file with sample IDs and trait values:

```
sample_id,yield,height,env
Sample1,45.2,120,field_A
Sample2,52.1,115,field_A
Sample3,48.7,125,field_B
```

- First column: `sample_id` (must match genotype column headers)
- Remaining columns: traits and covariates

## Step 2: Convert VCF (if needed)

If your genotypes are in VCF format, convert them first:

```bash
binx convert \
  --vcf your_data.vcf.gz \
  --format gwaspoly \
  --output genotypes.tsv
```

## Step 3: Run GWAS

Run a genome-wide association study:

```bash
binx gwas \
  --geno genotypes.tsv \
  --pheno phenotypes.csv \
  --trait yield \
  --ploidy 4 \
  --models additive \
  --out gwas_results.csv
```

### Understanding the Parameters

| Parameter | Description |
|-----------|-------------|
| `--geno` | Path to genotype file |
| `--pheno` | Path to phenotype file |
| `--trait` | Column name of the trait to analyze |
| `--ploidy` | Ploidy level (2, 4, 6, etc.) |
| `--models` | Genetic models to test |
| `--out` | Output file path |

## Step 4: Examine Results

The output CSV contains association results for each marker:

```csv
marker_id,chrom,pos,model,effect,stderr,pvalue,log10p
SNP001,1,1000,additive,0.52,0.12,3.2e-05,4.49
SNP002,1,2000,additive,0.08,0.11,0.47,0.33
...
```

### Key columns:
- `effect`: Effect size estimate
- `pvalue`: Association p-value
- `log10p`: -log10 transformed p-value (useful for plotting)

## Step 5: Visualize Results

Create a Manhattan plot:

```bash
binx plot \
  --input gwas_results.csv \
  --plot-type manhattan \
  --model additive \
  --threshold 5 \
  --output manhattan.svg
```

Create a QQ plot:

```bash
binx plot \
  --input gwas_results.csv \
  --plot-type qq \
  --model additive \
  --output qq.svg
```

## Step 6: Identify QTLs

Extract significant QTLs:

```bash
binx qtl \
  --input gwas_results.csv \
  --bp-window 10000000 \
  --output significant_qtls.csv
```

## Complete Example Script

Here's a complete analysis pipeline you can adapt:

```bash
#!/bin/bash

# Define input files
VCF="data/samples.vcf.gz"
PHENO="data/phenotypes.csv"
TRAIT="yield"
PLOIDY=4

# Create output directory
mkdir -p results

# Step 1: Convert VCF to Binx format
binx convert \
  --vcf $VCF \
  --format gwaspoly \
  --output results/genotypes.tsv

# Step 2: Run GWAS with multiple models
binx gwas \
  --geno results/genotypes.tsv \
  --pheno $PHENO \
  --trait $TRAIT \
  --ploidy $PLOIDY \
  --models additive,general,simplex,duplex \
  --out results/gwas_results.csv

# Step 3: Generate plots
binx plot \
  --input results/gwas_results.csv \
  --plot-type manhattan \
  --model additive \
  --threshold 5 \
  --output results/manhattan.svg

binx plot \
  --input results/gwas_results.csv \
  --plot-type qq \
  --model additive \
  --output results/qq.svg

# Step 4: Extract QTLs
binx qtl \
  --input results/gwas_results.csv \
  --bp-window 10000000 \
  --output results/qtls.csv

echo "Analysis complete! Results in results/"
```

## Next Steps

- Learn about [Input Formats](input-formats.md) in detail
- Explore [Genetic Models](../reference/genetic-models.md) available in Binx
- Follow the [First GWAS Tutorial](../tutorials/first-gwas.md) for a more detailed walkthrough
- Check the [Command Reference](../commands/overview.md) for all options
