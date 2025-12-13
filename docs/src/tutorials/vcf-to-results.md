# From VCF to Results

A complete pipeline tutorial showing how to go from raw VCF data to publication-ready GWAS results.

## Overview

This tutorial demonstrates the full Binx workflow:

1. Convert VCF to Binx format
2. Quality control and filtering
3. Compute kinship matrix
4. Run GWAS with multiple models
5. Generate visualizations
6. Extract and interpret QTLs

## Complete Pipeline Script

```bash
#!/bin/bash
set -e

# === Configuration ===
VCF="raw_data/variants.vcf.gz"
PHENO="raw_data/phenotypes.csv"
TRAIT="yield"
PLOIDY=4
OUTDIR="results"

mkdir -p $OUTDIR

# === Step 1: Convert VCF ===
echo "Converting VCF..."
binx convert \
  --vcf $VCF \
  --format gwaspoly \
  --min-maf 0.05 \
  --max-missing 0.20 \
  --output $OUTDIR/genotypes.tsv

# === Step 2: Compute Kinship ===
echo "Computing kinship..."
binx kinship \
  --geno $OUTDIR/genotypes.tsv \
  --ploidy $PLOIDY \
  --output $OUTDIR/kinship.tsv

# === Step 3: Run GWAS ===
echo "Running GWAS..."
binx gwas \
  --geno $OUTDIR/genotypes.tsv \
  --pheno $PHENO \
  --trait $TRAIT \
  --kinship $OUTDIR/kinship.tsv \
  --ploidy $PLOIDY \
  --models additive,general \
  --loco \
  --out $OUTDIR/gwas_results.csv

# === Step 4: Visualize ===
echo "Creating plots..."
binx plot \
  --input $OUTDIR/gwas_results.csv \
  --plot-type manhattan \
  --model additive \
  --threshold 5 \
  --output $OUTDIR/manhattan.svg

binx plot \
  --input $OUTDIR/gwas_results.csv \
  --plot-type qq \
  --model additive \
  --output $OUTDIR/qq.svg

# === Step 5: Extract QTLs ===
echo "Extracting QTLs..."
binx qtl \
  --input $OUTDIR/gwas_results.csv \
  --threshold 5 \
  --bp-window 5000000 \
  --output $OUTDIR/qtls.csv

echo "Done! Results in $OUTDIR/"
```

## See Also

- [Your First GWAS](first-gwas.md)
- [Command Reference](../commands/overview.md)
