# Your First GWAS Analysis

This tutorial walks you through a complete GWAS analysis from start to finish using Binx.

## What You'll Learn

- How to prepare your data for Binx
- Running a basic GWAS analysis
- Interpreting and visualizing results
- Extracting significant QTLs

## Prerequisites

- Binx installed ([Installation Guide](../getting-started/installation.md))
- Genotype data (VCF or dosage format)
- Phenotype data (CSV/TSV)

## Sample Data

For this tutorial, we'll use a simulated tetraploid potato dataset with:
- 200 samples
- 10,000 SNP markers
- 1 quantitative trait (tuber yield)

You can download sample data from the [Binx GitHub repository](https://github.com/alex-sandercock/Binx/tree/main/tests).

## Step 1: Examine Your Data

First, let's look at our input files:

```bash
# Check genotype file structure
head -3 genotypes.tsv
```

```
marker_id	chrom	pos	Sample001	Sample002	Sample003	...
SNP_1_1000	1	1000	0	2	4	...
SNP_1_2500	1	2500	1	1	3	...
```

```bash
# Check phenotype file
head -5 phenotypes.csv
```

```csv
sample_id,yield,environment
Sample001,45.2,field_A
Sample002,52.1,field_A
Sample003,48.7,field_B
...
```

Verify sample counts match:

```bash
# Count samples in genotype file (columns - 3)
head -1 genotypes.tsv | awk -F'\t' '{print NF-3, "samples"}'

# Count samples in phenotype file (lines - 1)
wc -l < phenotypes.csv | awk '{print $1-1, "samples"}'
```

## Step 2: Compute Kinship Matrix

The kinship matrix captures genetic relationships. Computing it separately allows reuse across multiple traits:

```bash
binx kinship \
  --geno genotypes.tsv \
  --ploidy 4 \
  --output kinship.tsv
```

Check the kinship matrix:

```bash
# View corner of matrix
head -5 kinship.tsv | cut -f1-5
```

Diagonal values should be approximately 1.0. Off-diagonal values represent relatedness between samples.

## Step 3: Run GWAS

Now run the association analysis:

```bash
binx gwas \
  --geno genotypes.tsv \
  --pheno phenotypes.csv \
  --trait yield \
  --kinship kinship.tsv \
  --ploidy 4 \
  --models additive \
  --out gwas_results.csv
```

This will:
1. Load genotypes and phenotypes
2. Match samples between files
3. Fit a mixed model for each marker
4. Output association statistics

### Understanding the Output

```bash
head gwas_results.csv
```

```csv
marker_id,chrom,pos,model,effect,stderr,pvalue,log10p,maf,n
SNP_1_1000,1,1000,additive,0.123,0.089,0.167,0.78,0.32,195
SNP_1_2500,1,2500,additive,0.521,0.102,3.2e-07,6.49,0.28,198
...
```

Key columns:
- `effect`: How much the trait changes per allele dosage unit
- `pvalue`: Probability of seeing this effect by chance
- `log10p`: -log10 transformed p-value (higher = more significant)

## Step 4: Calculate Significance Threshold

Determine the significance threshold:

```bash
binx threshold \
  --results gwas_results.csv \
  --method bonferroni \
  --alpha 0.05
```

Output:
```
Method: Bonferroni
Number of tests: 10000
P-value threshold: 5.00e-06
-log10(p) threshold: 5.30
```

## Step 5: Create Visualizations

### Manhattan Plot

```bash
binx plot \
  --input gwas_results.csv \
  --plot-type manhattan \
  --model additive \
  --threshold 5.3 \
  --title "Tuber Yield GWAS" \
  --output manhattan.svg
```

![Manhattan Plot Example](../assets/manhattan_example.png)

The Manhattan plot shows:
- X-axis: Genomic position (by chromosome)
- Y-axis: -log10(p-value)
- Red line: Significance threshold
- Peaks above the line are significant associations

### QQ Plot

```bash
binx plot \
  --input gwas_results.csv \
  --plot-type qq \
  --model additive \
  --output qq.svg
```

![QQ Plot Example](../assets/qq_example.png)

A good QQ plot shows:
- Points following the diagonal line (no inflation)
- Deviation at the tail (true associations)
- Lambda (λ) close to 1.0

## Step 6: Extract Significant QTLs

Identify significant loci:

```bash
binx qtl \
  --input gwas_results.csv \
  --threshold 5.3 \
  --bp-window 5000000 \
  --best-per-window \
  --output significant_qtls.csv
```

```bash
cat significant_qtls.csv
```

```csv
marker_id,chrom,pos,model,effect,pvalue,log10p,qtl_id
SNP_3_15234000,3,15234000,additive,0.82,1.2e-08,7.92,QTL_1
SNP_7_8234000,7,8234000,additive,0.45,3.1e-06,5.51,QTL_2
```

## Step 7: Interpret Results

For each significant QTL:

1. **Effect size**: A positive effect means the alternate allele increases the trait
2. **Position**: Look up genes near the QTL position
3. **MAF**: Very rare variants may be false positives

### Candidate Gene Analysis

Once you have QTL positions, you can:
- Look up nearby genes in genome browsers
- Check if known candidate genes are in the region
- Examine the LD block around the peak marker

## Complete Script

Here's the full analysis as a script:

```bash
#!/bin/bash
set -e

# Configuration
GENO="genotypes.tsv"
PHENO="phenotypes.csv"
TRAIT="yield"
PLOIDY=4
OUTDIR="results"

# Create output directory
mkdir -p $OUTDIR

# Step 1: Compute kinship
echo "Computing kinship matrix..."
binx kinship --geno $GENO --ploidy $PLOIDY --output $OUTDIR/kinship.tsv

# Step 2: Run GWAS
echo "Running GWAS..."
binx gwas \
  --geno $GENO \
  --pheno $PHENO \
  --trait $TRAIT \
  --kinship $OUTDIR/kinship.tsv \
  --ploidy $PLOIDY \
  --models additive \
  --out $OUTDIR/gwas_results.csv

# Step 3: Calculate threshold
echo "Calculating threshold..."
binx threshold --results $OUTDIR/gwas_results.csv --method bonferroni --alpha 0.05

# Step 4: Generate plots
echo "Creating plots..."
binx plot --input $OUTDIR/gwas_results.csv --plot-type manhattan --threshold 5.3 --output $OUTDIR/manhattan.svg
binx plot --input $OUTDIR/gwas_results.csv --plot-type qq --output $OUTDIR/qq.svg

# Step 5: Extract QTLs
echo "Extracting QTLs..."
binx qtl --input $OUTDIR/gwas_results.csv --threshold 5.3 --bp-window 5000000 --output $OUTDIR/qtls.csv

echo "Done! Results in $OUTDIR/"
```

## Next Steps

- Try different [genetic models](../reference/genetic-models.md)
- Use [LOCO](gwas.md#loco-analysis) for better p-value calibration
- Analyze [multiple environments](multi-environment.md)
- Explore [polyploid-specific models](polyploid-analysis.md)

## Troubleshooting

### "Sample ID mismatch" error

Ensure sample IDs in phenotype file exactly match genotype column headers (case-sensitive).

### Inflated QQ plot (λ > 1.1)

- Try including principal components (`--n-pc 5`)
- Check for population structure in your data
- Use LOCO kinship (`--loco`)

### No significant results

- Check if trait is heritable
- Ensure sufficient sample size (>100 recommended)
- Try different genetic models
