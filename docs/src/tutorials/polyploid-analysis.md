# Working with Polyploids

This tutorial covers GWAS analysis in polyploid species using Binx's specialized genetic models.

## Introduction to Polyploid GWAS

Polyploid species (tetraploids, hexaploids, etc.) have more than two copies of each chromosome, which creates unique challenges and opportunities for GWAS:

- **More allele combinations**: A tetraploid has 5 possible genotypes per locus (0-4 copies)
- **Complex inheritance**: Dominance relationships are more nuanced
- **Higher genetic diversity**: More combinations can influence traits

Binx implements the GWASpoly framework, which models various forms of allele dosage effects.

## Genetic Models for Polyploids

### Understanding Dosage Effects

In a tetraploid, the five genotypes (AAAA, AAAB, AABB, ABBB, BBBB) can affect traits differently:

| Model | Assumption | Best For |
|-------|------------|----------|
| **Additive** | Linear dosage effect | Quantitative traits with dosage dependence |
| **General** | No assumption (4 df) | Unknown inheritance; hypothesis generation |
| **Simplex dominant** | One B allele is sufficient | Traits with low-dosage dominance |
| **Duplex dominant** | Two B alleles are sufficient | Intermediate dominance |

### Choosing Models

**Start with additive + general:**

```bash
binx gwas \
  --geno genotypes.tsv \
  --pheno phenotypes.csv \
  --trait yield \
  --ploidy 4 \
  --models additive,general \
  --out results.csv
```

- **Additive** captures dosage-dependent effects
- **General** captures any pattern (exploratory)

**Then investigate specific dominance patterns:**

```bash
binx gwas \
  --geno genotypes.tsv \
  --pheno phenotypes.csv \
  --trait disease_resistance \
  --ploidy 4 \
  --models additive,simplex-dom,duplex-dom \
  --out results.csv
```

## Example: Tetraploid Potato GWAS

Let's analyze a tetraploid potato dataset for tuber yield.

### Step 1: Verify Ploidy in Data

Check that dosage values are appropriate:

```bash
# Find max dosage value
awk -F'\t' 'NR>1 {for(i=4;i<=NF;i++) if($i>max) max=$i} END {print "Max dosage:", max}' genotypes.tsv
```

For tetraploid data, max should be 4.

### Step 2: Compute Polyploid Kinship

```bash
binx kinship \
  --geno genotypes.tsv \
  --ploidy 4 \
  --output kinship.tsv
```

### Step 3: Run Multi-Model GWAS

```bash
binx gwas \
  --geno genotypes.tsv \
  --pheno phenotypes.csv \
  --trait yield \
  --kinship kinship.tsv \
  --ploidy 4 \
  --models additive,general,simplex-dom,duplex-dom \
  --loco \
  --out gwas_results.csv
```

### Step 4: Compare Models

Extract results by model:

```bash
# Count significant hits per model (threshold -log10p > 5)
awk -F',' 'NR>1 && $8>5 {count[$4]++} END {for(m in count) print m, count[m]}' gwas_results.csv
```

Create model-specific Manhattan plots:

```bash
for model in additive general simplex-dom duplex-dom; do
  binx plot \
    --input gwas_results.csv \
    --plot-type manhattan \
    --model $model \
    --threshold 5 \
    --title "Yield GWAS - $model model" \
    --output manhattan_${model}.svg
done
```

### Step 5: Interpret Model-Specific Results

If a QTL is significant under:

- **Additive only**: Dosage-dependent effect (each additional allele adds to trait)
- **Simplex-dominant only**: Presence/absence effect (one copy is enough)
- **General but not additive**: Complex dominance pattern
- **Multiple models**: Robust association, exact inheritance unclear

## Hexaploid Analysis

For hexaploid species (ploidy=6), the same workflow applies:

```bash
binx gwas \
  --geno genotypes.tsv \
  --pheno phenotypes.csv \
  --trait yield \
  --ploidy 6 \
  --models additive,general \
  --out results.csv
```

Hexaploids have 7 possible dosage values (0-6) and even more complex dominance patterns.

## Tips for Polyploid GWAS

### Sample Size

Polyploids need larger sample sizes due to:
- More parameters in genetic models
- Lower power to detect effects
- **Recommendation**: 200+ samples for tetraploids

### MAF Filtering

Be careful with MAF filtering in polyploids:

```bash
# More lenient MAF for polyploids
binx gwas \
  --geno genotypes.tsv \
  --pheno phenotypes.csv \
  --trait yield \
  --ploidy 4 \
  --min-maf 0.02 \
  --out results.csv
```

Low-frequency variants in polyploids can still be informative.

### Interpreting Effect Sizes

Effect sizes in polyploid models:

- **Additive**: Effect per dosage unit
- **General model**: Reports multiple effects (one per genotype class)

To get standardized effects:

```bash
# Effects are in trait units per dosage
# Multiply by ploidy for total genetic range
```

## Diploidized Analysis

Sometimes you want to treat polyploid data as diploid-like:

```bash
binx gwas \
  --geno genotypes.tsv \
  --pheno phenotypes.csv \
  --trait yield \
  --ploidy 4 \
  --models diplo-add,diplo-dom \
  --out results.csv
```

This collapses dosage categories:
- 0, 1 → "AA-like"
- 2 → "AB-like"  
- 3, 4 → "BB-like"

Useful when expecting diploid-like inheritance in a polyploid.

## See Also

- [Genetic Models Reference](../reference/genetic-models.md) - Detailed model descriptions
- [Multi-Environment Trials](multi-environment.md) - Complex experimental designs
- [Validation](../reference/validation.md) - Accuracy compared to R/GWASpoly
