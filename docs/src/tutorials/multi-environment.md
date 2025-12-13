# Multi-Environment Trials

This tutorial covers GWAS analysis with repeated measurements across multiple environments.

## Overview

Multi-environment trials (MET) are common in plant breeding where the same genotypes are evaluated across multiple locations, years, or treatments. Binx handles these designs through:

- Repeated sample IDs in phenotype files
- Environment as a covariate
- Appropriate mixed model handling

## Data Structure

Your phenotype file can include repeated measurements:

```csv
sample_id,yield,environment,replicate
Sample001,45.2,field_A,1
Sample001,47.8,field_A,2
Sample001,44.1,field_B,1
Sample002,52.1,field_A,1
Sample002,50.3,field_A,2
```

## Running MET GWAS

```bash
binx gwas \
  --geno genotypes.tsv \
  --pheno phenotypes.csv \
  --trait yield \
  --covariates environment \
  --ploidy 4 \
  --out results.csv
```

## Analyzing G×E Interactions

*Coming soon: Detailed tutorial on G×E analysis*

## See Also

- [Your First GWAS](first-gwas.md)
- [binx gwas](../commands/gwas.md)
