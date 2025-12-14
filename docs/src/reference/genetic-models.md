# Genetic Models

This page provides detailed information about the genetic models available in Binx for GWAS analysis.

## Overview

Genetic models define how allele dosage relates to phenotype. Binx implements models from GWASpoly (Rosyara et al., 2016) that accommodate various inheritance patterns in diploids and polyploids.

## Diploid Models (ploidy=2)

### Additive

The standard additive model assumes each allele copy contributes equally to the trait.

| Genotype | AA | AB | BB |
|----------|-----|-----|-----|
| Dosage | 0 | 1 | 2 |
| Model value | 0 | 1 | 2 |

**Use when:** Trait value scales linearly with allele count.

```bash
binx gwas --ploidy 2 --models additive ...
```

### 1-dom-ref (Reference Dominant)

Tests if the reference allele (A) is dominant.

| Genotype | AA | AB | BB |
|----------|-----|-----|-----|
| Dosage | 0 | 1 | 2 |
| Model value | 0 | 1 | 1 |

**Use when:** One copy of B is sufficient to express the B phenotype.

```bash
binx gwas --ploidy 2 --models 1-dom-ref ...
```

### 1-dom-alt (Alternate Dominant)

Tests if the alternate allele (B) is dominant.

| Genotype | AA | AB | BB |
|----------|-----|-----|-----|
| Dosage | 0 | 1 | 2 |
| Model value | 0 | 0 | 1 |

**Use when:** Two copies of B are needed to express the B phenotype.

```bash
binx gwas --ploidy 2 --models 1-dom-alt ...
```

## Tetraploid Models (ploidy=4)

### Additive

Linear dosage effect across all five genotype classes.

| Genotype | AAAA | AAAB | AABB | ABBB | BBBB |
|----------|------|------|------|------|------|
| Dosage | 0 | 1 | 2 | 3 | 4 |
| Model value | 0 | 1 | 2 | 3 | 4 |

**Use when:** Each B allele adds equally to trait value.

### General (4 degrees of freedom)

No assumption about inheritance pattern. Estimates separate effects for each genotype class.

| Genotype | AAAA | AAAB | AABB | ABBB | BBBB |
|----------|------|------|------|------|------|
| Dosage | 0 | 1 | 2 | 3 | 4 |
| Dummy 1 | 0 | 1 | 0 | 0 | 0 |
| Dummy 2 | 0 | 0 | 1 | 0 | 0 |
| Dummy 3 | 0 | 0 | 0 | 1 | 0 |
| Dummy 4 | 0 | 0 | 0 | 0 | 1 |

**Use when:** Exploring inheritance pattern; hypothesis generation.

**Note:** Performs a single joint test of all 4 dummy variables simultaneously, using more degrees of freedom. This reduces power compared to additive but can detect complex inheritance patterns.

### 1-dom (Simplex Dominant)

One copy of B is sufficient for effect. Using `1-dom` tests both `1-dom-ref` and `1-dom-alt`.

| Genotype | AAAA | AAAB | AABB | ABBB | BBBB |
|----------|------|------|------|------|------|
| Dosage | 0 | 1 | 2 | 3 | 4 |
| 1-dom-ref | 0 | 1 | 1 | 1 | 1 |
| 1-dom-alt | 0 | 0 | 0 | 0 | 1 |

**Use when:** Trait exhibits complete dominance; presence/absence effect.

```bash
binx gwas --ploidy 4 --models 1-dom ...
```

### 2-dom (Duplex Dominant)

Two copies of B are sufficient for effect. Using `2-dom` tests both `2-dom-ref` and `2-dom-alt`.

| Genotype | AAAA | AAAB | AABB | ABBB | BBBB |
|----------|------|------|------|------|------|
| Dosage | 0 | 1 | 2 | 3 | 4 |
| 2-dom-ref | 0 | 0 | 1 | 1 | 1 |
| 2-dom-alt | 0 | 0 | 0 | 1 | 1 |

**Use when:** Partial dominance; two copies needed for effect.

```bash
binx gwas --ploidy 4 --models 2-dom ...
```

### diplo-additive (Diploidized Additive)

Treats the tetraploid as if it were diploid.

| Genotype | AAAA | AAAB | AABB | ABBB | BBBB |
|----------|------|------|------|------|------|
| Dosage | 0 | 1 | 2 | 3 | 4 |
| Model value | 0 | 0.5 | 0.5 | 0.5 | 1 |

**Use when:** Expecting diploid-like inheritance in autopolyploid.

```bash
binx gwas --ploidy 4 --models diplo-additive ...
```

### diplo-general (Diploidized General)

Diploid-style general model in tetraploid context (heterozygotes collapsed).

| Genotype | AAAA | AAAB | AABB | ABBB | BBBB |
|----------|------|------|------|------|------|
| Dosage | 0 | 1 | 2 | 3 | 4 |
| Group | AA | Het | Het | Het | BB |

**Use when:** Expected diploid-like inheritance with unknown dominance.

```bash
binx gwas --ploidy 4 --models diplo-general ...
```

## Hexaploid Models (ploidy=6)

Similar patterns extend to hexaploids:

| Model | Encoding (dosage 0-6) |
|-------|----------------------|
| `additive` | 0, 1, 2, 3, 4, 5, 6 |
| `general` | 6 dummy variables |
| `1-dom` | 0, 1, 1, 1, 1, 1, 1 (ref) / 0, 0, 0, 0, 0, 0, 1 (alt) |
| `2-dom` | 0, 0, 1, 1, 1, 1, 1 (ref) / 0, 0, 0, 0, 0, 1, 1 (alt) |
| `3-dom` | 0, 0, 0, 1, 1, 1, 1 (ref) / 0, 0, 0, 0, 1, 1, 1 (alt) |

## Choosing Models

### Recommended Strategy

1. **Start broad**: Run additive + general models
2. **Compare results**: Look for QTLs significant in one but not other
3. **Refine hypotheses**: Test specific dominance models
4. **Validate**: Check if model assumptions match biology

### Model Selection Guide

| Scenario | Recommended Models |
|----------|-------------------|
| Unknown inheritance | `additive,general` |
| Quantitative trait | `additive` |
| Disease resistance | `additive,1-dom,2-dom` |
| Exploratory analysis | `additive,general` |
| Confirmation study | Model from prior evidence |

### Multiple Testing Considerations

Running multiple models increases false positive rate:

- Apply correction across all tests
- Or use Bonferroni within each model separately
- Consider the general model as a single 4-df test

## Statistical Details

### Effect Estimation

For each model, Binx estimates:

```
y = μ + Xβ + u + ε
```

Where:
- `y` = phenotype
- `μ` = intercept
- `X` = design matrix (model-specific)
- `β` = fixed marker effect
- `u` = random polygenic effect
- `ε` = residual

### Testing

The null hypothesis (H₀: β = 0) is tested using a Wald test:

```
W = β² / Var(β)
```

Which follows a χ² distribution with degrees of freedom depending on the model.

## Examples

### Compare Additive vs Dominant

```bash
# Run both models
binx gwas \
  --geno geno.tsv \
  --pheno pheno.csv \
  --trait yield \
  --ploidy 4 \
  --models additive,1-dom \
  --out results.csv

# Find markers significant in one but not other
awk -F',' 'NR==1 {print; next}
  {key=$1","$2","$3; if(key in seen) {
    if(($4=="additive" && $5>5 && seen[key]<5) ||
       ($4!="additive" && $5<5 && seen[key]>5))
      print key, "differs"
  }
  seen[key]=$5}' results.csv
```

### All Tetraploid Models

```bash
binx gwas \
  --geno geno.tsv \
  --pheno pheno.csv \
  --trait yield \
  --ploidy 4 \
  --models additive,general,1-dom,2-dom,diplo-additive,diplo-general \
  --out all_models.csv
```

## References

1. Rosyara, U.R., De Jong, W.S., Douches, D.S., & Endelman, J.B. (2016). Software for genome-wide association studies in autopolyploids and its application to potato. *The Plant Genome* 9(2).

2. Endelman, J.B. (2011). Ridge regression and other kernels for genomic selection with R package rrBLUP. *The Plant Genome* 4:250-255.

## See Also

- [binx gwas](../commands/gwas.md) - GWAS command reference
- [Working with Polyploids](../tutorials/polyploid-analysis.md) - Polyploid tutorial
