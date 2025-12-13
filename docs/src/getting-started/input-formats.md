# Input Formats

This page describes the input file formats that Binx accepts.

## Genotype File

The genotype file contains marker information and dosage values for each sample.

### Format Specification

- **File type**: TSV (tab-separated) or CSV (comma-separated)
- **Header**: Required (first row)
- **Columns**:
  - Column 1: `marker_id` - unique marker identifier
  - Column 2: `chrom` - chromosome name/number
  - Column 3: `pos` - base pair position (integer)
  - Columns 4+: Sample dosage values

### Example

```
marker_id	chrom	pos	Sample1	Sample2	Sample3	Sample4
SNP_1_1000	1	1000	0	2	4	1
SNP_1_2000	1	2000	1	1	3	2
SNP_1_3500	1	3500	4	0	2	2
SNP_2_500	2	500	2	2	2	3
SNP_2_1200	2	1200	0	1	1	0
```

### Dosage Values

Dosage values represent the count of the alternate allele:

| Ploidy | Valid Values | Meaning |
|--------|-------------|---------|
| Diploid (2) | 0, 1, 2 | 0=AA, 1=AB, 2=BB |
| Tetraploid (4) | 0, 1, 2, 3, 4 | 0=AAAA ... 4=BBBB |
| Hexaploid (6) | 0, 1, 2, 3, 4, 5, 6 | 0=AAAAAA ... 6=BBBBBB |

### Missing Values

Missing genotypes can be encoded as:
- `NA`
- Empty cell
- `.`

```
marker_id	chrom	pos	Sample1	Sample2	Sample3
SNP001	1	1000	0	NA	4
SNP002	1	2000	1	.	3
SNP003	1	3500		2	2
```

### Converting from VCF

Use `binx convert` to create a genotype file from VCF:

```bash
binx convert \
  --vcf input.vcf.gz \
  --format gwaspoly \
  --output genotypes.tsv
```

---

## Phenotype File

The phenotype file contains trait values and optional covariates for each sample.

### Format Specification

- **File type**: TSV or CSV
- **Header**: Required (first row)
- **Columns**:
  - Column 1: `sample_id` - must match genotype file column headers
  - Remaining columns: traits and/or covariates

### Example

```csv
sample_id,yield,height,flowering_date,environment,block
Sample1,45.2,120,156,field_A,1
Sample2,52.1,115,148,field_A,2
Sample3,48.7,125,152,field_B,1
Sample4,51.3,118,150,field_B,2
```

### Trait Values

- Numeric values for quantitative traits
- Can include `NA` or empty cells for missing values

### Covariates

Binx automatically detects covariate types:

| Type | Detection | Example |
|------|-----------|---------|
| Numeric | All values are numbers | `height`, `age` |
| Factor | Contains non-numeric values | `environment`, `block` |

Factor covariates are automatically dummy-coded during analysis.

### Multi-Environment Trials

For repeated measurements (same sample in multiple environments), repeat the sample ID:

```csv
sample_id,yield,environment,replicate
Sample1,45.2,field_A,1
Sample1,47.8,field_A,2
Sample1,44.1,field_B,1
Sample2,52.1,field_A,1
Sample2,50.3,field_A,2
```

---

## Kinship Matrix

The kinship matrix represents genetic relationships between samples.

### Format Specification

- **File type**: TSV or CSV
- **Structure**: Square symmetric matrix
- **Header**: Sample IDs
- **Row names**: Sample IDs (first column)

### Example

```
sample_id	Sample1	Sample2	Sample3	Sample4
Sample1	1.000	0.250	0.125	0.150
Sample2	0.250	1.000	0.200	0.180
Sample3	0.125	0.200	1.000	0.220
Sample4	0.150	0.180	0.220	1.000
```

### Computing a Kinship Matrix

Use `binx kinship` to compute from genotypes:

```bash
binx kinship \
  --geno genotypes.tsv \
  --ploidy 4 \
  --method vanraden \
  --out kinship.tsv
```

### Kinship Methods

| Method | Description |
|--------|-------------|
| `vanraden` | VanRaden (2008) method 1 (default) |
| `gwaspoly` | GWASpoly-style kinship |

---

## VCF Files

Binx can import VCF (Variant Call Format) files via `binx convert`.

### Supported Features

- Gzipped (`.vcf.gz`) and uncompressed (`.vcf`) files
- Diploid and polyploid genotypes
- GT (genotype) and AD (allelic depth) fields

### Example Conversion

```bash
# Convert VCF to GWASpoly format (using GT field)
binx convert \
  --vcf input.vcf.gz \
  --format gwaspoly \
  --output genotypes.tsv

# Convert VCF to allele depths for dosage estimation (using AD field)
binx convert \
  --vcf input.vcf.gz \
  --format csv \
  --output allele_depths.csv
```

---

## File Tips

### Sample ID Matching

Sample IDs in the phenotype file must exactly match the column headers in the genotype file:

```
# Genotype file header:
marker_id  chrom  pos  Sample_001  Sample_002  Sample_003

# Phenotype file:
sample_id,yield
Sample_001,45.2    ✓ matches
Sample_002,52.1    ✓ matches
sample_003,48.7    ✗ case mismatch!
```

### Chromosome Naming

Chromosome names can be:
- Numeric: `1`, `2`, `3`, ...
- String: `chr1`, `Chr1`, `chromosome1`

Be consistent within your genotype file.

### Large Files

For very large datasets:

1. **Use compressed VCF**: Keep VCF files gzipped
2. **Filter early**: Apply MAF and missing data filters during conversion
3. **Subset chromosomes**: Analyze one chromosome at a time if needed

### Validation

Check your files before analysis:

```bash
# Check genotype file structure
head -5 genotypes.tsv

# Count samples and markers
awk 'NR==1 {print "Samples:", NF-3} END {print "Markers:", NR-1}' genotypes.tsv

# Check phenotype file
head phenotypes.csv
wc -l phenotypes.csv
```
