# binx dosage

Estimate genotype dosages from sequencing read count data.

## Synopsis

```bash
binx dosage --ploidy <INT> [INPUT OPTIONS] [OPTIONS]
```

## Description

The `dosage` command estimates genotype dosages from read count data using algorithms based on the R/Updog package (Gerard et al., 2018). This is useful when working with genotyping-by-sequencing (GBS) or similar data where discrete genotype calls may be uncertain.

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--ploidy <INT>` | Ploidy level (e.g., 2, 4, 6) |

## Input Options

Choose **one** of the following input modes:

### VCF Mode (Recommended)

```bash
binx dosage --vcf <FILE> --ploidy 4 --output dosages.tsv
```

| Option | Description |
|--------|-------------|
| `--vcf <FILE>` | VCF file (plain or gzipped) with FORMAT/AD allele depths |
| `--chunk-size <INT>` | Chunk size for streaming VCF markers (default: stream one by one) |

### Two-Line CSV Mode

```bash
binx dosage --csv <FILE> --ploidy 4 --output dosages.tsv
```

| Option | Description |
|--------|-------------|
| `--csv <FILE>` | CSV file with alternating lines of Ref and Total counts per locus |

### Matrix Mode

```bash
binx dosage --counts --ref-path ref.tsv --total-path total.tsv --ploidy 4 --output dosages.tsv
```

| Option | Description |
|--------|-------------|
| `--counts` | Enable matrix mode |
| `--ref-path <FILE>` | Ref count matrix (markers in rows, samples in columns; first column marker ID) |
| `--total-path <FILE>` | Total count matrix (markers in rows, samples in columns; first column marker ID) |

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--output <FILE>` | stdout | Output file path |
| `--mode <MODE>` | auto | Optimization mode (see below) |
| `--format <FMT>` | matrix | Output format (see below) |
| `--compress <MODE>` | none | Compression: `none` or `gzip` |
| `--threads <INT>` | num_cpus | Number of threads for parallel processing |
| `--verbose` | false | Enable verbose output |

### Optimization Modes

| Mode | Description |
|------|-------------|
| `auto` | Automatically select best mode based on data |
| `updog` | Standard Updog algorithm |
| `updog-fast` | Faster Updog with approximations |
| `updog-exact` | Exact Updog (slower, more accurate) |
| `fast` | Fast estimation |
| `turbo` | Fastest estimation |
| `turboauto` | Turbo with automatic parameter selection |
| `turboauto-safe` | Turboauto with additional safety checks |

### Output Formats

| Format | Description |
|--------|-------------|
| `matrix` | Simple dosage matrix (markers x samples) |
| `stats` | Detailed statistics per marker |
| `beagle` | BEAGLE format for imputation |
| `vcf` | VCF format with dosage annotations |
| `plink` | PLINK raw format |
| `gwaspoly` | GWASpoly-compatible format (marker, chrom, pos, samples...) |

## Examples

### Basic Dosage Estimation from VCF

```bash
binx dosage \
  --vcf variants.vcf.gz \
  --ploidy 4 \
  --output dosages.tsv
```

### Output in GWASpoly Format

```bash
binx dosage \
  --vcf variants.vcf.gz \
  --ploidy 4 \
  --format gwaspoly \
  --output genotypes.tsv
```

### Parallel Processing with Chunks

```bash
binx dosage \
  --vcf variants.vcf.gz \
  --ploidy 4 \
  --chunk-size 1000 \
  --threads 8 \
  --output dosages.tsv
```

### From Two-Line CSV

```bash
binx dosage \
  --csv read_counts.csv \
  --ploidy 4 \
  --mode updog \
  --output dosages.tsv
```

### From Separate Ref/Total Matrices

```bash
binx dosage \
  --counts \
  --ref-path ref_counts.tsv \
  --total-path total_counts.tsv \
  --ploidy 4 \
  --output dosages.tsv
```

### Compressed Output

```bash
binx dosage \
  --vcf variants.vcf.gz \
  --ploidy 4 \
  --compress gzip \
  --output dosages.tsv.gz
```

## Input Formats

### VCF Format

The VCF file should contain the `AD` (Allelic Depths) field in the FORMAT column:

```
#CHROM  POS     ID      REF  ALT  QUAL  FILTER  INFO  FORMAT      Sample1   Sample2
chr1    1000    SNP001  A    T    .     .       .     GT:AD:DP    0/1:10,5:15   1/1:2,18:20
```

### Two-Line CSV Format

Alternating lines of reference and total counts:

```
locus,Sample1,Sample2,Sample3
SNP001,10,2,15
SNP001,15,20,18
SNP002,8,12,5
SNP002,16,25,10
```

### Matrix Format

Two separate files with matching structure:

**ref_counts.tsv:**
```
marker_id	Sample1	Sample2	Sample3
SNP001	10	2	15
SNP002	8	12	5
```

**total_counts.tsv:**
```
marker_id	Sample1	Sample2	Sample3
SNP001	15	20	18
SNP002	16	25	10
```

## Output Format

The default `matrix` format:

```
marker_id	Sample1	Sample2	Sample3
SNP001	1	4	2
SNP002	2	2	1
```

The `gwaspoly` format (suitable for `binx gwas`):

```
Marker	Chrom	Position	Sample1	Sample2	Sample3
SNP001	chr1	1000	1	4	2
SNP002	chr1	2000	2	2	1
```

## See Also

- [binx convert](convert.md) - Convert VCF to other formats (GT-based)
- [binx gwas](gwas.md) - Run GWAS with dosage data
- [Input Formats](../getting-started/input-formats.md)
