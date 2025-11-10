# Pangosity

Pangenome-based zygosity matrices.

<p align="center">
  <img src="pangosity-logo.png" alt="Pangosity logo" width="300"/>
</p>

## Installation

```bash
git clone --recursive https://github.com/AndreaGuarracino/pangosity.git
cd pangosity
cargo build --release
```

## Usage

```bash
# With input PACK format (coverage files) - no GFA needed
pangosity -s samples.txt -g genotypes.tsv

# With input GAF format - GFA required
pangosity -s samples.txt --gfa graph.gfa -g genotypes.tsv

# Save coverage (for GAF inputs)
pangosity -s samples.txt --gfa graph.gfa -g genotypes.tsv --save-coverage coverage_dir/

# Multiple output formats
pangosity -s samples.txt --gfa graph.gfa -g genotypes.tsv -d dosages.tsv -b dosages.bimbam

# Mixed input formats (automatically detected)
pangosity -s mixed_samples.txt --gfa graph.gfa -d dosages.tsv
```

### Input

Pangosity supports multiple input formats (uncompressed or gzipped) with automatic detection:

- **PACK** (coverage files) - Direct input, no conversion needed
- **GAF** (Graph Alignment Format) - Converted to coverage via `gafpack`; it requires the GFA graph file

Sample table file (tab-separated):
```
sample1	sample1.pack.gz
sample2	sample2.gaf.gz
```

#### PACK format

Coverage files in tall format:
```
##sample: sample_name
#coverage
123.45
234.56
...
```

### Parameters

**Input:**
- `-s, --sample-table`: Sample table file: sample_name<tab>coverage_file (supports .gz)

**Calling parameters:**
- `-p, --ploidy`: Ploidy level (1 or 2) [default: 2]
- `-m, --norm-method`: Normalization method (mean or median) [default: median]
- `-c, --calling-thresholds`: Genotype calling thresholds [default: 0.5 if ploidy=1; 0.25,0.75 if ploidy=2]

**Output:**
- `-g, --genotype-matrix`: Output genotype matrix file (0,1 if ploidy=1; 0/0,0/1,1/1 if ploidy=2)
- `-d, --dosage-matrix`: Output dosage matrix file (0,1 if ploidy=1; 0,1,2 if ploidy=2)
- `-b, --dosage-bimbam`: Output dosage matrix file in BIMBAM format (feature,ref,alt,dosages...)
- `--feature-cov-mask`: Output feature coverage mask file (1=keep, 0=filter outliers)
- `--save-coverage`: Save computed coverage to directory as compressed PACK files (for GAF inputs)

**General:**
- `-t, --threads`: Number of threads for parallel processing [default: 4]
- `-v, --verbose`: Verbosity level (0=error, 1=info, 2=debug) [default: 1]
- `-h, --help`: Print help
- `-V, --version`: Print version

### Output

**Genotype matrix** (tab-separated):

```
#feature   sample1     sample2     sample3
1          0/0         0/1         1/1
2          0/1         1/1         0/0
```

- **Haploid** (ploidy=1): `0` (absent), `1` (present), `.` (missing)
- **Diploid** (ploidy=2): `0/0` (absent), `0/1` (heterozygous), `1/1` (present), `./.` (missing)

**Dosage matrix** (tab-separated):

```
#feature   sample1     sample2     sample3
1          0           1           2
2          1           2           0
```

- **Haploid** (ploidy=1): `0` (absent), `1` (present), `NA` (missing)
- **Diploid** (ploidy=2): `0` (0/0), `1` (0/1), `2` (1/1), `NA` (missing)

**BIMBAM dosage matrix** (comma-separated):

```
N1,A,T,0,1,2
N2,A,T,1,2,0
N3,A,T,0,0,1
```

Format: `feature_id,ref_allele,alt_allele,dosage1,dosage2,dosage3,...`

- Features as rows, samples as columns
- Alleles: A (reference), T (alternate)
- Missing dosages written as `NA`

## Genotype calling

Genotypes are called using coverage thresholds relative to each sample's reference coverage (`ref`). Thresholds must be positive and strictly ascending. `ref` is defined as each sample's mean or median coverage across all features (computed using non-zero values only).

### Simple mode (default)

**Ploidy 1**: specify a `threshold` (default `threshold=0.5`):
- `coverage < threshold×ref → 0`
- `coverage ≥ threshold×ref → 1`

**Ploidy 2**: specify `low` and `high` thresholds (default `low=0.25, high=0.75`):
- `coverage < low×ref → 0/0`
- `low×ref ≤ coverage < high×ref → 0/1`
- `coverage ≥ high×ref → 1/1`

### No-call mode

To mark uncertain genotypes as missing:

**Ploidy 1**: specify 2 values to create a no-call zone:
```bash
pangosity -s samples.txt -g genotypes.tsv -p 1 -c 0.4,0.6
```
- `coverage < 0.4×ref → 0`
- `0.4×ref ≤ coverage < 0.6×ref → .` (no-call)
- `coverage ≥ 0.6×ref → 1`

**Ploidy 2**: specify 4 values to create no-call zones:
```bash
pangosity -s samples.txt -g genotypes.tsv -p 2 -c 0.15,0.35,0.65,0.85
```
- `coverage < 0.15×ref → 0/0`
- `0.15×ref ≤ coverage < 0.35×ref → ./.` (no-call)
- `0.35×ref ≤ coverage < 0.65×ref → 0/1`
- `0.65×ref ≤ coverage < 0.85×ref → ./.` (no-call)
- `coverage ≥ 0.85×ref → 1/1`

## Feature coverage filtering

Outlier detection using IQR method:

```bash
pangosity -s samples.txt -g genotypes.tsv --feature-cov-mask coverage_mask.txt
```

The filter:
- Computes mean coverage for each feature across all samples
- Calculates Q1 (25th percentile) and Q3 (75th percentile)
- Computes the Interquartile range (IQR) = Q3 - Q1
- Sets bounds at Q1 - 1.5×IQR and Q3 + 1.5×IQR
- Outputs mask file: `1` (keep feature) or `0` (filter outlier), one value per line

## Example

```bash
# Generate coverage with gafpack
gafpack --gfa graph.gfa --gaf sample1.gaf --coverage-column --len-scale | gzip > sample1.coverage.txt.gz
gafpack --gfa graph.gfa --gaf sample2.gaf --coverage-column --len-scale | gzip > sample2.coverage.txt.gz

# Create sample table
echo -e "sample1\tsample1.coverage.txt.gz\nsample2\tsample2.coverage.txt.gz" > samples.txt

# Call genotypes
pangosity -s samples.txt -g genotypes.tsv
```

## GWAS with GEMMA

For genome-wide association studies, pangosity can output dosage matrices in BIMBAM format:

```bash
# Generate BIMBAM dosage matrix for GEMMA
pangosity -s samples.txt -b dosages.bimbam

# Create phenotype file (one value per line, matching sample order in samples.txt)
echo -e "1.5\n2.3\n0.8" > phenotypes.txt

# Run GEMMA
# 1. Calculate kinship matrix
gemma -g dosages.bimbam -p phenotypes.txt -gk 1 -o kinship

# 2. Run linear mixed model association
gemma -g dosages.bimbam -p phenotypes.txt -k output/kinship.cXX.txt -lmm 1 -o results
```

**BIMBAM format** is recognized automatically by GEMMA. The phenotype file must have one value per line, in the same order as samples appear in your sample table.

## License

MIT
