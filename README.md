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

# With copy-number matrix (continuous values for QC)
pangosity -s samples.txt -g genotypes.tsv -n copynumbers.tsv

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
- `--min-coverage`: Minimum coverage for including features in normalization calculations [default: 0.0]

**Output:**
- `-g, --genotype-matrix`: Output genotype matrix file (0,1 if ploidy=1; 0/0,0/1,1/1 if ploidy=2)
- `-d, --dosage-matrix`: Output dosage matrix file (0,1 if ploidy=1; 0,1,2 if ploidy=2)
- `-b, --dosage-bimbam`: Output dosage matrix file in BIMBAM format (feature,ref,alt,dosages...)
- `-n, --copynumber-matrix`: Output copy-number matrix file (continuous relative coverage values)
- `--feature-cov-mask`: Output feature coverage mask file (1=keep, 0=filter outliers)
- `--save-coverage`: Save computed coverage to directory as compressed PACK files (for GAF inputs)

**General:**
- `-t, --threads`: Number of threads for parallel processing [default: 4]
- `-v, --verbose`: Verbosity level (0=error, 1=info, 2=debug) [default: 1]
- `-h, --help`: Print help
- `-V, --version`: Print version

### Output

**Which output format to use:**
- **GWAS with GEMMA** → `-b, --dosage-bimbam` (CSV format required by GEMMA)
- **VCF-style variant tools** → `-g, --genotype-matrix` (0/0, 0/1, 1/1 notation)
- **General quantitative analysis** → `-d, --dosage-matrix` (numeric copy counts)
- **QC, visualization, custom thresholds** → `-n, --copynumber-matrix` (continuous values)
- **Filtering outlier features** → `--feature-cov-mask` (binary keep/filter mask)

---

**Genotype matrix** (`-g`): VCF-style genotypes
```
#feature   sample1   sample2   sample3
1          0/0       0/1       1/1
```
Ploidy 1: `0`, `1`, `.` | Ploidy 2: `0/0`, `0/1`, `1/1`, `./.`

**Dosage matrix** (`-d`): Numeric copy counts
```
#feature   sample1   sample2   sample3
1          0         1         2
```
Ploidy 1: `0`, `1`, `NA` | Ploidy 2: `0`, `1`, `2`, `NA`

**BIMBAM dosage** (`-b`): CSV format for GEMMA
```
N1,A,T,0,1,2
N2,A,T,1,2,0
```
Format: `feature_id,ref,alt,dosage1,dosage2,...`

**Copy-number matrix** (`-n`): Continuous relative coverage
```
#feature   sample1   sample2   sample3
1          0.021     0.987     1.943
```
Values = `coverage / haploid_coverage`. For QC, visualization, borderline call detection. Values > ploidy indicate duplications.

## Genotype calling

Genotypes are called using a **copy-number model** where coverage is normalized to haploid (single-copy) coverage. Thresholds must be positive and strictly ascending.

**Haploid coverage:** Each sample's mean or median coverage across features (with coverage above `--min-coverage`), divided by ploidy to estimate single-copy coverage.

**Copy-number interpretation:**
- **Ploidy 1**: 0 copies = 0×haploid, 1 copy = 1×haploid
- **Ploidy 2**: 0 copies = 0×haploid, 1 copy (0/1) = 1×haploid, 2 copies (1/1) = 2×haploid

### Simple mode (default)

**Ploidy 1**: specify a `threshold` (default `threshold=0.5`):
- `coverage < threshold×haploid → 0` (absent)
- `coverage ≥ threshold×haploid → 1` (present)

**Ploidy 2**: specify `low` and `high` thresholds (default `low=0.5, high=1.5`):
- `coverage < low×haploid → 0/0` (0 copies, ~0×haploid)
- `low×haploid ≤ coverage < high×haploid → 0/1` (1 copy, ~1×haploid)
- `coverage ≥ high×haploid → 1/1` (2 copies, ~2×haploid)

### No-call mode

To mark uncertain genotypes as missing:

**Ploidy 1**: specify 2 values to create a no-call zone:
```bash
pangosity -s samples.txt -g genotypes.tsv -p 1 -c 0.3,0.7
```
- `coverage < 0.3×haploid → 0`
- `0.3×haploid ≤ coverage < 0.7×haploid → .` (no-call)
- `coverage ≥ 0.7×haploid → 1`

**Ploidy 2**: specify 4 values to create no-call zones:
```bash
pangosity -s samples.txt -g genotypes.tsv -p 2 -c 0.3,0.7,1.3,1.7
```
- `coverage < 0.3×haploid → 0/0` (confident absent)
- `0.3×haploid ≤ coverage < 0.7×haploid → ./.` (no-call, uncertain 0/0 vs 0/1)
- `0.7×haploid ≤ coverage < 1.3×haploid → 0/1` (confident heterozygous, ~1 copy)
- `1.3×haploid ≤ coverage < 1.7×haploid → ./.` (no-call, uncertain 0/1 vs 1/1)
- `coverage ≥ 1.7×haploid → 1/1` (confident homozygous, ~2 copies)

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
