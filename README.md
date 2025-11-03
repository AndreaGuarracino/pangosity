# Pangosity

Convert node coverage vectors from pangenome alignments into genotype matrices.

## Installation

```bash
git clone --recursive https://github.com/AndreaGuarracino/pangosity.git
cd pangosity
cargo build --release
```

## Usage

```bash
pangosity -s samples.txt -g genotypes.tsv
```

### Input

Sample table file (tab-separated):
```
sample1	sample1.coverage.txt.gz
sample2	sample2.coverage.txt.gz
```

Coverage files in tall format from `gafpack --coverage-column` (supports .gz):
```
##sample: sample_name
#coverage
123.45
234.56
...
```

### Parameters

- `-s, --sample-table`: Sample table (sample_name<tab>coverage_file)
- `-g, --genotype-matrix`: Output genotype matrix
- `-p, --ploidy`: Ploidy: 1 or 2 [default: 2]
- `-m, --norm-method`: Normalization: mean or median [default: median]
- `--min-coverage`: Minimum coverage [default: 0.0]
- `-t, --calling-thresholds`: Calling thresholds (ploidy 1: value; ploidy 2: lower,upper) [default: 0.5 / 0.25,0.75]
- `-v, --verbose`: Verbosity level: 0 (error), 1 (info), 2 (debug) [default: 1]
- `--node-filter-mask`: Output file for node coverage filter mask (1 = keep, 0 = filter outliers)

### Output

Tab-separated genotype matrix:
```
#node   sample1     sample2     sample3
1       0/0         0/1         1/1
2       0/1         1/1         0/0
```

**Haploid** (ploidy=1): `0` (absent), `1` (present), `.` (missing)

**Diploid** (ploidy=2): `0/0` (absent), `0/1` (heterozygous), `1/1` (present), `./.` (missing)

## Genotype calling

Genotypes are called using coverage thresholds relative to each sample's mean or median coverage (non-zero values only).

- **Haploid (ploidy 1)**: `coverage < t×ref → 0`, `coverage ≥ t×ref → 1` (default `t=0.5`)
- **Diploid (ploidy 2)**: `coverage < L×ref → 0/0`, `L×ref ≤ coverage < U×ref → 0/1`, `coverage ≥ U×ref → 1/1` (default `L=0.25, U=0.75`)

## Node coverage filtering

Outlier detection for node coverage filtering:

```bash
pangosity -s samples.txt -g genotypes.tsv --node-filter-mask filter.txt
```

The filter:
- Computes mean coverage for each node across all samples
- Calculates Q1 (25th percentile) and Q3 (75th percentile)
- Computes the Interquartile range (IQR) = Q3 - Q1
- Sets bounds at Q1 - 1.5×IQR and Q3 + 1.5×IQR
- Outputs mask file: `1` (keep node) or `0` (filter outlier), one value per line

## Example

```bash
# Generate coverage with gafpack
gafpack --gfa graph.gfa --gaf sample1.gaf --coverage-column --len-scale | gzip > sample1.coverage.txt.gz
gafpack --gfa graph.gfa --gaf sample2.gaf --coverage-column --len-scale | gzip > sample2.coverage.txt.gz

# Create sample table
echo -e "sample1\tsample1.coverage.txt.gz\nsample2\tsample2.coverage.txt.gz" > samples.txt

# Call genotypes
pangosity -s samples.txt -g genotypes.tsv -p 2 -m median
```

## License

MIT
