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

### Output

Tab-separated genotype matrix:
```
#node   sample1     sample2     sample3
1       0/0         0/1         1/1
2       0/1         1/1         0/0
```

**Haploid** (ploidy=1): `0` (absent), `1` (present), `.` (missing)

**Diploid** (ploidy=2): `0/0` (absent), `0/1` (heterozygous), `1/1` (present), `./.` (missing)

## Genotype Calling

Genotypes are called using coverage thresholds relative to each sample's mean or median coverage (non-zero values only).

- **Haploid (ploidy 1)**: `coverage < t×ref → 0`, `coverage ≥ t×ref → 1` (default `t=0.5`)

- **Diploid (ploidy 2)**: `coverage < L×ref → 0/0`, `L×ref ≤ coverage < U×ref → 0/1`, `coverage ≥ U×ref → 1/1` (default `L=0.25, U=0.75`)

Customize with `-t`: `--calling-thresholds 0.6` (ploidy 1) or `--calling-thresholds 0.3,0.7` (ploidy 2).

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
