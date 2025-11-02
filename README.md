# Pangosity

A tool to convert node coverage vectors from pangenome alignments into zygosity matrices.

## Overview

Pangosity takes node coverage vectors (e.g., from [gafpack](https://github.com/pangenome/gafpack)) and calls zygosity for each node based on coverage depth. It supports both haploid and diploid genomes.

## Installation

```bash
git clone --recursive https://github.com/AndreaGuarracino/pangosity.git
cd pangosity
cargo build --release
```

The binary will be available at `target/release/pangosity`.

## Usage

```bash
pangosity --input samples.txt --output output_prefix --ploidy 2 --method median
```

### Input Format

The input file should be a tab-separated file with two columns:
```
sample1	sample1.coverage.txt
sample2	sample2.coverage.txt
sample3	sample3.coverage.txt
```

Coverage files should be in the format produced by `gafpack`:

**Column format** (from `gafpack --coverage-column`):
```
##sample: sample_name
#coverage
123.45
234.56
...
```

**Row format** (from `gafpack` default output):
```
#sample	node.1	node.2	node.3	...
sample_name	123.45	234.56	345.67	...
```

### Parameters

- `--input, -i`: Input file list (sample_name<tab>coverage_file, one per line)
- `--output, -o`: Output prefix for the zygosity matrix
- `--ploidy, -p`: Ploidy level (1 or 2) [default: 2]
- `--method, -m`: Normalization method (mean or median) [default: median]
- `--min-coverage`: Minimum coverage threshold for calling genotypes [default: 0.0]

### Output Format

The tool produces a tab-separated matrix file (`{output_prefix}.tsv`):

```
#node	sample1	sample2	sample3
0	0/0	0/1	1/1
1	0/1	1/1	0/0
2	1/1	0/0	0/1
...
```

For haploid genomes (ploidy=1), genotypes are:
- `0`: Absent
- `1`: Present
- `.`: Missing/low coverage

For diploid genomes (ploidy=2), genotypes are:
- `0/0`: Homozygous reference (absent)
- `0/1`: Heterozygous (one copy)
- `1/1`: Homozygous alternate (two copies)
- `./.`: Missing/low coverage

## Zygosity Calling Logic

The tool uses coverage-based thresholds relative to the mean or median coverage:

### Ploidy 1 (Haploid)
- Coverage < 0.5 × reference → `0` (absent)
- Coverage ≥ 0.5 × reference → `1` (present)

### Ploidy 2 (Diploid)
- Coverage < 0.25 × reference → `0/0` (homozygous absent)
- 0.25 × reference ≤ coverage < 0.75 × reference → `0/1` (heterozygous)
- Coverage ≥ 0.75 × reference → `1/1` (homozygous present)

The reference coverage is computed as either the mean or median of all node coverages for each sample.

## Example Workflow

1. Generate coverage vectors with gafpack:
```bash
gafpack --gfa graph.gfa --gaf sample1.gaf --coverage-column --len-scale > sample1.coverage.txt
gafpack --gfa graph.gfa --gaf sample2.gaf --coverage-column --len-scale > sample2.coverage.txt
```

2. Create input file list:
```bash
echo -e "sample1\tsample1.coverage.txt\nsample2\tsample2.coverage.txt" > samples.txt
```

3. Run pangosity:
```bash
pangosity -i samples.txt -o zygosity -p 2 -m median
```

4. Output will be written to `zygosity.tsv`

## Related Tools

- [gafpack](https://github.com/ekg/gafpack): Generate node coverage vectors from GAF alignments
- [gfa2bin](https://github.com/MoinSebi/gfa2bin): Convert pangenome data to binary matrices

## License

MIT
