use clap::Parser;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::path::Path;

/// Pangosity: Call zygosity from node coverage vectors
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Sample table file: sample_name<tab>coverage_file (supports .gz)
    #[arg(short, long)]
    sample_table: String,

    /// Ploidy level: 1 or 2
    #[arg(short, long, default_value = "2")]
    ploidy: u8,

    /// Normalization method: mean or median
    #[arg(short = 'm', long, default_value = "median")]
    norm_method: String,

    /// Output genotype matrix file
    #[arg(short, long)]
    genotype_matrix: String,

    /// Minimum coverage threshold (below: missing)
    #[arg(long, default_value = "0.0")]
    min_coverage: f64,

    /// Calling thresholds (ploidy 1: value; ploidy 2: lower,upper)
    #[arg(short = 't', long)]
    calling_thresholds: Option<String>,
}

/// Sample coverage data
#[derive(Debug)]
struct Sample {
    name: String,
    coverage: Vec<f64>,
    mean_coverage: f64,
    median_coverage: f64,
}

/// Zygosity genotype for diploid
#[derive(Debug, Clone, Copy, PartialEq)]
enum Zygosity {
    HomRef,      // 0/0
    Het,         // 0/1
    HomAlt,      // 1/1
    Missing,     // ./.
}

impl Zygosity {
    fn to_string(&self) -> String {
        match self {
            Zygosity::HomRef => "0/0".to_string(),
            Zygosity::Het => "0/1".to_string(),
            Zygosity::HomAlt => "1/1".to_string(),
            Zygosity::Missing => "./.".to_string(),
        }
    }

    fn to_string_haploid(&self) -> String {
        match self {
            Zygosity::HomRef => "0".to_string(),
            Zygosity::HomAlt => "1".to_string(),
            _ => ".".to_string(),
        }
    }
}

/// Read a file with optional gzip decompression
fn create_reader(path: &Path) -> std::io::Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    let (reader, _compression) = niffler::get_reader(Box::new(file)).unwrap();
    Ok(Box::new(BufReader::new(reader)))
}

/// Parse a coverage file in tall/column format (from gafpack --coverage-column)
/// Supports compressed (.gz) and uncompressed files
/// Format:
/// ##sample: sample_name
/// #coverage
/// value1
/// value2
/// ...
fn parse_coverage_tall_format(path: &str) -> std::io::Result<(String, Vec<f64>)> {
    let file_path = Path::new(path);
    let mut reader = create_reader(file_path)?;
    let mut line = String::new();
    let mut sample_name = String::new();
    let mut coverage = Vec::new();

    let mut has_sample_line = false;

    loop {
        line.clear();
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 {
            break;
        }

        let line_str = line.trim();

        if line_str.starts_with("##sample:") {
            has_sample_line = true;
            sample_name = line_str
                .strip_prefix("##sample:")
                .unwrap()
                .trim()
                .to_string();
        } else if line_str.starts_with("#") {
            // Skip header lines
            continue;
        } else if !line_str.is_empty() {
            if let Ok(value) = line_str.parse::<f64>() {
                coverage.push(value);
            }
        }
    }

    // Validate that we got data
    if !has_sample_line || coverage.is_empty() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("Invalid tall format file: {}. Expected ##sample: header and coverage values", path),
        ));
    }

    Ok((sample_name, coverage))
}


/// Compute mean coverage (only non-zero values)
fn compute_mean(values: &[f64]) -> f64 {
    let non_zero: Vec<f64> = values.iter().filter(|&&x| x > 0.0).copied().collect();
    if non_zero.is_empty() {
        return 0.0;
    }
    non_zero.iter().sum::<f64>() / non_zero.len() as f64
}

/// Compute median coverage (only non-zero values)
fn compute_median(values: &[f64]) -> f64 {
    let mut non_zero: Vec<f64> = values.iter().filter(|&&x| x > 0.0).copied().collect();
    if non_zero.is_empty() {
        return 0.0;
    }
    non_zero.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mid = non_zero.len() / 2;
    if non_zero.len() % 2 == 0 {
        (non_zero[mid - 1] + non_zero[mid]) / 2.0
    } else {
        non_zero[mid]
    }
}

/// Call zygosity for a single node based on coverage and ploidy
///
/// For ploidy 1:
/// - coverage < haploid_threshold * ref_coverage = 0
/// - coverage >= haploid_threshold * ref_coverage = 1
///
/// For ploidy 2:
/// - coverage < het_lower * ref_coverage = 0/0
/// - het_lower * ref_coverage <= coverage < het_upper * ref_coverage = 0/1
/// - coverage >= het_upper * ref_coverage = 1/1
fn call_zygosity(
    coverage: f64,
    ref_coverage: f64,
    ploidy: u8,
    min_cov: f64,
    het_lower: f64,
    het_upper: f64,
    haploid_threshold: f64,
) -> Zygosity {
    if coverage < min_cov {
        return Zygosity::Missing;
    }

    if ploidy == 1 {
        if coverage >= haploid_threshold * ref_coverage {
            Zygosity::HomAlt  // Present (1)
        } else {
            Zygosity::HomRef  // Absent (0)
        }
    } else if ploidy == 2 {
        let threshold_lower = het_lower * ref_coverage;
        let threshold_upper = het_upper * ref_coverage;

        if coverage < threshold_lower {
            Zygosity::HomRef  // 0/0
        } else if coverage < threshold_upper {
            Zygosity::Het     // 0/1
        } else {
            Zygosity::HomAlt  // 1/1
        }
    } else {
        panic!("Unsupported ploidy: {}", ploidy);
    }
}

/// Load all samples from input file list
/// Each line should be: sample_name<tab>coverage_file_path
/// Coverage files must be in tall format (one value per line) and can be gzip compressed
fn load_samples(input_file: &str, _method: &str) -> std::io::Result<Vec<Sample>> {
    let file = File::open(input_file)?;
    let reader = BufReader::new(file);
    let mut samples = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 2 {
            eprintln!("Warning: Skipping malformed line: {}", line);
            continue;
        }

        let sample_name = fields[0].to_string();
        let coverage_file = fields[1];

        // Parse tall format coverage file (supports .gz compression)
        let coverage = match parse_coverage_tall_format(coverage_file) {
            Ok((_pack_name, cov)) => cov,
            Err(e) => {
                eprintln!("Error parsing {}: {}", coverage_file, e);
                continue;
            }
        };

        let mean_coverage = compute_mean(&coverage);
        let median_coverage = compute_median(&coverage);
        let non_zero_count = coverage.iter().filter(|&&x| x > 0.0).count();

        samples.push(Sample {
            name: sample_name,
            coverage,
            mean_coverage,
            median_coverage,
        });

        eprintln!(
            "Loaded sample: {} ({} total nodes, {} non-zero, mean={:.2}, median={:.2})",
            samples.last().unwrap().name,
            samples.last().unwrap().coverage.len(),
            non_zero_count,
            mean_coverage,
            median_coverage
        );
    }

    Ok(samples)
}

/// Write zygosity matrix to output file
fn write_zygosity_matrix(
    samples: &[Sample],
    output_path: &str,
    ploidy: u8,
    method: &str,
    min_coverage: f64,
    het_lower: f64,
    het_upper: f64,
    haploid_threshold: f64,
) -> std::io::Result<()> {
    if samples.is_empty() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "No samples to write",
        ));
    }

    let num_nodes = samples[0].coverage.len();

    // Check all samples have the same number of nodes
    for sample in samples {
        if sample.coverage.len() != num_nodes {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "Sample {} has {} nodes, expected {}",
                    sample.name,
                    sample.coverage.len(),
                    num_nodes
                ),
            ));
        }
    }

    let mut file = File::create(output_path)?;

    // Write header
    write!(file, "#node")?;
    for sample in samples {
        write!(file, "\t{}", sample.name)?;
    }
    writeln!(file)?;

    // Write zygosity for each node (1-based to match gafpack node.1, node.2, etc.)
    for node_idx in 0..num_nodes {
        write!(file, "{}", node_idx + 1)?;

        for sample in samples {
            let ref_coverage = if method == "mean" {
                sample.mean_coverage
            } else {
                sample.median_coverage
            };

            let coverage = sample.coverage[node_idx];
            let zygosity = call_zygosity(
                coverage,
                ref_coverage,
                ploidy,
                min_coverage,
                het_lower,
                het_upper,
                haploid_threshold,
            );

            let genotype = if ploidy == 1 {
                zygosity.to_string_haploid()
            } else {
                zygosity.to_string()
            };

            write!(file, "\t{}", genotype)?;
        }
        writeln!(file)?;
    }

    eprintln!("Wrote zygosity matrix to {}", output_path);
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    if args.ploidy != 1 && args.ploidy != 2 {
        eprintln!("Error: Ploidy must be 1 or 2");
        std::process::exit(1);
    }

    if args.norm_method != "mean" && args.norm_method != "median" {
        eprintln!("Error: Method must be 'mean' or 'median'");
        std::process::exit(1);
    }

    // Parse thresholds based on ploidy
    let (haploid_threshold, het_lower, het_upper) = if let Some(ref t) = args.calling_thresholds {
        let values: Vec<&str> = t.split(',').collect();
        if args.ploidy == 1 {
            if values.len() != 1 {
                eprintln!("Error: Ploidy 1 requires a single threshold value (e.g., --calling-thresholds 0.5)");
                std::process::exit(1);
            }
            let threshold = values[0].parse::<f64>().unwrap_or_else(|_| {
                eprintln!("Error: Invalid threshold value: {}", values[0]);
                std::process::exit(1);
            });
            (threshold, 0.25, 0.75) // het values don't matter for ploidy 1
        } else {
            if values.len() != 2 {
                eprintln!("Error: Ploidy 2 requires two comma-separated threshold values (e.g., --calling-thresholds 0.25,0.75)");
                std::process::exit(1);
            }
            let lower = values[0].parse::<f64>().unwrap_or_else(|_| {
                eprintln!("Error: Invalid lower threshold value: {}", values[0]);
                std::process::exit(1);
            });
            let upper = values[1].parse::<f64>().unwrap_or_else(|_| {
                eprintln!("Error: Invalid upper threshold value: {}", values[1]);
                std::process::exit(1);
            });
            (0.5, lower, upper) // haploid value doesn't matter for ploidy 2
        }
    } else {
        // Default thresholds
        (0.5, 0.25, 0.75)
    };

    eprintln!("Pangosity - Node coverage to zygosity matrix");
    eprintln!("Input: {}", args.sample_table);
    eprintln!("Ploidy: {}", args.ploidy);
    eprintln!("Method: {}", args.norm_method);
    eprintln!("Min coverage: {}", args.min_coverage);
    if args.ploidy == 1 {
        eprintln!("Haploid threshold: {}", haploid_threshold);
    } else {
        eprintln!("Diploid thresholds: {}, {}", het_lower, het_upper);
    }
    eprintln!();

    // Load samples
    let samples = load_samples(&args.sample_table, &args.norm_method)?;
    eprintln!("\nLoaded {} samples", samples.len());

    // Write zygosity matrix
    write_zygosity_matrix(
        &samples,
        &args.genotype_matrix,
        args.ploidy,
        &args.norm_method,
        args.min_coverage,
        het_lower,
        het_upper,
        haploid_threshold,
    )?;

    eprintln!("\nDone!");
    Ok(())
}
