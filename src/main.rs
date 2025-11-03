use clap::Parser;
use log::{debug, error, info, warn};
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufReader, prelude::*};
use std::path::Path;

/// Pangenome-based zygosity matrices from graph coverage data.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None, disable_help_flag = true, disable_version_flag = true)]
struct Args {
    /// Sample table file: sample_name<tab>coverage_file (supports .gz)
    #[arg(help_heading = "Input", short, long)]
    sample_table: String,

    /// Ploidy level (1 or 2)
    #[arg(help_heading = "Calling parameters", short, long, default_value = "2")]
    ploidy: u8,

    /// Normalization method (mean or median)
    #[arg(help_heading = "Calling parameters", short = 'm', long, default_value = "median")]
    norm_method: String,

    /// Calling thresholds (value if ploidy=1; lower,upper if ploidy=2) [default: 0.5 if ploidy=1; 0.25,0.75 if ploidy=2]
    #[arg(short, help_heading = "Calling parameters", long)]
    calling_thresholds: Option<String>,

    /// Minimum coverage threshold (below: missing)
    #[arg(help_heading = "Calling parameters", long, default_value = "0.0")]
    min_coverage: f64,

    /// Output genotype matrix file (0,1 if ploidy=1; 0/0,0/1,1/1 if ploidy=2)
    #[arg(help_heading = "Output", short, long)]
    genotype_matrix: Option<String>,

    /// Output dosage matrix file (0,1 if ploidy=1; 0,1,2 if ploidy=2)
    #[arg(help_heading = "Output", short, long)]
    dosage_matrix: Option<String>,

    /// Output dosage matrix file in BIMBAM format (variant,ref,alt,dosages...)
    #[arg(help_heading = "Output", short = 'b', long)]
    dosage_bimbam: Option<String>,

    /// Output node coverage mask file (1=keep, 0=filter outliers)
    #[arg(help_heading = "Output", long)]
    node_cov_mask: Option<String>,

    /// Number of threads for parallel processing
    #[arg(help_heading = "General", short, long, default_value = "4")]
    threads: usize,

    /// Verbosity level (0=error, 1=info, 2=debug)
    #[arg(help_heading = "General", short, long, default_value = "1")]
    verbose: u8,

    /// Print help
    #[arg(help_heading = "General", short, long, action = clap::ArgAction::Help)]
    help: Option<bool>,

    /// Print version
    #[arg(help_heading = "General", short = 'V', long, action = clap::ArgAction::Version)]
    version: Option<bool>,
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
    HomRef,  // 0/0
    Het,     // 0/1
    HomAlt,  // 1/1
    Missing, // ./.
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

    fn to_dosage(&self, ploidy: u8) -> String {
        if ploidy == 1 {
            match self {
                Zygosity::HomRef => "0".to_string(),
                Zygosity::HomAlt => "1".to_string(),
                _ => "NA".to_string(),
            }
        } else {
            match self {
                Zygosity::HomRef => "0".to_string(),
                Zygosity::Het => "1".to_string(),
                Zygosity::HomAlt => "2".to_string(),
                Zygosity::Missing => "NA".to_string(),
            }
        }
    }
}

/// Read a file with optional gzip decompression
fn create_reader(path: &Path) -> std::io::Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    let (reader, _compression) = niffler::get_reader(Box::new(file)).unwrap();
    Ok(Box::new(BufReader::new(reader)))
}

/// Parse a coverage file (from gafpack --coverage-column)
fn parse_coverage_file(path: &str) -> std::io::Result<(String, Vec<f64>)> {
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
            format!(
                "Invalid coverage file: {}. Expected ##sample: header and coverage values",
                path
            ),
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
            Zygosity::HomAlt // Present (1)
        } else {
            Zygosity::HomRef // Absent (0)
        }
    } else if ploidy == 2 {
        let threshold_lower = het_lower * ref_coverage;
        let threshold_upper = het_upper * ref_coverage;

        if coverage < threshold_lower {
            Zygosity::HomRef // 0/0
        } else if coverage < threshold_upper {
            Zygosity::Het // 0/1
        } else {
            Zygosity::HomAlt // 1/1
        }
    } else {
        panic!("Unsupported ploidy: {}", ploidy);
    }
}

/// Load all samples from input file list
/// Each line should be: sample_name<tab>coverage_file_path
/// Coverage files must be in tall format (one value per line) and can be gzip compressed
fn load_samples(input_file: &str) -> std::io::Result<Vec<Sample>> {
    let file = File::open(input_file)?;
    let reader = BufReader::new(file);
    let lines: Vec<String> = reader.lines().collect::<Result<_, _>>()?;

    let samples: Vec<Sample> = lines
        .par_iter()
        .filter_map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() != 2 {
                warn!("Skipping malformed line: {}", line);
                return None;
            }

            let sample_name = fields[0].to_string();
            let coverage_file = fields[1];

            let coverage = match parse_coverage_file(coverage_file) {
                Ok((_pack_name, cov)) => cov,
                Err(e) => {
                    error!("Error parsing {}: {}", coverage_file, e);
                    return None;
                }
            };

            let mean_coverage = compute_mean(&coverage);
            let median_coverage = compute_median(&coverage);
            let non_zero_count = coverage.iter().filter(|&&x| x > 0.0).count();

            debug!(
                "Loaded sample {} ({} nodes, mean={:.2} and median={:.2} coverage on {} non-zero nodes)",
                sample_name,
                coverage.len(),
                mean_coverage,
                median_coverage,
                non_zero_count
            );

            Some(Sample {
                name: sample_name,
                coverage,
                mean_coverage,
                median_coverage,
            })
        })
        .collect();

    Ok(samples)
}

/// Compute node-level coverage filter mask using IQR method
/// Returns a vector where 1 = keep node, 0 = filter node
fn compute_node_filter_mask(samples: &[Sample]) -> Vec<u8> {
    if samples.is_empty() {
        return Vec::new();
    }

    let num_nodes = samples[0].coverage.len();

    // Compute mean coverage for each node across all samples (parallel)
    let node_mean_coverages: Vec<f64> = (0..num_nodes)
        .into_par_iter()
        .map(|node_idx| {
            let sum: f64 = samples.iter().map(|s| s.coverage[node_idx]).sum();
            sum / samples.len() as f64
        })
        .collect();

    // Calculate quartiles and IQR
    let mut sorted_coverages = node_mean_coverages.clone();
    sorted_coverages.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let q1_idx = (sorted_coverages.len() as f64 * 0.25) as usize;
    let q3_idx = (sorted_coverages.len() as f64 * 0.75) as usize;
    let q1 = sorted_coverages[q1_idx];
    let q3 = sorted_coverages[q3_idx];
    let iqr = q3 - q1;

    let lower_bound = q1 - 1.5 * iqr;
    let upper_bound = q3 + 1.5 * iqr;

    // Create mask
    let mask: Vec<u8> = node_mean_coverages
        .par_iter()
        .map(|&cov| {
            if cov >= lower_bound && cov <= upper_bound {
                1
            } else {
                0
            }
        })
        .collect();

    let filtered_count = mask.iter().filter(|&&m| m == 0).count();
    debug!(
        "Node filter IQR: Q1={:.2}, Q3={:.2}, IQR={:.2}, bounds=[{:.2}, {:.2}]",
        q1, q3, iqr, lower_bound, upper_bound
    );
    debug!(
        "Node filtering: {} nodes kept, {} filtered",
        num_nodes - filtered_count,
        filtered_count
    );

    mask
}

/// Write node filter mask to output file
fn write_node_filter_mask(mask: &[u8], output_path: &str) -> std::io::Result<()> {
    let mut file = File::create(output_path)?;
    for &m in mask {
        writeln!(file, "{}", m)?;
    }
    Ok(())
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

    Ok(())
}

/// Write dosage matrix to output file
fn write_dosage_matrix(
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
    let mut file = File::create(output_path)?;

    // Write header
    write!(file, "#node")?;
    for sample in samples {
        write!(file, "\t{}", sample.name)?;
    }
    writeln!(file)?;

    // Write dosage for each node
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

            write!(file, "\t{}", zygosity.to_dosage(ploidy))?;
        }
        writeln!(file)?;
    }

    Ok(())
}

/// Write dosage matrix in BIMBAM format for GEMMA
/// Format: variant_id, ref_allele, alt_allele, dosage1, dosage2, dosage3, ...
/// Missing dosages are written as "NA"
fn write_dosage_bimbam(
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
    let mut file = File::create(output_path)?;

    // Write BIMBAM format (no header)
    // Each line: variant_id, ref_allele, alt_allele, dosage1, dosage2, ...
    for node_idx in 0..num_nodes {
        write!(file, "N{},A,T", node_idx + 1)?;

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

            write!(file, ",{}", zygosity.to_dosage(ploidy))?;
        }
        writeln!(file)?;
    }

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // Set log level based on verbosity
    env_logger::Builder::new()
        .filter_level(match args.verbose {
            0 => log::LevelFilter::Error,
            1 => log::LevelFilter::Info,
            _ => log::LevelFilter::Debug,
        })
        .init();

    // Configure rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    // Validate at least one output is specified
    if args.genotype_matrix.is_none() && args.dosage_matrix.is_none() && args.dosage_bimbam.is_none() {
        error!("At least one output must be specified: --genotype-matrix, --dosage-matrix, or --dosage-bimbam");
        std::process::exit(1);
    }

    if args.ploidy != 1 && args.ploidy != 2 {
        error!("Ploidy must be 1 or 2");
        std::process::exit(1);
    }

    if args.norm_method != "mean" && args.norm_method != "median" {
        error!("Method must be 'mean' or 'median'");
        std::process::exit(1);
    }

    // Parse thresholds based on ploidy
    let (haploid_threshold, het_lower, het_upper) = if let Some(ref t) = args.calling_thresholds {
        let values: Vec<&str> = t.split(',').collect();
        if args.ploidy == 1 {
            if values.len() != 1 {
                error!(
                    "Ploidy 1 requires a single threshold value (e.g., --calling-thresholds 0.5)"
                );
                std::process::exit(1);
            }
            let threshold = values[0].parse::<f64>().unwrap_or_else(|_| {
                error!("Invalid threshold value: {}", values[0]);
                std::process::exit(1);
            });
            (threshold, 0.25, 0.75) // het values don't matter for ploidy 1
        } else {
            if values.len() != 2 {
                error!(
                    "Ploidy 2 requires two comma-separated threshold values (e.g., --calling-thresholds 0.25,0.75)"
                );
                std::process::exit(1);
            }
            let lower = values[0].parse::<f64>().unwrap_or_else(|_| {
                error!("Invalid lower threshold value: {}", values[0]);
                std::process::exit(1);
            });
            let upper = values[1].parse::<f64>().unwrap_or_else(|_| {
                error!("Invalid upper threshold value: {}", values[1]);
                std::process::exit(1);
            });
            (0.5, lower, upper) // haploid value doesn't matter for ploidy 2
        }
    } else {
        // Default thresholds
        (0.5, 0.25, 0.75)
    };

    debug!("Input: {}", args.sample_table);
    if let Some(ref path) = args.genotype_matrix {
        debug!("Genotype matrix output: {}", path);
    }
    if let Some(ref path) = args.dosage_matrix {
        debug!("Dosage matrix output: {}", path);
    }
    debug!("Ploidy: {}", args.ploidy);
    debug!("Normalization method: {}", args.norm_method);
    if args.ploidy == 1 {
        debug!("Calling threshold: {}", haploid_threshold);
    } else {
        debug!("Calling thresholds: {}, {}", het_lower, het_upper);
    }
    if args.min_coverage > 0.0 {
        debug!("Minimum coverage: {}", args.min_coverage);
    }

    // Load samples
    info!("Loading samples from {}", args.sample_table);
    let samples = load_samples(&args.sample_table)?;
    debug!("Loaded {} samples", samples.len());

    // Generate and write node coverage mask if requested
    if let Some(ref mask_path) = args.node_cov_mask {
        info!("Writing node coverage mask to {}", mask_path);
        let mask = compute_node_filter_mask(&samples);
        write_node_filter_mask(&mask, mask_path)?;
    }

    // Write genotype matrix if requested
    if let Some(ref genotype_path) = args.genotype_matrix {
        info!("Writing genotype matrix to {}", genotype_path);
        write_zygosity_matrix(
            &samples,
            genotype_path,
            args.ploidy,
            &args.norm_method,
            args.min_coverage,
            het_lower,
            het_upper,
            haploid_threshold,
        )?;
    }

    // Write dosage matrix if requested
    if let Some(ref dosage_path) = args.dosage_matrix {
        info!("Writing dosage matrix to {}", dosage_path);
        write_dosage_matrix(
            &samples,
            dosage_path,
            args.ploidy,
            &args.norm_method,
            args.min_coverage,
            het_lower,
            het_upper,
            haploid_threshold,
        )?;
    }

    // Write BIMBAM dosage matrix if requested
    if let Some(ref bimbam_path) = args.dosage_bimbam {
        info!("Writing dosage BIMBAM to {}", bimbam_path);
        write_dosage_bimbam(
            &samples,
            bimbam_path,
            args.ploidy,
            &args.norm_method,
            args.min_coverage,
            het_lower,
            het_upper,
            haploid_threshold,
        )?;
    }

    info!("Done!");
    Ok(())
}
