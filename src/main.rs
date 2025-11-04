use clap::Parser;
use log::{debug, error, info, warn};
use rayon::prelude::*;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::sync::Arc;

/// Input file format
#[derive(Debug, PartialEq)]
enum InputFormat {
    Pack, // Coverage file (PACK format)
    Gaf,  // GAF alignment file
}

/// Detect input format based on file extension and content
fn detect_format(path: &str) -> std::io::Result<InputFormat> {
    let path_lower = path.to_lowercase();

    // Check extension first
    if path_lower.ends_with(".gaf") || path_lower.ends_with(".gaf.gz") {
        return Ok(InputFormat::Gaf);
    }

    // Check content for ambiguous cases
    let file = File::open(path)?;
    let (reader, _) = niffler::get_reader(Box::new(file))
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
    let mut buf_reader = BufReader::new(reader);
    let mut line = String::new();

    buf_reader.read_line(&mut line)?;

    if line.starts_with("##sample:") {
        Ok(InputFormat::Pack)
    } else if !line.is_empty() {
        // Check if it's GAF by looking for graph path with > or < characters
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 6 && (fields[5].contains('>') || fields[5].contains('<')) {
            Ok(InputFormat::Gaf)
        } else {
            Ok(InputFormat::Pack) // Default fallback
        }
    } else {
        Ok(InputFormat::Pack)
    }
}

/// Process input file using pre-parsed GFA segments (efficient for parallel processing)
fn process_input_to_coverage(
    input_path: &str,
    gfa_data: &Option<Arc<(Vec<usize>, usize)>>,
) -> std::io::Result<Vec<f64>> {
    let format = detect_format(input_path)?;
    debug!("Detected format for {}: {:?}", input_path, format);

    match format {
        InputFormat::Pack => {
            // Read PACK format directly
            parse_coverage_file(input_path).map(|(_, cov)| cov)
        }
        InputFormat::Gaf => {
            // Convert GAF to coverage using pre-parsed segments
            let gfa = gfa_data.as_ref().expect("GFA data should be available for GAF input");
            gafpack::compute_coverage_with_segments(&gfa.0, gfa.1, input_path, true, false)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))
        }
    }
}

/// Parse a coverage file (from gafpack --coverage-column)
fn parse_coverage_file(path: &str) -> std::io::Result<(String, Vec<f64>)> {
    let file_path = std::path::Path::new(path);
    let file = File::open(file_path)?;
    let (reader, _compression) = niffler::get_reader(Box::new(file))
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
    let mut reader = BufReader::new(reader);
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
        } else if !line_str.is_empty()
            && let Ok(value) = line_str.parse::<f64>()
        {
            coverage.push(value);
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

/// Pangenome-based zygosity matrices from graph coverage data.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None, disable_help_flag = true, disable_version_flag = true)]
struct Args {
    /// Sample table file: sample_name<tab>input_file (supports PACK and GAF formats)
    #[arg(help_heading = "Input", short, long)]
    sample_table: String,

    /// GFA pangenome graph file (required for GAF inputs)
    #[arg(help_heading = "Input", short, long)]
    gfa: Option<String>,

    /// Ploidy level (1 or 2)
    #[arg(help_heading = "Calling parameters", short, long, default_value = "2")]
    ploidy: u8,

    /// Normalization method (mean or median)
    #[arg(
        help_heading = "Calling parameters",
        short = 'm',
        long,
        default_value = "median"
    )]
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

    /// Output feature coverage mask file (1=keep, 0=filter outliers)
    #[arg(help_heading = "Output", long)]
    feature_cov_mask: Option<String>,

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
    if args.genotype_matrix.is_none()
        && args.dosage_matrix.is_none()
        && args.dosage_bimbam.is_none()
    {
        error!(
            "At least one output must be specified: --genotype-matrix, --dosage-matrix, or --dosage-bimbam"
        );
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

    // Log parameters
    if args.verbose > 1 {
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
    }

    // Parse GFA only if needed (i.e., if there are GAF files in the sample table)
    let gfa_data = if let Some(ref gfa_path) = args.gfa {
        if has_gaf_files(&args.sample_table)? {
            info!("Parsing GFA graph from {}", gfa_path);
            let gfa = gafpack::parse_gfa(gfa_path)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
            debug!(
                "Parsed GFA: {} segments, min_id={}",
                gfa.0.len(),
                gfa.1
            );
            Some(Arc::new(gfa))
        } else {
            warn!("No GAF files detected in sample table, skipping GFA parsing");
            None
        }
    } else {
        None
    };

    // Load samples
    info!("Loading samples from {}", args.sample_table);
    let samples = load_samples(&args.sample_table, &gfa_data)?;
    debug!("Loaded {} samples", samples.len());

    // Generate and write feature coverage mask if requested
    if let Some(ref mask_path) = args.feature_cov_mask {
        info!("Writing feature coverage mask to {}", mask_path);
        let mask = compute_feature_filter_mask(&samples);
        write_feature_filter_mask(&mask, mask_path)?;
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

/// Check if sample table contains any GAF files using format detection
fn has_gaf_files(sample_table: &str) -> std::io::Result<bool> {
    let file = File::open(sample_table)?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 2 {
            let input_path = fields[1];
            // Use detect_format to properly detect GAF files
            if let Ok(InputFormat::Gaf) = detect_format(input_path) {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// Load all samples from input file list
/// Each line should be: sample_name<tab>input_file_path
/// Supports PACK and GAF formats (automatically detected)
fn load_samples(
    input_file: &str,
    gfa_data: &Option<Arc<(Vec<usize>, usize)>>,
) -> std::io::Result<Vec<Sample>> {
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
            let input_file = fields[1];

            let coverage = match process_input_to_coverage(input_file, gfa_data) {
                Ok(cov) => cov,
                Err(e) => {
                    error!("Error processing {}: {}", input_file, e);
                    return None;
                }
            };

            let mean_coverage = compute_mean(&coverage);
            let median_coverage = compute_median(&coverage);
            let non_zero_count = coverage.iter().filter(|&&x| x > 0.0).count();

            debug!(
                "Loaded sample {} ({} features, mean={:.2} and median={:.2} coverage on {} non-zero features)",
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

/// Compute feature-level coverage filter mask using IQR method
/// Returns a vector where 1 = keep feature, 0 = filter feature
fn compute_feature_filter_mask(samples: &[Sample]) -> Vec<u8> {
    if samples.is_empty() {
        return Vec::new();
    }

    let num_features = samples[0].coverage.len();

    // Compute mean coverage for each feature across all samples (parallel)
    let feature_mean_coverages: Vec<f64> = (0..num_features)
        .into_par_iter()
        .map(|feature_idx| {
            let sum: f64 = samples.iter().map(|s| s.coverage[feature_idx]).sum();
            sum / samples.len() as f64
        })
        .collect();

    // Calculate quartiles and IQR
    let mut sorted_coverages = feature_mean_coverages.clone();
    sorted_coverages.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let q1_idx = (sorted_coverages.len() as f64 * 0.25) as usize;
    let q3_idx = (sorted_coverages.len() as f64 * 0.75) as usize;
    let q1 = sorted_coverages[q1_idx];
    let q3 = sorted_coverages[q3_idx];
    let iqr = q3 - q1;

    let lower_bound = q1 - 1.5 * iqr;
    let upper_bound = q3 + 1.5 * iqr;

    // Create mask
    let mask: Vec<u8> = feature_mean_coverages
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
        "Feature filter IQR: Q1={:.2}, Q3={:.2}, IQR={:.2}, bounds=[{:.2}, {:.2}]",
        q1, q3, iqr, lower_bound, upper_bound
    );
    debug!(
        "Feature filtering: {} features kept, {} filtered",
        num_features - filtered_count,
        filtered_count
    );

    mask
}

/// Call zygosity for a single feature based on coverage and ploidy
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

/// Write feature filter mask to output file
fn write_feature_filter_mask(mask: &[u8], output_path: &str) -> std::io::Result<()> {
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

    let num_features = samples[0].coverage.len();

    // Check all samples have the same number of features
    for sample in samples {
        if sample.coverage.len() != num_features {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "Sample {} has {} features, expected {}",
                    sample.name,
                    sample.coverage.len(),
                    num_features
                ),
            ));
        }
    }

    let mut file = File::create(output_path)?;

    // Write header
    write!(file, "#feature")?;
    for sample in samples {
        write!(file, "\t{}", sample.name)?;
    }
    writeln!(file)?;

    // Write zygosity for each feature (1-based indexing)
    for feature_idx in 0..num_features {
        write!(file, "{}", feature_idx + 1)?;

        for sample in samples {
            let ref_coverage = if method == "mean" {
                sample.mean_coverage
            } else {
                sample.median_coverage
            };

            let coverage = sample.coverage[feature_idx];
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

    let num_features = samples[0].coverage.len();
    let mut file = File::create(output_path)?;

    // Write header
    write!(file, "#feature")?;
    for sample in samples {
        write!(file, "\t{}", sample.name)?;
    }
    writeln!(file)?;

    // Write dosage for each feature
    for feature_idx in 0..num_features {
        write!(file, "{}", feature_idx + 1)?;

        for sample in samples {
            let ref_coverage = if method == "mean" {
                sample.mean_coverage
            } else {
                sample.median_coverage
            };

            let coverage = sample.coverage[feature_idx];
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

    let num_features = samples[0].coverage.len();
    let mut file = File::create(output_path)?;

    // Write BIMBAM format (no header)
    // Each line: feature_id, ref_allele, alt_allele, dosage1, dosage2, ...
    for feature_idx in 0..num_features {
        write!(file, "N{},A,T", feature_idx + 1)?;

        for sample in samples {
            let ref_coverage = if method == "mean" {
                sample.mean_coverage
            } else {
                sample.median_coverage
            };

            let coverage = sample.coverage[feature_idx];
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
    if non_zero.len().is_multiple_of(2) {
        (non_zero[mid - 1] + non_zero[mid]) / 2.0
    } else {
        non_zero[mid]
    }
}