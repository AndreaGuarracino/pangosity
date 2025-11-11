use clap::Parser;
use flate2::Compression;
use flate2::write::GzEncoder;
use log::{debug, error, info, warn};
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufReader, prelude::*};
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
        .map_err(|e| std::io::Error::other(e))?;
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
            let gfa = gfa_data.as_ref().ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    format!(
                        "Detected GAF input '{}' but no GFA graph was provided. Use --gfa to specify the graph file.",
                        input_path
                    ),
                )
            })?;
            gafpack::compute_coverage_with_segments(&gfa.0, gfa.1, input_path, true, false)
                .map_err(|e| std::io::Error::other(e))
        }
    }
}

/// Parse a coverage file (from gafpack --coverage-column)
fn parse_coverage_file(path: &str) -> std::io::Result<(String, Vec<f64>)> {
    let file_path = std::path::Path::new(path);
    let file = File::open(file_path)?;
    let (reader, _compression) = niffler::get_reader(Box::new(file))
        .map_err(|e| std::io::Error::other(e))?;
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
    #[arg(help_heading = "Input", long)]
    gfa: Option<String>,

    /// Normalization method (mean or median)
    #[arg(
        help_heading = "Calling parameters",
        short = 'm',
        long,
        default_value = "median"
    )]
    norm_method: String,

    /// Ploidy level (1 or 2)
    #[arg(help_heading = "Calling parameters", short, long, default_value = "2")]
    ploidy: u8,

    /// Genotype calling thresholds [default: 0.5 if ploidy=1; 0.5,1.5 if ploidy=2]
    #[arg(short, help_heading = "Calling parameters", long)]
    calling_thresholds: Option<String>,

    /// Minimum coverage threshold for including features in normalization calculations
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

    /// Output copy-number matrix file (continuous relative coverage values)
    #[arg(help_heading = "Output", short = 'n', long)]
    copynumber_matrix: Option<String>,

    /// Output feature coverage mask file (1=keep, 0=filter outliers)
    #[arg(help_heading = "Output", long)]
    feature_cov_mask: Option<String>,

    /// Save computed coverage to directory (for GAF inputs)
    #[arg(help_heading = "Output", long)]
    save_coverage: Option<String>,

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
    haploid_coverage: f64,
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
        && args.copynumber_matrix.is_none()
    {
        error!(
            "At least one output must be specified: --genotype-matrix, --dosage-matrix, --dosage-bimbam, or --copynumber-matrix"
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

    if args.min_coverage < 0.0 {
        error!("Minimum coverage threshold must be non-negative");
        std::process::exit(1);
    }

    // Parse thresholds based on ploidy
    let calling_thresholds: Vec<f64> = if let Some(ref t) = args.calling_thresholds {
        let values: Vec<&str> = t.split(',').collect();

        // Validate threshold count based on ploidy
        let valid_count = match (args.ploidy, values.len()) {
            (1, 1) => true,  // Simple mode: single threshold
            (1, 2) => true,  // No-call mode: lower,upper
            (2, 2) => true,  // Simple mode: het_lower,het_upper
            (2, 4) => true,  // No-call mode: t1,t2,t3,t4
            (1, n) => {
                error!(
                    "Ploidy 1 requires 1 threshold (simple) or 2 thresholds (with no-calls), got {}",
                    n
                );
                false
            }
            (2, n) => {
                error!(
                    "Ploidy 2 requires 2 thresholds (simple) or 4 thresholds (with no-calls), got {}",
                    n
                );
                false
            }
            _ => false,
        };

        if !valid_count {
            std::process::exit(1);
        }

        // Parse threshold values
        values
            .iter()
            .enumerate()
            .map(|(i, &v)| {
                v.parse::<f64>().unwrap_or_else(|_| {
                    error!("Invalid threshold value at position {}: {}", i + 1, v);
                    std::process::exit(1);
                })
            })
            .collect()
    } else {
        // Default thresholds (copy-number model)
        if args.ploidy == 1 {
            vec![0.5]
        } else {
            vec![0.5, 1.5]
        }
    };

    // Validate thresholds are positive and strictly ascending
    if calling_thresholds[0] <= 0.0 {
        error!("Thresholds must be positive, got {}", calling_thresholds[0]);
        std::process::exit(1);
    }
    for i in 1..calling_thresholds.len() {
        if calling_thresholds[i] <= calling_thresholds[i - 1] {
            error!(
                "Thresholds must be strictly ascending: {} ≤ {}",
                calling_thresholds[i - 1], calling_thresholds[i]
            );
            std::process::exit(1);
        }
    }

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
        debug!("Minimum coverage threshold: {}", args.min_coverage);
        match (args.ploidy, calling_thresholds.len()) {
            (1, 1) => debug!("Calling threshold: {}", calling_thresholds[0]),
            (1, 2) => debug!(
                "Calling thresholds with no-calls: <{} → 0, [{},{}] → missing, ≥{} → 1",
                calling_thresholds[0],
                calling_thresholds[0],
                calling_thresholds[1],
                calling_thresholds[1]
            ),
            (2, 2) => debug!(
                "Calling thresholds: <{} → 0/0, [{},{}) → 0/1, ≥{} → 1/1",
                calling_thresholds[0], calling_thresholds[0], calling_thresholds[1], calling_thresholds[1]
            ),
            (2, 4) => debug!(
                "Calling thresholds with no-calls: <{} → 0/0, [{},{}) → missing, [{},{}) → 0/1, [{},{}) → missing, ≥{} → 1/1",
                calling_thresholds[0],
                calling_thresholds[0], calling_thresholds[1],
                calling_thresholds[1], calling_thresholds[2],
                calling_thresholds[2], calling_thresholds[3],
                calling_thresholds[3]
            ),
            _ => {}
        }
    }

    // Validate sample table and check if GAF files are present
    let has_gaf = validate_sample_table(&args.sample_table)?;
    if has_gaf {
        info!("GAF input detected: coverage will be computed on-the-fly (slower than PACK)");
        if args.gfa.is_none() {
            error!(
                "GAF files detected in sample table but no GFA graph provided. Use --gfa to specify the graph file."
            );
            std::process::exit(1);
        }
    }

    // Create coverage output directory if requested
    if let Some(ref dir) = args.save_coverage {
        if let Err(e) = std::fs::create_dir_all(dir) {
            error!("Failed to create coverage directory '{}': {}", dir, e);
            std::process::exit(1);
        }
        if !has_gaf {
            warn!("--save-coverage specified but no GAF files detected.");
        }
    }

    // Parse GFA graph if provided (needed for GAF inputs)
    let gfa_data = if let Some(ref gfa_path) = args.gfa {
        if has_gaf {
            info!("Parsing GFA graph from {}", gfa_path);
        } else {
            warn!("--gfa provided but no GAF files were detected by extension in the sample table; parsing graph anyway");
        }
        let gfa = gafpack::parse_gfa(gfa_path)
            .map_err(|e| std::io::Error::other(e))?;
        debug!("Parsed GFA: {} segments, min_id={}", gfa.0.len(), gfa.1);
        Some(Arc::new(gfa))
    } else {
        None
    };

    // Load samples
    info!("Loading samples from {}", args.sample_table);
    let samples = load_samples(
        &args.sample_table,
        &gfa_data,
        &args.norm_method,
        args.ploidy,
        args.min_coverage,
        args.save_coverage.as_deref(),
    )?;
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
            &calling_thresholds,
        )?;
    }

    // Write dosage matrix if requested
    if let Some(ref dosage_path) = args.dosage_matrix {
        info!("Writing dosage matrix to {}", dosage_path);
        write_dosage_matrix(
            &samples,
            dosage_path,
            args.ploidy,
            &calling_thresholds,
        )?;
    }

    // Write BIMBAM dosage matrix if requested
    if let Some(ref bimbam_path) = args.dosage_bimbam {
        info!("Writing dosage BIMBAM to {}", bimbam_path);
        write_dosage_bimbam(
            &samples,
            bimbam_path,
            args.ploidy,
            &calling_thresholds,
        )?;
    }

    // Write copy-number matrix if requested
    if let Some(ref copynumber_path) = args.copynumber_matrix {
        info!("Writing copy-number matrix to {}", copynumber_path);
        write_copynumber_matrix(&samples, copynumber_path)?;
    }

    info!("Done!");
    Ok(())
}

/// Validate sample table format and check if it contains GAF files
fn validate_sample_table(sample_table: &str) -> std::io::Result<bool> {
    let file = File::open(sample_table)?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 2 {
            error!("Malformed line in sample table: {}", line);
            std::process::exit(1);
        }
        let input_path = fields[1];
        if let Ok(InputFormat::Gaf) = detect_format(input_path) {
            return Ok(true);
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
    norm_method: &str,
    ploidy: u8,
    min_cov_threshold: f64,
    save_coverage_dir: Option<&str>,
) -> std::io::Result<Vec<Sample>> {
    let file = File::open(input_file)?;
    let reader = BufReader::new(file);
    let lines: Vec<String> = reader.lines().collect::<Result<_, _>>()?;

    let samples: Vec<Sample> = lines
        .par_iter()
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let sample_name = fields[0].to_string();
            let input_file = fields[1];

            let coverage = process_input_to_coverage(input_file, gfa_data).unwrap_or_else(|e| {
                error!("Error processing {}: {}", input_file, e);
                std::process::exit(1);
            });

            // Compute typical coverage across features
            let typical_coverage = if norm_method == "mean" {
                compute_mean(&coverage, min_cov_threshold)
            } else {
                compute_median(&coverage, min_cov_threshold)
            };

            if typical_coverage <= 0.0 {
                error!(
                    "Sample '{}' has no coverage values above --min-coverage ({}). Cannot estimate haploid coverage.",
                    sample_name, min_cov_threshold
                );
                std::process::exit(1);
            }

            // Convert to haploid coverage estimate
            let haploid_coverage = typical_coverage / ploidy as f64;

            debug!(
                "Loaded sample {} ({} features, {}={:.2}, haploid={:.2} on {} features > {})",
                sample_name,
                coverage.len(),
                norm_method,
                typical_coverage,
                haploid_coverage,
                coverage.iter().filter(|&&x| x > min_cov_threshold).count(),
                min_cov_threshold
            );

            if let Some(dir) = save_coverage_dir {
                let output_path = format!("{}/{}.pack.gz", dir, sample_name);
                let content = gafpack::format_coverage_column(&sample_name, &coverage);
                if let Ok(file) = File::create(&output_path) {
                    let mut encoder = GzEncoder::new(file, Compression::default());
                    let _ = encoder
                        .write_all(content.as_bytes())
                        .and_then(|_| encoder.finish().map(|_| ()));
                }
            }

            Sample {
                name: sample_name,
                coverage,
                haploid_coverage,
            }
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
/// Uses copy-number model where haploid_coverage = single-copy coverage
///
/// Copy-number interpretation:
/// - Ploidy 1: 0 copies = 0x haploid, 1 copy = 1x haploid
/// - Ploidy 2: 0 copies = 0x haploid, 1 copy = 1x haploid, 2 copies = 2x haploid
///
/// For ploidy 1 with 1 threshold [t]:
/// - rel_cov < t → 0 (absent)
/// - rel_cov >= t → 1 (present)
///
/// For ploidy 1 with 2 thresholds [lower, upper]:
/// - rel_cov < lower → 0
/// - lower <= rel_cov < upper → missing (no-call)
/// - rel_cov >= upper → 1
///
/// For ploidy 2 with 2 thresholds [t1, t2] (default: 0.5, 1.5):
/// - rel_cov < t1 → 0/0 (0 copies, coverage ~0x haploid)
/// - t1 <= rel_cov < t2 → 0/1 (1 copy, coverage ~1x haploid)
/// - rel_cov >= t2 → 1/1 (2 copies, coverage ~2x haploid)
///
/// For ploidy 2 with 4 thresholds [t1, t2, t3, t4]:
/// - rel_cov < t1 → 0/0
/// - t1 <= rel_cov < t2 → missing (no-call)
/// - t2 <= rel_cov < t3 → 0/1
/// - t3 <= rel_cov < t4 → missing (no-call)
/// - rel_cov >= t4 → 1/1
///
/// where rel_cov = coverage / haploid_coverage
fn call_zygosity(
    coverage: f64,
    haploid_coverage: f64,
    ploidy: u8,
    thresholds: &[f64],
) -> Zygosity {
    let rel_cov = coverage / haploid_coverage;

    match (ploidy, thresholds.len()) {
        // Ploidy 1, simple mode: single threshold
        (1, 1) => {
            if rel_cov >= thresholds[0] {
                Zygosity::HomAlt // 1
            } else {
                Zygosity::HomRef // 0
            }
        }
        // Ploidy 1, no-call mode: lower and upper thresholds
        (1, 2) => {
            if rel_cov < thresholds[0] {
                Zygosity::HomRef // 0
            } else if rel_cov < thresholds[1] {
                Zygosity::Missing // uncertain, no-call
            } else {
                Zygosity::HomAlt // 1
            }
        }
        // Ploidy 2, simple mode: het_lower and het_upper
        (2, 2) => {
            if rel_cov < thresholds[0] {
                Zygosity::HomRef // 0/0
            } else if rel_cov < thresholds[1] {
                Zygosity::Het // 0/1
            } else {
                Zygosity::HomAlt // 1/1
            }
        }
        // Ploidy 2, no-call mode: 4 boundaries
        (2, 4) => {
            if rel_cov < thresholds[0] {
                Zygosity::HomRef // 0/0
            } else if rel_cov < thresholds[1] {
                Zygosity::Missing // uncertain between 0/0 and 0/1
            } else if rel_cov < thresholds[2] {
                Zygosity::Het // 0/1
            } else if rel_cov < thresholds[3] {
                Zygosity::Missing // uncertain between 0/1 and 1/1
            } else {
                Zygosity::HomAlt // 1/1
            }
        }
        _ => panic!(
            "Invalid combination: ploidy={}, threshold_count={}",
            ploidy,
            thresholds.len()
        ),
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
    thresholds: &[f64],
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
            let coverage = sample.coverage[feature_idx];
            let zygosity = call_zygosity(
                coverage,
                sample.haploid_coverage,
                ploidy,
                thresholds,
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
    thresholds: &[f64],
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
            let coverage = sample.coverage[feature_idx];
            let zygosity = call_zygosity(
                coverage,
                sample.haploid_coverage,
                ploidy,
                thresholds,
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
    thresholds: &[f64],
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
            let coverage = sample.coverage[feature_idx];
            let zygosity = call_zygosity(
                coverage,
                sample.haploid_coverage,
                ploidy,
                thresholds,
            );

            write!(file, ",{}", zygosity.to_dosage(ploidy))?;
        }
        writeln!(file)?;
    }

    Ok(())
}

/// Write copy-number matrix (continuous relative coverage values)
fn write_copynumber_matrix(
    samples: &[Sample],
    output_path: &str,
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

    // Write relative coverage for each feature
    for feature_idx in 0..num_features {
        write!(file, "{}", feature_idx + 1)?;

        for sample in samples {
            let coverage = sample.coverage[feature_idx];
            let rel_cov = coverage / sample.haploid_coverage;
            write!(file, "\t{:.3}", rel_cov)?;
        }
        writeln!(file)?;
    }

    Ok(())
}

/// Compute mean coverage (only values above threshold)
fn compute_mean(values: &[f64], threshold: f64) -> f64 {
    let (sum, count) = values
        .iter()
        .filter(|&&x| x > threshold)
        .fold((0.0, 0), |(sum, count), &x| (sum + x, count + 1));
    if count == 0 { 0.0 } else { sum / count as f64 }
}

/// Compute median coverage (only values above threshold)
fn compute_median(values: &[f64], threshold: f64) -> f64 {
    let mut filtered: Vec<f64> = values.iter().filter(|&&x| x > threshold).copied().collect();
    if filtered.is_empty() {
        return 0.0;
    }
    filtered.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mid = filtered.len() / 2;
    if filtered.len().is_multiple_of(2) {
        (filtered[mid - 1] + filtered[mid]) / 2.0
    } else {
        filtered[mid]
    }
}
