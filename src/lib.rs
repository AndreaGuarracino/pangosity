use log::debug;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// Input file format
#[derive(Debug, PartialEq)]
pub enum InputFormat {
    Pack, // Coverage file (PACK format)
    Gaf,  // GAF alignment file
}

/// Detect input format based on file extension and content
pub fn detect_format(path: &str) -> std::io::Result<InputFormat> {
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

/// Convert GAF to coverage (PACK format)
pub fn gaf_to_coverage(
    gaf_path: &str,
    gfa_path: &str,
    _sample_name: &str,
) -> std::io::Result<Vec<f64>> {
    let (coverage, _min_id) = gafpack::compute_coverage(gfa_path, gaf_path, true, false)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
    Ok(coverage)
}

/// Process input file and return coverage vector
/// Automatically detects format and applies necessary conversions
pub fn process_input_to_coverage(
    input_path: &str,
    gfa_path: Option<&str>,
    sample_name: &str,
) -> std::io::Result<Vec<f64>> {
    let format = detect_format(input_path)?;
    debug!("Detected format for {}: {:?}", input_path, format);

    match format {
        InputFormat::Pack => {
            // Read PACK format directly
            crate::parse_coverage_file(input_path).map(|(_, cov)| cov)
        }
        InputFormat::Gaf => {
            // Convert GAF to coverage
            let gfa = gfa_path.ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "GFA file required for GAF input",
                )
            })?;
            gaf_to_coverage(input_path, gfa, sample_name)
        }
    }
}

/// Parse a coverage file (from gafpack --coverage-column)
pub fn parse_coverage_file(path: &str) -> std::io::Result<(String, Vec<f64>)> {
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
