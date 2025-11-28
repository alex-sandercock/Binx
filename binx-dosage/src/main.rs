use binx_dosage::{run_dosage, FitMode, InputSource};
use std::env;

fn main() -> anyhow::Result<()> {
    let argv: Vec<String> = env::args().skip(1).collect();
    if argv.is_empty() {
        eprintln!("Usage:");
        eprintln!("  binx-dosage <csv_file> <ploidy> [--mode <auto|updog|updog-fast|updog-exact|fast>] [--format <matrix|stats|beagle|vcf|plink|gwaspoly>] [--compress <none|gzip>] [--threads N] [--output <path>] [--verbose]");
        eprintln!("    CSV Format: Alternating lines of Ref counts and Total counts per locus.");
        eprintln!("  binx-dosage --counts --ref-path <ref_matrix.csv> --total-path <total_matrix.csv> <ploidy> [--mode ...] [--format ...] [--compress ...] [--threads N] [--output <path>] [--verbose]");
        eprintln!("    Matrix Format: markers in rows, samples in columns; first column is marker ID, header row has samples.");
        eprintln!("  binx-dosage --vcf <file.vcf[.gz]> <ploidy> [--chunk-size N] [--mode ...] [--format ...] [--compress ...] [--threads N] [--output <path>] [--verbose]");
        eprintln!("    VCF Format: uses FORMAT/AD (allele depths) per sample; streaming and gzipped supported.");
        eprintln!();
        eprintln!("Options:");
        eprintln!("  --mode auto     Hybrid sprint (3 starts: 0.5, 1.0, 2.0) - Default");
        eprintln!("  --mode updog    Full validation (5 starts like R package) with Binx bounds");
        eprintln!("  --mode updog-fast Hybrid sprint with 5 Updog starts (faster, still thorough)");
        eprintln!("  --mode updog-exact Full validation with relaxed bounds matching Updog");
        eprintln!("  --mode fast     Single start at bias=1.0 (fastest)");
        eprintln!("  --format matrix | stats | beagle | vcf | plink | gwaspoly (default: matrix)");
        eprintln!("  --compress none | gzip (default: none)");
        eprintln!("  --chunk-size N  Process VCF records in batches of N (default: 256; set 0 to stream one by one)");
        eprintln!("  --threads N     Limit Rayon worker threads (default: Rayon chooses)");
        eprintln!("  --output PATH   Write results to PATH instead of stdout");
        eprintln!("  --verbose       Print detailed progress");
        std::process::exit(1);
    }

    let mut csv_file: Option<String> = None;
    let mut ref_path: Option<String> = None;
    let mut total_path: Option<String> = None;
    let mut vcf_path: Option<String> = None;
    let mut ploidy: Option<usize> = None;
    let mut mode: FitMode = FitMode::Auto;
    let mut format = binx_dosage::OutputFormat::Matrix;
    let mut compress = binx_dosage::CompressMode::None;
    let mut verbose = false;
    let mut counts_mode = false;
    let mut chunk_size: Option<usize> = None;
    let mut threads: Option<usize> = None;
    let mut output: Option<String> = None;

    let mut idx = 0;
    while idx < argv.len() {
        match argv[idx].as_str() {
            "--counts" => {
                counts_mode = true;
                idx += 1;
            }
            "--ref" | "--ref-path" => {
                if idx + 1 >= argv.len() {
                    eprintln!("--ref-path requires a path argument");
                    std::process::exit(1);
                }
                ref_path = Some(argv[idx + 1].clone());
                idx += 2;
            }
            "--total" | "--total-path" => {
                if idx + 1 >= argv.len() {
                    eprintln!("--total-path requires a path argument");
                    std::process::exit(1);
                }
                total_path = Some(argv[idx + 1].clone());
                idx += 2;
            }
            "--vcf" => {
                if idx + 1 >= argv.len() {
                    eprintln!("--vcf requires a path argument");
                    std::process::exit(1);
                }
                vcf_path = Some(argv[idx + 1].clone());
                idx += 2;
            }
            "--mode" => {
                if idx + 1 >= argv.len() {
                    eprintln!("--mode requires an argument");
                    std::process::exit(1);
                }
                let mode_str = argv[idx + 1].as_str();
                mode = match mode_str {
                    "auto" => FitMode::Auto,
                    "updog" => FitMode::Updog,
                    "updog-fast" => FitMode::UpdogFast,
                    "fast" => FitMode::Fast,
                    "updog-exact" => FitMode::UpdogExact,
                    _ => {
                        eprintln!("Invalid mode: {}. Use 'auto', 'updog', 'updog-fast', 'updog-exact', or 'fast'", mode_str);
                        std::process::exit(1);
                    }
                };
                idx += 2;
            }
            "--verbose" => {
                verbose = true;
                idx += 1;
            }
            "--chunk-size" => {
                if idx + 1 >= argv.len() {
                    eprintln!("--chunk-size requires an integer argument");
                    std::process::exit(1);
                }
                chunk_size = Some(argv[idx + 1].parse().unwrap_or_else(|_| {
                    eprintln!("Invalid chunk-size value: {}", argv[idx + 1]);
                    std::process::exit(1);
                }));
                idx += 2;
            }
            "--threads" => {
                if idx + 1 >= argv.len() {
                    eprintln!("--threads requires an integer argument");
                    std::process::exit(1);
                }
                threads = Some(argv[idx + 1].parse().unwrap_or_else(|_| {
                    eprintln!("Invalid threads value: {}", argv[idx + 1]);
                    std::process::exit(1);
                }));
                idx += 2;
            }
            "--format" => {
                if idx + 1 >= argv.len() {
                    eprintln!("--format requires an argument");
                    std::process::exit(1);
                }
                format = match argv[idx + 1].as_str() {
                    "matrix" => binx_dosage::OutputFormat::Matrix,
                    "stats" => binx_dosage::OutputFormat::Stats,
                    "beagle" => binx_dosage::OutputFormat::Beagle,
                    "vcf" => binx_dosage::OutputFormat::Vcf,
                    "plink" => binx_dosage::OutputFormat::PlinkRaw,
                    "gwaspoly" => binx_dosage::OutputFormat::GwasPoly,
                    other => {
                        eprintln!("Invalid format: {}", other);
                        std::process::exit(1);
                    }
                };
                idx += 2;
            }
            "--compress" => {
                if idx + 1 >= argv.len() {
                    eprintln!("--compress requires an argument");
                    std::process::exit(1);
                }
                compress = match argv[idx + 1].as_str() {
                    "none" => binx_dosage::CompressMode::None,
                    "gzip" => binx_dosage::CompressMode::Gzip,
                    other => {
                        eprintln!("Invalid compress: {}", other);
                        std::process::exit(1);
                    }
                };
                idx += 2;
            }
            "--output" => {
                if idx + 1 >= argv.len() {
                    eprintln!("--output requires a path argument");
                    std::process::exit(1);
                }
                output = Some(argv[idx + 1].clone());
                idx += 2;
            }
            other => {
                // Positional arguments: first non-flag is CSV (legacy) or ploidy
                if csv_file.is_none() && !counts_mode {
                    csv_file = Some(other.to_string());
                } else if ploidy.is_none() {
                    ploidy = Some(other.parse().unwrap_or_else(|_| {
                        eprintln!("Invalid ploidy value: {}", other);
                        std::process::exit(1);
                    }));
                } else {
                    eprintln!("Unrecognized argument: {}", other);
                    std::process::exit(1);
                }
                idx += 1;
            }
        }
    }

    let ploidy = ploidy.unwrap_or_else(|| {
        eprintln!("Ploidy must be specified");
        std::process::exit(1);
    });

    // Enforce mutually exclusive input modes
    let modes_selected = counts_mode as u8 + (vcf_path.is_some() as u8) + (csv_file.is_some() as u8);
    if modes_selected == 0 {
        eprintln!("Provide one of: CSV positional path, --counts with ref/total matrices, or --vcf <path>");
        std::process::exit(1);
    }
    if modes_selected > 1 {
        eprintln!("Inputs are mutually exclusive: choose only one of CSV, --counts, or --vcf");
        std::process::exit(1);
    }

    let input = if let Some(vcf) = vcf_path {
        InputSource::Vcf { path: vcf, chunk_size }
    } else if counts_mode {
        let ref_path = ref_path.unwrap_or_else(|| {
            eprintln!("--counts requires --ref-path <path>");
            std::process::exit(1);
        });
        let total_path = total_path.unwrap_or_else(|| {
            eprintln!("--counts requires --total-path <path>");
            std::process::exit(1);
        });
        InputSource::RefTotalMatrices { ref_path, total_path }
    } else {
        let csv = csv_file.unwrap_or_else(|| {
            eprintln!("CSV file path is required when --counts and --vcf are not set");
            std::process::exit(1);
        });
        InputSource::TwoLineCsv(csv)
    };

    run_dosage(input, ploidy, mode, verbose, threads, output.as_deref(), format, compress)
}
