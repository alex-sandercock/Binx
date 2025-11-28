use binx_dosage::{run_dosage, FitMode, InputSource};
use std::env;

fn main() -> anyhow::Result<()> {
    let argv: Vec<String> = env::args().skip(1).collect();
    if argv.is_empty() {
        eprintln!("Usage:");
        eprintln!("  binx-dosage <csv_file> <ploidy> [--mode <auto|updog|updog-fast|updog-exact|fast>] [--verbose]");
        eprintln!("    CSV Format: Alternating lines of Ref counts and Total counts per locus.");
        eprintln!("  binx-dosage --counts --ref-path <ref_matrix.csv> --total-path <total_matrix.csv> <ploidy> [--mode ...] [--verbose]");
        eprintln!("    Matrix Format: markers in rows, samples in columns; first column is marker ID, header row has samples.");
        eprintln!();
        eprintln!("Options:");
        eprintln!("  --mode auto     Hybrid sprint (3 starts: 0.5, 1.0, 2.0) - Default");
        eprintln!("  --mode updog    Full validation (5 starts like R package) with Binx bounds");
        eprintln!("  --mode updog-fast Hybrid sprint with 5 Updog starts (faster, still thorough)");
        eprintln!("  --mode updog-exact Full validation with relaxed bounds matching Updog");
        eprintln!("  --mode fast     Single start at bias=1.0 (fastest)");
        eprintln!("  --verbose       Print detailed progress");
        std::process::exit(1);
    }

    let mut csv_file: Option<String> = None;
    let mut ref_path: Option<String> = None;
    let mut total_path: Option<String> = None;
    let mut ploidy: Option<usize> = None;
    let mut mode: FitMode = FitMode::Auto;
    let mut verbose = false;
    let mut counts_mode = false;

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

    let input = if counts_mode {
        let ref_path = ref_path.unwrap_or_else(|| {
            eprintln!("--counts requires --ref-path <path>");
            std::process::exit(1);
        });
        let total_path = total_path.unwrap_or_else(|| {
            eprintln!("--counts requires --total-path <path>");
            std::process::exit(1);
        });
        if csv_file.is_some() {
            eprintln!("Do not provide a CSV positional argument when using --counts");
            std::process::exit(1);
        }
        InputSource::RefTotalMatrices { ref_path, total_path }
    } else {
        let csv = csv_file.unwrap_or_else(|| {
            eprintln!("CSV file path is required when --counts is not set");
            std::process::exit(1);
        });
        InputSource::TwoLineCsv(csv)
    };

    run_dosage(input, ploidy, mode, verbose)
}
