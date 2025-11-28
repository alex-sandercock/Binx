use binx_dosage::{run_dosage, FitMode};
use std::env;

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: binx-dosage <csv_file> <ploidy> [--mode <auto|updog|updog-exact|fast>] [--verbose]");
        eprintln!("CSV Format: Alternating lines of Ref counts and Total counts per locus.");
        eprintln!();
        eprintln!("Options:");
        eprintln!("  --mode auto     Hybrid sprint (3 starts: 0.5, 1.0, 2.0) - Default");
        eprintln!("  --mode updog    Full validation (5 starts like R package) with Binx bounds");
        eprintln!("  --mode updog-exact Full validation with relaxed bounds matching Updog");
        eprintln!("  --mode fast     Single start at bias=1.0 (fastest)");
        eprintln!("  --verbose       Print detailed progress");
        std::process::exit(1);
    }

    let filepath = &args[1];
    let ploidy: usize = args[2].parse().expect("Ploidy must be an integer");

    // Parse --mode flag
    let mode = if let Some(mode_idx) = args.iter().position(|arg| arg == "--mode") {
        let mode_str = args.get(mode_idx + 1).expect("--mode requires an argument");
        match mode_str.as_str() {
            "auto" => FitMode::Auto,
            "updog" => FitMode::Updog,
            "fast" => FitMode::Fast,
            "updog-exact" => FitMode::UpdogExact,
            _ => {
                eprintln!("Invalid mode: {}. Use 'auto', 'updog', 'updog-exact', or 'fast'", mode_str);
                std::process::exit(1);
            }
        }
    } else {
        FitMode::Auto  // Default
    };

    let verbose = args.contains(&"--verbose".to_string());

    run_dosage(filepath, ploidy, mode, verbose)
}
