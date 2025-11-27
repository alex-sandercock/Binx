use binx_dosage::run_dosage;
use std::env;

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: binx-dosage <csv_file> <ploidy> [--verbose]");
        eprintln!("CSV Format: Alternating lines of Ref counts and Total counts per locus.");
        std::process::exit(1);
    }

    let filepath = &args[1];
    let ploidy: usize = args[2].parse().expect("Ploidy must be an integer");
    let verbose = args.contains(&"--verbose".to_string());

    run_dosage(filepath, ploidy, verbose)
}
