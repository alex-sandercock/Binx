use anyhow::{anyhow, bail, Context, Result};
use csv::StringRecord;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::NamedTempFile;

#[derive(Debug, Clone)]
struct MarkerRow {
    effect: Option<f64>,
    p_value: Option<f64>,
}

fn project_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("binx-cli has a parent dir")
        .to_path_buf()
}

fn repo_path(rel: &str) -> PathBuf {
    project_root().join(rel)
}

fn parse_f64_opt(raw: &str) -> Option<f64> {
    let trimmed = raw.trim();
    if trimmed.is_empty() || trimmed.eq_ignore_ascii_case("na") {
        None
    } else {
        trimmed.parse::<f64>().ok()
    }
}

fn find_col(headers: &StringRecord, names: &[&str]) -> Option<usize> {
    headers.iter().position(|h| names.iter().any(|n| h.eq_ignore_ascii_case(n)))
}

fn load_table(path: &Path) -> Result<HashMap<String, MarkerRow>> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("opening {}", path.display()))?;
    let headers = rdr
        .headers()
        .with_context(|| format!("reading headers from {}", path.display()))?
        .clone();

    let marker_idx = find_col(
        &headers,
        &["marker_id", "marker", "Marker", "SNP", "Snp", "ID", "id"],
    )
    .ok_or_else(|| anyhow!("marker column missing in {}", path.display()))?;
    let effect_idx = find_col(
        &headers,
        &["effect", "Effect", "beta", "Beta", "EffB", "additive", "Additive"],
    );
    let p_idx = find_col(&headers, &["p_value", "pValue", "P", "p", "pval", "p-val", "pvalue"])
        .ok_or_else(|| anyhow!("p-value column missing in {}", path.display()))?;

    let mut map = HashMap::new();
    for rec in rdr.records() {
        let rec = rec?;
        let marker = rec
            .get(marker_idx)
            .map(|s| s.to_string())
            .unwrap_or_default();
        if marker.is_empty() {
            continue;
        }
        let effect = effect_idx.and_then(|idx| rec.get(idx)).and_then(parse_f64_opt);
        let p_value = rec.get(p_idx).and_then(parse_f64_opt);
        map.insert(marker, MarkerRow { effect, p_value });
    }
    Ok(map)
}

fn compare_tables(
    binx_path: &Path,
    ref_path: &Path,
    epsilon_p: f64,
    epsilon_beta: f64,
) -> Result<()> {
    let binx = load_table(binx_path)?;
    let gwaspoly = load_table(ref_path)?;

    if binx.is_empty() || gwaspoly.is_empty() {
        bail!(
            "Empty table(s): binx rows={}, gwaspoly rows={}",
            binx.len(),
            gwaspoly.len()
        );
    }

    let binx_markers: std::collections::HashSet<_> = binx.keys().cloned().collect();
    let ref_markers: std::collections::HashSet<_> = gwaspoly.keys().cloned().collect();

    let missing_in_binx: Vec<_> = ref_markers.difference(&binx_markers).cloned().collect();
    let missing_in_ref: Vec<_> = binx_markers.difference(&ref_markers).cloned().collect();

    if !missing_in_binx.is_empty() || !missing_in_ref.is_empty() {
        bail!(
            "Marker mismatch: missing_in_binx={} missing_in_ref={}",
            missing_in_binx.len(),
            missing_in_ref.len()
        );
    }

    let mut max_delta_p = 0.0;
    let mut max_delta_beta = 0.0;
    for (marker, b) in binx.iter() {
        let r = gwaspoly
            .get(marker)
            .ok_or_else(|| anyhow::anyhow!("missing marker in ref: {}", marker))?;

        if let (Some(pb), Some(pr)) = (b.p_value, r.p_value) {
            if pb > 0.0 && pr > 0.0 {
                let dp = (pb.log10() - pr.log10()).abs();
                if dp > max_delta_p {
                    max_delta_p = dp;
                }
            }
        }
        if let (Some(bb), Some(br)) = (b.effect, r.effect) {
            let db = (bb - br).abs();
            if db > max_delta_beta {
                max_delta_beta = db;
            }
        }
    }

    if max_delta_p > epsilon_p {
        bail!(
            "p-value delta too large: max |Δ -log10(p)| = {} (eps = {})",
            max_delta_p,
            epsilon_p
        );
    }
    if max_delta_beta > epsilon_beta {
        bail!(
            "effect delta too large: max |Δ beta| = {} (eps = {})",
            max_delta_beta,
            epsilon_beta
        );
    }
    Ok(())
}

fn should_run_parity() -> bool {
    std::env::var("BINX_PARITY").map(|v| v == "1" || v == "true").unwrap_or(false)
}

fn ensure_fixture(path: &Path) -> bool {
    if path.exists() {
        true
    } else {
        eprintln!("Skipping parity: missing fixture {}", path.display());
        false
    }
}

fn run_binx_gwas(args: &[&str], out_path: &Path) -> Result<()> {
    let mut full_args = vec!["run", "-p", "binx-cli", "--", "gwas"];
    full_args.extend_from_slice(args);
    full_args.extend_from_slice(&["--out", out_path.to_str().unwrap()]);

    let status = Command::new("cargo")
        .current_dir(project_root())
        .args(&full_args)
        .status()
        .context("running binx gwas via cargo")?;
    if !status.success() {
        bail!("binx gwas command failed with status {}", status);
    }
    Ok(())
}

#[test]
fn parity_toy_additive() -> Result<()> {
    if !should_run_parity() {
        eprintln!("Skipping parity tests (set BINX_PARITY=1 to run)");
        return Ok(());
    }

    let geno = repo_path("tests/parity/data/toy/toy.geno.tsv");
    let pheno = repo_path("tests/parity/data/toy/toy.pheno.tsv");
    let reference = repo_path("tests/parity/fixtures/toy_additive.tsv");
    if !ensure_fixture(&reference) {
        return Ok(());
    }
    if !geno.exists() || !pheno.exists() {
        eprintln!("Skipping toy parity: missing geno/pheno fixtures");
        return Ok(());
    }

    let outfile = NamedTempFile::new().context("creating temp output")?;
    run_binx_gwas(
        &[
            "--geno",
            geno.to_str().unwrap(),
            "--pheno",
            pheno.to_str().unwrap(),
            "--ploidy",
            "4",
            "--trait-name",
            "Trait1",
            "--model",
            "additive",
        ],
        outfile.path(),
    )?;

    compare_tables(outfile.path(), &reference, 1e-4, 1e-6)
}

#[test]
fn parity_potato_additive_env() -> Result<()> {
    if !should_run_parity() {
        eprintln!("Skipping parity tests (set BINX_PARITY=1 to run)");
        return Ok(());
    }

    // Skip until binx-gwas supports repeated IDs / env fixed effects parity.
    if std::env::var("BINX_PARITY_POTATO")
        .map(|v| v == "1" || v == "true")
        .unwrap_or(false)
        == false
    {
        eprintln!("Skipping potato parity (set BINX_PARITY_POTATO=1 to enable; requires repeated-ID support)");
        return Ok(());
    }

    let geno = repo_path("tests/parity/data/potato/new_potato_geno.csv");
    let pheno = repo_path("tests/parity/data/potato/new_potato_pheno.csv");
    let reference = repo_path("tests/parity/fixtures/potato_additive_env.tsv");
    if !ensure_fixture(&reference) {
        return Ok(());
    }
    if !geno.exists() || !pheno.exists() {
        eprintln!("Skipping potato parity: missing geno/pheno fixtures");
        return Ok(());
    }

    let outfile = NamedTempFile::new().context("creating temp output")?;
    run_binx_gwas(
        &[
            "--geno",
            geno.to_str().unwrap(),
            "--pheno",
            pheno.to_str().unwrap(),
            "--ploidy",
            "4",
            "--trait-name",
            "vine.maturity",
            "--model",
            "additive",
        ],
        outfile.path(),
    )?;

    compare_tables(outfile.path(), &reference, 1e-4, 1e-6)
}
