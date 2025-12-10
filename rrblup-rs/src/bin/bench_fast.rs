//! Benchmark comparing mixed_solve vs mixed_solve_fast
//!
//! Run with: cargo run --bin bench_fast --release

use nalgebra::DMatrix;
use rand::Rng;
use rrblup_rs::{mixed_solve_reml, mixed_solve_fast, MixedSolveOptions};
use std::time::Instant;

fn generate_test_data(n: usize, seed: u64) -> (Vec<f64>, DMatrix<f64>) {
    let mut rng = rand::thread_rng();

    // Generate y with some random values
    let y: Vec<f64> = (0..n).map(|i| (i as f64) + rng.gen::<f64>()).collect();

    // Generate positive definite K (kinship-like)
    let mut k = DMatrix::<f64>::zeros(n, n);
    for i in 0..n {
        k[(i, i)] = 1.0;
        for j in (i+1)..n {
            let cov = 0.3 * (-0.01 * ((i as f64 - j as f64).powi(2))).exp();
            k[(i, j)] = cov;
            k[(j, i)] = cov;
        }
    }

    (y, k)
}

fn benchmark_size(n: usize, runs: usize) {
    println!("\nn = {}:", n);
    println!("{:-<40}", "");

    let (y, k) = generate_test_data(n, 42);
    let opts = MixedSolveOptions::default();

    // Warmup
    let _ = mixed_solve_reml(&y, None, Some(&k), None, Some(opts.clone()));
    let _ = mixed_solve_fast(&y, None, Some(&k), None, Some(opts.clone()));

    // Benchmark original
    let mut times_orig = Vec::with_capacity(runs);
    for _ in 0..runs {
        let start = Instant::now();
        let _ = mixed_solve_reml(&y, None, Some(&k), None, Some(opts.clone()));
        times_orig.push(start.elapsed().as_secs_f64() * 1000.0);
    }
    let avg_orig: f64 = times_orig.iter().sum::<f64>() / (runs as f64);

    // Benchmark fast
    let mut times_fast = Vec::with_capacity(runs);
    for _ in 0..runs {
        let start = Instant::now();
        let _ = mixed_solve_fast(&y, None, Some(&k), None, Some(opts.clone()));
        times_fast.push(start.elapsed().as_secs_f64() * 1000.0);
    }
    let avg_fast: f64 = times_fast.iter().sum::<f64>() / (runs as f64);

    let speedup = avg_orig / avg_fast;

    println!("  mixed_solve (nalgebra): {:.2} ms", avg_orig);
    println!("  mixed_solve_fast (faer): {:.2} ms", avg_fast);
    println!("  Speedup: {:.2}x", speedup);

    // Verify results match
    let result_orig = mixed_solve_reml(&y, None, Some(&k), None, Some(opts.clone())).unwrap();
    let result_fast = mixed_solve_fast(&y, None, Some(&k), None, Some(opts.clone())).unwrap();

    let vu_diff = (result_orig.vu - result_fast.vu).abs();
    let ve_diff = (result_orig.ve - result_fast.ve).abs();
    let beta_diff = (result_orig.beta[0] - result_fast.beta[0]).abs();

    if vu_diff < 1e-5 && ve_diff < 1e-5 && beta_diff < 1e-5 {
        println!("  Results match: OK");
    } else {
        println!("  Results differ! vu_diff={:.2e}, ve_diff={:.2e}, beta_diff={:.2e}",
                 vu_diff, ve_diff, beta_diff);
    }
}

fn main() {
    println!("mixed_solve vs mixed_solve_fast Benchmark");
    println!("==========================================");
    println!("\nComparing nalgebra::SymmetricEigen vs faer eigendecomposition");

    let runs = 3;

    for &n in &[50, 100, 200, 500, 1000] {
        benchmark_size(n, runs);
    }

    println!("\nDone.");
}
