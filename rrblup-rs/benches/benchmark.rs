//! Benchmark for rrblup-rs mixed_solve and kin_blup

use nalgebra::DMatrix;
use rrblup_rs::{
    kin_blup, mixed_solve_new, KinBlupData, KinBlupOptions,
};
use std::time::Instant;

fn generate_kinship(n: usize, seed: u64) -> DMatrix<f64> {
    use rand::{Rng, SeedableRng};
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

    let ncol = 30;
    let g: DMatrix<f64> = DMatrix::from_fn(n, ncol, |_, _| rng.gen::<f64>() * 2.0 - 1.0);
    let mut k = &g * g.transpose() / (ncol as f64);

    // Add small diagonal for positive definiteness
    for i in 0..n {
        k[(i, i)] += 0.1;
    }
    k
}

fn benchmark<F>(name: &str, mut f: F, times: usize)
where
    F: FnMut(),
{
    let mut timings = Vec::with_capacity(times);

    for _ in 0..times {
        let start = Instant::now();
        f();
        let elapsed = start.elapsed().as_secs_f64() * 1000.0;
        timings.push(elapsed);
    }

    let mean: f64 = timings.iter().sum::<f64>() / times as f64;
    let min = timings.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = timings.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    println!("{}: {:.2} ms (mean of {} runs, range: {:.2} - {:.2} ms)",
             name, mean, times, min, max);
}

fn main() {
    println!("rrblup-rs (Rust) Benchmark");
    println!("==========================\n");

    // ============================================================
    // mixed_solve with Z (marker matrix) - the common genomic case
    // ============================================================
    println!("mixed_solve with marker matrix Z (n genotypes x m markers):");
    println!("------------------------------------------------------------");

    // Varying markers with fixed genotypes
    let n = 500usize;
    for m in [1000, 5000, 10000, 20000] {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        // Marker matrix (centered)
        let z = DMatrix::from_fn(n, m, |_, _| {
            let val = match rng.gen_range(0..3) {
                0 => -1.0,
                1 => 0.0,
                _ => 1.0,
            };
            val
        });
        // Center columns
        let col_means: Vec<f64> = (0..m).map(|j| {
            (0..n).map(|i| z[(i, j)]).sum::<f64>() / n as f64
        }).collect();
        let z_centered = DMatrix::from_fn(n, m, |i, j| z[(i, j)] - col_means[j]);

        // Phenotype
        let beta: Vec<f64> = (0..m).map(|_| rng.gen::<f64>() * 0.02 - 0.01).collect();
        let y: Vec<f64> = (0..n).map(|i| {
            let zb: f64 = (0..m).map(|j| z_centered[(i, j)] * beta[j]).sum();
            zb + rng.gen::<f64>() * 0.5
        }).collect();

        let z_clone = z_centered.clone();
        let y_clone = y.clone();

        benchmark(&format!("n={}, m={} markers", n, m), move || {
            let _ = mixed_solve_new(
                &y_clone,
                Some(&z_clone),
                None,  // No K - use Z directly
                None,
                None,
            );
        }, 3);
    }

    println!();

    // Varying genotypes with fixed markers
    let m = 10000usize;
    for n in [100, 200, 500, 1000] {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        let z = DMatrix::from_fn(n, m, |_, _| {
            match rng.gen_range(0..3) {
                0 => -1.0,
                1 => 0.0,
                _ => 1.0,
            }
        });
        let col_means: Vec<f64> = (0..m).map(|j| {
            (0..n).map(|i| z[(i, j)]).sum::<f64>() / n as f64
        }).collect();
        let z_centered = DMatrix::from_fn(n, m, |i, j| z[(i, j)] - col_means[j]);

        let beta: Vec<f64> = (0..m).map(|_| rng.gen::<f64>() * 0.02 - 0.01).collect();
        let y: Vec<f64> = (0..n).map(|i| {
            let zb: f64 = (0..m).map(|j| z_centered[(i, j)] * beta[j]).sum();
            zb + rng.gen::<f64>() * 0.5
        }).collect();

        let z_clone = z_centered.clone();
        let y_clone = y.clone();

        benchmark(&format!("n={} genotypes, m={}", n, m), move || {
            let _ = mixed_solve_new(
                &y_clone,
                Some(&z_clone),
                None,
                None,
                None,
            );
        }, 3);
    }

    println!();

    // (K-based benchmarks skipped for speed - see earlier results)

    // ============================================================
    // mixed_solve with K (kinship matrix) benchmarks - reduced
    // ============================================================
    println!("mixed_solve with kinship matrix K:");
    println!("-----------------------------------");

    for n in [100, 500] {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        let p = 50usize;

        // Generate X matrix (intercept + p-1 random columns)
        let x = DMatrix::from_fn(n, p, |i, j| {
            if j == 0 { 1.0 } else { rng.gen::<f64>() * 2.0 - 1.0 }
        });

        // Generate K
        let k = generate_kinship(n, 42);

        // Generate y
        let beta: Vec<f64> = (0..p).map(|_| rng.gen::<f64>() * 2.0 - 1.0).collect();
        let u: Vec<f64> = (0..n).map(|_| rng.gen::<f64>() * 0.5).collect();
        let y: Vec<f64> = (0..n)
            .map(|i| {
                let xb: f64 = (0..p).map(|j| x[(i, j)] * beta[j]).sum();
                xb + u[i] + rng.gen::<f64>() * 0.3
            })
            .collect();

        let z = DMatrix::identity(n, n);

        let y_clone = y.clone();
        let z_clone = z.clone();
        let k_clone = k.clone();
        let x_clone = x.clone();

        benchmark(&format!("n={}", n), move || {
            let _ = mixed_solve_new(
                &y_clone,
                Some(&z_clone),
                Some(&k_clone),
                Some(&x_clone),
                None,
            );
        }, 3);
    }

    println!();

    // ============================================================
    // kin_blup benchmarks
    // ============================================================
    println!("kin_blup benchmarks:");
    println!("--------------------");

    for n_geno in [100, 200] {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        let n_obs = n_geno * 3;

        // Generate kinship
        let k = generate_kinship(n_geno, 42);
        let geno_names: Vec<String> = (0..n_geno).map(|i| format!("G{}", i + 1)).collect();

        // Generate phenotypes
        let geno_effects: Vec<f64> = (0..n_geno).map(|_| rng.gen::<f64>() * 2.0 - 1.0).collect();
        let obs_geno: Vec<String> = (0..n_obs)
            .map(|_| geno_names[rng.gen_range(0..n_geno)].clone())
            .collect();
        let y: Vec<f64> = obs_geno
            .iter()
            .map(|g| {
                let idx = geno_names.iter().position(|x| x == g).unwrap();
                5.0 + geno_effects[idx] * 0.5 + rng.gen::<f64>() * 0.3
            })
            .collect();

        // Create data
        let data = KinBlupData {
            geno_ids: obs_geno.clone(),
            pheno: y.clone(),
            fixed: None,
            covariates: None,
        };

        let options = KinBlupOptions {
            k: Some(k.clone()),
            k_ids: Some(geno_names.clone()),
            ..Default::default()
        };

        benchmark(&format!("n_geno={}, n_obs={}", n_geno, n_obs), move || {
            let _ = kin_blup(&data, Some(options.clone()));
        }, 3);
    }

    println!("\nDone.");
}
