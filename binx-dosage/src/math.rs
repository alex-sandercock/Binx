use statrs::function::gamma::ln_gamma;

/// Calculate probability of observing a reference read given a genotype, error, and bias.
/// Based on Gerard et al. (2018)
///
/// Params:
/// - k: The dosage (number of ref alleles)
/// - ploidy: The organism ploidy
/// - epsilon: Sequencing error rate
/// - bias: Bias parameter (h).
pub fn prob_ref_read(k: usize, ploidy: usize, epsilon: f64, bias: f64) -> f64 {
    // 1. Dosage proportion
    let p = k as f64 / ploidy as f64;
    
    // 2. Apply Sequencing Error (f)
    // f = p(1-e) + (1-p)e
    let f = p * (1.0 - epsilon) + (1.0 - p) * epsilon;

    // 3. Apply Bias (xi) via Odds Ratio
    // xi = (f * h) / (f * h + (1 - f))
    let term = f * bias;
    term / (term + (1.0 - f))
}

/// Log PDF of the Beta-Binomial Distribution
/// 
/// Params:
/// - k: Number of successes (ref reads)
/// - n: Total trials (total reads)
/// - p: Probability of success (xi)
/// - rho: Overdispersion parameter (0 <= rho < 1)
pub fn log_beta_binomial_pdf(k: u32, n: u32, p: f64, rho: f64) -> f64 {
    let k = k as f64;
    let n = n as f64;

    // Handle edge cases
    if rho <= 1e-9 {
        // Converge to Binomial
        return log_binomial_pdf(k, n, p);
    }
    
    // alpha = p * (1/rho - 1)
    // beta = (1-p) * (1/rho - 1)
    let m = (1.0 / rho) - 1.0;
    let alpha = p * m;
    let beta = (1.0 - p) * m;

    let binom_const = ln_gamma(n + 1.0) - ln_gamma(k + 1.0) - ln_gamma(n - k + 1.0);
    
    let numer = ln_gamma(k + alpha) + ln_gamma(n - k + beta) - ln_gamma(n + m);
    let denom = ln_gamma(alpha) + ln_gamma(beta) - ln_gamma(m);

    binom_const + numer - denom
}

pub fn log_binomial_pdf(k: f64, n: f64, p: f64) -> f64 {
    if p == 0.0 {
        return if k == 0.0 { 0.0 } else { f64::NEG_INFINITY };
    }
    if p == 1.0 {
        return if k == n { 0.0 } else { f64::NEG_INFINITY };
    }
    let binom_const = ln_gamma(n + 1.0) - ln_gamma(k + 1.0) - ln_gamma(n - k + 1.0);
    binom_const + k * p.ln() + (n - k) * (1.0 - p).ln()
}

/// Log-Sum-Exp trick for numerical stability
pub fn log_sum_exp(vals: &[f64]) -> f64 {
    if vals.is_empty() { return f64::NEG_INFINITY; }
    let max_val = vals.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    if max_val == f64::NEG_INFINITY { return f64::NEG_INFINITY; }
    let sum = vals.iter().map(|v| (v - max_val).exp()).sum::<f64>();
    max_val + sum.ln()
}

/// Discretized Normal Distribution Probability
pub fn discretized_normal_probs(mu: f64, sigma: f64, ploidy: usize) -> Vec<f64> {
    let mut log_probs = Vec::with_capacity(ploidy + 1);
    let two_sigma_sq = 2.0 * sigma * sigma;
    
    for k in 0..=ploidy {
        let diff = k as f64 - mu;
        log_probs.push(-(diff * diff) / two_sigma_sq);
    }

    let denom = log_sum_exp(&log_probs);
    log_probs.iter().map(|lp| (lp - denom).exp()).collect()
}
