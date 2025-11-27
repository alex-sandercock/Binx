use crate::math::{prob_ref_read, log_sum_exp, discretized_normal_probs, log_beta_binomial_pdf};
use crate::GenotypeResult;
use ndarray::Array1;

const MAX_ITER: usize = 200;
const TOLERANCE: f64 = 1e-6;

pub fn fit_norm(
    ref_counts: &Array1<u32>,
    total_counts: &Array1<u32>,
    ploidy: usize,
) -> anyhow::Result<GenotypeResult> {
    let n_samples = ref_counts.len();
    
    // --- 1. Initialization ---
    let mut bias = 1.0; 
    let mut rho = 0.01; // Initial overdispersion
    let seq_error = 0.005; // Still fixed for Phase 1
    let mut mu = (ploidy as f64) / 2.0;
    let mut sigma = 1.0; 

    let mut posteriors = vec![vec![0.0; ploidy + 1]; n_samples];
    let mut current_log_lik = f64::NEG_INFINITY;
    
    for _iter in 0..MAX_ITER {
        // --- 2. E-Step: Calculate Genotype Posteriors & LogLik ---
        let prior_probs = discretized_normal_probs(mu, sigma, ploidy);
        let mut new_log_lik = 0.0;

        for i in 0..n_samples {
            let r = ref_counts[i];
            let n = total_counts[i];
            
            let mut joint_log_probs = Vec::with_capacity(ploidy + 1);
            
            for k in 0..=ploidy {
                let p_xi = prob_ref_read(k, ploidy, seq_error, bias);
                let log_pmf = log_beta_binomial_pdf(r, n, p_xi, rho);
                let prior = prior_probs[k];
                let log_prior = if prior > 0.0 { prior.ln() } else { -1e9 };
                
                joint_log_probs.push(log_pmf + log_prior);
            }

            let log_marginal_i = log_sum_exp(&joint_log_probs);
            new_log_lik += log_marginal_i;

            for k in 0..=ploidy {
                posteriors[i][k] = (joint_log_probs[k] - log_marginal_i).exp();
            }
        }

        // --- Check Convergence ---
        if (new_log_lik - current_log_lik).abs() < TOLERANCE {
            current_log_lik = new_log_lik;
            break;
        }
        current_log_lik = new_log_lik;

        // --- 3. M-Step: Update Parameters ---
        
        // 3a. Update Mu/Sigma (Weighted Moments)
        let mut sum_w = 0.0;
        let mut sum_w_k = 0.0;

        for i in 0..n_samples {
            for k in 0..=ploidy {
                let w = posteriors[i][k];
                sum_w += w;
                sum_w_k += w * (k as f64);
            }
        }
        if sum_w < 1e-9 { sum_w = 1.0; } 

        mu = sum_w_k / sum_w;
        mu = mu.max(-5.0).min((ploidy as f64) + 5.0);

        let mut sum_w_sq_diff = 0.0;
        for i in 0..n_samples {
            for k in 0..=ploidy {
                let w = posteriors[i][k];
                let diff = (k as f64) - mu;
                sum_w_sq_diff += w * diff * diff;
            }
        }
        sigma = (sum_w_sq_diff / sum_w).sqrt();
        if sigma < 0.01 { sigma = 0.01; }

        // 3b. Update Bias and Rho (Coordinate Descent)
        for _sub_iter in 0..3 {
            bias = optimize_param_gss(
                bias, 0.1, 10.0, 
                |b| q_func(b, rho, ref_counts, total_counts, &posteriors, ploidy, seq_error)
            );

            rho = optimize_param_gss(
                rho, 0.0001, 0.5, 
                |r| q_func(bias, r, ref_counts, total_counts, &posteriors, ploidy, seq_error)
            );
        }
    }

    // --- Finalize ---
    let mut best_genotypes = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let mut best_k = 0;
        let mut max_p = -1.0;
        for k in 0..=ploidy {
            if posteriors[i][k] > max_p {
                max_p = posteriors[i][k];
                best_k = k;
            }
        }
        best_genotypes.push(best_k);
    }

    Ok(GenotypeResult {
        locus_id: "unknown".to_string(),
        genotype_probs: posteriors,
        best_genotypes,
        bias,
        seq_error,
        overdispersion: rho,
        model_mu: mu,
        model_sigma: sigma,
        final_log_lik: current_log_lik
    })
}

fn q_func(
    bias: f64,
    rho: f64,
    ref_counts: &Array1<u32>,
    total_counts: &Array1<u32>,
    posteriors: &[Vec<f64>], 
    ploidy: usize, 
    seq_error: f64
) -> f64 {
    let mut q = 0.0;
    for i in 0..ref_counts.len() {
        let r = ref_counts[i];
        let n = total_counts[i];
        for k in 0..=ploidy {
            let w = posteriors[i][k];
            if w < 1e-12 { continue; }

            let p_xi = prob_ref_read(k, ploidy, seq_error, bias);
            let log_lik = log_beta_binomial_pdf(r, n, p_xi, rho);
            
            q += w * log_lik;
        }
    }
    q
}

fn optimize_param_gss<F>(
    current_val: f64,
    min_bound: f64,
    max_bound: f64,
    mut func: F
) -> f64 
where F: FnMut(f64) -> f64 
{
    let gr = (5.0_f64.sqrt() - 1.0) / 2.0;
    let width = (max_bound - min_bound) * 0.5;
    let mut a = (current_val - width).max(min_bound);
    let mut b = (current_val + width).min(max_bound);

    let mut c = b - gr * (b - a);
    let mut d = a + gr * (b - a);
    
    let mut fc = func(c);
    let mut fd = func(d);
    
    for _ in 0..15 {
        if fc > fd {
            b = d;
            d = c;
            fd = fc;
            c = b - gr * (b - a);
            fc = func(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + gr * (b - a);
            fd = func(d);
        }
    }
    (b + a) / 2.0
}
