use crate::math::{prob_ref_read, log_sum_exp, discretized_normal_probs, log_beta_binomial_pdf};
use crate::GenotypeResult;
use ndarray::Array1;

const MAX_ITER: usize = 200;
const TOLERANCE: f64 = 1e-6;
const SPRINT_ITERS: usize = 8; // Short run for multi-start comparison
const PRIOR_MEAN_LOG_BIAS: f64 = 0.0;
const PRIOR_VAR_LOG_BIAS: f64 = 0.7 * 0.7; // Updog default 0.7^2
const PRIOR_MEAN_LOGIT_SEQ: f64 = -4.7;
const PRIOR_VAR_LOGIT_SEQ: f64 = 1.0;
const PRIOR_MEAN_LOGIT_RHO: f64 = -5.5;
const PRIOR_VAR_LOGIT_RHO: f64 = 0.5 * 0.5; // Updog default 0.25
const MIN_SIGMA: f64 = 1e-3;
const NEG_INF_FALLBACK: f64 = -1.0e30;
const LN_2PI: f64 = 1.8378770664093453; // ln(2*pi)

/// Mode for multi-start strategy
#[derive(Debug, Clone, Copy)]
pub enum FitMode {
    /// Hybrid sprint: 3 starts [0.5, 1.0, 2.0], run 8 iters, pick best (Default)
    Auto,
    /// Full Updog validation: 5 starts [0.37, 0.61, 1.0, 1.65, 2.72]
    Updog,
    /// Fast mode: single start at bias=1.0
    Fast,
}

/// Encapsulates the state of the EM algorithm for norm model fitting.
/// Allows pausing/resuming EM iterations for multi-start strategies.
pub struct EMState<'a> {
    // Zero-copy references to data
    ref_counts: &'a Array1<u32>,
    total_counts: &'a Array1<u32>,
    ploidy: usize,

    // Model parameters
    pub bias: f64,
    pub rho: f64,
    pub seq_error: f64,
    pub mu: f64,
    pub sigma: f64,

    // Posterior cache (allocated once, reused each iteration)
    pub posteriors: Vec<Vec<f64>>,

    // Current state
    pub current_log_lik: f64,
    pub iter_count: usize,
}

impl<'a> EMState<'a> {
    /// Create a new EM state with data-driven mu/sigma initialization and explicit bias
    pub fn new_with_bias(
        ref_counts: &'a Array1<u32>,
        total_counts: &'a Array1<u32>,
        ploidy: usize,
        bias: f64,
    ) -> Self {
        let n_samples = ref_counts.len();

        // Data-driven initialization for mu and sigma
        let total_reads: f64 = total_counts.iter().map(|&n| n as f64).sum();
        let total_refs: f64 = ref_counts.iter().map(|&r| r as f64).sum();
        let ref_prop = if total_reads > 0.0 {
            (total_refs / total_reads).clamp(0.0, 1.0)
        } else {
            0.5
        };

        let mu = (ref_prop * ploidy as f64).clamp(0.0, ploidy as f64);
        let sigma = (ploidy as f64 * ref_prop * (1.0 - ref_prop)).sqrt().max(0.1);

        // Allocate posterior cache
        let posteriors = vec![vec![0.0; ploidy + 1]; n_samples];

        EMState {
            ref_counts,
            total_counts,
            ploidy,
            bias,
            rho: 0.01,
            seq_error: 0.005,
            mu,
            sigma,
            posteriors,
            current_log_lik: f64::NEG_INFINITY,
            iter_count: 0,
        }
    }

    /// Run a single EM iteration: E-step + M-step
    /// Updates self.posteriors, self.current_log_lik, and all parameters
    pub fn em_step(&mut self) {
        // --- E-Step: Calculate Genotype Posteriors & LogLik ---
        let prior_probs = discretized_normal_probs(self.mu, self.sigma, self.ploidy);
        let mut new_log_lik = 0.0;

        for i in 0..self.ref_counts.len() {
            let r = self.ref_counts[i];
            let n = self.total_counts[i];

            let mut joint_log_probs = Vec::with_capacity(self.ploidy + 1);

            for k in 0..=self.ploidy {
                let p_xi = prob_ref_read(k, self.ploidy, self.seq_error, self.bias);
                let log_pmf = log_beta_binomial_pdf(r, n, p_xi, self.rho);
                let prior = prior_probs[k];
                let log_prior = if prior > 0.0 { prior.ln() } else { -1e9 };

                joint_log_probs.push(log_pmf + log_prior);
            }

            let log_marginal_i = log_sum_exp(&joint_log_probs);
            new_log_lik += log_marginal_i;

            for k in 0..=self.ploidy {
                self.posteriors[i][k] = (joint_log_probs[k] - log_marginal_i).exp();
            }
        }

        self.current_log_lik = new_log_lik;

        // --- M-Step: Update Parameters ---

        // 3a. Update Mu/Sigma via discrete normal objective
        let mut weights = vec![0.0; self.ploidy + 1];
        for i in 0..self.ref_counts.len() {
            for k in 0..=self.ploidy {
                weights[k] += self.posteriors[i][k];
            }
        }
        let (new_mu, new_sigma) = optimize_mu_sigma(self.mu, self.sigma, &weights, self.ploidy);
        self.mu = new_mu;
        self.sigma = new_sigma;

        // 3b. Update Bias and Rho (Coordinate Descent)
        for _sub_iter in 0..3 {
            self.bias = optimize_param_gss(
                self.bias,
                1e-3,
                10.0,
                |b| {
                    q_func(
                        b,
                        self.rho,
                        self.seq_error,
                        self.ref_counts,
                        self.total_counts,
                        &self.posteriors,
                        self.ploidy,
                    )
                },
            );

            self.rho = optimize_param_gss(
                self.rho,
                1e-6,
                0.999,
                |r| {
                    q_func(
                        self.bias,
                        r,
                        self.seq_error,
                        self.ref_counts,
                        self.total_counts,
                        &self.posteriors,
                        self.ploidy,
                    )
                },
            );

            self.seq_error = optimize_param_gss(
                self.seq_error,
                1e-5,
                0.05,
                |s| {
                    q_func(
                        self.bias,
                        self.rho,
                        s,
                        self.ref_counts,
                        self.total_counts,
                        &self.posteriors,
                        self.ploidy,
                    )
                },
            );
        }

        self.iter_count += 1;
    }

    /// Check if EM has converged
    pub fn has_converged(&self, prev_log_lik: f64) -> bool {
        (self.current_log_lik - prev_log_lik).abs() < TOLERANCE
    }

    /// Convert to final GenotypeResult
    pub fn into_result(self, locus_id: String) -> GenotypeResult {
        let n_samples = self.ref_counts.len();
        let mut best_genotypes = Vec::with_capacity(n_samples);

        for i in 0..n_samples {
            let mut best_k = 0;
            let mut max_p = -1.0;
            for k in 0..=self.ploidy {
                if self.posteriors[i][k] > max_p {
                    max_p = self.posteriors[i][k];
                    best_k = k;
                }
            }
            best_genotypes.push(best_k);
        }

        GenotypeResult {
            locus_id,
            genotype_probs: self.posteriors,
            best_genotypes,
            bias: self.bias,
            seq_error: self.seq_error,
            overdispersion: self.rho,
            model_mu: self.mu,
            model_sigma: self.sigma,
            final_log_lik: self.current_log_lik,
        }
    }
}

/// Fit norm model with hybrid sprint multi-start strategy (default)
pub fn fit_norm(
    ref_counts: &Array1<u32>,
    total_counts: &Array1<u32>,
    ploidy: usize,
) -> anyhow::Result<GenotypeResult> {
    fit_norm_with_mode(ref_counts, total_counts, ploidy, FitMode::Auto)
}

/// Fit norm model with configurable multi-start strategy
pub fn fit_norm_with_mode(
    ref_counts: &Array1<u32>,
    total_counts: &Array1<u32>,
    ploidy: usize,
    mode: FitMode,
) -> anyhow::Result<GenotypeResult> {
    match mode {
        FitMode::Auto => {
            // Hybrid sprint: 3 starts, 8 iters each, pick winner, run to convergence
            fit_norm_hybrid_sprint(ref_counts, total_counts, ploidy, vec![0.5, 1.0, 2.0])
        }
        FitMode::Updog => {
            // Full multi-start: ALL 5 starts to convergence, pick best (matches R package)
            fit_norm_full_multistart(ref_counts, total_counts, ploidy, vec![0.368, 0.607, 1.0, 1.649, 2.718])
        }
        FitMode::Fast => {
            // Single start, no sprint
            fit_norm_single_start(ref_counts, total_counts, ploidy, 1.0)
        }
    }
}

/// Run a single EM chain to convergence
fn fit_norm_single_start(
    ref_counts: &Array1<u32>,
    total_counts: &Array1<u32>,
    ploidy: usize,
    bias: f64,
) -> anyhow::Result<GenotypeResult> {
    let mut state = EMState::new_with_bias(ref_counts, total_counts, ploidy, bias);

    for _ in 0..MAX_ITER {
        let prev_log_lik = state.current_log_lik;
        state.em_step();

        if state.has_converged(prev_log_lik) {
            break;
        }
    }

    Ok(state.into_result("unknown".to_string()))
}

/// Hybrid sprint: Short sprint on all starts, pick winner, run to convergence
fn fit_norm_hybrid_sprint(
    ref_counts: &Array1<u32>,
    total_counts: &Array1<u32>,
    ploidy: usize,
    bias_starts: Vec<f64>,
) -> anyhow::Result<GenotypeResult> {
    let mut best_state: Option<EMState> = None;
    let mut max_log_lik = f64::NEG_INFINITY;

    // Sprint phase: run limited iterations on all starts
    for bias in &bias_starts {
        let mut state = EMState::new_with_bias(ref_counts, total_counts, ploidy, *bias);

        for _ in 0..SPRINT_ITERS {
            state.em_step();
        }

        if state.current_log_lik > max_log_lik {
            max_log_lik = state.current_log_lik;
            best_state = Some(state);
        }
    }

    let mut winner = best_state.ok_or_else(|| anyhow::anyhow!("No valid EM state found"))?;

    // Marathon: run winner to convergence
    for _ in 0..MAX_ITER {
        let prev_log_lik = winner.current_log_lik;
        winner.em_step();

        if winner.has_converged(prev_log_lik) {
            break;
        }
    }

    Ok(winner.into_result("unknown".to_string()))
}

/// Full multi-start: Run ALL starts to full convergence, pick best (True Updog behavior)
/// This matches the R package exactly: each bias init runs to convergence independently
fn fit_norm_full_multistart(
    ref_counts: &Array1<u32>,
    total_counts: &Array1<u32>,
    ploidy: usize,
    bias_starts: Vec<f64>,
) -> anyhow::Result<GenotypeResult> {
    let mut best_result: Option<GenotypeResult> = None;
    let mut max_log_lik = f64::NEG_INFINITY;

    // Run each start to full convergence (expensive but thorough)
    for bias in &bias_starts {
        let mut state = EMState::new_with_bias(ref_counts, total_counts, ploidy, *bias);

        // Run to convergence
        for _ in 0..MAX_ITER {
            let prev_log_lik = state.current_log_lik;
            state.em_step();

            if state.has_converged(prev_log_lik) {
                break;
            }
        }

        // Keep best converged result
        if state.current_log_lik > max_log_lik {
            max_log_lik = state.current_log_lik;
            best_result = Some(state.into_result("unknown".to_string()));
        }
    }

    best_result.ok_or_else(|| anyhow::anyhow!("No valid solution found from any start"))
}

fn q_func(
    bias: f64,
    rho: f64,
    seq_error: f64,
    ref_counts: &Array1<u32>,
    total_counts: &Array1<u32>,
    posteriors: &[Vec<f64>], 
    ploidy: usize
) -> f64 {
    if !is_valid_prob(seq_error) || !is_valid_prob(rho) || bias <= 0.0 {
        return NEG_INF_FALLBACK;
    }

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
    q + prior_penalty(bias, seq_error, rho)
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

fn optimize_mu_sigma(
    current_mu: f64,
    current_sigma: f64,
    weights: &[f64],
    ploidy: usize,
) -> (f64, f64) {
    let mut mu = current_mu;
    let mut sigma = current_sigma.max(MIN_SIGMA);
    for _ in 0..5 {
        mu = optimize_param_gss(
            mu,
            -1.0,
            (ploidy as f64) + 1.0,
            |m| lnorm_obj(m, sigma, weights, ploidy),
        );
        sigma = optimize_param_gss(
            sigma,
            MIN_SIGMA,
            (ploidy as f64) + 2.0,
            |s| lnorm_obj(mu, s, weights, ploidy),
        );
    }
    (mu, sigma)
}

fn lnorm_obj(mu: f64, sigma: f64, weights: &[f64], ploidy: usize) -> f64 {
    if sigma <= 0.0 {
        return NEG_INF_FALLBACK;
    }
    let mut log_probs = Vec::with_capacity(ploidy + 1);
    for k in 0..=ploidy {
        let diff = (k as f64) - mu;
        let logp = -0.5 * (LN_2PI + 2.0 * sigma.ln()) - (diff * diff) / (2.0 * sigma * sigma);
        log_probs.push(logp);
    }
    let lse = log_sum_exp(&log_probs);
    let sum_w: f64 = weights.iter().sum();
    if sum_w <= 0.0 {
        return NEG_INF_FALLBACK;
    }
    let mut obj = 0.0;
    for k in 0..=ploidy {
        obj += weights[k] * (log_probs[k] - lse);
    }
    obj
}

fn prior_penalty(bias: f64, seq_error: f64, rho: f64) -> f64 {
    pen_bias(bias, PRIOR_MEAN_LOG_BIAS, PRIOR_VAR_LOG_BIAS)
        + pen_seq_error(seq_error, PRIOR_MEAN_LOGIT_SEQ, PRIOR_VAR_LOGIT_SEQ)
        + pen_seq_error(rho, PRIOR_MEAN_LOGIT_RHO, PRIOR_VAR_LOGIT_RHO)
}

fn pen_bias(h: f64, mu_h: f64, var_h: f64) -> f64 {
    if h <= 0.0 {
        return NEG_INF_FALLBACK;
    }
    if !var_h.is_finite() || var_h <= 0.0 {
        return 0.0;
    }
    let diff = h.ln() - mu_h;
    -diff * diff / (2.0 * var_h)
}

fn pen_seq_error(eps: f64, mu_eps: f64, var_eps: f64) -> f64 {
    if !is_valid_prob(eps) {
        return NEG_INF_FALLBACK;
    }
    if !var_eps.is_finite() || var_eps <= 0.0 {
        return 0.0;
    }
    let logit = logit(eps);
    - (eps * (1.0 - eps)).ln() - (logit - mu_eps).powi(2) / (2.0 * var_eps)
}

fn is_valid_prob(x: f64) -> bool {
    x > 0.0 && x < 1.0
}

fn logit(x: f64) -> f64 {
    (x / (1.0 - x)).ln()
}
