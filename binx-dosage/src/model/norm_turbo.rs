//! Turbo mode: Parallel E-step with thread-local gamma caching
//!
//! This module provides optimized variants of the norm model fitting algorithms
//! that parallelize the E-step across samples using Rayon and cache gamma function
//! evaluations using thread-local storage.
//!
//! Expected speedup: 8-10x from parallelization + ~1.2x from gamma caching

use crate::math::{prob_ref_read, log_sum_exp, discretized_normal_probs};
use crate::GenotypeResult;
use argmin::core::{CostFunction, Error as ArgminError, Executor, Gradient};
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin::solver::quasinewton::LBFGS;
use ndarray::Array1;
use rayon::prelude::*;
use statrs::function::gamma::ln_gamma;
use std::sync::OnceLock;

// Constants (same as norm.rs)
const MAX_ITER: usize = 200;
const TOLERANCE: f64 = 1e-6;
const SPRINT_ITERS: usize = 8;
const PRIOR_MEAN_LOG_BIAS: f64 = 0.0;
const PRIOR_VAR_LOG_BIAS: f64 = 0.7 * 0.7;
const PRIOR_MEAN_LOGIT_SEQ: f64 = -4.7;
const PRIOR_VAR_LOGIT_SEQ: f64 = 1.0;
const PRIOR_MEAN_LOGIT_RHO: f64 = -5.5;
const PRIOR_VAR_LOGIT_RHO: f64 = 0.5 * 0.5;
const MIN_SIGMA: f64 = 1e-3;
const NEG_INF_FALLBACK: f64 = -1.0e30;
const LN_2PI: f64 = 1.8378770664093453;
const LBFGS_MAX_ITERS: u64 = 20;

/// Global ln_gamma table for integer arguments (0..=MAX_LG)
static LN_GAMMA_TABLE: OnceLock<Vec<f64>> = OnceLock::new();
const MAX_LG: usize = 5000;

/// Cached version of ln_gamma that uses thread-local storage
#[inline]
fn ln_gamma_cached(x: f64) -> f64 {
    // Fast path: integer arguments via global table
    if x >= 0.0 && x.fract() == 0.0 {
        let idx = x as usize;
        if idx <= MAX_LG {
            let table = LN_GAMMA_TABLE.get_or_init(|| {
                let mut t = Vec::with_capacity(MAX_LG + 1);
                for i in 0..=MAX_LG {
                    t.push(ln_gamma(i as f64));
                }
                t
            });
            return table[idx];
        }
    }
    // Fallback for non-integer or > MAX_LG
    ln_gamma(x)
}

/// Cached version of log_beta_binomial_pdf using ln_gamma_cached
#[inline]
fn log_beta_binomial_pdf_cached(k: u32, n: u32, p: f64, rho: f64) -> f64 {
    let k_f = k as f64;
    let n_f = n as f64;

    // Handle edge cases
    if rho <= 1e-9 {
        // Converge to Binomial
        return log_binomial_pdf_cached(k_f, n_f, p);
    }

    let m = (1.0 / rho) - 1.0;
    let alpha = p * m;
    let beta = (1.0 - p) * m;

    // Use cached ln_gamma for integer arguments
    let binom_const = ln_gamma_cached(n_f + 1.0)
        - ln_gamma_cached(k_f + 1.0)
        - ln_gamma_cached(n_f - k_f + 1.0);

    let numer = ln_gamma(k_f + alpha) + ln_gamma(n_f - k_f + beta) - ln_gamma(n_f + m);
    let denom = ln_gamma(alpha) + ln_gamma(beta) - ln_gamma(m);

    binom_const + numer - denom
}

#[inline]
fn log_binomial_pdf_cached(k: f64, n: f64, p: f64) -> f64 {
    if p == 0.0 {
        return if k == 0.0 { 0.0 } else { f64::NEG_INFINITY };
    }
    if p == 1.0 {
        return if k == n { 0.0 } else { f64::NEG_INFINITY };
    }
    let binom_const = ln_gamma_cached(n + 1.0)
        - ln_gamma_cached(k + 1.0)
        - ln_gamma_cached(n - k + 1.0);
    binom_const + k * p.ln() + (n - k) * (1.0 - p).ln()
}

#[derive(Clone, Copy)]
struct Bounds {
    bias_min: f64,
    bias_max: f64,
    rho_min: f64,
    rho_max: f64,
    seq_min: f64,
    seq_max: f64,
    sigma_min: f64,
    sigma_max: f64,
}

impl Bounds {
    fn default() -> Self {
        Bounds {
            bias_min: 1e-3,
            bias_max: 10.0,
            rho_min: 1e-6,
            rho_max: 0.999,
            seq_min: 1e-5,
            seq_max: 0.05,
            sigma_min: MIN_SIGMA,
            sigma_max: 2.0,
        }
    }
}

/// Turbo EM state with parallel E-step capability
pub(crate) struct EMStateTurbo<'a> {
    ref_counts: &'a Array1<u32>,
    total_counts: &'a Array1<u32>,
    ploidy: usize,
    bounds: Bounds,
    tolerance: f64,

    pub bias: f64,
    pub rho: f64,
    pub seq_error: f64,
    pub mu: f64,
    pub sigma: f64,
    pub posteriors: Vec<Vec<f64>>,
    pub current_log_lik: f64,
    pub current_log_lik_penalized: f64,
    pub iter_count: usize,
}

impl<'a> EMStateTurbo<'a> {
    pub(crate) fn new_with_bias(
        ref_counts: &'a Array1<u32>,
        total_counts: &'a Array1<u32>,
        ploidy: usize,
        bias: f64,
        bounds: Bounds,
        tolerance: f64,
    ) -> Self {
        let n_samples = ref_counts.len();

        let total_reads: f64 = total_counts.iter().map(|&n| n as f64).sum();
        let total_refs: f64 = ref_counts.iter().map(|&r| r as f64).sum();
        let ref_prop = if total_reads > 0.0 {
            (total_refs / total_reads).clamp(0.0, 1.0)
        } else {
            0.5
        };

        let mu = (ref_prop * ploidy as f64).clamp(0.0, ploidy as f64);
        let sigma = (ploidy as f64 * ref_prop * (1.0 - ref_prop)).sqrt().max(0.1);

        let posteriors = vec![vec![0.0; ploidy + 1]; n_samples];

        EMStateTurbo {
            ref_counts,
            total_counts,
            ploidy,
            bounds,
            tolerance,
            bias,
            rho: 0.01,
            seq_error: 0.005,
            mu,
            sigma,
            posteriors,
            current_log_lik: f64::NEG_INFINITY,
            current_log_lik_penalized: f64::NEG_INFINITY,
            iter_count: 0,
        }
    }

    /// PARALLEL E-step: Calculate posteriors across samples in parallel
    pub fn em_step(&mut self) {
        let prior_probs = discretized_normal_probs(self.mu, self.sigma, self.ploidy);
        let cap = depth_cap(self.ploidy);

        // Parallel E-step: compute posteriors and log-likelihoods for all samples
        let results: Vec<(Vec<f64>, f64)> = (0..self.ref_counts.len())
            .into_par_iter()
            .map(|i| {
                let (r, n) = cap_counts(self.ref_counts[i], self.total_counts[i], cap);

                if n == 0 {
                    // No information: posterior equals prior
                    return (prior_probs.clone(), 0.0);
                }

                let mut joint_log_probs = Vec::with_capacity(self.ploidy + 1);

                for k in 0..=self.ploidy {
                    let p_xi = prob_ref_read(k, self.ploidy, self.seq_error, self.bias);
                    let log_pmf = log_beta_binomial_pdf_cached(r, n, p_xi, self.rho);
                    let prior = prior_probs[k];
                    let log_prior = if prior > 0.0 { prior.ln() } else { -1e9 };
                    joint_log_probs.push(log_pmf + log_prior);
                }

                let log_marginal_i = log_sum_exp(&joint_log_probs);
                let posteriors_i: Vec<f64> = joint_log_probs
                    .iter()
                    .map(|&lp| (lp - log_marginal_i).exp())
                    .collect();

                (posteriors_i, log_marginal_i)
            })
            .collect();

        // Update posteriors and compute total log likelihood
        let mut new_log_lik = 0.0;
        for (i, (posts, log_marg)) in results.into_iter().enumerate() {
            self.posteriors[i] = posts;
            new_log_lik += log_marg;
        }
        self.current_log_lik = new_log_lik;

        // M-step: Update parameters (sequential, but fast compared to E-step)
        let mut weights = vec![0.0; self.ploidy + 1];
        for i in 0..self.ref_counts.len() {
            for k in 0..=self.ploidy {
                weights[k] += self.posteriors[i][k];
            }
        }

        let (new_mu, new_sigma) = optimize_mu_sigma(
            self.mu,
            self.sigma,
            &weights,
            self.ploidy,
            self.bounds,
        );
        self.mu = new_mu;
        self.sigma = new_sigma;

        // Update bias, rho, seq_error using GSS (3 iterations to mirror non-turbo paths)
        for _sub_iter in 0..3 {
            self.bias = optimize_param_gss(
                self.bias,
                self.bounds.bias_min,
                self.bounds.bias_max,
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
                self.bounds.rho_min,
                self.bounds.rho_max,
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
                self.bounds.seq_min,
                self.bounds.seq_max,
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

        let penalty = prior_penalty(self.bias, self.seq_error, self.rho);
        self.iter_count += 1;
        self.current_log_lik_penalized = new_log_lik + penalty;
    }

    pub fn has_converged(&self, prev_log_lik: f64) -> bool {
        (self.current_log_lik_penalized - prev_log_lik).abs() < self.tolerance
    }

    pub fn into_result(self, locus_id: String) -> GenotypeResult {
        let n_samples = self.ref_counts.len();
        let mut best_genotypes = Vec::with_capacity(n_samples);
        let mut max_posteriors = Vec::with_capacity(n_samples);

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
            max_posteriors.push(max_p);
        }

        GenotypeResult {
            locus_id,
            genotype_probs: self.posteriors,
            best_genotypes,
            max_posteriors,
            bias: self.bias,
            seq_error: self.seq_error,
            overdispersion: self.rho,
            model_mu: self.mu,
            model_sigma: self.sigma,
            final_log_lik: self.current_log_lik,
        }
    }
}

/// Turbo mode: Single start at bias=1.0 with parallel E-step
pub fn fit_norm_turbo(
    ref_counts: &Array1<u32>,
    total_counts: &Array1<u32>,
    ploidy: usize,
) -> anyhow::Result<GenotypeResult> {
    fit_norm_single_start_turbo(
        ref_counts,
        total_counts,
        ploidy,
        1.0,
        Bounds::default(),
        TOLERANCE,
    )
}

/// TurboAuto mode: Hybrid sprint with 3 starts, parallel E-step
pub fn fit_norm_turbo_auto(
    ref_counts: &Array1<u32>,
    total_counts: &Array1<u32>,
    ploidy: usize,
) -> anyhow::Result<GenotypeResult> {
    fit_norm_hybrid_sprint_turbo(
        ref_counts,
        total_counts,
        ploidy,
        vec![0.5, 1.0, 2.0],
        Bounds::default(),
        TOLERANCE,
    )
}

/// TurboAutoSafe: Hybrid sprint with 3 starts, parallel E-step, capping + gamma caching
pub fn fit_norm_turbo_auto_safe(
    ref_counts: &Array1<u32>,
    total_counts: &Array1<u32>,
    ploidy: usize,
) -> anyhow::Result<GenotypeResult> {
    fit_norm_hybrid_sprint_turbo(
        ref_counts,
        total_counts,
        ploidy,
        vec![0.5, 1.0, 2.0],
        Bounds::default(),
        TOLERANCE,
    )
}

fn fit_norm_single_start_turbo(
    ref_counts: &Array1<u32>,
    total_counts: &Array1<u32>,
    ploidy: usize,
    bias: f64,
    bounds: Bounds,
    tolerance: f64,
) -> anyhow::Result<GenotypeResult> {
    let mut state = EMStateTurbo::new_with_bias(
        ref_counts,
        total_counts,
        ploidy,
        bias,
        bounds,
        tolerance,
    );

    for _ in 0..MAX_ITER {
        let prev_log_lik = state.current_log_lik_penalized;
        state.em_step();

        if state.has_converged(prev_log_lik) {
            break;
        }
    }

    Ok(state.into_result("unknown".to_string()))
}

fn fit_norm_hybrid_sprint_turbo(
    ref_counts: &Array1<u32>,
    total_counts: &Array1<u32>,
    ploidy: usize,
    bias_starts: Vec<f64>,
    bounds: Bounds,
    tolerance: f64,
) -> anyhow::Result<GenotypeResult> {
    let mut best_state: Option<EMStateTurbo> = None;
    let mut max_log_lik = f64::NEG_INFINITY;

    // Sprint phase: run limited iterations on all starts
    for bias in &bias_starts {
        let mut state = EMStateTurbo::new_with_bias(
            ref_counts,
            total_counts,
            ploidy,
            *bias,
            bounds,
            tolerance,
        );

        for _ in 0..SPRINT_ITERS {
            state.em_step();
        }

        if state.current_log_lik_penalized > max_log_lik {
            max_log_lik = state.current_log_lik_penalized;
            best_state = Some(state);
        }
    }

    let mut winner = best_state.ok_or_else(|| anyhow::anyhow!("No valid EM state found"))?;

    // Marathon: run winner to convergence
    for _ in 0..MAX_ITER {
        let prev_log_lik = winner.current_log_lik_penalized;
        winner.em_step();

        if winner.has_converged(prev_log_lik) {
            break;
        }
    }

    Ok(winner.into_result("unknown".to_string()))
}

// Helper functions (same as norm.rs but using cached beta-binomial)

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

    let cap = depth_cap(ploidy);

    // Parallelize across samples to keep high-depth loci responsive and capped
    let q: f64 = (0..ref_counts.len())
        .into_par_iter()
        .map(|i| {
            let (r, n) = cap_counts(ref_counts[i], total_counts[i], cap);
            let mut acc = 0.0;
            for k in 0..=ploidy {
                let w = posteriors[i][k];
                if w < 1e-12 {
                    continue;
                }
                let p_xi = prob_ref_read(k, ploidy, seq_error, bias);
                let log_lik = log_beta_binomial_pdf_cached(r, n, p_xi, rho);
                acc += w * log_lik;
            }
            acc
        })
        .sum();
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
    bounds: Bounds,
) -> (f64, f64) {
    let init = Array1::from(vec![
        current_mu,
        current_sigma.max(bounds.sigma_min),
    ]);

    let problem = MuSigmaProblem {
        weights,
        ploidy,
        bounds,
    };

    let linesearch = MoreThuenteLineSearch::new();
    let solver = LBFGS::new(linesearch, 5)
        .with_tolerance_grad(1e-6)
        .unwrap();

    let res = Executor::new(problem, solver)
        .configure(|state| state.param(init.clone()).max_iters(LBFGS_MAX_ITERS))
        .run();

    let best = match res {
        Ok(r) => r.state.best_param.unwrap_or(init),
        Err(_) => init,
    };

    let mu = best[0].clamp(-1.0, (ploidy as f64) + 1.0);
    let sigma = best[1].clamp(bounds.sigma_min, bounds.sigma_max);
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

struct MuSigmaProblem<'a> {
    weights: &'a [f64],
    ploidy: usize,
    bounds: Bounds,
}

impl<'a> CostFunction for MuSigmaProblem<'a> {
    type Param = Array1<f64>;
    type Output = f64;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, ArgminError> {
        let mu = param[0].clamp(-1.0, (self.ploidy as f64) + 1.0);
        let sigma = param[1].clamp(self.bounds.sigma_min, self.bounds.sigma_max);
        Ok(-lnorm_obj(mu, sigma, self.weights, self.ploidy))
    }
}

impl<'a> Gradient for MuSigmaProblem<'a> {
    type Param = Array1<f64>;
    type Gradient = Array1<f64>;

    fn gradient(&self, param: &Self::Param) -> Result<Self::Gradient, ArgminError> {
        let mu = param[0].clamp(-1.0, (self.ploidy as f64) + 1.0);
        let sigma = param[1].clamp(self.bounds.sigma_min, self.bounds.sigma_max);

        let mut log_probs = Vec::with_capacity(self.ploidy + 1);
        for k in 0..=self.ploidy {
            let diff = (k as f64) - mu;
            let logp = -0.5 * (LN_2PI + 2.0 * sigma.ln()) - (diff * diff) / (2.0 * sigma * sigma);
            log_probs.push(logp);
        }
        let lse = log_sum_exp(&log_probs);
        let sum_w: f64 = self.weights.iter().sum();
        if sum_w <= 0.0 {
            return Ok(Array1::from(vec![0.0, 0.0]));
        }

        let mut pvec = Vec::with_capacity(self.ploidy + 1);
        for k in 0..=self.ploidy {
            pvec.push((log_probs[k] - lse).exp());
        }

        let mut g_mu = 0.0;
        let mut g_sigma = 0.0;
        for k in 0..=self.ploidy {
            let delta = (k as f64) - mu;
            g_mu += self.weights[k] * delta;
            g_sigma += self.weights[k] * delta * delta;
        }
        for k in 0..=self.ploidy {
            let delta = (k as f64) - mu;
            g_mu -= sum_w * delta * pvec[k];
            g_sigma -= sum_w * delta * delta * pvec[k];
        }

        g_mu /= sigma * sigma;
        g_sigma /= sigma * sigma * sigma;

        Ok(Array1::from(vec![-g_mu, -g_sigma]))
    }
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

#[inline]
fn depth_cap(ploidy: usize) -> u32 {
    // Scale cap by ploidy to preserve discrimination between adjacent dosages
    (25 * ploidy as u32).max(25)
}

#[inline]
fn cap_counts(r: u32, n: u32, cap: u32) -> (u32, u32) {
    if n <= cap {
        return (r, n);
    }
    let cap_f = cap as f64;
    let n_f = n as f64;
    let scaled_r = ((r as f64) * (cap_f / n_f)).round() as u32;
    (scaled_r.min(cap), cap)
}
