//! Survival analysis — Kaplan-Meier estimator, log-rank test, Cox proportional hazards.
//!
//! # Overview
//!
//! This module provides the three core methods of survival analysis:
//!
//! - [`kaplan_meier`] — non-parametric survival curve estimation with Greenwood
//!   standard errors and log-transform confidence intervals
//! - [`log_rank_test`] — chi-squared test comparing survival between groups
//! - [`cox_ph`] — semi-parametric Cox proportional hazards regression via
//!   Newton-Raphson on Breslow partial likelihood
//!
//! # Example
//!
//! ```
//! use cyanea_stats::survival::kaplan_meier;
//!
//! let times  = [1.0, 2.0, 3.0, 4.0, 5.0];
//! let status = [true, false, true, false, true]; // true = event
//! let km = kaplan_meier(&times, &status).unwrap();
//! assert!(km.steps[0].survival < 1.0);
//! ```

#![allow(non_snake_case)]

use cyanea_core::{CyaneaError, Result};

use crate::distribution::{ChiSquared, Distribution, Normal};

// ── Result types ────────────────────────────────────────────────────────────

/// A single step in the Kaplan-Meier survival curve.
#[derive(Debug, Clone)]
pub struct KmStep {
    /// Event time.
    pub time: f64,
    /// Number at risk just before this time.
    pub n_risk: usize,
    /// Number of events at this time.
    pub n_events: usize,
    /// Number censored in the interval ending at this time (exclusive of events).
    pub n_censored: usize,
    /// Product-limit survival estimate S(t).
    pub survival: f64,
    /// Greenwood standard error of S(t).
    pub std_err: f64,
    /// Lower bound of 95% confidence interval (log-transform method).
    pub ci_lower: f64,
    /// Upper bound of 95% confidence interval (log-transform method).
    pub ci_upper: f64,
}

/// Result of the Kaplan-Meier estimator.
#[derive(Debug, Clone)]
pub struct KmResult {
    /// Steps at each unique event time, in ascending order.
    pub steps: Vec<KmStep>,
    /// Median survival time: first time where S(t) ≤ 0.5, or `None` if S(t)
    /// never reaches 0.5.
    pub median_survival: Option<f64>,
    /// Total number of observations.
    pub n_total: usize,
    /// Total number of events.
    pub n_events: usize,
}

/// Result of the log-rank test.
#[derive(Debug, Clone)]
pub struct LogRankResult {
    /// Chi-squared test statistic.
    pub statistic: f64,
    /// Degrees of freedom (number of groups − 1).
    pub degrees_of_freedom: usize,
    /// p-value from the chi-squared distribution.
    pub p_value: f64,
    /// Number of groups compared.
    pub n_groups: usize,
    /// Observed number of events per group.
    pub observed: Vec<f64>,
    /// Expected number of events per group under H₀.
    pub expected: Vec<f64>,
}

/// Result of Cox proportional hazards regression.
#[derive(Debug, Clone)]
pub struct CoxPhResult {
    /// Estimated regression coefficients (β).
    pub coefficients: Vec<f64>,
    /// Standard errors of coefficients.
    pub std_errors: Vec<f64>,
    /// Wald z-statistics (β / SE).
    pub z_values: Vec<f64>,
    /// Two-sided p-values from standard normal.
    pub p_values: Vec<f64>,
    /// Hazard ratios exp(β).
    pub hazard_ratios: Vec<f64>,
    /// Lower 95% CI for hazard ratios.
    pub hr_ci_lower: Vec<f64>,
    /// Upper 95% CI for hazard ratios.
    pub hr_ci_upper: Vec<f64>,
    /// Number of observations.
    pub n_observations: usize,
    /// Number of events.
    pub n_events: usize,
    /// Number of covariates.
    pub n_covariates: usize,
    /// Log partial likelihood at convergence.
    pub log_likelihood: f64,
    /// Number of Newton-Raphson iterations performed.
    pub n_iterations: usize,
    /// Whether the optimizer converged.
    pub converged: bool,
}

// ── Helpers ─────────────────────────────────────────────────────────────────

/// Validate times/status arrays and return indices sorted by time (ascending).
fn validate_and_sort(times: &[f64], status: &[bool]) -> Result<Vec<usize>> {
    if times.is_empty() {
        return Err(CyaneaError::InvalidInput("empty input arrays".into()));
    }
    if times.len() != status.len() {
        return Err(CyaneaError::InvalidInput(format!(
            "times length ({}) != status length ({})",
            times.len(),
            status.len()
        )));
    }
    for (i, &t) in times.iter().enumerate() {
        if t < 0.0 {
            return Err(CyaneaError::InvalidInput(format!(
                "negative time {} at index {}",
                t, i
            )));
        }
        if t.is_nan() {
            return Err(CyaneaError::InvalidInput(format!(
                "NaN time at index {}",
                i
            )));
        }
    }
    let mut indices: Vec<usize> = (0..times.len()).collect();
    indices.sort_by(|&a, &b| times[a].partial_cmp(&times[b]).unwrap());
    Ok(indices)
}

/// Solve a linear system Ax = b via Gaussian elimination with partial pivoting.
/// `a` is row-major n×n, `b` is length n. Returns x.
fn solve_linear_system(a: &[f64], b: &[f64], n: usize) -> Result<Vec<f64>> {
    let mut aug = vec![0.0; n * (n + 1)];
    for i in 0..n {
        for j in 0..n {
            aug[i * (n + 1) + j] = a[i * n + j];
        }
        aug[i * (n + 1) + n] = b[i];
    }

    for col in 0..n {
        // Partial pivoting
        let mut max_row = col;
        let mut max_val = aug[col * (n + 1) + col].abs();
        for row in (col + 1)..n {
            let val = aug[row * (n + 1) + col].abs();
            if val > max_val {
                max_val = val;
                max_row = row;
            }
        }
        if max_val < 1e-15 {
            return Err(CyaneaError::InvalidInput(
                "singular matrix in linear solve".into(),
            ));
        }
        if max_row != col {
            for j in 0..=n {
                aug.swap(col * (n + 1) + j, max_row * (n + 1) + j);
            }
        }
        // Eliminate below
        let pivot = aug[col * (n + 1) + col];
        for row in (col + 1)..n {
            let factor = aug[row * (n + 1) + col] / pivot;
            for j in col..=n {
                let above = aug[col * (n + 1) + j];
                aug[row * (n + 1) + j] -= factor * above;
            }
        }
    }

    // Back-substitute
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = aug[i * (n + 1) + n];
        for j in (i + 1)..n {
            sum -= aug[i * (n + 1) + j] * x[j];
        }
        x[i] = sum / aug[i * (n + 1) + i];
    }
    Ok(x)
}

/// Invert a symmetric positive-definite matrix (row-major n×n) via Gaussian
/// elimination with partial pivoting. Returns the inverse in row-major order.
fn invert_matrix(a: &[f64], n: usize) -> Result<Vec<f64>> {
    // Augment with identity
    let mut aug = vec![0.0; n * 2 * n];
    for i in 0..n {
        for j in 0..n {
            aug[i * 2 * n + j] = a[i * n + j];
        }
        aug[i * 2 * n + n + i] = 1.0;
    }

    let cols = 2 * n;
    for col in 0..n {
        // Partial pivoting
        let mut max_row = col;
        let mut max_val = aug[col * cols + col].abs();
        for row in (col + 1)..n {
            let val = aug[row * cols + col].abs();
            if val > max_val {
                max_val = val;
                max_row = row;
            }
        }
        if max_val < 1e-15 {
            return Err(CyaneaError::InvalidInput(
                "singular information matrix — covariates may be collinear".into(),
            ));
        }
        if max_row != col {
            for j in 0..cols {
                aug.swap(col * cols + j, max_row * cols + j);
            }
        }
        // Scale pivot row
        let pivot = aug[col * cols + col];
        for j in 0..cols {
            aug[col * cols + j] /= pivot;
        }
        // Eliminate all other rows
        for row in 0..n {
            if row == col {
                continue;
            }
            let factor = aug[row * cols + col];
            for j in 0..cols {
                let above = aug[col * cols + j];
                aug[row * cols + j] -= factor * above;
            }
        }
    }

    // Extract inverse
    let mut inv = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            inv[i * n + j] = aug[i * cols + n + j];
        }
    }
    Ok(inv)
}

// ── Kaplan-Meier ────────────────────────────────────────────────────────────

/// Compute the Kaplan-Meier survival curve.
///
/// Returns survival estimates at each unique event time with Greenwood standard
/// errors and 95% log-transform confidence intervals.
///
/// # Arguments
///
/// * `times` — observed times (non-negative)
/// * `status` — event indicators (`true` = event, `false` = censored)
///
/// # Errors
///
/// Returns an error if inputs are empty, have mismatched lengths, or contain
/// negative times.
///
/// # Example
///
/// ```
/// use cyanea_stats::survival::kaplan_meier;
///
/// let times  = [1.0, 2.0, 3.0, 4.0, 5.0];
/// let status = [true, true, true, true, true];
/// let km = kaplan_meier(&times, &status).unwrap();
/// assert_eq!(km.steps.len(), 5);
/// assert_eq!(km.n_events, 5);
/// ```
pub fn kaplan_meier(times: &[f64], status: &[bool]) -> Result<KmResult> {
    let order = validate_and_sort(times, status)?;
    let n = times.len();

    // Collect unique event times and their event/censored counts.
    // We process sorted observations, tracking the risk set.
    struct TimeInfo {
        time: f64,
        events: usize,
        censored_before: usize, // censored strictly between previous event time and this one
    }

    let mut time_infos: Vec<TimeInfo> = Vec::new();
    let mut censored_accum = 0usize;

    let mut i = 0;
    while i < n {
        let t = times[order[i]];
        let is_event = status[order[i]];

        if !is_event {
            censored_accum += 1;
            i += 1;
            continue;
        }

        // Count all events (and censored) at this exact time
        let mut events_at_t = 0usize;
        let mut censored_at_t = 0usize;
        let j_start = i;
        while i < n && times[order[i]] == t {
            if status[order[i]] {
                events_at_t += 1;
            } else {
                censored_at_t += 1;
            }
            i += 1;
        }
        // Censored at exact event time count as censored before the next event
        // (convention: events processed before censoring at the same time)
        let _ = j_start;
        time_infos.push(TimeInfo {
            time: t,
            events: events_at_t,
            censored_before: censored_accum,
        });
        censored_accum = censored_at_t;
    }

    let total_events: usize = time_infos.iter().map(|ti| ti.events).sum();

    // Build KM steps
    let mut steps = Vec::with_capacity(time_infos.len());
    let mut n_risk = n;
    let mut survival = 1.0;
    let mut greenwood_sum = 0.0;
    let mut median = None;

    for ti in &time_infos {
        // Remove censored subjects who left before this event time
        n_risk -= ti.censored_before;

        let d = ti.events;
        survival *= 1.0 - (d as f64 / n_risk as f64);

        // Greenwood variance term
        if n_risk > d {
            greenwood_sum += d as f64 / (n_risk as f64 * (n_risk as f64 - d as f64));
        }

        let std_err = survival * greenwood_sum.sqrt();

        // Log-transform CI: transform to θ = log(-log(S)), compute CI, back-transform.
        let (ci_lower, ci_upper) = if survival > 0.0 && survival < 1.0 {
            let log_neg_log_s = (-survival.ln()).ln();
            // SE on log(-log(S)) scale: SE(S) / (S * |log(S)|)
            let se_theta = if greenwood_sum > 0.0 {
                greenwood_sum.sqrt() / survival.ln().abs()
            } else {
                0.0
            };
            let z = 1.96;
            let theta_lo = log_neg_log_s - z * se_theta;
            let theta_hi = log_neg_log_s + z * se_theta;
            let lo = (-theta_hi.exp()).exp(); // note: swapped because of double-log
            let hi = (-theta_lo.exp()).exp();
            (lo.max(0.0), hi.min(1.0))
        } else if survival <= 0.0 {
            (0.0, 0.0)
        } else {
            // survival == 1.0 (shouldn't happen at event times, but be safe)
            (1.0, 1.0)
        };

        let n_censored = ti.censored_before;

        steps.push(KmStep {
            time: ti.time,
            n_risk,
            n_events: d,
            n_censored,
            survival,
            std_err,
            ci_lower,
            ci_upper,
        });

        // Median: first time S(t) ≤ 0.5
        if median.is_none() && survival <= 0.5 {
            median = Some(ti.time);
        }

        // Remove event subjects from risk set
        n_risk -= d;
    }

    Ok(KmResult {
        steps,
        median_survival: median,
        n_total: n,
        n_events: total_events,
    })
}

// ── Log-rank test ───────────────────────────────────────────────────────────

/// Compare survival curves between groups using the log-rank test.
///
/// # Arguments
///
/// * `times` — observed times (non-negative)
/// * `status` — event indicators (`true` = event, `false` = censored)
/// * `groups` — group labels (0-indexed; must contain at least 2 distinct values)
///
/// # Errors
///
/// Returns an error if inputs are empty, have mismatched lengths, contain
/// negative times, or have fewer than 2 groups.
///
/// # Example
///
/// ```
/// use cyanea_stats::survival::log_rank_test;
///
/// let times  = [1.0, 2.0, 3.0, 4.0, 1.5, 2.5, 3.5, 4.5];
/// let status = [true, true, true, true, true, true, true, true];
/// let groups = [0, 0, 0, 0, 1, 1, 1, 1];
/// let lr = log_rank_test(&times, &status, &groups).unwrap();
/// assert_eq!(lr.degrees_of_freedom, 1);
/// ```
pub fn log_rank_test(
    times: &[f64],
    status: &[bool],
    groups: &[usize],
) -> Result<LogRankResult> {
    let order = validate_and_sort(times, status)?;
    let n = times.len();

    if groups.len() != n {
        return Err(CyaneaError::InvalidInput(format!(
            "times length ({}) != groups length ({})",
            n,
            groups.len()
        )));
    }

    let n_groups = *groups.iter().max().unwrap() + 1;
    if n_groups < 2 {
        return Err(CyaneaError::InvalidInput(
            "log-rank test requires at least 2 groups".into(),
        ));
    }

    // Check all group labels are within range
    for (i, &g) in groups.iter().enumerate() {
        if g >= n_groups {
            return Err(CyaneaError::InvalidInput(format!(
                "group label {} at index {} exceeds max group {}",
                g,
                i,
                n_groups - 1
            )));
        }
    }

    // Per-group risk counts
    let mut risk_per_group = vec![0usize; n_groups];
    for &idx in &order {
        risk_per_group[groups[idx]] += 1;
    }

    let mut observed = vec![0.0; n_groups];
    let mut expected = vec![0.0; n_groups];

    let mut i = 0;
    while i < n {
        let t = times[order[i]];

        // Gather all observations at this time
        let mut d_total = 0usize;
        let mut d_per_group = vec![0usize; n_groups];
        let mut c_per_group = vec![0usize; n_groups];
        let j_start = i;
        while i < n && times[order[i]] == t {
            let g = groups[order[i]];
            if status[order[i]] {
                d_total += 1;
                d_per_group[g] += 1;
            } else {
                c_per_group[g] += 1;
            }
            i += 1;
        }

        let n_risk_total: usize = risk_per_group.iter().sum();
        if n_risk_total > 0 && d_total > 0 {
            for g in 0..n_groups {
                observed[g] += d_per_group[g] as f64;
                expected[g] +=
                    risk_per_group[g] as f64 * d_total as f64 / n_risk_total as f64;
            }
        }

        // Remove subjects at this time from risk sets
        let _ = j_start;
        for g in 0..n_groups {
            risk_per_group[g] -= d_per_group[g] + c_per_group[g];
        }
    }

    // Chi-squared statistic
    let mut chi2 = 0.0;
    for g in 0..n_groups {
        if expected[g] > 0.0 {
            chi2 += (observed[g] - expected[g]).powi(2) / expected[g];
        }
    }

    let df = n_groups - 1;
    let chi_dist = ChiSquared::new(df as f64)?;
    let p_value = 1.0 - chi_dist.cdf(chi2);

    Ok(LogRankResult {
        statistic: chi2,
        degrees_of_freedom: df,
        p_value,
        n_groups,
        observed,
        expected,
    })
}

// ── Cox proportional hazards ────────────────────────────────────────────────

/// Fit a Cox proportional hazards model via Breslow partial likelihood.
///
/// Uses Newton-Raphson optimization (max 100 iterations, tolerance 1e-9).
///
/// # Arguments
///
/// * `times` — observed times (non-negative)
/// * `status` — event indicators (`true` = event, `false` = censored)
/// * `covariates` — covariate matrix in row-major order (n × p), where n is
///   the number of observations and p is `n_covariates`
/// * `n_covariates` — number of covariates (columns)
///
/// # Errors
///
/// Returns an error if inputs are empty, have mismatched lengths, `n_covariates`
/// is zero, or the covariate matrix has wrong dimensions.
///
/// # Example
///
/// ```
/// use cyanea_stats::survival::cox_ph;
///
/// let times  = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
/// let status = [true, true, false, true, false, true];
/// let covs   = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]; // single binary covariate
/// let result = cox_ph(&times, &status, &covs, 1).unwrap();
/// assert_eq!(result.n_covariates, 1);
/// assert!(result.converged);
/// ```
pub fn cox_ph(
    times: &[f64],
    status: &[bool],
    covariates: &[f64],
    n_covariates: usize,
) -> Result<CoxPhResult> {
    let order = validate_and_sort(times, status)?;
    let n = times.len();
    let p = n_covariates;

    if p == 0 {
        return Err(CyaneaError::InvalidInput(
            "n_covariates must be at least 1".into(),
        ));
    }
    if covariates.len() != n * p {
        return Err(CyaneaError::InvalidInput(format!(
            "covariates length ({}) != n ({}) * p ({})",
            covariates.len(),
            n,
            p
        )));
    }

    let total_events = status.iter().filter(|&&s| s).count();
    if total_events == 0 {
        return Err(CyaneaError::InvalidInput(
            "no events in the data".into(),
        ));
    }

    // Accessor for covariate j of subject i
    let x = |i: usize, j: usize| covariates[i * p + j];

    // Precompute unique event times and the subjects who have events at each time
    let mut event_times: Vec<f64> = Vec::new();
    let mut event_at_time: Vec<Vec<usize>> = Vec::new();
    // For each unique time (event or not), record sorted-order index ranges
    // so we can efficiently build risk sets via backward scan.
    {
        let mut i = 0;
        while i < n {
            let t = times[order[i]];
            let mut events_here = Vec::new();
            while i < n && times[order[i]] == t {
                if status[order[i]] {
                    events_here.push(order[i]);
                }
                i += 1;
            }
            if !events_here.is_empty() {
                event_times.push(t);
                event_at_time.push(events_here);
            }
        }
    }

    // Compute Breslow partial log-likelihood, score, and information matrix
    // for a given beta.
    let compute_ll = |beta: &[f64]| -> (f64, Vec<f64>, Vec<f64>) {
        let mut score = vec![0.0; p];
        let mut info = vec![0.0; p * p];
        let mut ll = 0.0;
        let mut s0 = 0.0;
        let mut s1 = vec![0.0; p];
        let mut s2 = vec![0.0; p * p];
        let mut right = n;

        for ti in (0..event_times.len()).rev() {
            let t = event_times[ti];

            // Add all subjects with time >= t to the risk set
            while right > 0 && times[order[right - 1]] >= t {
                right -= 1;
                let subj = order[right];
                let eta: f64 = (0..p).map(|j| beta[j] * x(subj, j)).sum();
                let w = eta.exp();
                s0 += w;
                for j in 0..p {
                    s1[j] += x(subj, j) * w;
                    for k in 0..p {
                        s2[j * p + k] += x(subj, j) * x(subj, k) * w;
                    }
                }
            }

            let d = event_at_time[ti].len() as f64;
            if s0 <= 0.0 {
                continue;
            }

            for &subj in &event_at_time[ti] {
                let eta: f64 = (0..p).map(|j| beta[j] * x(subj, j)).sum();
                ll += eta;
                for j in 0..p {
                    score[j] += x(subj, j);
                }
            }

            ll -= d * s0.ln();
            for j in 0..p {
                score[j] -= d * s1[j] / s0;
                for k in 0..p {
                    info[j * p + k] +=
                        d * (s2[j * p + k] / s0 - (s1[j] * s1[k]) / (s0 * s0));
                }
            }
        }

        (ll, score, info)
    };

    // Newton-Raphson with step-halving
    let max_iter = 100;
    let tol = 1e-9;
    let mut beta = vec![0.0; p];
    let mut converged = false;
    let mut n_iter = 0;

    let (mut ll, mut score, mut info_matrix) = compute_ll(&beta);

    for iter in 0..max_iter {
        n_iter = iter + 1;

        // Newton-Raphson direction: delta = I^{-1} * U
        let delta = match solve_linear_system(&info_matrix, &score, p) {
            Ok(d) => d,
            Err(_) => {
                // Near-singular information matrix — stop iteration
                break;
            }
        };

        let mut max_delta = 0.0_f64;
        for j in 0..p {
            max_delta = max_delta.max(delta[j].abs());
        }

        if max_delta < tol {
            converged = true;
            break;
        }

        // Step-halving: ensure log-likelihood improves
        let mut step = 1.0;
        let mut new_beta = vec![0.0; p];
        let mut new_ll;
        let mut new_score;
        let mut new_info;
        loop {
            for j in 0..p {
                new_beta[j] = beta[j] + step * delta[j];
            }
            let result = compute_ll(&new_beta);
            new_ll = result.0;
            new_score = result.1;
            new_info = result.2;
            if new_ll >= ll - 1e-10 || step < 1e-4 {
                break;
            }
            step *= 0.5;
        }

        beta = new_beta;
        ll = new_ll;
        score = new_score;
        info_matrix = new_info;
    }

    // Standard errors from inverse of information matrix
    let inv_info = invert_matrix(&info_matrix, p)?;
    let normal = Normal::standard();

    let mut std_errors = Vec::with_capacity(p);
    let mut z_values = Vec::with_capacity(p);
    let mut p_values = Vec::with_capacity(p);
    let mut hazard_ratios = Vec::with_capacity(p);
    let mut hr_ci_lower = Vec::with_capacity(p);
    let mut hr_ci_upper = Vec::with_capacity(p);

    for j in 0..p {
        let var = inv_info[j * p + j];
        let se = if var > 0.0 { var.sqrt() } else { 0.0 };
        let z = if se > 0.0 { beta[j] / se } else { 0.0 };
        let pval = 2.0 * (1.0 - normal.cdf(z.abs()));
        let hr = beta[j].exp();

        std_errors.push(se);
        z_values.push(z);
        p_values.push(pval);
        hazard_ratios.push(hr);
        hr_ci_lower.push((beta[j] - 1.96 * se).exp());
        hr_ci_upper.push((beta[j] + 1.96 * se).exp());
    }

    Ok(CoxPhResult {
        coefficients: beta,
        std_errors,
        z_values,
        p_values,
        hazard_ratios,
        hr_ci_lower,
        hr_ci_upper,
        n_observations: n,
        n_events: total_events,
        n_covariates: p,
        log_likelihood: ll,
        n_iterations: n_iter,
        converged,
    })
}

// ── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-6;

    // ── Kaplan-Meier ────────────────────────────────────────────────────

    #[test]
    fn km_no_censoring() {
        // 5 subjects, all events at t = 1..5
        let times = [1.0, 2.0, 3.0, 4.0, 5.0];
        let status = [true, true, true, true, true];
        let km = kaplan_meier(&times, &status).unwrap();

        assert_eq!(km.steps.len(), 5);
        assert_eq!(km.n_total, 5);
        assert_eq!(km.n_events, 5);

        // S(1) = 4/5 = 0.8
        assert!((km.steps[0].survival - 0.8).abs() < TOL);
        // S(2) = 0.8 * 3/4 = 0.6
        assert!((km.steps[1].survival - 0.6).abs() < TOL);
        // S(3) = 0.6 * 2/3 = 0.4
        assert!((km.steps[2].survival - 0.4).abs() < TOL);
        // S(4) = 0.4 * 1/2 = 0.2
        assert!((km.steps[3].survival - 0.2).abs() < TOL);
        // S(5) = 0.0
        assert!((km.steps[4].survival - 0.0).abs() < TOL);
    }

    #[test]
    fn km_with_censoring() {
        // Classic textbook example
        let times = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let status = [true, false, true, false, true, false];
        let km = kaplan_meier(&times, &status).unwrap();

        assert_eq!(km.n_events, 3);
        assert_eq!(km.steps.len(), 3); // only event times

        // At t=1: n_risk=6, d=1, S = 5/6
        assert_eq!(km.steps[0].n_risk, 6);
        assert_eq!(km.steps[0].n_events, 1);
        assert!((km.steps[0].survival - 5.0 / 6.0).abs() < TOL);

        // At t=3: n_risk=4 (lost t=2 censored), d=1, S = 5/6 * 3/4 = 5/8
        assert_eq!(km.steps[1].n_risk, 4);
        assert!((km.steps[1].survival - 5.0 / 8.0).abs() < TOL);

        // At t=5: n_risk=2 (lost t=4 censored), d=1, S = 5/8 * 1/2 = 5/16
        assert_eq!(km.steps[2].n_risk, 2);
        assert!((km.steps[2].survival - 5.0 / 16.0).abs() < TOL);
    }

    #[test]
    fn km_all_censored() {
        let times = [1.0, 2.0, 3.0];
        let status = [false, false, false];
        let km = kaplan_meier(&times, &status).unwrap();

        assert_eq!(km.n_events, 0);
        assert!(km.steps.is_empty());
        assert!(km.median_survival.is_none());
    }

    #[test]
    fn km_single_event() {
        let times = [5.0];
        let status = [true];
        let km = kaplan_meier(&times, &status).unwrap();

        assert_eq!(km.steps.len(), 1);
        assert_eq!(km.steps[0].n_risk, 1);
        assert!((km.steps[0].survival - 0.0).abs() < TOL);
    }

    #[test]
    fn km_single_censored() {
        let times = [5.0];
        let status = [false];
        let km = kaplan_meier(&times, &status).unwrap();

        assert!(km.steps.is_empty());
        assert_eq!(km.n_events, 0);
    }

    #[test]
    fn km_tied_times() {
        // Two events and one censored at the same time
        let times = [3.0, 3.0, 3.0, 5.0];
        let status = [true, true, false, true];
        let km = kaplan_meier(&times, &status).unwrap();

        assert_eq!(km.steps.len(), 2);
        // At t=3: n_risk=4, d=2, S = 2/4 = 0.5
        assert_eq!(km.steps[0].n_risk, 4);
        assert_eq!(km.steps[0].n_events, 2);
        assert!((km.steps[0].survival - 0.5).abs() < TOL);
    }

    #[test]
    fn km_median_survival() {
        // Median at t=3 where S first ≤ 0.5
        let times = [1.0, 2.0, 3.0, 4.0, 5.0];
        let status = [true, true, true, true, true];
        let km = kaplan_meier(&times, &status).unwrap();

        // S(3) = 0.4 ≤ 0.5, so median = 3.0
        assert_eq!(km.median_survival, Some(3.0));
    }

    #[test]
    fn km_median_not_reached() {
        // Only one event out of 5, S never reaches 0.5
        let times = [1.0, 2.0, 3.0, 4.0, 5.0];
        let status = [true, false, false, false, false];
        let km = kaplan_meier(&times, &status).unwrap();

        // S(1) = 4/5 = 0.8 > 0.5, no further event times
        assert!(km.median_survival.is_none());
    }

    #[test]
    fn km_ci_bounds() {
        let times = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let status = [true, false, true, false, true, false, true, false];
        let km = kaplan_meier(&times, &status).unwrap();

        for step in &km.steps {
            assert!(step.ci_lower >= 0.0, "CI lower < 0: {}", step.ci_lower);
            assert!(step.ci_upper <= 1.0, "CI upper > 1: {}", step.ci_upper);
            assert!(
                step.ci_lower <= step.survival + TOL,
                "CI lower > survival"
            );
            assert!(
                step.ci_upper >= step.survival - TOL,
                "CI upper < survival"
            );
        }
    }

    #[test]
    fn km_se_non_negative() {
        let times = [1.0, 2.0, 3.0, 4.0, 5.0];
        let status = [true, true, false, true, true];
        let km = kaplan_meier(&times, &status).unwrap();

        for step in &km.steps {
            assert!(step.std_err >= 0.0, "negative SE: {}", step.std_err);
        }
    }

    #[test]
    fn km_monotone_decreasing() {
        let times = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
        let status = [true, true, true, false, true, false, true, true, false, true];
        let km = kaplan_meier(&times, &status).unwrap();

        for i in 1..km.steps.len() {
            assert!(
                km.steps[i].survival <= km.steps[i - 1].survival + TOL,
                "survival increased at step {}",
                i
            );
        }
    }

    // ── Log-rank test ───────────────────────────────────────────────────

    #[test]
    fn logrank_significant_difference() {
        // Group 0: early events; Group 1: late events
        let times = [1.0, 2.0, 3.0, 4.0, 10.0, 11.0, 12.0, 13.0];
        let status = [true, true, true, true, true, true, true, true];
        let groups = [0, 0, 0, 0, 1, 1, 1, 1];

        let lr = log_rank_test(&times, &status, &groups).unwrap();
        assert_eq!(lr.degrees_of_freedom, 1);
        assert!(lr.p_value < 0.05, "expected significant, got p={}", lr.p_value);
    }

    #[test]
    fn logrank_no_difference() {
        // Interleaved times, same pattern
        let times = [1.0, 2.0, 3.0, 4.0, 1.5, 2.5, 3.5, 4.5];
        let status = [true, true, true, true, true, true, true, true];
        let groups = [0, 0, 0, 0, 1, 1, 1, 1];

        let lr = log_rank_test(&times, &status, &groups).unwrap();
        assert!(
            lr.p_value > 0.05,
            "expected non-significant, got p={}",
            lr.p_value
        );
    }

    #[test]
    fn logrank_tied_times() {
        let times = [1.0, 1.0, 2.0, 2.0, 3.0, 3.0];
        let status = [true, true, true, true, true, true];
        let groups = [0, 1, 0, 1, 0, 1];

        let lr = log_rank_test(&times, &status, &groups).unwrap();
        // Symmetric groups should have equal observed and expected
        assert!((lr.observed[0] - lr.expected[0]).abs() < TOL);
        assert!((lr.observed[1] - lr.expected[1]).abs() < TOL);
    }

    #[test]
    fn logrank_three_groups() {
        let times = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let status = [true, true, true, true, true, true, true, true, true];
        let groups = [0, 0, 0, 1, 1, 1, 2, 2, 2];

        let lr = log_rank_test(&times, &status, &groups).unwrap();
        assert_eq!(lr.n_groups, 3);
        assert_eq!(lr.degrees_of_freedom, 2);
    }

    #[test]
    fn logrank_observed_expected_sum() {
        let times = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let status = [true, true, true, true, true, true];
        let groups = [0, 0, 0, 1, 1, 1];

        let lr = log_rank_test(&times, &status, &groups).unwrap();

        let obs_sum: f64 = lr.observed.iter().sum();
        let exp_sum: f64 = lr.expected.iter().sum();
        assert!(
            (obs_sum - exp_sum).abs() < TOL,
            "O sum ({}) != E sum ({})",
            obs_sum,
            exp_sum
        );
    }

    // ── Cox PH ──────────────────────────────────────────────────────────

    #[test]
    fn cox_single_binary_covariate() {
        // Group 0 (x=0) tends to have earlier events, group 1 (x=1) later,
        // but with overlap to avoid perfect separation.
        let times = [1.0, 2.0, 4.0, 5.0, 3.0, 6.0, 7.0, 8.0];
        let status = [true, true, true, false, true, true, true, false];
        let covs = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];

        let result = cox_ph(&times, &status, &covs, 1).unwrap();
        assert!(result.converged);
        assert_eq!(result.n_covariates, 1);
        // x=1 has later events → lower hazard → beta < 0
        assert!(
            result.coefficients[0] < 0.0,
            "expected negative beta, got {}",
            result.coefficients[0]
        );
    }

    #[test]
    fn cox_known_hr_recovery() {
        // Overlapping times: x=0 tends early, x=1 tends late, but ranges overlap
        let times = [1.0, 2.0, 4.0, 6.0, 8.0, 3.0, 5.0, 7.0, 9.0, 10.0];
        let status = [true, true, true, true, false, true, true, true, true, true];
        let covs = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0];

        let result = cox_ph(&times, &status, &covs, 1).unwrap();
        assert!(result.converged);
        // HR should be < 1 (x=1 subjects survive longer on average)
        assert!(
            result.hazard_ratios[0] < 1.0,
            "expected HR < 1, got {}",
            result.hazard_ratios[0]
        );
    }

    #[test]
    fn cox_multiple_covariates() {
        let times = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let status = [true, true, true, true, true, true, true, true];
        // Two covariates
        #[rustfmt::skip]
        let covs = [
            0.0, 1.0,
            0.0, 2.0,
            1.0, 1.0,
            1.0, 3.0,
            0.0, 4.0,
            1.0, 2.0,
            0.0, 5.0,
            1.0, 4.0,
        ];

        let result = cox_ph(&times, &status, &covs, 2).unwrap();
        assert!(result.converged);
        assert_eq!(result.coefficients.len(), 2);
        assert_eq!(result.std_errors.len(), 2);
        assert_eq!(result.hazard_ratios.len(), 2);
    }

    #[test]
    fn cox_hr_ci_contains_point() {
        // Overlapping groups to avoid separation
        let times = [1.0, 3.0, 5.0, 2.0, 4.0, 6.0, 7.0, 8.0];
        let status = [true, true, true, true, true, true, true, true];
        let covs = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0];

        let result = cox_ph(&times, &status, &covs, 1).unwrap();

        assert!(
            result.hr_ci_lower[0] <= result.hazard_ratios[0],
            "HR CI lower ({}) > HR ({})",
            result.hr_ci_lower[0],
            result.hazard_ratios[0]
        );
        assert!(
            result.hr_ci_upper[0] >= result.hazard_ratios[0],
            "HR CI upper ({}) < HR ({})",
            result.hr_ci_upper[0],
            result.hazard_ratios[0]
        );
    }

    #[test]
    fn cox_null_effect() {
        // Random-ish covariate unrelated to survival
        let times = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let status = [true, true, true, true, true, true, true, true];
        let covs = [1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0]; // alternating, uncorrelated

        let result = cox_ph(&times, &status, &covs, 1).unwrap();
        assert!(result.converged);
        // With no real effect, p-value should be large
        assert!(
            result.p_values[0] > 0.05,
            "expected non-significant, got p={}",
            result.p_values[0]
        );
    }

    // ── Error cases ─────────────────────────────────────────────────────

    #[test]
    fn error_empty_input() {
        assert!(kaplan_meier(&[], &[]).is_err());
        assert!(log_rank_test(&[], &[], &[]).is_err());
        assert!(cox_ph(&[], &[], &[], 1).is_err());
    }

    #[test]
    fn error_mismatched_lengths() {
        assert!(kaplan_meier(&[1.0, 2.0], &[true]).is_err());
        assert!(log_rank_test(&[1.0], &[true], &[0, 1]).is_err());
        assert!(cox_ph(&[1.0, 2.0], &[true, true], &[0.0], 1).is_err());
    }

    #[test]
    fn error_negative_time() {
        assert!(kaplan_meier(&[-1.0, 2.0], &[true, true]).is_err());
    }

    #[test]
    fn error_single_group_logrank() {
        let times = [1.0, 2.0, 3.0];
        let status = [true, true, true];
        let groups = [0, 0, 0];
        assert!(log_rank_test(&times, &status, &groups).is_err());
    }

    #[test]
    fn error_zero_covariates_cox() {
        assert!(cox_ph(&[1.0], &[true], &[], 0).is_err());
    }
}
