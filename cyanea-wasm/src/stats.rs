//! Statistical wrappers with JSON input/output.
//!
//! Types like [`DescriptiveStats`] and [`TestResult`] do not derive `Serialize`,
//! so we provide thin wrappers (`Js*` types) that mirror their fields and add
//! `#[derive(Serialize)]`.

use serde::Serialize;

use cyanea_stats::correlation;
use cyanea_stats::correction;
use cyanea_stats::descriptive;
use cyanea_stats::diversity;
use cyanea_stats::null_model;
use cyanea_stats::popgen;
use cyanea_stats::survival;
use cyanea_stats::testing;

use crate::error::{wasm_err, wasm_ok};

#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::*;

// ── Wrapper types ────────────────────────────────────────────────────────

/// Serializable mirror of [`cyanea_stats::DescriptiveStats`].
#[derive(Debug, Serialize)]
pub struct JsDescriptiveStats {
    pub count: usize,
    pub mean: f64,
    pub median: f64,
    pub variance: f64,
    pub sample_variance: f64,
    pub std_dev: f64,
    pub sample_std_dev: f64,
    pub min: f64,
    pub max: f64,
    pub range: f64,
    pub q1: f64,
    pub q3: f64,
    pub iqr: f64,
    pub skewness: f64,
    pub kurtosis: f64,
}

impl From<descriptive::DescriptiveStats> for JsDescriptiveStats {
    fn from(s: descriptive::DescriptiveStats) -> Self {
        Self {
            count: s.count,
            mean: s.mean,
            median: s.median,
            variance: s.variance,
            sample_variance: s.sample_variance,
            std_dev: s.std_dev,
            sample_std_dev: s.sample_std_dev,
            min: s.min,
            max: s.max,
            range: s.range,
            q1: s.q1,
            q3: s.q3,
            iqr: s.iqr,
            skewness: s.skewness,
            kurtosis: s.kurtosis,
        }
    }
}

/// Serializable mirror of [`cyanea_stats::TestResult`].
#[derive(Debug, Serialize)]
pub struct JsTestResult {
    pub statistic: f64,
    pub p_value: f64,
    pub degrees_of_freedom: Option<f64>,
    pub method: String,
}

impl From<testing::TestResult> for JsTestResult {
    fn from(r: testing::TestResult) -> Self {
        Self {
            statistic: r.statistic,
            p_value: r.p_value,
            degrees_of_freedom: r.degrees_of_freedom,
            method: r.method,
        }
    }
}

/// Serializable mirror of [`cyanea_stats::survival::KmStep`].
#[derive(Debug, Serialize)]
pub struct JsKmStep {
    pub time: f64,
    pub n_risk: usize,
    pub n_events: usize,
    pub survival: f64,
    pub std_err: f64,
    pub ci_lower: f64,
    pub ci_upper: f64,
}

/// Serializable mirror of [`cyanea_stats::survival::KmResult`].
#[derive(Debug, Serialize)]
pub struct JsKmResult {
    pub steps: Vec<JsKmStep>,
    pub median_survival: Option<f64>,
    pub n_total: usize,
    pub n_events: usize,
}

/// Serializable mirror of [`cyanea_stats::survival::LogRankResult`].
#[derive(Debug, Serialize)]
pub struct JsLogRankResult {
    pub statistic: f64,
    pub p_value: f64,
    pub n_groups: usize,
}

/// Serializable mirror of [`cyanea_stats::survival::CoxPhResult`].
#[derive(Debug, Serialize)]
pub struct JsCoxPhResult {
    pub coefficients: Vec<f64>,
    pub std_errors: Vec<f64>,
    pub hazard_ratios: Vec<f64>,
    pub p_values: Vec<f64>,
    pub n_observations: usize,
    pub n_events: usize,
    pub converged: bool,
}

/// Serializable mirror of [`cyanea_stats::null_model::WrightFisherResult`].
#[derive(Debug, Serialize)]
pub struct JsWrightFisherResult {
    pub frequencies: Vec<f64>,
    pub fixation_gen: Option<usize>,
    pub final_freq: f64,
}

/// Serializable Tajima's D result.
#[derive(Debug, Serialize)]
pub struct JsTajimaD {
    pub d: f64,
    pub pi: f64,
    pub theta_w: f64,
    pub segregating_sites: usize,
    pub n_sequences: usize,
}

/// Serializable Fst result.
#[derive(Debug, Serialize)]
pub struct JsFstResult {
    pub fst: f64,
}

// ── JSON boundary functions ──────────────────────────────────────────────

fn parse_f64_array(json: &str) -> Result<Vec<f64>, String> {
    serde_json::from_str::<Vec<f64>>(json).map_err(|e| format!("invalid JSON array: {e}"))
}

fn parse_usize_array(json: &str) -> Result<Vec<usize>, String> {
    serde_json::from_str::<Vec<usize>>(json).map_err(|e| format!("invalid JSON array: {e}"))
}

fn parse_bool_array(json: &str) -> Result<Vec<bool>, String> {
    serde_json::from_str::<Vec<bool>>(json).map_err(|e| format!("invalid JSON array: {e}"))
}

/// Compute descriptive statistics from a JSON array of numbers.
///
/// Input: `"[1.0, 2.0, 3.0]"` — Output: JSON `JsDescriptiveStats`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn describe(data_json: &str) -> String {
    let data = match parse_f64_array(data_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    match descriptive::describe(&data) {
        Ok(stats) => wasm_ok(&JsDescriptiveStats::from(stats)),
        Err(e) => wasm_err(e),
    }
}

/// Pearson correlation between two JSON arrays.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn pearson(x_json: &str, y_json: &str) -> String {
    let x = match parse_f64_array(x_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let y = match parse_f64_array(y_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    match correlation::pearson(&x, &y) {
        Ok(r) => wasm_ok(&r),
        Err(e) => wasm_err(e),
    }
}

/// One-sample t-test on a JSON array against hypothesised mean `mu`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn t_test(data_json: &str, mu: f64) -> String {
    let data = match parse_f64_array(data_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    match testing::t_test_one_sample(&data, mu) {
        Ok(r) => wasm_ok(&JsTestResult::from(r)),
        Err(e) => wasm_err(e),
    }
}

/// Spearman rank correlation between two JSON arrays.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn spearman(x_json: &str, y_json: &str) -> String {
    let x = match parse_f64_array(x_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let y = match parse_f64_array(y_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    match correlation::spearman(&x, &y) {
        Ok(r) => wasm_ok(&r),
        Err(e) => wasm_err(e),
    }
}

/// Two-sample t-test (Student's or Welch's) on two JSON arrays.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn t_test_two_sample(x_json: &str, y_json: &str, equal_var: bool) -> String {
    let x = match parse_f64_array(x_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let y = match parse_f64_array(y_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    match testing::t_test_two_sample(&x, &y, equal_var) {
        Ok(r) => wasm_ok(&JsTestResult::from(r)),
        Err(e) => wasm_err(e),
    }
}

/// Mann-Whitney U test (non-parametric) on two JSON arrays.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn mann_whitney_u(x_json: &str, y_json: &str) -> String {
    let x = match parse_f64_array(x_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let y = match parse_f64_array(y_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    match testing::mann_whitney_u(&x, &y) {
        Ok(r) => wasm_ok(&JsTestResult::from(r)),
        Err(e) => wasm_err(e),
    }
}

/// Bonferroni p-value correction on a JSON array of p-values.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn bonferroni(p_json: &str) -> String {
    let p = match parse_f64_array(p_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    match correction::bonferroni(&p) {
        Ok(corrected) => wasm_ok(&corrected),
        Err(e) => wasm_err(e),
    }
}

/// Benjamini-Hochberg FDR correction on a JSON array of p-values.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn benjamini_hochberg(p_json: &str) -> String {
    let p = match parse_f64_array(p_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    match correction::benjamini_hochberg(&p) {
        Ok(corrected) => wasm_ok(&corrected),
        Err(e) => wasm_err(e),
    }
}

/// Kaplan-Meier survival curve from JSON arrays of times and event status.
///
/// Input: times as `"[1.0, 2.0, 3.0]"`, status as `"[true, false, true]"`.
/// Output: JSON `JsKmResult`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn kaplan_meier(times_json: &str, status_json: &str) -> String {
    let times = match parse_f64_array(times_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let status = match parse_bool_array(status_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    match survival::kaplan_meier(&times, &status) {
        Ok(km) => {
            let js = JsKmResult {
                steps: km
                    .steps
                    .iter()
                    .map(|s| JsKmStep {
                        time: s.time,
                        n_risk: s.n_risk,
                        n_events: s.n_events,
                        survival: s.survival,
                        std_err: s.std_err,
                        ci_lower: s.ci_lower,
                        ci_upper: s.ci_upper,
                    })
                    .collect(),
                median_survival: km.median_survival,
                n_total: km.n_total,
                n_events: km.n_events,
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Log-rank test comparing survival between two groups.
///
/// Input: times and status for each group as JSON arrays.
/// Output: JSON `JsLogRankResult`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn log_rank_test(
    times1_json: &str,
    status1_json: &str,
    times2_json: &str,
    status2_json: &str,
) -> String {
    let times1 = match parse_f64_array(times1_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let status1 = match parse_bool_array(status1_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let times2 = match parse_f64_array(times2_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let status2 = match parse_bool_array(status2_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };

    let n1 = times1.len();
    let n2 = times2.len();
    let mut times = times1;
    times.extend_from_slice(&times2);
    let mut status = status1;
    status.extend_from_slice(&status2);
    let mut groups: Vec<usize> = vec![0; n1];
    groups.extend(vec![1; n2]);

    match survival::log_rank_test(&times, &status, &groups) {
        Ok(lr) => {
            let js = JsLogRankResult {
                statistic: lr.statistic,
                p_value: lr.p_value,
                n_groups: lr.n_groups,
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Cox proportional hazards regression.
///
/// Input: times, status, and flattened covariate matrix as JSON arrays.
/// Output: JSON `JsCoxPhResult`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn cox_ph(
    times_json: &str,
    status_json: &str,
    covariates_json: &str,
    n_covariates: usize,
) -> String {
    let times = match parse_f64_array(times_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let status = match parse_bool_array(status_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let covariates = match parse_f64_array(covariates_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    match survival::cox_ph(&times, &status, &covariates, n_covariates) {
        Ok(result) => {
            let js = JsCoxPhResult {
                coefficients: result.coefficients,
                std_errors: result.std_errors,
                hazard_ratios: result.hazard_ratios,
                p_values: result.p_values,
                n_observations: result.n_observations,
                n_events: result.n_events,
                converged: result.converged,
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Wright-Fisher allele frequency drift simulation.
///
/// Output: JSON `JsWrightFisherResult`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn wright_fisher(
    pop_size: usize,
    initial_freq: f64,
    n_generations: usize,
    seed: u64,
) -> String {
    match null_model::wright_fisher(pop_size, initial_freq, n_generations, seed) {
        Ok(result) => {
            let js = JsWrightFisherResult {
                frequencies: result.frequencies,
                fixation_gen: result.fixation_gen,
                final_freq: result.final_freq,
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Permutation test: generates a null distribution of mean differences
/// between groups.
///
/// Input: pooled values and group sizes as JSON arrays.
/// Output: JSON array of permutation statistic values.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn permutation_test(
    values_json: &str,
    group_sizes_json: &str,
    n_permutations: usize,
    seed: u64,
) -> String {
    let values = match parse_f64_array(values_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let group_sizes = match parse_usize_array(group_sizes_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let diff_means = |groups: &[&[f64]]| {
        if groups.len() < 2 || groups[0].is_empty() || groups[1].is_empty() {
            return 0.0;
        }
        let m0: f64 = groups[0].iter().sum::<f64>() / groups[0].len() as f64;
        let m1: f64 = groups[1].iter().sum::<f64>() / groups[1].len() as f64;
        (m0 - m1).abs()
    };
    match null_model::permutation_null(&values, &group_sizes, &diff_means, n_permutations, seed) {
        Ok(dist) => wasm_ok(&dist),
        Err(e) => wasm_err(e),
    }
}

/// Bootstrap confidence interval: generates a null distribution of means
/// via bootstrap resampling.
///
/// Input: data as a JSON array.
/// Output: JSON array of bootstrap statistic values.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn bootstrap_ci(data_json: &str, n_bootstrap: usize, seed: u64) -> String {
    let data = match parse_f64_array(data_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let mean_fn = |d: &[f64]| {
        if d.is_empty() {
            return 0.0;
        }
        d.iter().sum::<f64>() / d.len() as f64
    };
    let dist = null_model::bootstrap_null(&data, &mean_fn, n_bootstrap, seed);
    wasm_ok(&dist)
}

/// Shannon diversity index from a JSON array of species counts.
///
/// Output: JSON f64 value.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn shannon_index(counts_json: &str) -> String {
    let counts = match parse_usize_array(counts_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    match diversity::shannon_index(&counts) {
        Ok(h) => wasm_ok(&h),
        Err(e) => wasm_err(e),
    }
}

/// Simpson diversity index from a JSON array of species counts.
///
/// Output: JSON f64 value.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn simpson_index(counts_json: &str) -> String {
    let counts = match parse_usize_array(counts_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    match diversity::simpson_index(&counts) {
        Ok(d) => wasm_ok(&d),
        Err(e) => wasm_err(e),
    }
}

/// Bray-Curtis dissimilarity between two samples (JSON arrays of counts).
///
/// Output: JSON f64 value (0 = identical, 1 = completely different).
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn bray_curtis(a_json: &str, b_json: &str) -> String {
    let a = match parse_usize_array(a_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let b = match parse_usize_array(b_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    match diversity::bray_curtis(&a, &b) {
        Ok(bc) => wasm_ok(&bc),
        Err(e) => wasm_err(e),
    }
}

/// Hudson's Fst estimator between two populations.
///
/// Input: two JSON arrays of genotype matrices, each `Vec<Vec<Option<u8>>>`.
/// Each inner Vec is the genotype vector for one locus.
/// Output: JSON `JsFstResult`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn fst_hudson(pop1_json: &str, pop2_json: &str) -> String {
    let pop1: Vec<Vec<Option<u8>>> = match serde_json::from_str(pop1_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(format!("invalid JSON for pop1: {e}")),
    };
    let pop2: Vec<Vec<Option<u8>>> = match serde_json::from_str(pop2_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(format!("invalid JSON for pop2: {e}")),
    };
    let pop1_refs: Vec<&[Option<u8>]> = pop1.iter().map(|v| v.as_slice()).collect();
    let pop2_refs: Vec<&[Option<u8>]> = pop2.iter().map(|v| v.as_slice()).collect();
    match popgen::fst_hudson(&pop1_refs, &pop2_refs) {
        Ok(result) => {
            let js = JsFstResult { fst: result.fst };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Tajima's D statistic from pre-computed summary values.
///
/// Output: JSON `JsTajimaD`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn tajimas_d(
    segregating_sites: usize,
    n_sequences: usize,
    avg_pairwise_diff: f64,
) -> String {
    match popgen::tajimas_d(segregating_sites, n_sequences, avg_pairwise_diff) {
        Ok(result) => {
            let js = JsTajimaD {
                d: result.d,
                pi: result.pi,
                theta_w: result.theta_w,
                segregating_sites: result.segregating_sites,
                n_sequences: result.n_sequences,
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn describe_known_data() {
        let json = describe("[2,4,4,4,5,5,7,9]");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let stats = &v["ok"];
        assert_eq!(stats["count"], 8);
        assert!((stats["mean"].as_f64().unwrap() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn describe_empty_error() {
        let json = describe("[]");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    #[test]
    fn describe_invalid_json() {
        let json = describe("not json");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].as_str().unwrap().contains("invalid JSON"));
    }

    #[test]
    fn pearson_positive() {
        let json = pearson("[1,2,3,4,5]", "[2,4,6,8,10]");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!((v["ok"].as_f64().unwrap() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn pearson_mismatched_lengths() {
        let json = pearson("[1,2,3]", "[1,2]");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    #[test]
    fn t_test_known() {
        // Sample with known mean=3.0, testing against mu=0
        let json = t_test("[1,2,3,4,5]", 0.0);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let result = &v["ok"];
        assert!(result["statistic"].as_f64().unwrap() > 0.0);
        assert!(result["p_value"].as_f64().unwrap() < 0.05);
        assert_eq!(result["method"].as_str().unwrap(), "One-sample t-test");
    }

    #[test]
    fn spearman_perfect_rank() {
        let json = spearman("[1,2,3,4,5]", "[2,4,6,8,10]");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!((v["ok"].as_f64().unwrap() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn t_test_two_sample_known() {
        let json = t_test_two_sample("[1,2,3,4,5]", "[6,7,8,9,10]", false);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let result = &v["ok"];
        assert!(result["statistic"].as_f64().unwrap().abs() > 0.0);
        assert!(result["p_value"].as_f64().unwrap() < 0.05);
    }

    #[test]
    fn mann_whitney_known() {
        let json = mann_whitney_u("[1,2,3,4,5]", "[6,7,8,9,10]");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let result = &v["ok"];
        assert!(result["statistic"].is_number());
        assert!(result["p_value"].is_number());
    }

    #[test]
    fn bonferroni_known() {
        let json = bonferroni("[0.01, 0.04, 0.03]");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let corrected = v["ok"].as_array().unwrap();
        assert_eq!(corrected.len(), 3);
        // 0.01 * 3 = 0.03
        assert!((corrected[0].as_f64().unwrap() - 0.03).abs() < 1e-10);
    }

    #[test]
    fn benjamini_hochberg_known() {
        let json = benjamini_hochberg("[0.01, 0.04, 0.03]");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let corrected = v["ok"].as_array().unwrap();
        assert_eq!(corrected.len(), 3);
    }

    // ── Survival analysis ────────────────────────────────────────────────

    #[test]
    fn kaplan_meier_basic() {
        let json = kaplan_meier("[1,2,3,4,5]", "[true,false,true,false,true]");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let result = &v["ok"];
        let steps = result["steps"].as_array().unwrap();
        assert!(!steps.is_empty(), "expected non-empty KM steps");
        assert_eq!(result["n_total"], 5);
        assert_eq!(result["n_events"], 3);
    }

    #[test]
    fn kaplan_meier_all_events() {
        let json = kaplan_meier("[1,2,3,4,5]", "[true,true,true,true,true]");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let result = &v["ok"];
        assert_eq!(result["n_events"], 5);
        let steps = result["steps"].as_array().unwrap();
        assert_eq!(steps.len(), 5);
        // Last step survival should be 0.0
        let last = steps.last().unwrap();
        assert!((last["survival"].as_f64().unwrap() - 0.0).abs() < 1e-10);
    }

    #[test]
    fn log_rank_test_basic() {
        // Group 1: early events; Group 2: late events
        let json = log_rank_test(
            "[1,2,3,4]",
            "[true,true,true,true]",
            "[10,11,12,13]",
            "[true,true,true,true]",
        );
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let result = &v["ok"];
        assert_eq!(result["n_groups"], 2);
        assert!(result["statistic"].as_f64().unwrap() >= 0.0);
        assert!(result["p_value"].as_f64().unwrap() >= 0.0);
    }

    #[test]
    fn cox_ph_basic() {
        let json = cox_ph(
            "[1,2,3,4,5,6]",
            "[true,true,false,true,false,true]",
            "[0,1,0,1,0,1]",
            1,
        );
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let result = &v["ok"];
        assert!(result["converged"].as_bool().unwrap());
        assert_eq!(result["n_observations"], 6);
        assert_eq!(result["n_events"], 4);
        assert_eq!(result["coefficients"].as_array().unwrap().len(), 1);
    }

    // ── Null models ──────────────────────────────────────────────────────

    #[test]
    fn wright_fisher_basic() {
        let json = wright_fisher(100, 0.5, 10, 42);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let result = &v["ok"];
        let freqs = result["frequencies"].as_array().unwrap();
        // Should have n_generations + 1 entries (gen 0..=10)
        assert_eq!(freqs.len(), 11);
        let final_freq = result["final_freq"].as_f64().unwrap();
        assert!(final_freq >= 0.0 && final_freq <= 1.0);
    }

    #[test]
    fn permutation_test_basic() {
        let json = permutation_test("[1,2,3,4,5,6]", "[3,3]", 100, 42);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let dist = v["ok"].as_array().unwrap();
        assert_eq!(dist.len(), 100);
        for val in dist {
            assert!(val.as_f64().unwrap().is_finite());
        }
    }

    #[test]
    fn bootstrap_ci_basic() {
        let json = bootstrap_ci("[1,2,3,4,5]", 200, 42);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let dist = v["ok"].as_array().unwrap();
        assert_eq!(dist.len(), 200);
        for val in dist {
            let x = val.as_f64().unwrap();
            assert!(x >= 1.0 && x <= 5.0, "bootstrap mean {} out of range", x);
        }
    }

    // ── Diversity ────────────────────────────────────────────────────────

    #[test]
    fn shannon_index_basic() {
        // Uniform distribution: H = ln(3) for 3 equal-abundance species
        let json = shannon_index("[10,10,10]");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let h = v["ok"].as_f64().unwrap();
        let expected = (3.0f64).ln();
        assert!((h - expected).abs() < 1e-10, "H={}, expected={}", h, expected);
    }

    #[test]
    fn simpson_index_basic() {
        // Uniform distribution with 3 species, each count 10
        // D = 3 * 10*9 / (30*29) = 270/870
        let json = simpson_index("[10,10,10]");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let d = v["ok"].as_f64().unwrap();
        let expected = 270.0 / 870.0;
        assert!((d - expected).abs() < 1e-10, "D={}, expected={}", d, expected);
    }

    #[test]
    fn bray_curtis_basic() {
        // Identical samples should give BC = 0
        let json = bray_curtis("[1,2,3]", "[1,2,3]");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let bc = v["ok"].as_f64().unwrap();
        assert!(bc.abs() < 1e-10, "BC={}, expected 0", bc);
    }

    // ── Population genetics ──────────────────────────────────────────────

    #[test]
    fn fst_hudson_basic() {
        // Two differentiated populations: pop1 mostly ref, pop2 mostly alt
        let pop1 = "[[0,0,0,1,1]]"; // one locus, 5 individuals
        let pop2 = "[[1,1,2,2,2]]";
        let json = fst_hudson(pop1, pop2);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let result = &v["ok"];
        let fst = result["fst"].as_f64().unwrap();
        assert!(fst > 0.0 && fst <= 1.0, "fst={}", fst);
    }

    #[test]
    fn tajimas_d_basic() {
        // n=20, S=10, pi = theta_w → D ≈ 0
        let json = tajimas_d(10, 20, 3.0);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let result = &v["ok"];
        assert!(result["d"].as_f64().unwrap().is_finite());
        assert_eq!(result["segregating_sites"], 10);
        assert_eq!(result["n_sequences"], 20);
        assert!(result["pi"].as_f64().unwrap() > 0.0);
        assert!(result["theta_w"].as_f64().unwrap() > 0.0);
    }
}
