//! Statistical wrappers with JSON input/output.
//!
//! Types like [`DescriptiveStats`] and [`TestResult`] do not derive `Serialize`,
//! so we provide thin wrappers (`Js*` types) that mirror their fields and add
//! `#[derive(Serialize)]`.

use serde::Serialize;

use cyanea_stats::correlation;
use cyanea_stats::descriptive;
use cyanea_stats::testing;

use crate::error::{wasm_err, wasm_ok};

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

// ── JSON boundary functions ──────────────────────────────────────────────

fn parse_f64_array(json: &str) -> Result<Vec<f64>, String> {
    serde_json::from_str::<Vec<f64>>(json).map_err(|e| format!("invalid JSON array: {e}"))
}

/// Compute descriptive statistics from a JSON array of numbers.
///
/// Input: `"[1.0, 2.0, 3.0]"` — Output: JSON `JsDescriptiveStats`.
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
}
