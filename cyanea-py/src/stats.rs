//! Python bindings for cyanea-stats: statistics and hypothesis testing.

use pyo3::prelude::*;

use crate::error::IntoPyResult;

// ---------------------------------------------------------------------------
// Result classes
// ---------------------------------------------------------------------------

/// Descriptive statistics (15 fields).
#[pyclass(frozen, get_all)]
pub struct DescriptiveStats {
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

impl From<cyanea_stats::DescriptiveStats> for DescriptiveStats {
    fn from(s: cyanea_stats::DescriptiveStats) -> Self {
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

/// Hypothesis test result.
#[pyclass(frozen, get_all)]
pub struct TestResult {
    pub statistic: f64,
    pub p_value: f64,
    pub degrees_of_freedom: Option<f64>,
    pub method: String,
}

impl From<cyanea_stats::TestResult> for TestResult {
    fn from(r: cyanea_stats::TestResult) -> Self {
        Self {
            statistic: r.statistic,
            p_value: r.p_value,
            degrees_of_freedom: r.degrees_of_freedom,
            method: r.method,
        }
    }
}

// ---------------------------------------------------------------------------
// Module functions
// ---------------------------------------------------------------------------

/// Compute descriptive statistics for a list of floats.
#[pyfunction]
fn describe(data: Vec<f64>) -> PyResult<DescriptiveStats> {
    cyanea_stats::descriptive::describe(&data)
        .map(DescriptiveStats::from)
        .into_pyresult()
}

/// Pearson product-moment correlation coefficient.
#[pyfunction]
fn pearson(x: Vec<f64>, y: Vec<f64>) -> PyResult<f64> {
    cyanea_stats::correlation::pearson(&x, &y).into_pyresult()
}

/// Spearman rank correlation coefficient.
#[pyfunction]
fn spearman(x: Vec<f64>, y: Vec<f64>) -> PyResult<f64> {
    cyanea_stats::correlation::spearman(&x, &y).into_pyresult()
}

/// One-sample t-test.
#[pyfunction]
#[pyo3(signature = (data, *, mu=0.0))]
fn t_test(data: Vec<f64>, mu: f64) -> PyResult<TestResult> {
    cyanea_stats::testing::t_test_one_sample(&data, mu)
        .map(TestResult::from)
        .into_pyresult()
}

/// Two-sample t-test (Student's or Welch's).
#[pyfunction]
#[pyo3(signature = (x, y, *, equal_var=false))]
fn t_test_two_sample(x: Vec<f64>, y: Vec<f64>, equal_var: bool) -> PyResult<TestResult> {
    cyanea_stats::testing::t_test_two_sample(&x, &y, equal_var)
        .map(TestResult::from)
        .into_pyresult()
}

/// Mann-Whitney U test (non-parametric).
#[pyfunction]
fn mann_whitney_u(x: Vec<f64>, y: Vec<f64>) -> PyResult<TestResult> {
    cyanea_stats::testing::mann_whitney_u(&x, &y)
        .map(TestResult::from)
        .into_pyresult()
}

/// Bonferroni p-value correction.
#[pyfunction]
fn bonferroni(p_values: Vec<f64>) -> PyResult<Vec<f64>> {
    cyanea_stats::correction::bonferroni(&p_values).into_pyresult()
}

/// Benjamini-Hochberg FDR correction.
#[pyfunction]
fn benjamini_hochberg(p_values: Vec<f64>) -> PyResult<Vec<f64>> {
    cyanea_stats::correction::benjamini_hochberg(&p_values).into_pyresult()
}

// ---------------------------------------------------------------------------
// Submodule registration
// ---------------------------------------------------------------------------

pub fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent.py(), "stats")?;
    m.add_class::<DescriptiveStats>()?;
    m.add_class::<TestResult>()?;
    m.add_function(wrap_pyfunction!(describe, &m)?)?;
    m.add_function(wrap_pyfunction!(pearson, &m)?)?;
    m.add_function(wrap_pyfunction!(spearman, &m)?)?;
    m.add_function(wrap_pyfunction!(t_test, &m)?)?;
    m.add_function(wrap_pyfunction!(t_test_two_sample, &m)?)?;
    m.add_function(wrap_pyfunction!(mann_whitney_u, &m)?)?;
    m.add_function(wrap_pyfunction!(bonferroni, &m)?)?;
    m.add_function(wrap_pyfunction!(benjamini_hochberg, &m)?)?;
    parent.add_submodule(&m)?;
    Ok(())
}
