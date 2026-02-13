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
// Effect sizes
// ---------------------------------------------------------------------------

/// Cohen's d standardized mean difference between two groups.
#[pyfunction]
fn cohens_d(group1: Vec<f64>, group2: Vec<f64>) -> PyResult<f64> {
    cyanea_stats::effect_size::cohens_d(&group1, &group2).into_pyresult()
}

/// Eta-squared (proportion of variance explained) for one-way ANOVA.
#[pyfunction]
fn eta_squared(groups: Vec<Vec<f64>>) -> PyResult<f64> {
    let refs: Vec<&[f64]> = groups.iter().map(|g| g.as_slice()).collect();
    cyanea_stats::effect_size::eta_squared(&refs).into_pyresult()
}

/// Odds ratio from a 2×2 contingency table [[a, b], [c, d]].
#[pyfunction]
fn odds_ratio(table: [[usize; 2]; 2]) -> PyResult<f64> {
    cyanea_stats::effect_size::odds_ratio(&table).into_pyresult()
}

/// Relative risk from a 2×2 contingency table [[a, b], [c, d]].
#[pyfunction]
fn relative_risk(table: [[usize; 2]; 2]) -> PyResult<f64> {
    cyanea_stats::effect_size::relative_risk(&table).into_pyresult()
}

// ---------------------------------------------------------------------------
// Distributions
// ---------------------------------------------------------------------------

/// Normal distribution CDF at x with parameters mu and sigma.
#[pyfunction]
#[pyo3(signature = (x, mu=0.0, sigma=1.0))]
fn normal_cdf(x: f64, mu: f64, sigma: f64) -> PyResult<f64> {
    use cyanea_stats::distribution::Distribution;
    let dist = cyanea_stats::distribution::Normal::new(mu, sigma).into_pyresult()?;
    Ok(dist.cdf(x))
}

/// Normal distribution PDF at x with parameters mu and sigma.
#[pyfunction]
#[pyo3(signature = (x, mu=0.0, sigma=1.0))]
fn normal_pdf(x: f64, mu: f64, sigma: f64) -> PyResult<f64> {
    use cyanea_stats::distribution::Distribution;
    let dist = cyanea_stats::distribution::Normal::new(mu, sigma).into_pyresult()?;
    Ok(dist.pdf(x))
}

/// Error function (erf).
#[pyfunction]
fn erf(x: f64) -> f64 {
    cyanea_stats::distribution::erf(x)
}

/// Natural log of the gamma function.
#[pyfunction]
fn ln_gamma(x: f64) -> f64 {
    cyanea_stats::distribution::ln_gamma(x)
}

// ---------------------------------------------------------------------------
// Bayesian conjugate priors
// ---------------------------------------------------------------------------

/// Beta-binomial conjugate update.
///
/// Given prior (alpha, beta) and observed (successes, trials),
/// returns posterior (alpha', beta').
#[pyfunction]
fn bayesian_beta_update(alpha: f64, beta: f64, successes: u64, trials: u64) -> PyResult<(f64, f64)> {
    let prior = cyanea_stats::bayesian::Beta::new(alpha, beta).into_pyresult()?;
    let posterior = prior.update_binomial(successes, trials);
    Ok((posterior.alpha(), posterior.beta()))
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
    m.add_function(wrap_pyfunction!(cohens_d, &m)?)?;
    m.add_function(wrap_pyfunction!(eta_squared, &m)?)?;
    m.add_function(wrap_pyfunction!(odds_ratio, &m)?)?;
    m.add_function(wrap_pyfunction!(relative_risk, &m)?)?;
    m.add_function(wrap_pyfunction!(normal_cdf, &m)?)?;
    m.add_function(wrap_pyfunction!(normal_pdf, &m)?)?;
    m.add_function(wrap_pyfunction!(erf, &m)?)?;
    m.add_function(wrap_pyfunction!(ln_gamma, &m)?)?;
    m.add_function(wrap_pyfunction!(bayesian_beta_update, &m)?)?;
    parent.add_submodule(&m)?;
    Ok(())
}
