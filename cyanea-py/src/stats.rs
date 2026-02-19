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
// Survival analysis
// ---------------------------------------------------------------------------

/// A single step in the Kaplan-Meier survival curve.
#[pyclass(frozen, get_all)]
pub struct PyKmStep {
    pub time: f64,
    pub survival: f64,
    pub n_at_risk: usize,
    pub n_events: usize,
}

/// Result of Kaplan-Meier survival estimation.
#[pyclass(frozen, get_all)]
pub struct PyKmResult {
    pub steps: Vec<pyo3::Py<PyKmStep>>,
    pub median_survival: Option<f64>,
    pub total_events: usize,
    pub total_censored: usize,
}

/// Result of the log-rank test comparing two survival curves.
#[pyclass(frozen, get_all)]
pub struct PyLogRankResult {
    pub statistic: f64,
    pub p_value: f64,
    pub degrees_of_freedom: usize,
}

/// Result of Cox proportional hazards regression.
#[pyclass(frozen, get_all)]
pub struct PyCoxPhResult {
    pub coefficients: Vec<f64>,
    pub standard_errors: Vec<f64>,
    pub hazard_ratios: Vec<f64>,
    pub log_likelihood: f64,
}

/// Kaplan-Meier survival estimator.
#[pyfunction]
fn kaplan_meier(py: Python<'_>, times: Vec<f64>, status: Vec<bool>) -> PyResult<PyKmResult> {
    let km = cyanea_stats::survival::kaplan_meier(&times, &status).into_pyresult()?;
    let steps: Vec<pyo3::Py<PyKmStep>> = km
        .steps
        .iter()
        .map(|s| {
            pyo3::Py::new(
                py,
                PyKmStep {
                    time: s.time,
                    survival: s.survival,
                    n_at_risk: s.n_risk,
                    n_events: s.n_events,
                },
            )
        })
        .collect::<PyResult<Vec<_>>>()?;
    let total_censored = km.n_total.saturating_sub(km.n_events);
    Ok(PyKmResult {
        steps,
        median_survival: km.median_survival,
        total_events: km.n_events,
        total_censored,
    })
}

/// Log-rank test comparing two survival curves.
#[pyfunction]
fn log_rank_test(
    t1: Vec<f64>,
    s1: Vec<bool>,
    t2: Vec<f64>,
    s2: Vec<bool>,
) -> PyResult<PyLogRankResult> {
    let mut times = t1;
    let mut status = s1;
    let n1 = times.len();
    times.extend(t2);
    status.extend(s2);
    let mut groups: Vec<usize> = vec![0; n1];
    groups.extend(vec![1; times.len() - n1]);
    let lr =
        cyanea_stats::survival::log_rank_test(&times, &status, &groups).into_pyresult()?;
    Ok(PyLogRankResult {
        statistic: lr.statistic,
        p_value: lr.p_value,
        degrees_of_freedom: lr.degrees_of_freedom,
    })
}

/// Cox proportional hazards regression.
#[pyfunction]
fn cox_ph(
    times: Vec<f64>,
    status: Vec<bool>,
    covariates: Vec<f64>,
    n_features: usize,
) -> PyResult<PyCoxPhResult> {
    let result =
        cyanea_stats::survival::cox_ph(&times, &status, &covariates, n_features)
            .into_pyresult()?;
    Ok(PyCoxPhResult {
        coefficients: result.coefficients,
        standard_errors: result.std_errors,
        hazard_ratios: result.hazard_ratios,
        log_likelihood: result.log_likelihood,
    })
}

// ---------------------------------------------------------------------------
// Null models
// ---------------------------------------------------------------------------

/// Result of a Wright-Fisher drift simulation.
#[pyclass(frozen, get_all)]
pub struct PyWrightFisherResult {
    pub trajectory: Vec<f64>,
    pub fixation_generation: Option<usize>,
    pub final_frequency: f64,
}

/// Wright-Fisher allele frequency drift simulation.
#[pyfunction]
fn wright_fisher(
    pop_size: usize,
    init_freq: f64,
    n_gens: usize,
    seed: u64,
) -> PyResult<PyWrightFisherResult> {
    let result =
        cyanea_stats::null_model::wright_fisher(pop_size, init_freq, n_gens, seed)
            .into_pyresult()?;
    Ok(PyWrightFisherResult {
        trajectory: result.frequencies,
        fixation_generation: result.fixation_gen,
        final_frequency: result.final_freq,
    })
}

/// Permutation test: computes the difference-of-means statistic under permutation.
///
/// Returns a TestResult with the observed statistic, a two-sided p-value,
/// and degrees_of_freedom = None.
#[pyfunction]
fn permutation_test(
    values: Vec<f64>,
    group_sizes: Vec<usize>,
    n_perms: usize,
    seed: u64,
) -> PyResult<TestResult> {
    // Compute observed statistic: difference of group means
    let diff_means = |groups: &[&[f64]]| -> f64 {
        if groups.len() < 2 {
            return 0.0;
        }
        let m0: f64 = groups[0].iter().sum::<f64>() / groups[0].len().max(1) as f64;
        let m1: f64 = groups[1].iter().sum::<f64>() / groups[1].len().max(1) as f64;
        (m0 - m1).abs()
    };

    // Compute observed statistic from the original grouping
    let mut offset = 0;
    let mut groups: Vec<&[f64]> = Vec::with_capacity(group_sizes.len());
    for &sz in &group_sizes {
        groups.push(&values[offset..offset + sz]);
        offset += sz;
    }
    let observed = diff_means(&groups);

    let null_dist = cyanea_stats::null_model::permutation_null(
        &values,
        &group_sizes,
        &diff_means,
        n_perms,
        seed,
    )
    .into_pyresult()?;

    let count = null_dist.iter().filter(|&&x| x >= observed).count();
    let p_value = (count as f64 / n_perms as f64).max(1.0 / n_perms as f64);

    Ok(TestResult {
        statistic: observed,
        p_value,
        degrees_of_freedom: None,
        method: "permutation_test".to_string(),
    })
}

/// Bootstrap confidence interval: returns (lower, mean, upper) at 95% CI.
#[pyfunction]
fn bootstrap_ci(data: Vec<f64>, n_bootstrap: usize, seed: u64) -> PyResult<(f64, f64, f64)> {
    let mean_fn = |d: &[f64]| d.iter().sum::<f64>() / d.len().max(1) as f64;
    let mut null_dist =
        cyanea_stats::null_model::bootstrap_null(&data, &mean_fn, n_bootstrap, seed);
    null_dist.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = null_dist.len();
    if n == 0 {
        return Ok((f64::NAN, f64::NAN, f64::NAN));
    }
    let lower_idx = ((0.025 * n as f64) as usize).min(n - 1);
    let upper_idx = ((0.975 * n as f64) as usize).min(n - 1);
    let mean = null_dist.iter().sum::<f64>() / n as f64;
    Ok((null_dist[lower_idx], mean, null_dist[upper_idx]))
}

// ---------------------------------------------------------------------------
// Diversity indices
// ---------------------------------------------------------------------------

/// Shannon diversity index H = -sum(p_i * ln(p_i)).
#[pyfunction]
fn shannon_index(counts: Vec<f64>) -> PyResult<f64> {
    let int_counts: Vec<usize> = counts.iter().map(|&c| c as usize).collect();
    cyanea_stats::diversity::shannon_index(&int_counts).into_pyresult()
}

/// Simpson's diversity index D = sum(n_i*(n_i-1)) / (N*(N-1)).
#[pyfunction]
fn simpson_index(counts: Vec<f64>) -> PyResult<f64> {
    let int_counts: Vec<usize> = counts.iter().map(|&c| c as usize).collect();
    cyanea_stats::diversity::simpson_index(&int_counts).into_pyresult()
}

/// Bray-Curtis dissimilarity between two samples.
#[pyfunction]
fn bray_curtis(a: Vec<f64>, b: Vec<f64>) -> PyResult<f64> {
    let a_int: Vec<usize> = a.iter().map(|&c| c as usize).collect();
    let b_int: Vec<usize> = b.iter().map(|&c| c as usize).collect();
    cyanea_stats::diversity::bray_curtis(&a_int, &b_int).into_pyresult()
}

// ---------------------------------------------------------------------------
// Population genetics
// ---------------------------------------------------------------------------

/// Hudson's Fst estimator between two populations.
///
/// Each element in `pop1` and `pop2` is an allele frequency (alternate allele)
/// at one locus. Internally converts to genotype representation for the
/// Hudson estimator.
#[pyfunction]
fn fst_hudson(pop1: Vec<f64>, pop2: Vec<f64>) -> PyResult<f64> {
    // Convert allele frequency vectors to Option<u8> genotype slices.
    // Each frequency is treated as a single-locus summary.
    // We create synthetic genotype arrays from the frequency:
    // one locus per entry, using 100 synthetic diploid individuals.
    let n_synthetic = 100usize;
    let mut loci_pop1: Vec<Vec<Option<u8>>> = Vec::with_capacity(pop1.len());
    let mut loci_pop2: Vec<Vec<Option<u8>>> = Vec::with_capacity(pop2.len());

    for &freq in &pop1 {
        let n_alt = (freq * 2.0 * n_synthetic as f64).round() as usize;
        let mut genotypes = Vec::with_capacity(n_synthetic);
        let mut alt_remaining = n_alt;
        for _ in 0..n_synthetic {
            let g = if alt_remaining >= 2 {
                alt_remaining -= 2;
                2u8
            } else if alt_remaining == 1 {
                alt_remaining -= 1;
                1u8
            } else {
                0u8
            };
            genotypes.push(Some(g));
        }
        loci_pop1.push(genotypes);
    }

    for &freq in &pop2 {
        let n_alt = (freq * 2.0 * n_synthetic as f64).round() as usize;
        let mut genotypes = Vec::with_capacity(n_synthetic);
        let mut alt_remaining = n_alt;
        for _ in 0..n_synthetic {
            let g = if alt_remaining >= 2 {
                alt_remaining -= 2;
                2u8
            } else if alt_remaining == 1 {
                alt_remaining -= 1;
                1u8
            } else {
                0u8
            };
            genotypes.push(Some(g));
        }
        loci_pop2.push(genotypes);
    }

    let refs1: Vec<&[Option<u8>]> = loci_pop1.iter().map(|v| v.as_slice()).collect();
    let refs2: Vec<&[Option<u8>]> = loci_pop2.iter().map(|v| v.as_slice()).collect();

    let result = cyanea_stats::popgen::fst_hudson(&refs1, &refs2).into_pyresult()?;
    Ok(result.fst)
}

/// Tajima's D statistic.
///
/// Takes a genotype matrix where each inner Vec<u8> is a sequence of
/// 0/1 (ancestral/derived) allele calls at each site.
/// Computes segregating sites, average pairwise differences, and D.
#[pyfunction]
fn tajimas_d(genotypes: Vec<Vec<u8>>) -> PyResult<f64> {
    if genotypes.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "genotypes must be non-empty",
        ));
    }
    let n_sequences = genotypes.len();
    let n_sites = genotypes[0].len();

    // Compute segregating sites and average pairwise differences
    let refs: Vec<&[u8]> = genotypes.iter().map(|v| v.as_slice()).collect();
    let mut seg_sites = 0usize;
    let mut total_diffs = 0.0_f64;
    let n_pairs = (n_sequences * (n_sequences - 1)) / 2;

    for site in 0..n_sites {
        let mut has_variation = false;
        let first = refs[0][site];
        for seq in &refs {
            if seq[site] != first {
                has_variation = true;
                break;
            }
        }
        if has_variation {
            seg_sites += 1;
            // Count pairwise differences at this site
            let mut diffs_at_site = 0usize;
            for i in 0..n_sequences {
                for j in (i + 1)..n_sequences {
                    if refs[i][site] != refs[j][site] {
                        diffs_at_site += 1;
                    }
                }
            }
            total_diffs += diffs_at_site as f64;
        }
    }

    let avg_pairwise_diff = if n_pairs > 0 {
        total_diffs / n_pairs as f64
    } else {
        0.0
    };

    let result =
        cyanea_stats::popgen::tajimas_d(seg_sites, n_sequences, avg_pairwise_diff)
            .into_pyresult()?;
    Ok(result.d)
}

// ---------------------------------------------------------------------------
// Enrichment analysis
// ---------------------------------------------------------------------------

/// Result of GSEA preranked analysis for one gene set.
#[pyclass(frozen, get_all)]
pub struct PyGseaResult {
    pub name: String,
    pub es: f64,
    pub nes: f64,
    pub p_value: f64,
}

/// Result of over-representation analysis for one gene set.
#[pyclass(frozen, get_all)]
pub struct PyOraResult {
    pub name: String,
    pub p_value: f64,
    pub odds_ratio: f64,
    pub overlap: usize,
}

/// GSEA preranked enrichment analysis.
///
/// `ranked` is a list of (gene_name, score) tuples in rank order.
/// `gene_sets` is a dict mapping set names to lists of gene names.
/// `n_perms` is the number of permutations for significance estimation.
#[pyfunction]
fn gsea_preranked(
    ranked: Vec<(String, f64)>,
    gene_sets: std::collections::HashMap<String, Vec<String>>,
    n_perms: usize,
) -> PyResult<Vec<PyGseaResult>> {
    // Map gene names to integer indices
    let gene_indices: Vec<usize> = (0..ranked.len()).collect();
    let scores: Vec<f64> = ranked.iter().map(|(_, s)| *s).collect();
    let name_to_idx: std::collections::HashMap<&str, usize> = ranked
        .iter()
        .enumerate()
        .map(|(i, (name, _))| (name.as_str(), i))
        .collect();

    let gs: Vec<cyanea_stats::enrichment::GeneSet> = gene_sets
        .iter()
        .map(|(name, genes)| cyanea_stats::enrichment::GeneSet {
            name: name.clone(),
            genes: genes
                .iter()
                .filter_map(|g| name_to_idx.get(g.as_str()).copied())
                .collect(),
        })
        .collect();

    let results = cyanea_stats::enrichment::gsea_preranked(
        &gene_indices,
        &scores,
        &gs,
        1.0,
        n_perms,
    )
    .into_pyresult()?;

    Ok(results
        .into_iter()
        .map(|r| PyGseaResult {
            name: r.gene_set,
            es: r.enrichment_score,
            nes: r.normalized_es,
            p_value: r.p_value,
        })
        .collect())
}

/// Over-representation analysis (hypergeometric test).
///
/// `genes` is a list of significant gene names.
/// `gene_sets` is a dict mapping set names to lists of gene names.
/// `background` is the total number of genes in the universe.
#[pyfunction]
fn ora(
    genes: Vec<String>,
    gene_sets: std::collections::HashMap<String, Vec<String>>,
    background: usize,
) -> PyResult<Vec<PyOraResult>> {
    // Build a mapping from gene name to index (assigning sequential indices).
    let mut name_to_idx: std::collections::HashMap<String, usize> =
        std::collections::HashMap::new();
    let mut next_idx = 0usize;

    // Register significant genes first
    for g in &genes {
        if !name_to_idx.contains_key(g) {
            name_to_idx.insert(g.clone(), next_idx);
            next_idx += 1;
        }
    }
    // Register gene set members
    for members in gene_sets.values() {
        for g in members {
            if !name_to_idx.contains_key(g) {
                name_to_idx.insert(g.clone(), next_idx);
                next_idx += 1;
            }
        }
    }

    let sig_indices: Vec<usize> = genes
        .iter()
        .filter_map(|g| name_to_idx.get(g).copied())
        .collect();

    let gs: Vec<cyanea_stats::enrichment::GeneSet> = gene_sets
        .iter()
        .map(|(name, members)| cyanea_stats::enrichment::GeneSet {
            name: name.clone(),
            genes: members
                .iter()
                .filter_map(|g| name_to_idx.get(g).copied())
                .collect(),
        })
        .collect();

    let results =
        cyanea_stats::enrichment::ora(&sig_indices, &gs, background).into_pyresult()?;

    Ok(results
        .into_iter()
        .map(|r| {
            // Compute a simple odds ratio from overlap counts
            let k = r.overlap;
            let n = sig_indices.len();
            let big_k = r.gene_set_size;
            let big_n = background;
            let a = k;
            let b = n.saturating_sub(k);
            let c = big_k.saturating_sub(k);
            let d = big_n.saturating_sub(n).saturating_sub(c);
            let or = if b > 0 && c > 0 {
                (a as f64 * d as f64) / (b as f64 * c as f64)
            } else {
                f64::INFINITY
            };
            PyOraResult {
                name: r.gene_set,
                p_value: r.p_value,
                odds_ratio: or,
                overlap: r.overlap,
            }
        })
        .collect())
}

// ---------------------------------------------------------------------------
// Submodule registration
// ---------------------------------------------------------------------------

pub fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent.py(), "stats")?;
    m.add_class::<DescriptiveStats>()?;
    m.add_class::<TestResult>()?;
    m.add_class::<PyKmStep>()?;
    m.add_class::<PyKmResult>()?;
    m.add_class::<PyLogRankResult>()?;
    m.add_class::<PyCoxPhResult>()?;
    m.add_class::<PyWrightFisherResult>()?;
    m.add_class::<PyGseaResult>()?;
    m.add_class::<PyOraResult>()?;
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
    // Survival analysis
    m.add_function(wrap_pyfunction!(kaplan_meier, &m)?)?;
    m.add_function(wrap_pyfunction!(log_rank_test, &m)?)?;
    m.add_function(wrap_pyfunction!(cox_ph, &m)?)?;
    // Null models
    m.add_function(wrap_pyfunction!(wright_fisher, &m)?)?;
    m.add_function(wrap_pyfunction!(permutation_test, &m)?)?;
    m.add_function(wrap_pyfunction!(bootstrap_ci, &m)?)?;
    // Diversity indices
    m.add_function(wrap_pyfunction!(shannon_index, &m)?)?;
    m.add_function(wrap_pyfunction!(simpson_index, &m)?)?;
    m.add_function(wrap_pyfunction!(bray_curtis, &m)?)?;
    // Population genetics
    m.add_function(wrap_pyfunction!(fst_hudson, &m)?)?;
    m.add_function(wrap_pyfunction!(tajimas_d, &m)?)?;
    // Enrichment analysis
    m.add_function(wrap_pyfunction!(gsea_preranked, &m)?)?;
    m.add_function(wrap_pyfunction!(ora, &m)?)?;
    parent.add_submodule(&m)?;
    Ok(())
}
