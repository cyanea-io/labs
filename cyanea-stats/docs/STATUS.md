# cyanea-stats

Statistical methods for life sciences: descriptive statistics, correlation, hypothesis testing, probability distributions, multiple testing correction, and effect sizes.

## Status: Complete

All statistical functionality is implemented including descriptive statistics, correlation, hypothesis testing, distributions, multiple testing correction, effect sizes, and PCA dimensionality reduction.

## Public API

### Descriptive statistics (`descriptive.rs`)

| Type/Function | Description |
|---------------|-------------|
| `DescriptiveStats` | 15 fields: `count`, `mean`, `median`, `variance`, `sample_variance`, `std_dev`, `sample_std_dev`, `min`, `max`, `range`, `q1`, `q3`, `iqr`, `skewness`, `kurtosis` |
| `describe(data) -> Result<DescriptiveStats>` | Compute all 15 statistics |
| `mean(data) -> Result<f64>` | Arithmetic mean |
| `median(data) -> Result<f64>` | Median (interpolated) |
| `variance(data, ddof) -> Result<f64>` | Population or sample variance |
| `std_dev(data, ddof) -> Result<f64>` | Standard deviation |
| `quantile(data, q) -> Result<f64>` | Arbitrary quantile (0.0-1.0) |
| `iqr(data) -> Result<f64>` | Interquartile range |
| `mad(data) -> Result<f64>` | Median absolute deviation |

### Correlation (`correlation.rs`)

| Function | Description |
|----------|-------------|
| `pearson(x, y) -> Result<f64>` | Pearson product-moment correlation |
| `spearman(x, y) -> Result<f64>` | Spearman rank correlation |
| `CorrelationMatrix` | Pairwise correlation/covariance matrix computation |

### Ranking (`rank.rs`)

| Type/Function | Description |
|---------------|-------------|
| `RankMethod` | Enum: `Average`, `Min`, `Max`, `First` |
| `rank(data, method) -> Vec<f64>` | Compute ranks with tie-breaking |

### Hypothesis testing (`testing.rs`)

| Type/Function | Description |
|---------------|-------------|
| `TestResult` | `statistic`, `p_value`, `degrees_of_freedom`, `method` |
| `t_test_one_sample(data, mu) -> Result<TestResult>` | One-sample t-test |
| `t_test_two_sample(x, y, equal_var) -> Result<TestResult>` | Student's (equal_var=true) or Welch's (equal_var=false) |
| `mann_whitney_u(x, y) -> Result<TestResult>` | Non-parametric rank-sum test |

### Distributions (`distribution.rs`)

| Type/Function | Description |
|---------------|-------------|
| `trait Distribution` | `cdf(x)`, `pdf(x)` methods |
| `Normal` | Gaussian distribution (mean, std_dev) |
| `Poisson` | Poisson distribution (lambda) |
| `Binomial` | Binomial distribution (n, p) |
| `ChiSquared` | Chi-squared distribution (df) |
| `FDistribution` | F distribution (df1, df2) |
| `erf(x) -> f64` | Error function |
| `ln_gamma(x) -> f64` | Log gamma function (Lanczos approximation) |
| `betai(a, b, x) -> Result<f64>` | Regularized incomplete beta function |

### Multiple testing correction (`correction.rs`)

| Type/Function | Description |
|---------------|-------------|
| `CorrectionMethod` | Enum: `Bonferroni`, `BenjaminiHochberg` |
| `correct(p_values, method) -> Result<Vec<f64>>` | Apply correction by method |
| `bonferroni(p_values) -> Result<Vec<f64>>` | Bonferroni correction |
| `benjamini_hochberg(p_values) -> Result<Vec<f64>>` | BH FDR correction |

### Effect sizes (`effect_size.rs`)

| Function | Description |
|----------|-------------|
| `cohens_d(group1, group2) -> Result<f64>` | Cohen's d: standardized mean difference (pooled SD). Each group needs >= 2 observations |
| `eta_squared(groups) -> Result<f64>` | Eta-squared: proportion of variance explained by group membership (ANOVA). Returns value in [0, 1] |
| `odds_ratio(table) -> Result<f64>` | Odds ratio for a 2x2 contingency table `[[a, b], [c, d]]`. OR = (a*d) / (b*c) |
| `relative_risk(table) -> Result<f64>` | Relative risk (risk ratio) for a 2x2 contingency table. RR = [a/(a+b)] / [c/(c+d)] |

### Dimensionality reduction (`reduction.rs`)

| Type/Function | Description |
|---------------|-------------|
| `PcaResult` | `components`, `explained_variance`, `explained_variance_ratio`, `transformed`, `mean` |
| `pca(data, n_components) -> Result<PcaResult>` | Principal Component Analysis |

When the `blas` feature is enabled, `pca` dispatches to an ndarray-based implementation (`blas_pca.rs`) that uses matrix multiply for covariance computation and power iteration, which can leverage BLAS for acceleration.

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support |
| `wasm` | No | WASM target |
| `serde` | No | Serialization support |
| `parallel` | No | Rayon parallelism |
| `blas` | No | BLAS-accelerated PCA via ndarray |

## Dependencies

- `cyanea-core` -- error types
- `ndarray` (optional, `blas` feature) -- matrix operations for BLAS-accelerated PCA
- `rayon` (optional, `parallel` feature) -- data parallelism

## Tests

127 tests across 10 source files.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 45 | Module declarations, re-exports |
| `descriptive.rs` | 354 | Descriptive statistics |
| `correlation.rs` | 244 | Pearson/Spearman correlation |
| `testing.rs` | 605 | t-tests, Mann-Whitney U, Fisher's exact, chi-squared, ANOVA |
| `distribution.rs` | 781 | Normal, Poisson, Binomial, ChiSquared, F, erf, gamma, beta |
| `correction.rs` | 161 | Bonferroni, Benjamini-Hochberg |
| `rank.rs` | 141 | Rank computation with tie-breaking |
| `effect_size.rs` | 253 | Cohen's d, eta-squared, odds ratio, relative risk |
| `reduction.rs` | 284 | PCA dimensionality reduction |
| `blas_pca.rs` | 117 | BLAS-accelerated PCA via ndarray (feature-gated, internal to reduction) |
