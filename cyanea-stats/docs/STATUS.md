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

### Bayesian statistics (`bayesian.rs`)

Conjugate prior distributions for Bayesian updating.

| Type | Description |
|------|-------------|
| `Beta` | Beta distribution (conjugate prior for binomial likelihood) |
| `Gamma` | Gamma distribution (conjugate prior for Poisson likelihood) |
| `NormalConjugate` | Normal-Normal conjugate model (known observation variance) |
| `Dirichlet` | Dirichlet distribution (conjugate prior for multinomial likelihood) |

**Beta methods:**

| Method | Description |
|--------|-------------|
| `new(alpha, beta) -> Result<Self>` | Create with shape parameters (alpha, beta > 0) |
| `update_binomial(&self, successes, trials) -> Self` | Posterior after observing binomial data |

Implements `Distribution` trait (`pdf`, `cdf`, `mean`, `variance`).

**Gamma methods:**

| Method | Description |
|--------|-------------|
| `new(shape, rate) -> Result<Self>` | Create with shape and rate (both > 0) |
| `update_poisson(&self, count) -> Self` | Posterior after observing a Poisson count |
| `update_poisson_batch(&self, counts) -> Self` | Posterior after observing multiple Poisson counts |

Implements `Distribution` trait (`pdf`, `cdf`, `mean`, `variance`).

**NormalConjugate methods:**

| Method | Description |
|--------|-------------|
| `new(prior_mu, prior_var, obs_var) -> Result<Self>` | Create with prior mean, prior variance, and observation variance |
| `update(&self, observation) -> Self` | Posterior after one observation |
| `update_batch(&self, observations) -> Self` | Posterior after multiple observations |
| `posterior_mean(&self) -> f64` | Posterior mean |
| `posterior_variance(&self) -> f64` | Posterior variance |

**Dirichlet methods:**

| Method | Description |
|--------|-------------|
| `new(alpha: Vec<f64>) -> Result<Self>` | Create from concentration parameters (all > 0) |
| `symmetric(k, alpha) -> Result<Self>` | Create symmetric Dirichlet with k categories |
| `update_multinomial(&self, counts) -> Self` | Posterior after observing multinomial counts |
| `mean(&self) -> Vec<f64>` | Expected category probabilities |
| `variance(&self) -> Vec<f64>` | Variance for each category |
| `ln_pdf(&self, x) -> Result<f64>` | Log probability density |

### Combinatorics (`combinatorics.rs`)

Exact and log-space combinatorial functions with overflow protection.

| Function | Description |
|----------|-------------|
| `factorial(n) -> Option<u64>` | Exact factorial (None if n > 20) |
| `ln_factorial(n) -> f64` | Log-space factorial via ln_gamma(n+1) |
| `binomial(n, k) -> Option<u64>` | Exact binomial coefficient with overflow check |
| `ln_binomial(n, k) -> Result<f64>` | Log-space binomial coefficient |
| `permutations(n, k) -> Option<u64>` | Exact k-permutations of n |
| `ln_permutations(n, k) -> Result<f64>` | Log-space permutations |
| `multinomial(n, counts) -> Option<u64>` | Exact multinomial coefficient |
| `ln_multinomial(n, counts) -> Result<f64>` | Log-space multinomial coefficient |
| `combinations(n, k) -> Combinations` | Iterator over all k-element subsets of [0, n) |

| Type | Description |
|------|-------------|
| `Combinations` | Iterator yielding `Vec<usize>` subsets in lexicographic order |

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

### Population genetics (`popgen.rs`)

| Type/Function | Description |
|---------------|-------------|
| `AlleleFrequencies` | Per-locus allele frequency summary (major/minor allele, MAF, observed/expected heterozygosity) |
| `allele_frequencies(genotypes, n_samples) -> Result<Vec<AlleleFrequencies>>` | Compute allele frequencies from 0/1/2 genotype matrix |
| `HweResult` | Hardy-Weinberg equilibrium test result (chi-squared statistic, p-value) |
| `hwe_test(genotypes, n_samples) -> Result<Vec<HweResult>>` | Per-locus HWE exact test |
| `FstMethod` | Enum: `WeirCockerham`, `Hudson` |
| `FstResult` | Fst estimate with per-locus values |
| `fst(genotypes, n_samples, populations, method) -> Result<FstResult>` | Fixation index between populations |
| `DiversityStats` | Nucleotide diversity: pi, theta (Watterson), S (segregating sites) |
| `nucleotide_diversity(genotypes, n_samples, seq_length) -> Result<DiversityStats>` | Pi, theta, S |
| `TajimaD` | Tajima's D statistic with variance components |
| `tajimas_d(genotypes, n_samples, seq_length) -> Result<TajimaD>` | Tajima's D neutrality test |
| `LdResult` | Linkage disequilibrium: r², D, D' |
| `ld(genotypes, n_samples, locus_a, locus_b) -> Result<LdResult>` | Pairwise LD via EM haplotype estimation |
| `genotype_pca(genotypes, n_samples, n_components) -> Result<PcaResult>` | PCA on centered genotype matrix |

### Normalization (`normalization.rs`)

| Function | Description |
|----------|-------------|
| `tpm(counts, lengths) -> Result<Vec<f64>>` | Transcripts Per Million |
| `fpkm(counts, lengths, total_reads) -> Result<Vec<f64>>` | Fragments Per Kilobase per Million |
| `cpm(counts) -> Result<Vec<f64>>` | Counts Per Million |
| `deseq2_size_factors(matrix, n_genes, n_samples) -> Result<Vec<f64>>` | DESeq2-style median-of-ratios normalization |

### Differential expression (`diffexpr.rs`)

| Type/Function | Description |
|---------------|-------------|
| `DeMethod` | Enum: `NegBinomialWald`, `Wilcoxon` |
| `DeGeneResult` | Per-gene result: log2FC, p-value, adjusted p-value, base mean, test statistic |
| `DeResults` | Collection of `DeGeneResult` with metadata |
| `VolcanoPoint` | Data point for volcano plots (log2FC, −log10 p-value, significance flag) |
| `de_test(matrix, n_genes, n_samples, group, method) -> Result<DeResults>` | Run differential expression analysis |
| `volcano_data(results, fc_threshold, p_threshold) -> Vec<VolcanoPoint>` | Generate volcano plot data |

### Enrichment analysis (`enrichment.rs`)

| Type/Function | Description |
|---------------|-------------|
| `GeneSet` | Named gene set: `name` (String) + `genes` (Vec<usize>, 0-based indices) |
| `OraResult` | Per-set ORA result: `gene_set`, `overlap`, `expected`, `gene_set_size`, `p_value`, `p_adjusted` |
| `ora(significant, gene_sets, n_total) -> Result<Vec<OraResult>>` | Over-representation analysis via hypergeometric upper-tail test with BH correction |
| `GseaResult` | Per-set GSEA result: `gene_set`, `enrichment_score`, `normalized_es`, `p_value`, `p_adjusted`, `leading_edge_size`, `gene_set_size` |
| `gsea_preranked(genes, scores, gene_sets, weight, n_permutations) -> Result<Vec<GseaResult>>` | Preranked GSEA (Subramanian et al. 2005) with score weighting, permutation p-values, NES, and BH correction |
| `GoNamespace` | Enum: `BiologicalProcess`, `MolecularFunction`, `CellularComponent` |
| `GoTerm` | GO term: `id`, `name`, `namespace`, `genes` (Vec<usize>) |
| `GoAnnotation` | Collection of GO term annotations with `new()`, `from_entries()`, `n_terms()`, `n_genes()`, `terms_for_gene()`, `filter_namespace()` |
| `GoEnrichmentConfig` | Config: `min_genes` (default 5), `max_genes` (default 500), `namespace` (optional filter) |
| `GoEnrichmentResult` | Per-term result: `term_id`, `term_name`, `namespace`, `overlap`, `expected`, `gene_set_size`, `p_value`, `p_adjusted` |
| `go_enrichment(significant, annotation, n_total, config) -> Result<Vec<GoEnrichmentResult>>` | GO enrichment via ORA — filters by namespace/size, delegates to `ora()`, returns BH-corrected results |

### Survival analysis (`survival.rs`)

| Type/Function | Description |
|---------------|-------------|
| `KmStep` | Single step: `time`, `n_risk`, `n_events`, `n_censored`, `survival`, `std_err`, `ci_lower`, `ci_upper` |
| `KmResult` | Full result: `steps`, `median_survival`, `n_total`, `n_events` |
| `kaplan_meier(times, status) -> Result<KmResult>` | Kaplan-Meier survival curve with Greenwood SE and log-transform 95% CI |
| `LogRankResult` | Chi-squared statistic, df, p-value, observed/expected per group |
| `log_rank_test(times, status, groups) -> Result<LogRankResult>` | Log-rank test comparing survival between groups |
| `CoxPhResult` | Coefficients, SE, z-values, p-values, hazard ratios with 95% CI, convergence info |
| `cox_ph(times, status, covariates, n_covariates) -> Result<CoxPhResult>` | Cox proportional hazards via Breslow partial likelihood (Newton-Raphson) |

## Tests

319 unit tests + 9 doc tests across 17 source files.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 62 | Module declarations, re-exports |
| `descriptive.rs` | 412 | Descriptive statistics |
| `correlation.rs` | 272 | Pearson/Spearman correlation |
| `testing.rs` | 604 | t-tests, Mann-Whitney U, Fisher's exact, chi-squared, ANOVA |
| `distribution.rs` | 974 | Normal, Poisson, Binomial, ChiSquared, F, NegativeBinomial, erf, gamma, beta |
| `correction.rs` | 161 | Bonferroni, Benjamini-Hochberg |
| `rank.rs` | 141 | Rank computation with tie-breaking |
| `effect_size.rs` | 252 | Cohen's d, eta-squared, odds ratio, relative risk |
| `reduction.rs` | 283 | PCA dimensionality reduction |
| `blas_pca.rs` | 117 | BLAS-accelerated PCA via ndarray (feature-gated, internal to reduction) |
| `bayesian.rs` | 568 | Bayesian conjugate priors (Beta, Gamma, NormalConjugate, Dirichlet) |
| `combinatorics.rs` | 311 | Combinatorial functions (factorial, binomial, permutations, multinomial, combinations) |
| `popgen.rs` | 1445 | Population genetics (allele freq, HWE, Fst, diversity, Tajima's D, LD, PCA) |
| `normalization.rs` | 357 | Expression normalization (TPM, FPKM, CPM, DESeq2 size factors) |
| `diffexpr.rs` | 546 | Differential expression (NB Wald test, Wilcoxon, volcano plot) |
| `enrichment.rs` | ~760 | Gene set enrichment (ORA, GSEA preranked, GO enrichment) |
| `survival.rs` | 1164 | Survival analysis (Kaplan-Meier, log-rank test, Cox PH) |
