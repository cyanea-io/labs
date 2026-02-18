# cyanea-stats

> Statistical methods for life sciences: descriptive statistics, hypothesis testing, distributions, and domain-specific analysis.

## What's Inside

- **Descriptive statistics** -- mean, median, variance, std dev, quantiles, IQR, MAD, skewness, kurtosis
- **Correlation** -- Pearson, Spearman rank, pairwise correlation matrices
- **Hypothesis testing** -- t-tests (one-sample, two-sample/Welch's), Mann-Whitney U, Fisher's exact, chi-squared, ANOVA
- **Distributions** -- Normal, Poisson, Binomial, ChiSquared, F, NegativeBinomial with PDF/CDF
- **Multiple testing correction** -- Bonferroni, Benjamini-Hochberg FDR
- **Effect sizes** -- Cohen's d, eta-squared, odds ratio, relative risk
- **PCA** -- principal component analysis (with optional BLAS acceleration)
- **Bayesian conjugate priors** -- Beta, Gamma, NormalConjugate, Dirichlet
- **Combinatorics** -- factorial, binomial coefficients, permutations, multinomial (exact and log-space)
- **Population genetics** -- allele frequencies, HWE, Fst (Weir-Cockerham/Hudson), nucleotide diversity, Tajima's D, linkage disequilibrium, genotype PCA
- **Expression normalization** -- TPM, FPKM, CPM, DESeq2 size factors
- **Differential expression** -- negative binomial Wald test, Wilcoxon, volcano plot data
- **Enrichment analysis** -- ORA (hypergeometric), GSEA (preranked with permutations), GO enrichment
- **Ecological diversity** -- Shannon, Simpson, Chao1, Bray-Curtis, rarefaction curves
- **Survival analysis** -- Kaplan-Meier, log-rank test, Cox proportional hazards

## Quick Start

```toml
[dependencies]
cyanea-stats = "0.1"
```

```rust
use cyanea_stats::{describe, pearson, t_test_two_sample};

let stats = describe(&[1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
println!("Mean: {}, Median: {}", stats.mean, stats.median);

let r = pearson(&[1.0, 2.0, 3.0], &[2.0, 4.0, 6.0]).unwrap();
println!("Pearson r: {}", r);

let test = t_test_two_sample(&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0], false).unwrap();
println!("p = {}", test.p_value);
```

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support |
| `wasm` | No | WASM target marker |
| `serde` | No | Serialize/Deserialize derives |
| `parallel` | No | Rayon parallelism |
| `blas` | No | BLAS-accelerated PCA via ndarray |

## Modules

| Module | Description |
|--------|-------------|
| `descriptive` | Mean, median, variance, quantiles, IQR, MAD |
| `correlation` | Pearson, Spearman, correlation matrices |
| `rank` | Rank computation with tie-breaking methods |
| `testing` | t-tests, Mann-Whitney U, Fisher's exact, chi-squared, ANOVA |
| `distribution` | Normal, Poisson, Binomial, ChiSquared, F, erf, ln_gamma, betai |
| `correction` | Bonferroni, Benjamini-Hochberg |
| `effect_size` | Cohen's d, eta-squared, odds ratio, relative risk |
| `reduction` | PCA dimensionality reduction |
| `bayesian` | Beta, Gamma, NormalConjugate, Dirichlet conjugate priors |
| `combinatorics` | Factorial, binomial coefficients, permutations, combinations |
| `popgen` | Allele frequencies, HWE, Fst, diversity, Tajima's D, LD |
| `normalization` | TPM, FPKM, CPM, DESeq2 size factors |
| `diffexpr` | Differential expression testing |
| `enrichment` | ORA, GSEA, GO enrichment |
| `diversity` | Alpha/beta diversity, rarefaction |
| `survival` | Kaplan-Meier, log-rank, Cox PH |

## See Also

- [API Reference (STATUS.md)](docs/STATUS.md)
- [Architecture](../ARCHITECTURE.md)
- [Build Guide](../docs/BUILDING.md)
