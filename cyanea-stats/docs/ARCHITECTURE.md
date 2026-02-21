# cyanea-stats Architecture

Internal design of the statistical analysis modules.

## Distribution Trait Hierarchy

The `Distribution` trait provides a common interface for all probability distributions:

```
trait Distribution {
    fn pdf(&self, x: f64) -> f64;
    fn cdf(&self, x: f64) -> f64;
}
```

Implemented by: `Normal`, `Poisson`, `Binomial`, `ChiSquared`, `FDistribution`, `NegativeBinomial`, `Beta` (in bayesian module), `Gamma` (in bayesian module).

CDF implementations use numerical approximations:
- Normal: error function (`erf`) via Abramowitz and Stegun polynomial approximation
- Chi-squared, F, Binomial: regularized incomplete beta function (`betai`) via continued fraction expansion
- Poisson: regularized incomplete gamma function
- Log-gamma: Lanczos approximation with 7 coefficients

## Bayesian Conjugate Priors

Four conjugate families, each with analytical posterior update:

| Prior | Likelihood | Update Rule |
|-------|-----------|-------------|
| `Beta(a, b)` | Binomial(n, p) | `Beta(a + successes, b + failures)` |
| `Gamma(shape, rate)` | Poisson(lambda) | `Gamma(shape + sum(counts), rate + n)` |
| `NormalConjugate(mu, var)` | Normal(mu, sigma^2_known) | Precision-weighted combination |
| `Dirichlet(alpha)` | Multinomial | `Dirichlet(alpha_i + count_i)` |

All updates are O(1) per observation (or O(k) for Dirichlet with k categories). The `Beta` and `Gamma` types implement the `Distribution` trait for PDF/CDF queries on the posterior.

## Population Genetics

Genotype data is stored as a flat `&[u8]` array in locus-major order, with values 0/1/2 representing the number of alternate alleles. This encoding is standard for biallelic SNP data (PLINK format).

- **Allele frequencies**: O(n) per locus, counting alleles directly.
- **HWE test**: Compares observed genotype counts to expected under Hardy-Weinberg equilibrium using a chi-squared test.
- **Fst**: Weir-Cockerham estimator uses variance components (among populations, within populations) computed per locus, then averaged. Hudson estimator uses the ratio of between-population to total nucleotide diversity.
- **Tajima's D**: Computes S (segregating sites) and pi (pairwise nucleotide differences), then forms D = (pi - theta_W) / sqrt(Var), where variance components use Tajima's exact formulae based on sample size.
- **LD (r-squared, D')**: Uses EM algorithm to estimate two-locus haplotype frequencies from unphased genotypes, then computes D = p_AB - p_A * p_B.
- **Genotype PCA**: Centers the genotype matrix by per-locus mean, then runs PCA (power iteration) on the covariance matrix.

## PERMANOVA

The PERMANOVA framework partitions the total sum of squared distances into among-group and within-group components:

1. Compute the pseudo-F statistic: F = (SS_among / (k-1)) / (SS_within / (n-k))
2. Generate a null distribution by permuting group labels `n_permutations` times
3. p-value = proportion of permuted F values >= observed F

The distance matrix is provided externally (typically Bray-Curtis from the diversity module). Group labels are 0-indexed integers.

ANOSIM, Mantel test, and AMOVA follow the same permutation framework with different test statistics:
- ANOSIM: R statistic based on rank dissimilarities
- Mantel: Pearson or Spearman correlation between two distance matrices
- AMOVA: Phi_ST from hierarchical variance partitioning

## Ordination

- **PCoA**: Double-centers the squared distance matrix to form the Gower matrix G, then extracts eigenvalues/eigenvectors via power iteration with deflation. Coordinates = eigenvector * sqrt(max(eigenvalue, 0)).
- **NMDS**: Iterative stress minimization (Kruskal's stress-1). Starts from PCoA initialization, then gradient descent to minimize the mismatch between distances in the embedding and the rank-order of the original distances.
- **RDA**: Constrained PCA. Regresses the response matrix on the predictor matrix, then runs PCA on the fitted values.
- **CCA**: Chi-squared-weighted RDA on the species matrix.
- **Procrustes**: Optimal rotation/scaling via SVD of the cross-product matrix, with permutation-based M-squared significance test.

## Survival Analysis

### Kaplan-Meier

The product-limit estimator processes events in time order. At each event time t_i:
- S(t_i) = S(t_{i-1}) * (1 - d_i / n_i)
- where d_i = events and n_i = at-risk count

Standard errors use Greenwood's formula. 95% confidence intervals use the log-transform method: CI = exp(log(S) +/- 1.96 * SE / S), which guarantees intervals in [0, 1].

### Log-rank test

Compares observed vs expected events across groups at each distinct event time. Uses the chi-squared approximation with k-1 degrees of freedom.

### Cox Proportional Hazards

Fits the partial likelihood using Newton-Raphson iteration:
1. Compute log-partial likelihood and its first/second derivatives (score and Hessian)
2. Update coefficients: beta_new = beta_old + Hessian^{-1} * score
3. Iterate until convergence (gradient norm < threshold or max iterations)

Handles ties using the Breslow approximation. Convergence typically requires 5-15 iterations.

## GSEA

Gene Set Enrichment Analysis (preranked) computes a running sum statistic:
1. Sort genes by score (decreasing)
2. Walk down the list; increment the running sum when hitting a gene in the set (weighted by |score|^p), decrement when missing
3. Enrichment score = maximum deviation from zero
4. Normalized ES = ES / mean(|ES_null|) from gene-set-size-matched permutations
5. P-value from the permutation null distribution, BH-corrected across gene sets

## Null Models

- **Permutation null**: Shuffles group labels (Fisher-Yates) and recomputes the user-supplied statistic function for each permutation.
- **Bootstrap null**: Resamples with replacement from the data, computes the statistic on each resample.
- **Wright-Fisher**: Simulates allele frequency drift over discrete generations using binomial sampling: each generation draws from Binomial(2N, p) where p is the current frequency.
