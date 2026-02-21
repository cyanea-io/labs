# cyanea-stats Usage Guide

Practical examples for statistical analysis in life sciences: descriptive stats, hypothesis testing, population genetics, survival analysis, enrichment, ordination, and multivariate testing.

## Descriptive Statistics

```rust
use cyanea_stats::{describe, mean, median, quantile, iqr, mad};

let data = [2.3, 4.1, 3.7, 5.2, 1.9, 4.8, 3.3];

let stats = describe(&data).unwrap();
println!("Mean: {:.2}", stats.mean);
println!("Median: {:.2}", stats.median);
println!("Std dev: {:.2}", stats.sample_std_dev);
println!("Skewness: {:.2}", stats.skewness);
println!("Kurtosis: {:.2}", stats.kurtosis);

// Individual functions
let q75 = quantile(&data, 0.75).unwrap();
let range = iqr(&data).unwrap();
let robust_spread = mad(&data).unwrap();
```

## Hypothesis Testing

```rust
use cyanea_stats::testing::*;

// Two-sample t-test (Welch's)
let control = [5.1, 4.8, 5.3, 4.9, 5.0];
let treatment = [6.2, 5.9, 6.5, 6.1, 6.3];
let result = t_test_two_sample(&control, &treatment, false).unwrap();
println!("t = {:.3}, p = {:.4}", result.statistic, result.p_value);

// Mann-Whitney U (non-parametric)
let u_result = mann_whitney_u(&control, &treatment).unwrap();
println!("U statistic = {:.1}, p = {:.4}", u_result.statistic, u_result.p_value);

// Fisher's exact test on 2x2 table
let table = [[10, 5], [3, 12]];
let fisher = fisher_exact(&table).unwrap();
println!("p = {:.4}", fisher.p_value);

// One-way ANOVA
let groups = vec![
    vec![1.0, 2.0, 3.0],
    vec![4.0, 5.0, 6.0],
    vec![7.0, 8.0, 9.0],
];
let anova_result = anova(&groups).unwrap();
```

## Correlation

```rust
use cyanea_stats::correlation::*;

let x = [1.0, 2.0, 3.0, 4.0, 5.0];
let y = [2.1, 3.9, 6.2, 7.8, 10.1];

// Pearson (linear)
let r = pearson(&x, &y).unwrap();
println!("Pearson r = {:.3}", r);

// Spearman (rank-based)
let rho = spearman(&x, &y).unwrap();
println!("Spearman rho = {:.3}", rho);

// Pairwise correlation matrix
let matrix_data = vec![x.to_vec(), y.to_vec()];
let cor_matrix = CorrelationMatrix::pearson(&matrix_data).unwrap();
```

## PCA

```rust
use cyanea_stats::reduction::{pca, PcaResult};

// 5 samples, 3 features each (row-major flat array)
let data = vec![
    1.0, 2.0, 3.0,
    4.0, 5.0, 6.0,
    7.0, 8.0, 9.0,
    2.0, 3.0, 4.0,
    5.0, 6.0, 7.0,
];
let result = pca(&data, 2).unwrap();
println!("Variance explained: {:?}", result.explained_variance_ratio);
```

## Population Genetics

```rust
use cyanea_stats::popgen::*;

// Genotype matrix: 5 loci x 10 individuals, values 0/1/2
let genotypes: Vec<u8> = vec![
    0, 1, 2, 0, 1, 1, 0, 2, 1, 0,  // locus 1
    1, 1, 0, 0, 2, 1, 1, 0, 0, 1,  // locus 2
    0, 0, 1, 0, 0, 1, 0, 0, 0, 1,  // locus 3
    2, 1, 1, 2, 1, 0, 1, 1, 2, 1,  // locus 4
    0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  // locus 5
];
let n_samples = 10;

// Allele frequencies
let freqs = allele_frequencies(&genotypes, n_samples).unwrap();
for f in &freqs {
    println!("MAF = {:.3}, Het_obs = {:.3}", f.minor_allele_freq, f.observed_het);
}

// Hardy-Weinberg equilibrium test
let hwe = hwe_test(&genotypes, n_samples).unwrap();

// Fst between two populations
let populations = vec![0, 0, 0, 0, 0, 1, 1, 1, 1, 1]; // 5 per pop
let fst_result = fst(&genotypes, n_samples, &populations, FstMethod::WeirCockerham).unwrap();
println!("Fst = {:.4}", fst_result.fst);

// Tajima's D
let tajima = tajimas_d(&genotypes, n_samples, 1000).unwrap();
println!("Tajima's D = {:.3}", tajima.d);

// Linkage disequilibrium
let ld_result = ld(&genotypes, n_samples, 0, 1).unwrap();
println!("r^2 = {:.3}, D' = {:.3}", ld_result.r_squared, ld_result.d_prime);

// Genotype PCA for population structure
let pca = genotype_pca(&genotypes, n_samples, 2).unwrap();
```

## Survival Analysis

```rust
use cyanea_stats::survival::*;

// Kaplan-Meier survival curve
let times = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
let status = vec![true, false, true, true, false, true, false, true]; // true = event

let km = kaplan_meier(&times, &status).unwrap();
println!("Median survival: {:?}", km.median_survival);
for step in &km.steps {
    println!("t={:.1}: S={:.3} ({:.3}-{:.3})",
        step.time, step.survival, step.ci_lower, step.ci_upper);
}

// Log-rank test comparing two groups
let groups = vec![0, 0, 0, 0, 1, 1, 1, 1];
let lr = log_rank_test(&times, &status, &groups).unwrap();
println!("Log-rank chi2 = {:.3}, p = {:.4}", lr.chi_squared, lr.p_value);

// Cox proportional hazards
let covariates = vec![
    0.5, 1.2,  // patient 1: age_scaled, treatment
    0.3, 0.0,  // patient 2
    0.8, 1.0,  // ...
    0.6, 0.0,
    0.4, 1.0,
    0.9, 0.0,
    0.7, 1.0,
    0.2, 0.0,
];
let cox = cox_ph(&times, &status, &covariates, 2).unwrap();
for i in 0..cox.coefficients.len() {
    println!("beta={:.3}, HR={:.3}, p={:.4}",
        cox.coefficients[i], cox.hazard_ratios[i], cox.p_values[i]);
}
```

## Enrichment Analysis

```rust
use cyanea_stats::enrichment::*;

// Over-representation analysis (ORA)
let significant_genes = vec![0, 5, 12, 23, 45]; // gene indices
let gene_sets = vec![
    GeneSet { name: "pathway_A".into(), genes: vec![0, 5, 10, 15, 20] },
    GeneSet { name: "pathway_B".into(), genes: vec![12, 23, 34, 45, 56] },
];
let ora_results = ora(&significant_genes, &gene_sets, 1000).unwrap();
for r in &ora_results {
    println!("{}: overlap={}, p_adj={:.4}", r.gene_set, r.overlap, r.p_adjusted);
}

// GSEA preranked
let all_genes: Vec<usize> = (0..1000).collect();
let scores: Vec<f64> = (0..1000).map(|i| 1.0 - (i as f64 / 500.0)).collect();
let gsea_results = gsea_preranked(&all_genes, &scores, &gene_sets, 1.0, 1000).unwrap();
for r in &gsea_results {
    println!("{}: NES={:.3}, p_adj={:.4}", r.gene_set, r.normalized_es, r.p_adjusted);
}
```

## Microbiome Ordination

```rust
use cyanea_stats::diversity::*;
use cyanea_stats::ordination::*;

// Alpha diversity
let sample_counts = vec![100, 50, 30, 20, 10, 5, 3, 1, 1];
let alpha = alpha_diversity(&sample_counts).unwrap();
println!("Shannon: {:.3}, Chao1: {:.1}", alpha.shannon, alpha.chao1);

// Beta diversity: Bray-Curtis distance matrix
let samples = vec![
    vec![100, 50, 30, 20],
    vec![80, 60, 40, 10],
    vec![10, 5, 90, 80],
    vec![15, 8, 85, 70],
];
let bc_matrix = bray_curtis_matrix(&samples).unwrap();

// PCoA ordination
let pcoa_result = pcoa(&bc_matrix, 2).unwrap();
println!("Variance explained: {:?}", pcoa_result.proportion_explained);

// NMDS (alternative to PCoA)
let nmds_result = nmds(&bc_matrix, 2, Default::default()).unwrap();
println!("Stress: {:.4}", nmds_result.stress);
```

## Multivariate Tests

```rust
use cyanea_stats::multivariate::*;

// PERMANOVA: test if microbial communities differ by treatment
let distances = vec![/* pairwise Bray-Curtis matrix */];
let groups = vec![0, 0, 0, 0, 1, 1, 1, 1]; // treatment groups

let result = permanova(&distances, &groups, 999, 42).unwrap();
println!("F = {:.3}, R^2 = {:.3}, p = {:.3}",
    result.f_statistic, result.r_squared, result.p_value);

// ANOSIM
let anosim_result = anosim(&distances, &groups, 999, 42).unwrap();
println!("R = {:.3}, p = {:.3}", anosim_result.r_statistic, anosim_result.p_value);

// Mantel test: correlation between two distance matrices
let env_distances = vec![/* environmental distance matrix */];
let mantel_result = mantel(&distances, &env_distances, 999, 42).unwrap();
println!("r = {:.3}, p = {:.3}", mantel_result.r_statistic, mantel_result.p_value);
```

## Bayesian Conjugate Priors

```rust
use cyanea_stats::bayesian::*;

// Beta-Binomial: estimate mutation rate
let prior = Beta::new(1.0, 1.0).unwrap(); // uniform prior
let posterior = prior.update_binomial(15, 100); // 15 mutations in 100 sites
println!("Posterior mean: {:.3}", posterior.mean());

// Gamma-Poisson: estimate expression count rate
let prior = Gamma::new(2.0, 1.0).unwrap();
let posterior = prior.update_poisson_batch(&[5, 8, 3, 7, 6]);

// Dirichlet-Multinomial: nucleotide composition
let prior = Dirichlet::symmetric(4, 1.0).unwrap(); // 4 bases, uniform
let counts = vec![30.0, 20.0, 25.0, 25.0]; // A, C, G, T counts
let posterior = prior.update_multinomial(&counts);
println!("Expected freqs: {:?}", posterior.mean());
```
