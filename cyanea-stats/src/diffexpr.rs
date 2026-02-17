//! Differential expression analysis for count data.
//!
//! Provides two methods for identifying differentially expressed genes between
//! two conditions:
//!
//! - **Negative binomial Wald test** ([`DeMethod::NegativeBinomial`]) — models
//!   count overdispersion (DESeq2-style pipeline with size-factor normalization,
//!   method-of-moments dispersion, and Wald z-test).
//! - **Wilcoxon rank-sum test** ([`DeMethod::Wilcoxon`]) — non-parametric
//!   alternative using the existing [`crate::testing::mann_whitney_u`].
//!
//! Both methods apply Benjamini-Hochberg correction via
//! [`crate::correction::benjamini_hochberg`].
//!
//! The [`volcano_plot`] function converts results into points suitable for
//! plotting (log2 fold-change vs. −log10 adjusted p-value).

use cyanea_core::{CyaneaError, Result};

use crate::correction;
use crate::distribution::{Distribution, Normal};
use crate::normalization;
use crate::testing;

// ── Result types ─────────────────────────────────────────────────────────────

/// Which statistical method to use for differential expression.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DeMethod {
    /// Negative binomial Wald test (DESeq2-style).
    NegativeBinomial,
    /// Wilcoxon rank-sum (Mann-Whitney U) test.
    Wilcoxon,
}

/// Per-gene differential expression result.
#[derive(Debug, Clone)]
pub struct DeGeneResult {
    /// Index of the gene in the input matrix.
    pub gene_index: usize,
    /// Log2 fold-change (condition / control).
    pub log2_fold_change: f64,
    /// Mean of normalized counts across all samples.
    pub base_mean: f64,
    /// Test statistic (Wald z for NB, U for Wilcoxon).
    pub statistic: f64,
    /// Raw p-value.
    pub p_value: f64,
    /// Benjamini-Hochberg adjusted p-value.
    pub p_adjusted: f64,
}

/// Aggregate results from a differential expression analysis.
#[derive(Debug, Clone)]
pub struct DeResults {
    /// Per-gene results, sorted by p_value ascending.
    pub genes: Vec<DeGeneResult>,
    /// Method used.
    pub method: DeMethod,
    /// Number of genes tested.
    pub n_genes: usize,
    /// Number of condition (treatment) samples.
    pub n_condition: usize,
    /// Number of control samples.
    pub n_control: usize,
}

/// A point for a volcano plot.
#[derive(Debug, Clone)]
pub struct VolcanoPoint {
    /// Gene index from the original matrix.
    pub gene_index: usize,
    /// Log2 fold-change.
    pub log2_fold_change: f64,
    /// −log10(adjusted p-value), clamped to 300.
    pub neg_log10_padj: f64,
    /// Whether this gene passes significance thresholds.
    pub significant: bool,
}

// ── Main entry point ─────────────────────────────────────────────────────────

/// Run differential expression analysis on a count matrix.
///
/// - `counts`: row-major `n_genes × n_samples` count matrix.
/// - `condition`: boolean mask (`true` = treatment, `false` = control), one
///   per sample.
/// - `method`: which test to apply.
///
/// Returns [`DeResults`] with genes sorted by ascending p-value.
pub fn differential_expression(
    counts: &[f64],
    n_genes: usize,
    n_samples: usize,
    condition: &[bool],
    method: DeMethod,
) -> Result<DeResults> {
    // Validate dimensions
    if n_genes == 0 || n_samples == 0 {
        return Err(CyaneaError::InvalidInput(
            "differential_expression: need at least 1 gene and 1 sample".into(),
        ));
    }
    if counts.len() != n_genes * n_samples {
        return Err(CyaneaError::InvalidInput(format!(
            "differential_expression: counts length ({}) != n_genes ({}) * n_samples ({})",
            counts.len(),
            n_genes,
            n_samples,
        )));
    }
    if condition.len() != n_samples {
        return Err(CyaneaError::InvalidInput(format!(
            "differential_expression: condition length ({}) != n_samples ({})",
            condition.len(),
            n_samples,
        )));
    }

    let n_cond = condition.iter().filter(|&&c| c).count();
    let n_ctrl = n_samples - n_cond;
    if n_cond < 2 || n_ctrl < 2 {
        return Err(CyaneaError::InvalidInput(
            "differential_expression: need at least 2 samples per group".into(),
        ));
    }

    // Build sample index lists
    let cond_idx: Vec<usize> = (0..n_samples).filter(|&j| condition[j]).collect();
    let ctrl_idx: Vec<usize> = (0..n_samples).filter(|&j| !condition[j]).collect();

    // Size-factor normalization
    let sf = normalization::size_factors(counts, n_genes, n_samples)?;
    let normed = normalization::normalize_by_size_factors(counts, n_genes, n_samples, &sf)?;

    let mut gene_results: Vec<DeGeneResult> = match method {
        DeMethod::NegativeBinomial => nb_wald(&normed, n_genes, n_samples, &cond_idx, &ctrl_idx)?,
        DeMethod::Wilcoxon => wilcoxon_de(&normed, n_genes, n_samples, &cond_idx, &ctrl_idx)?,
    };

    // BH correction
    let raw_p: Vec<f64> = gene_results.iter().map(|g| g.p_value).collect();
    let adj_p = correction::benjamini_hochberg(&raw_p)?;
    for (g, &padj) in gene_results.iter_mut().zip(adj_p.iter()) {
        g.p_adjusted = padj;
    }

    // Sort by p-value ascending
    gene_results.sort_by(|a, b| a.p_value.total_cmp(&b.p_value));

    Ok(DeResults {
        genes: gene_results,
        method,
        n_genes,
        n_condition: n_cond,
        n_control: n_ctrl,
    })
}

// ── Negative binomial Wald test ──────────────────────────────────────────────

fn nb_wald(
    normed: &[f64],
    n_genes: usize,
    n_samples: usize,
    cond_idx: &[usize],
    ctrl_idx: &[usize],
) -> Result<Vec<DeGeneResult>> {
    let normal = Normal::standard();
    let pseudo = 0.5;

    let mut results = Vec::with_capacity(n_genes);

    for i in 0..n_genes {
        let row = &normed[i * n_samples..(i + 1) * n_samples];

        // Group means
        let mu_cond: f64 = cond_idx.iter().map(|&j| row[j]).sum::<f64>() / cond_idx.len() as f64;
        let mu_ctrl: f64 = ctrl_idx.iter().map(|&j| row[j]).sum::<f64>() / ctrl_idx.len() as f64;
        let base_mean: f64 = row.iter().sum::<f64>() / n_samples as f64;

        // log2 fold-change with pseudocount
        let log2fc = ((mu_cond + pseudo) / (mu_ctrl + pseudo)).log2();

        // Overall mean and variance for dispersion estimation
        let overall_mean = base_mean;
        let overall_var = if n_samples > 1 {
            let ss: f64 = row.iter().map(|&x| (x - overall_mean).powi(2)).sum();
            ss / (n_samples - 1) as f64
        } else {
            0.0
        };

        // Method-of-moments dispersion: alpha = (var - mean) / mean^2
        let alpha = if overall_mean > 0.0 {
            ((overall_var - overall_mean) / (overall_mean * overall_mean)).clamp(1e-8, 1e8)
        } else {
            1e-8
        };

        // Standard error via delta method on NB variance
        // Var(X) = mu + alpha * mu^2 for NB
        // SE of group mean = sqrt(Var / n), SE of log2FC ≈ sqrt(SE_cond^2 + SE_ctrl^2) / ln(2)
        let var_cond = mu_cond + alpha * mu_cond * mu_cond;
        let var_ctrl = mu_ctrl + alpha * mu_ctrl * mu_ctrl;
        let se_cond = (var_cond / cond_idx.len() as f64).sqrt();
        let se_ctrl = (var_ctrl / ctrl_idx.len() as f64).sqrt();
        // log2FC = log2((mu_c + pc) / (mu_t + pc)), SE on log2 scale via delta method
        let se_log2fc = ((se_cond / (mu_cond + pseudo)).powi(2)
            + (se_ctrl / (mu_ctrl + pseudo)).powi(2))
        .sqrt()
            / 2.0_f64.ln();

        // Wald z-statistic and two-tailed p-value
        let (z, p_value) = if se_log2fc > 1e-15 {
            let z = log2fc / se_log2fc;
            let p = 2.0 * (1.0 - normal.cdf(z.abs()));
            (z, p.min(1.0))
        } else {
            (0.0, 1.0)
        };

        results.push(DeGeneResult {
            gene_index: i,
            log2_fold_change: log2fc,
            base_mean,
            statistic: z,
            p_value,
            p_adjusted: 1.0, // filled in later
        });
    }

    Ok(results)
}

// ── Wilcoxon rank-sum DE ─────────────────────────────────────────────────────

fn wilcoxon_de(
    normed: &[f64],
    n_genes: usize,
    n_samples: usize,
    cond_idx: &[usize],
    ctrl_idx: &[usize],
) -> Result<Vec<DeGeneResult>> {
    let pseudo = 0.5;
    let mut results = Vec::with_capacity(n_genes);

    for i in 0..n_genes {
        let row = &normed[i * n_samples..(i + 1) * n_samples];

        let cond_vals: Vec<f64> = cond_idx.iter().map(|&j| row[j]).collect();
        let ctrl_vals: Vec<f64> = ctrl_idx.iter().map(|&j| row[j]).collect();

        let mu_cond = cond_vals.iter().sum::<f64>() / cond_vals.len() as f64;
        let mu_ctrl = ctrl_vals.iter().sum::<f64>() / ctrl_vals.len() as f64;
        let base_mean: f64 = row.iter().sum::<f64>() / n_samples as f64;
        let log2fc = ((mu_cond + pseudo) / (mu_ctrl + pseudo)).log2();

        let test_result = testing::mann_whitney_u(&cond_vals, &ctrl_vals)?;

        results.push(DeGeneResult {
            gene_index: i,
            log2_fold_change: log2fc,
            base_mean,
            statistic: test_result.statistic,
            p_value: test_result.p_value,
            p_adjusted: 1.0,
        });
    }

    Ok(results)
}

// ── Volcano plot ─────────────────────────────────────────────────────────────

/// Convert DE results into volcano-plot points.
///
/// - `padj_threshold`: adjusted p-value cutoff (e.g. 0.05).
/// - `fc_threshold`: absolute log2 fold-change cutoff (e.g. 1.0).
///
/// A gene is marked `significant` if `p_adjusted < padj_threshold` **and**
/// `|log2_fold_change| > fc_threshold`.
pub fn volcano_plot(
    results: &DeResults,
    padj_threshold: f64,
    fc_threshold: f64,
) -> Vec<VolcanoPoint> {
    results
        .genes
        .iter()
        .map(|g| {
            let neg_log10 = if g.p_adjusted > 0.0 {
                (-g.p_adjusted.log10()).min(300.0)
            } else {
                300.0
            };
            VolcanoPoint {
                gene_index: g.gene_index,
                log2_fold_change: g.log2_fold_change,
                neg_log10_padj: neg_log10,
                significant: g.p_adjusted < padj_threshold
                    && g.log2_fold_change.abs() > fc_threshold,
            }
        })
        .collect()
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: build a simple count matrix with one strongly upregulated gene,
    /// one downregulated gene, and several unchanged genes.
    ///
    /// Layout: 5 genes × 6 samples (3 ctrl + 3 cond)
    fn test_counts() -> (Vec<f64>, usize, usize, Vec<bool>) {
        let n_genes = 5;
        let n_samples = 6;
        // condition: first 3 control, last 3 treatment
        let condition = vec![false, false, false, true, true, true];

        #[rustfmt::skip]
        let counts = vec![
            // gene 0: upregulated in treatment (ctrl ~10, cond ~200)
            10.0, 12.0, 11.0, 200.0, 210.0, 190.0,
            // gene 1: downregulated (ctrl ~200, cond ~10)
            200.0, 190.0, 210.0, 10.0, 12.0, 11.0,
            // gene 2: unchanged (~100)
            100.0, 105.0, 95.0, 98.0, 102.0, 100.0,
            // gene 3: unchanged (~50)
            50.0, 52.0, 48.0, 49.0, 51.0, 50.0,
            // gene 4: unchanged (~75)
            75.0, 78.0, 72.0, 74.0, 76.0, 75.0,
        ];
        (counts, n_genes, n_samples, condition)
    }

    #[test]
    fn nb_detects_upregulated() {
        let (counts, ng, ns, cond) = test_counts();
        let res = differential_expression(&counts, ng, ns, &cond, DeMethod::NegativeBinomial).unwrap();
        let gene0 = res.genes.iter().find(|g| g.gene_index == 0).unwrap();
        assert!(gene0.log2_fold_change > 2.0, "log2fc={}", gene0.log2_fold_change);
        assert!(gene0.p_adjusted < 0.05, "padj={}", gene0.p_adjusted);
    }

    #[test]
    fn nb_detects_downregulated() {
        let (counts, ng, ns, cond) = test_counts();
        let res = differential_expression(&counts, ng, ns, &cond, DeMethod::NegativeBinomial).unwrap();
        let gene1 = res.genes.iter().find(|g| g.gene_index == 1).unwrap();
        assert!(gene1.log2_fold_change < -2.0, "log2fc={}", gene1.log2_fold_change);
        assert!(gene1.p_adjusted < 0.05, "padj={}", gene1.p_adjusted);
    }

    #[test]
    fn nb_unchanged_genes_high_p() {
        let (counts, ng, ns, cond) = test_counts();
        let res = differential_expression(&counts, ng, ns, &cond, DeMethod::NegativeBinomial).unwrap();
        for idx in [2, 3, 4] {
            let gene = res.genes.iter().find(|g| g.gene_index == idx).unwrap();
            assert!(
                gene.p_value > 0.05,
                "gene {idx} should not be significant: p={}",
                gene.p_value
            );
        }
    }

    #[test]
    fn nb_log2fc_direction() {
        let (counts, ng, ns, cond) = test_counts();
        let res = differential_expression(&counts, ng, ns, &cond, DeMethod::NegativeBinomial).unwrap();
        let gene0 = res.genes.iter().find(|g| g.gene_index == 0).unwrap();
        let gene1 = res.genes.iter().find(|g| g.gene_index == 1).unwrap();
        assert!(gene0.log2_fold_change > 0.0);
        assert!(gene1.log2_fold_change < 0.0);
    }

    #[test]
    fn nb_padj_ge_pvalue() {
        let (counts, ng, ns, cond) = test_counts();
        let res = differential_expression(&counts, ng, ns, &cond, DeMethod::NegativeBinomial).unwrap();
        for g in &res.genes {
            assert!(
                g.p_adjusted >= g.p_value - 1e-15,
                "gene {}: padj={} < p={}",
                g.gene_index,
                g.p_adjusted,
                g.p_value
            );
        }
    }

    #[test]
    fn nb_results_sorted() {
        let (counts, ng, ns, cond) = test_counts();
        let res = differential_expression(&counts, ng, ns, &cond, DeMethod::NegativeBinomial).unwrap();
        for w in res.genes.windows(2) {
            assert!(
                w[0].p_value <= w[1].p_value + 1e-15,
                "not sorted: {} > {}",
                w[0].p_value,
                w[1].p_value
            );
        }
    }

    #[test]
    fn wilcoxon_detects_de() {
        let (counts, ng, ns, cond) = test_counts();
        let res = differential_expression(&counts, ng, ns, &cond, DeMethod::Wilcoxon).unwrap();
        let gene0 = res.genes.iter().find(|g| g.gene_index == 0).unwrap();
        assert!(gene0.log2_fold_change > 2.0);
        assert!(gene0.p_value < 0.1, "p={}", gene0.p_value);
    }

    #[test]
    fn wilcoxon_matches_direct_mwu() {
        // Verify that the Wilcoxon pathway produces the same p-value as
        // calling mann_whitney_u directly on the same normalized data.
        let (counts, ng, ns, cond) = test_counts();
        let sf = normalization::size_factors(&counts, ng, ns).unwrap();
        let normed = normalization::normalize_by_size_factors(&counts, ng, ns, &sf).unwrap();

        let res = differential_expression(&counts, ng, ns, &cond, DeMethod::Wilcoxon).unwrap();

        let cond_idx: Vec<usize> = (0..ns).filter(|&j| cond[j]).collect();
        let ctrl_idx: Vec<usize> = (0..ns).filter(|&j| !cond[j]).collect();

        for gene_res in &res.genes {
            let i = gene_res.gene_index;
            let row = &normed[i * ns..(i + 1) * ns];
            let cond_vals: Vec<f64> = cond_idx.iter().map(|&j| row[j]).collect();
            let ctrl_vals: Vec<f64> = ctrl_idx.iter().map(|&j| row[j]).collect();
            let direct = testing::mann_whitney_u(&cond_vals, &ctrl_vals).unwrap();
            assert!(
                (gene_res.p_value - direct.p_value).abs() < 1e-10,
                "gene {}: de_p={}, direct_p={}",
                i,
                gene_res.p_value,
                direct.p_value
            );
        }
    }

    #[test]
    fn dispersion_poisson_like() {
        // When data follows Poisson (variance ≈ mean), dispersion should be small
        let counts = vec![
            100.0, 101.0, 99.0, 100.0, 102.0, 98.0,
        ];
        let cond = vec![false, false, false, true, true, true];
        let res = differential_expression(&counts, 1, 6, &cond, DeMethod::NegativeBinomial).unwrap();
        // With nearly identical groups, p should be large
        assert!(res.genes[0].p_value > 0.5, "p={}", res.genes[0].p_value);
    }

    #[test]
    fn dispersion_overdispersed() {
        // Highly variable data should still work
        #[rustfmt::skip]
        let counts = vec![
            1.0, 50.0, 200.0, 500.0, 1000.0, 2000.0,
        ];
        let cond = vec![false, false, false, true, true, true];
        let res = differential_expression(&counts, 1, 6, &cond, DeMethod::NegativeBinomial);
        assert!(res.is_ok());
    }

    #[test]
    fn volcano_thresholds() {
        let (counts, ng, ns, cond) = test_counts();
        let res = differential_expression(&counts, ng, ns, &cond, DeMethod::NegativeBinomial).unwrap();
        let points = volcano_plot(&res, 0.05, 1.0);

        assert_eq!(points.len(), ng);
        // Gene 0 and 1 should be significant (large FC, low padj)
        let sig_genes: Vec<usize> = points.iter().filter(|p| p.significant).map(|p| p.gene_index).collect();
        assert!(sig_genes.contains(&0), "gene 0 should be significant");
        assert!(sig_genes.contains(&1), "gene 1 should be significant");

        // Unchanged genes should not be significant
        for idx in [2, 3, 4] {
            let pt = points.iter().find(|p| p.gene_index == idx).unwrap();
            assert!(!pt.significant, "gene {idx} should not be significant");
        }

        // neg_log10_padj should be non-negative
        for pt in &points {
            assert!(pt.neg_log10_padj >= 0.0);
        }
    }

    #[test]
    fn error_dimension_mismatch() {
        let cond = vec![false, true, false, true];
        assert!(differential_expression(&[1.0, 2.0], 2, 4, &cond, DeMethod::NegativeBinomial).is_err());
    }

    #[test]
    fn error_condition_length() {
        let counts = vec![1.0; 8];
        let cond = vec![false, true]; // too short
        assert!(differential_expression(&counts, 2, 4, &cond, DeMethod::NegativeBinomial).is_err());
    }

    #[test]
    fn error_too_few_per_group() {
        let counts = vec![10.0, 20.0, 30.0, 40.0];
        // Only 1 control
        let cond = vec![false, true, true, true];
        assert!(differential_expression(&counts, 1, 4, &cond, DeMethod::NegativeBinomial).is_err());
    }

    #[test]
    fn error_single_group() {
        let counts = vec![10.0, 20.0, 30.0, 40.0];
        let cond = vec![true, true, true, true];
        assert!(differential_expression(&counts, 1, 4, &cond, DeMethod::NegativeBinomial).is_err());
    }

    #[test]
    fn volcano_clamps_neg_log10() {
        // Create a result with p_adjusted = 0 to test clamping
        let results = DeResults {
            genes: vec![DeGeneResult {
                gene_index: 0,
                log2_fold_change: 5.0,
                base_mean: 100.0,
                statistic: 10.0,
                p_value: 0.0,
                p_adjusted: 0.0,
            }],
            method: DeMethod::NegativeBinomial,
            n_genes: 1,
            n_condition: 3,
            n_control: 3,
        };
        let points = volcano_plot(&results, 0.05, 1.0);
        assert_eq!(points[0].neg_log10_padj, 300.0);
        assert!(points[0].significant);
    }
}
