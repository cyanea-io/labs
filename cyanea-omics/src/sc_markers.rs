//! Single-cell marker gene identification: differential expression and filtering.

use std::collections::HashMap;

use cyanea_core::{CyaneaError, Result};

use crate::single_cell::{AnnData, ColumnData};

// ── Types ──────────────────────────────────────────────────────────────────

/// Statistical test method for marker gene identification.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MarkerMethod {
    /// Welch's t-test (unequal variance).
    TTest,
    /// Wilcoxon rank-sum test (Mann-Whitney U).
    Wilcoxon,
    /// Logistic regression via iteratively reweighted least squares.
    LogisticRegression,
}

/// Configuration for marker gene identification.
#[derive(Debug, Clone)]
pub struct MarkerConfig {
    /// Statistical test method.
    pub method: MarkerMethod,
    /// Key in `obs` containing cluster labels.
    pub cluster_key: String,
    /// Minimum log2 fold change threshold for filtering.
    pub log2fc_threshold: f64,
    /// Minimum fraction of cells expressing the gene in the cluster.
    pub min_pct: f64,
    /// Adjusted p-value threshold for filtering.
    pub padj_threshold: f64,
    /// Maximum number of genes to report per cluster (None = all).
    pub n_genes: Option<usize>,
}

impl Default for MarkerConfig {
    fn default() -> Self {
        Self {
            method: MarkerMethod::TTest,
            cluster_key: "leiden".into(),
            log2fc_threshold: 0.25,
            min_pct: 0.1,
            padj_threshold: 0.05,
            n_genes: None,
        }
    }
}

/// A single marker gene result.
#[derive(Debug, Clone)]
pub struct MarkerGene {
    /// Gene name.
    pub gene_name: String,
    /// Gene index in var_names.
    pub gene_index: usize,
    /// Log2 fold change (cluster vs rest).
    pub log2_fold_change: f64,
    /// Fraction of cells expressing the gene in the cluster.
    pub pct_in: f64,
    /// Fraction of cells expressing the gene outside the cluster.
    pub pct_out: f64,
    /// Test statistic.
    pub statistic: f64,
    /// Raw p-value.
    pub p_value: f64,
    /// BH-adjusted p-value.
    pub p_adjusted: f64,
}

/// Results of marker gene analysis.
#[derive(Debug, Clone)]
pub struct MarkerResults {
    /// Marker genes per cluster, keyed by cluster label.
    pub markers: HashMap<String, Vec<MarkerGene>>,
    /// Method used.
    pub method: MarkerMethod,
    /// Number of clusters.
    pub n_clusters: usize,
}

// ── Rank Genes Groups ──────────────────────────────────────────────────────

/// Identify marker genes for each cluster via one-vs-rest differential expression.
///
/// For each cluster, compares gene expression in the cluster vs all other cells
/// using the specified test method. Results are BH-corrected and sorted by
/// adjusted p-value.
pub fn rank_genes_groups(adata: &AnnData, config: &MarkerConfig) -> Result<MarkerResults> {
    let n_obs = adata.n_obs();
    let n_vars = adata.n_vars();

    // Get cluster labels
    let cluster_col = adata
        .get_obs(&config.cluster_key)
        .ok_or_else(|| {
            CyaneaError::InvalidInput(format!(
                "obs['{}'] not found; run clustering first",
                config.cluster_key
            ))
        })?;

    let labels: Vec<String> = match cluster_col {
        ColumnData::Strings(v) => v.clone(),
        ColumnData::Numeric(v) => v.iter().map(|x| x.to_string()).collect(),
        ColumnData::Categorical { codes, categories } => codes
            .iter()
            .map(|&c| {
                categories
                    .get(c as usize)
                    .cloned()
                    .unwrap_or_else(|| c.to_string())
            })
            .collect(),
    };

    let mut unique_labels: Vec<String> = labels.clone();
    unique_labels.sort();
    unique_labels.dedup();
    let n_clusters = unique_labels.len();

    let mut markers: HashMap<String, Vec<MarkerGene>> = HashMap::new();

    for cluster_label in &unique_labels {
        let in_mask: Vec<bool> = labels.iter().map(|l| l == cluster_label).collect();
        let in_indices: Vec<usize> = (0..n_obs).filter(|&i| in_mask[i]).collect();
        let out_indices: Vec<usize> = (0..n_obs).filter(|&i| !in_mask[i]).collect();

        if in_indices.is_empty() || out_indices.is_empty() {
            markers.insert(cluster_label.clone(), Vec::new());
            continue;
        }

        let mut gene_results: Vec<MarkerGene> = Vec::with_capacity(n_vars);
        let mut p_values: Vec<f64> = Vec::with_capacity(n_vars);

        for j in 0..n_vars {
            let in_vals: Vec<f64> = in_indices.iter().map(|&i| adata.x().get(i, j)).collect();
            let out_vals: Vec<f64> = out_indices.iter().map(|&i| adata.x().get(i, j)).collect();

            // Compute log2FC
            let mean_in = mean_or_zero(&in_vals);
            let mean_out = mean_or_zero(&out_vals);
            let log2fc = if mean_out.abs() < 1e-15 {
                if mean_in > 1e-15 {
                    10.0 // cap
                } else {
                    0.0
                }
            } else {
                (mean_in / mean_out).max(1e-15).log2()
            };

            // Fraction expressing
            let pct_in = in_vals.iter().filter(|&&v| v > 0.0).count() as f64 / in_vals.len() as f64;
            let pct_out = if out_vals.is_empty() {
                0.0
            } else {
                out_vals.iter().filter(|&&v| v > 0.0).count() as f64 / out_vals.len() as f64
            };

            // Statistical test
            let (statistic, p_value) = match config.method {
                MarkerMethod::TTest => {
                    match cyanea_stats::testing::t_test_two_sample(&in_vals, &out_vals, false) {
                        Ok(result) if result.p_value.is_finite() => (result.statistic, result.p_value),
                        _ => (0.0, 1.0), // degenerate case (zero variance)
                    }
                }
                MarkerMethod::Wilcoxon => {
                    match cyanea_stats::testing::mann_whitney_u(&in_vals, &out_vals) {
                        Ok(result) if result.p_value.is_finite() => (result.statistic, result.p_value),
                        _ => (0.0, 1.0),
                    }
                }
                MarkerMethod::LogisticRegression => {
                    match logistic_regression_test(&in_vals, &out_vals) {
                        Ok((s, p)) if p.is_finite() => (s, p),
                        _ => (0.0, 1.0),
                    }
                }
            };

            p_values.push(p_value);
            gene_results.push(MarkerGene {
                gene_name: adata.var_names()[j].clone(),
                gene_index: j,
                log2_fold_change: log2fc,
                pct_in,
                pct_out,
                statistic,
                p_value,
                p_adjusted: 0.0, // filled after BH correction
            });
        }

        // BH correction
        let adjusted = cyanea_stats::correction::benjamini_hochberg(&p_values)?;
        for (gene, &padj) in gene_results.iter_mut().zip(adjusted.iter()) {
            gene.p_adjusted = padj;
        }

        // Sort by adjusted p-value
        gene_results.sort_by(|a, b| {
            a.p_adjusted
                .partial_cmp(&b.p_adjusted)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        // Limit number of genes if requested
        if let Some(n) = config.n_genes {
            gene_results.truncate(n);
        }

        markers.insert(cluster_label.clone(), gene_results);
    }

    Ok(MarkerResults {
        markers,
        method: config.method,
        n_clusters,
    })
}

/// Filter marker gene results by log2FC, pct expressing, and adjusted p-value.
pub fn filter_markers(results: &MarkerResults, config: &MarkerConfig) -> MarkerResults {
    let mut filtered: HashMap<String, Vec<MarkerGene>> = HashMap::new();

    for (cluster, genes) in &results.markers {
        let filt: Vec<MarkerGene> = genes
            .iter()
            .filter(|g| {
                g.log2_fold_change >= config.log2fc_threshold
                    && g.pct_in >= config.min_pct
                    && g.p_adjusted <= config.padj_threshold
            })
            .cloned()
            .collect();
        filtered.insert(cluster.clone(), filt);
    }

    MarkerResults {
        markers: filtered,
        method: results.method,
        n_clusters: results.n_clusters,
    }
}

fn mean_or_zero(vals: &[f64]) -> f64 {
    if vals.is_empty() {
        return 0.0;
    }
    vals.iter().sum::<f64>() / vals.len() as f64
}

/// Logistic regression test: fit y~x by IRLS, return (z-statistic, p-value).
fn logistic_regression_test(in_vals: &[f64], out_vals: &[f64]) -> Result<(f64, f64)> {
    let n_in = in_vals.len();
    let n_out = out_vals.len();
    let n = n_in + n_out;

    // y = 1 for in-cluster, 0 for out-cluster
    // X = [1, expression_value] (intercept + gene expression)
    let mut x = vec![0.0; n * 2]; // n × 2 design matrix in row-major
    let mut y = vec![0.0; n];

    for (i, &v) in in_vals.iter().enumerate() {
        x[i * 2] = 1.0;
        x[i * 2 + 1] = v;
        y[i] = 1.0;
    }
    for (i, &v) in out_vals.iter().enumerate() {
        let idx = n_in + i;
        x[idx * 2] = 1.0;
        x[idx * 2 + 1] = v;
        y[idx] = 0.0;
    }

    // IRLS for logistic regression
    let mut beta = [0.0, 0.0]; // [intercept, coefficient]
    let max_iter = 25;

    for _ in 0..max_iter {
        // Compute predictions: p = sigmoid(X * beta)
        let mut p = vec![0.0; n];
        let mut converged = true;

        for i in 0..n {
            let eta = x[i * 2] * beta[0] + x[i * 2 + 1] * beta[1];
            p[i] = sigmoid(eta);
        }

        // W = diag(p * (1-p)), z = X*beta + W^{-1}(y - p)
        // Solve (X^T W X) beta_new = X^T W z
        let mut xtwx = [0.0; 4]; // 2×2 matrix
        let mut xtwy = [0.0; 2]; // 2-vector

        for i in 0..n {
            let w = (p[i] * (1.0 - p[i])).max(1e-10);
            let z = x[i * 2] * beta[0] + x[i * 2 + 1] * beta[1] + (y[i] - p[i]) / w;

            // X^T W X
            xtwx[0] += w * x[i * 2] * x[i * 2];
            xtwx[1] += w * x[i * 2] * x[i * 2 + 1];
            xtwx[2] += w * x[i * 2 + 1] * x[i * 2];
            xtwx[3] += w * x[i * 2 + 1] * x[i * 2 + 1];

            // X^T W z
            xtwy[0] += w * x[i * 2] * z;
            xtwy[1] += w * x[i * 2 + 1] * z;
        }

        // Solve 2×2 system
        let det = xtwx[0] * xtwx[3] - xtwx[1] * xtwx[2];
        if det.abs() < 1e-15 {
            break;
        }
        let new_beta = [
            (xtwx[3] * xtwy[0] - xtwx[1] * xtwy[1]) / det,
            (-xtwx[2] * xtwy[0] + xtwx[0] * xtwy[1]) / det,
        ];

        if (new_beta[0] - beta[0]).abs() > 1e-8 || (new_beta[1] - beta[1]).abs() > 1e-8 {
            converged = false;
        }
        beta = new_beta;
        if converged {
            break;
        }
    }

    // Standard error of beta[1] from Fisher information
    let mut info_11 = 0.0; // (1,1) element of X^T W X
    for i in 0..n {
        let eta = x[i * 2] * beta[0] + x[i * 2 + 1] * beta[1];
        let p = sigmoid(eta);
        let w = (p * (1.0 - p)).max(1e-10);
        info_11 += w * x[i * 2 + 1] * x[i * 2 + 1];
    }

    // Need to invert 2×2 Fisher information to get variance
    let mut info = [0.0; 4];
    for i in 0..n {
        let eta = x[i * 2] * beta[0] + x[i * 2 + 1] * beta[1];
        let p = sigmoid(eta);
        let w = (p * (1.0 - p)).max(1e-10);
        info[0] += w * x[i * 2] * x[i * 2];
        info[1] += w * x[i * 2] * x[i * 2 + 1];
        info[2] += w * x[i * 2 + 1] * x[i * 2];
        info[3] += w * x[i * 2 + 1] * x[i * 2 + 1];
    }
    let det = info[0] * info[3] - info[1] * info[2];
    let se = if det.abs() > 1e-15 {
        (info[0] / det).max(0.0).sqrt() // variance of beta[1] = info_inv[1][1]
    } else {
        1.0
    };

    // z-statistic and p-value (two-sided)
    let z = if se > 1e-15 {
        beta[1] / se
    } else {
        0.0
    };

    // Approximate two-sided p-value using normal CDF approximation
    let p_value = 2.0 * normal_cdf_complement(z.abs());

    Ok((z, p_value))
}

fn sigmoid(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

/// Approximate 1 - Φ(x) for x ≥ 0 using Abramowitz & Stegun.
fn normal_cdf_complement(x: f64) -> f64 {
    if x < 0.0 {
        return 1.0 - normal_cdf_complement(-x);
    }
    // Rational approximation (A&S 26.2.17)
    let b1 = 0.319381530;
    let b2 = -0.356563782;
    let b3 = 1.781477937;
    let b4 = -1.821255978;
    let b5 = 1.330274429;
    let p = 0.2316419;

    let t = 1.0 / (1.0 + p * x);
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let t5 = t4 * t;

    let phi = (-x * x / 2.0).exp() / (2.0 * std::f64::consts::PI).sqrt();
    phi * (b1 * t + b2 * t2 + b3 * t3 + b4 * t4 + b5 * t5)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::single_cell::MatrixData;

    fn make_marker_adata() -> AnnData {
        // 10 cells, 4 genes, 2 clusters
        // Cluster A (cells 0-4): high gene_0, gene_1; low gene_2, gene_3
        // Cluster B (cells 5-9): low gene_0, gene_1; high gene_2, gene_3
        let mut data = vec![vec![0.0; 4]; 10];
        for i in 0..5 {
            data[i][0] = 10.0 + i as f64;
            data[i][1] = 8.0 + i as f64;
            data[i][2] = 0.5;
            data[i][3] = 0.3;
        }
        for i in 5..10 {
            data[i][0] = 0.5;
            data[i][1] = 0.3;
            data[i][2] = 10.0 + (i - 5) as f64;
            data[i][3] = 8.0 + (i - 5) as f64;
        }

        let obs_names: Vec<String> = (0..10).map(|i| format!("cell_{}", i)).collect();
        let var_names = vec![
            "gene_0".into(),
            "gene_1".into(),
            "gene_2".into(),
            "gene_3".into(),
        ];
        let mut adata = AnnData::new(MatrixData::Dense(data), obs_names, var_names).unwrap();

        let labels: Vec<String> = (0..10)
            .map(|i| if i < 5 { "A".into() } else { "B".into() })
            .collect();
        adata
            .add_obs_column("leiden", ColumnData::Strings(labels))
            .unwrap();

        adata
    }

    // ── rank_genes_groups tests ──

    #[test]
    fn rank_genes_ttest() {
        let adata = make_marker_adata();
        let results = rank_genes_groups(
            &adata,
            &MarkerConfig {
                method: MarkerMethod::TTest,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(results.n_clusters, 2);
        assert!(results.markers.contains_key("A"));
        assert!(results.markers.contains_key("B"));

        // Cluster A should have gene_0, gene_1 as top markers
        let a_markers = &results.markers["A"];
        assert!(!a_markers.is_empty());
        // Top marker should be gene_0 or gene_1 (highest FC)
        assert!(
            a_markers[0].gene_name == "gene_0" || a_markers[0].gene_name == "gene_1",
            "top marker for A: {}",
            a_markers[0].gene_name
        );
    }

    #[test]
    fn rank_genes_wilcoxon() {
        let adata = make_marker_adata();
        let results = rank_genes_groups(
            &adata,
            &MarkerConfig {
                method: MarkerMethod::Wilcoxon,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(results.method, MarkerMethod::Wilcoxon);
        let b_markers = &results.markers["B"];
        assert!(!b_markers.is_empty());
    }

    #[test]
    fn rank_genes_logistic() {
        let adata = make_marker_adata();
        let results = rank_genes_groups(
            &adata,
            &MarkerConfig {
                method: MarkerMethod::LogisticRegression,
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(results.method, MarkerMethod::LogisticRegression);
        assert!(results.markers.contains_key("A"));
    }

    #[test]
    fn rank_genes_missing_cluster() {
        let x = MatrixData::Dense(vec![vec![1.0, 2.0]]);
        let adata = AnnData::new(x, vec!["c0".into()], vec!["g0".into(), "g1".into()]).unwrap();
        let result = rank_genes_groups(&adata, &MarkerConfig::default());
        assert!(result.is_err());
    }

    #[test]
    fn rank_genes_n_genes_limit() {
        let adata = make_marker_adata();
        let results = rank_genes_groups(
            &adata,
            &MarkerConfig {
                n_genes: Some(1),
                ..Default::default()
            },
        )
        .unwrap();

        for genes in results.markers.values() {
            assert!(genes.len() <= 1);
        }
    }

    #[test]
    fn rank_genes_sorted_by_padj() {
        let adata = make_marker_adata();
        let results = rank_genes_groups(&adata, &MarkerConfig::default()).unwrap();

        for genes in results.markers.values() {
            for i in 1..genes.len() {
                assert!(
                    genes[i - 1].p_adjusted <= genes[i].p_adjusted + 1e-10,
                    "markers not sorted by p_adjusted"
                );
            }
        }
    }

    #[test]
    fn rank_genes_log2fc_sign() {
        let adata = make_marker_adata();
        let results = rank_genes_groups(&adata, &MarkerConfig::default()).unwrap();

        // Cluster A: gene_0,1 should have positive log2FC; gene_2,3 negative
        let a_markers = &results.markers["A"];
        let gene0 = a_markers.iter().find(|g| g.gene_name == "gene_0").unwrap();
        let gene2 = a_markers.iter().find(|g| g.gene_name == "gene_2").unwrap();
        assert!(gene0.log2_fold_change > 0.0, "gene_0 l2fc = {}", gene0.log2_fold_change);
        assert!(gene2.log2_fold_change < 0.0, "gene_2 l2fc = {}", gene2.log2_fold_change);
    }

    #[test]
    fn rank_genes_pct_values() {
        let adata = make_marker_adata();
        let results = rank_genes_groups(&adata, &MarkerConfig::default()).unwrap();

        let a_markers = &results.markers["A"];
        for gene in a_markers {
            assert!((0.0..=1.0).contains(&gene.pct_in));
            assert!((0.0..=1.0).contains(&gene.pct_out));
        }
    }

    #[test]
    fn rank_genes_numeric_labels() {
        let mut data = vec![vec![0.0; 2]; 6];
        for i in 0..3 {
            data[i][0] = 10.0;
        }
        for i in 3..6 {
            data[i][1] = 10.0;
        }
        let obs_names: Vec<String> = (0..6).map(|i| format!("c{}", i)).collect();
        let var_names = vec!["g0".into(), "g1".into()];
        let mut adata = AnnData::new(MatrixData::Dense(data), obs_names, var_names).unwrap();
        adata
            .add_obs_column(
                "leiden",
                ColumnData::Numeric(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]),
            )
            .unwrap();

        let results = rank_genes_groups(&adata, &MarkerConfig::default()).unwrap();
        assert_eq!(results.n_clusters, 2);
    }

    // ── filter_markers tests ──

    #[test]
    fn filter_markers_by_log2fc() {
        let adata = make_marker_adata();
        let results = rank_genes_groups(&adata, &MarkerConfig::default()).unwrap();
        let filtered = filter_markers(
            &results,
            &MarkerConfig {
                log2fc_threshold: 2.0,
                min_pct: 0.0,
                padj_threshold: 1.0,
                ..Default::default()
            },
        );

        for genes in filtered.markers.values() {
            for gene in genes {
                assert!(gene.log2_fold_change >= 2.0);
            }
        }
    }

    #[test]
    fn filter_markers_by_pct() {
        let adata = make_marker_adata();
        let results = rank_genes_groups(&adata, &MarkerConfig::default()).unwrap();
        let filtered = filter_markers(
            &results,
            &MarkerConfig {
                log2fc_threshold: 0.0,
                min_pct: 0.5,
                padj_threshold: 1.0,
                ..Default::default()
            },
        );

        for genes in filtered.markers.values() {
            for gene in genes {
                assert!(gene.pct_in >= 0.5);
            }
        }
    }

    #[test]
    fn filter_markers_by_padj() {
        let adata = make_marker_adata();
        let results = rank_genes_groups(&adata, &MarkerConfig::default()).unwrap();
        let filtered = filter_markers(
            &results,
            &MarkerConfig {
                log2fc_threshold: 0.0,
                min_pct: 0.0,
                padj_threshold: 0.01,
                ..Default::default()
            },
        );

        for genes in filtered.markers.values() {
            for gene in genes {
                assert!(gene.p_adjusted <= 0.01);
            }
        }
    }

    #[test]
    fn filter_markers_preserves_method() {
        let adata = make_marker_adata();
        let results = rank_genes_groups(
            &adata,
            &MarkerConfig {
                method: MarkerMethod::Wilcoxon,
                ..Default::default()
            },
        )
        .unwrap();
        let filtered = filter_markers(&results, &MarkerConfig::default());
        assert_eq!(filtered.method, MarkerMethod::Wilcoxon);
    }

    // ── Logistic regression tests ──

    #[test]
    fn logistic_regression_separable() {
        // Perfectly separable data
        let in_vals = vec![10.0, 11.0, 12.0, 13.0, 14.0];
        let out_vals = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let (z, p) = logistic_regression_test(&in_vals, &out_vals).unwrap();
        assert!(z.abs() > 0.0, "z-statistic should be nonzero");
        assert!((0.0..=1.0).contains(&p), "p-value should be in [0,1]");
    }

    #[test]
    fn logistic_regression_identical() {
        // Same distribution → no significant difference
        let in_vals = vec![5.0, 5.0, 5.0, 5.0, 5.0];
        let out_vals = vec![5.0, 5.0, 5.0, 5.0, 5.0];
        let (z, p) = logistic_regression_test(&in_vals, &out_vals).unwrap();
        assert!(z.abs() < 1e-6 || p > 0.5, "identical groups should not be significant");
    }

    // ── Helper tests ──

    #[test]
    fn sigmoid_values() {
        assert!((sigmoid(0.0) - 0.5).abs() < 1e-10);
        assert!(sigmoid(10.0) > 0.99);
        assert!(sigmoid(-10.0) < 0.01);
    }

    #[test]
    fn normal_cdf_complement_values() {
        // Φ(0) = 0.5, so complement = 0.5
        let c0 = normal_cdf_complement(0.0);
        assert!((c0 - 0.5).abs() < 0.01, "complement at 0 = {}", c0);

        // At z=2, complement ≈ 0.0228
        let c2 = normal_cdf_complement(2.0);
        assert!((c2 - 0.0228).abs() < 0.005, "complement at 2 = {}", c2);

        // At z=3, complement ≈ 0.00135
        let c3 = normal_cdf_complement(3.0);
        assert!(c3 < 0.005, "complement at 3 = {}", c3);
    }

    #[test]
    fn mean_or_zero_empty() {
        assert_eq!(mean_or_zero(&[]), 0.0);
    }

    #[test]
    fn mean_or_zero_values() {
        assert!((mean_or_zero(&[2.0, 4.0, 6.0]) - 4.0).abs() < 1e-10);
    }

    // ── Pipeline test: full cluster → markers ──

    #[test]
    fn full_marker_pipeline() {
        let adata = make_marker_adata();
        let results = rank_genes_groups(&adata, &MarkerConfig::default()).unwrap();
        let filtered = filter_markers(
            &results,
            &MarkerConfig {
                log2fc_threshold: 1.0,
                padj_threshold: 0.05,
                min_pct: 0.5,
                ..Default::default()
            },
        );

        // Cluster A should have at least 1 significant marker (gene_0 or gene_1)
        let a_markers = &filtered.markers["A"];
        assert!(
            !a_markers.is_empty(),
            "cluster A should have significant markers"
        );
        // All should be upregulated in A
        for gene in a_markers {
            assert!(gene.log2_fold_change >= 1.0);
            assert!(gene.pct_in >= 0.5);
            assert!(gene.p_adjusted <= 0.05);
        }
    }

    #[test]
    fn three_cluster_markers() {
        let mut data = vec![vec![0.0; 3]; 15];
        for i in 0..5 {
            data[i][0] = 10.0 + i as f64;
        }
        for i in 5..10 {
            data[i][1] = 10.0 + (i - 5) as f64;
        }
        for i in 10..15 {
            data[i][2] = 10.0 + (i - 10) as f64;
        }

        let obs_names: Vec<String> = (0..15).map(|i| format!("c{}", i)).collect();
        let var_names = vec!["g0".into(), "g1".into(), "g2".into()];
        let mut adata = AnnData::new(MatrixData::Dense(data), obs_names, var_names).unwrap();

        let labels: Vec<String> = (0..15)
            .map(|i| match i / 5 {
                0 => "A",
                1 => "B",
                _ => "C",
            }
            .into())
            .collect();
        adata
            .add_obs_column("leiden", ColumnData::Strings(labels))
            .unwrap();

        let results = rank_genes_groups(&adata, &MarkerConfig::default()).unwrap();
        assert_eq!(results.n_clusters, 3);
        assert!(results.markers.contains_key("C"));
    }
}
