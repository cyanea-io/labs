//! Gene regulatory network (GRN) inference.
//!
//! Infer regulatory relationships between transcription factors and target genes
//! from expression data using correlation, mutual information, and CLR algorithms.

use crate::graph::{Edge, Graph, GraphType};
use cyanea_core::{CyaneaError, Result};
use std::collections::HashMap;

/// Expression matrix: genes × samples.
///
/// Each row is a gene's expression across samples.
#[derive(Debug, Clone)]
pub struct ExpressionMatrix {
    /// Gene identifiers (row labels).
    pub genes: Vec<String>,
    /// Sample identifiers (column labels).
    pub samples: Vec<String>,
    /// Expression values (genes × samples), row-major.
    pub data: Vec<Vec<f64>>,
}

impl ExpressionMatrix {
    pub fn new(genes: Vec<String>, samples: Vec<String>, data: Vec<Vec<f64>>) -> Result<Self> {
        if data.len() != genes.len() {
            return Err(CyaneaError::InvalidInput(
                "row count doesn't match gene count".into(),
            ));
        }
        for row in &data {
            if row.len() != samples.len() {
                return Err(CyaneaError::InvalidInput(
                    "column count doesn't match sample count".into(),
                ));
            }
        }
        Ok(Self {
            genes,
            samples,
            data,
        })
    }

    /// Number of genes.
    pub fn num_genes(&self) -> usize {
        self.genes.len()
    }

    /// Number of samples.
    pub fn num_samples(&self) -> usize {
        self.samples.len()
    }

    /// Get expression vector for a gene.
    pub fn gene_expression(&self, gene_idx: usize) -> &[f64] {
        &self.data[gene_idx]
    }
}

/// Result of GRN inference.
#[derive(Debug, Clone)]
pub struct GrnResult {
    /// The inferred network.
    pub network: Graph,
    /// Edge scores: (source, target) → score.
    pub scores: HashMap<(String, String), f64>,
    /// Method used for inference.
    pub method: String,
}

/// Pearson correlation coefficient between two vectors.
fn pearson(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    if n < 2.0 {
        return 0.0;
    }

    let mean_x: f64 = x.iter().sum::<f64>() / n;
    let mean_y: f64 = y.iter().sum::<f64>() / n;

    let mut cov = 0.0;
    let mut var_x = 0.0;
    let mut var_y = 0.0;

    for i in 0..x.len() {
        let dx = x[i] - mean_x;
        let dy = y[i] - mean_y;
        cov += dx * dy;
        var_x += dx * dx;
        var_y += dy * dy;
    }

    if var_x < 1e-15 || var_y < 1e-15 {
        return 0.0;
    }

    cov / (var_x.sqrt() * var_y.sqrt())
}

/// Spearman rank correlation between two vectors.
fn spearman(x: &[f64], y: &[f64]) -> f64 {
    let rank_x = ranks(x);
    let rank_y = ranks(y);
    pearson(&rank_x, &rank_y)
}

/// Compute ranks (average rank for ties).
fn ranks(values: &[f64]) -> Vec<f64> {
    let n = values.len();
    let mut indexed: Vec<(usize, f64)> = values.iter().copied().enumerate().collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut result = vec![0.0; n];
    let mut i = 0;
    while i < n {
        let mut j = i;
        while j < n && (indexed[j].1 - indexed[i].1).abs() < 1e-15 {
            j += 1;
        }
        let avg_rank = (i + j - 1) as f64 / 2.0 + 1.0;
        for k in i..j {
            result[indexed[k].0] = avg_rank;
        }
        i = j;
    }

    result
}

/// Mutual information between two continuous variables.
///
/// Uses histogram-based estimation with sqrt(n) bins.
fn mutual_information(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len();
    if n < 4 {
        return 0.0;
    }

    let num_bins = (n as f64).sqrt().ceil() as usize;
    let num_bins = num_bins.max(2);

    let x_min = x.iter().copied().fold(f64::INFINITY, f64::min);
    let x_max = x.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let y_min = y.iter().copied().fold(f64::INFINITY, f64::min);
    let y_max = y.iter().copied().fold(f64::NEG_INFINITY, f64::max);

    let x_range = (x_max - x_min).max(1e-15);
    let y_range = (y_max - y_min).max(1e-15);

    let mut joint = vec![vec![0usize; num_bins]; num_bins];
    let mut marginal_x = vec![0usize; num_bins];
    let mut marginal_y = vec![0usize; num_bins];

    for i in 0..n {
        let bx = (((x[i] - x_min) / x_range) * (num_bins - 1) as f64).round() as usize;
        let by = (((y[i] - y_min) / y_range) * (num_bins - 1) as f64).round() as usize;
        let bx = bx.min(num_bins - 1);
        let by = by.min(num_bins - 1);
        joint[bx][by] += 1;
        marginal_x[bx] += 1;
        marginal_y[by] += 1;
    }

    let nf = n as f64;
    let mut mi = 0.0;
    for bx in 0..num_bins {
        for by in 0..num_bins {
            if joint[bx][by] > 0 && marginal_x[bx] > 0 && marginal_y[by] > 0 {
                let pxy = joint[bx][by] as f64 / nf;
                let px = marginal_x[bx] as f64 / nf;
                let py = marginal_y[by] as f64 / nf;
                mi += pxy * (pxy / (px * py)).ln();
            }
        }
    }

    mi.max(0.0)
}

/// Correlation method for GRN inference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CorrelationMethod {
    Pearson,
    Spearman,
    MutualInformation,
}

/// Infer a GRN using pairwise correlation/MI.
///
/// Builds a regulatory network where edges represent regulatory relationships
/// between TFs (specified) and all genes. Edge weight = |correlation| or MI score.
///
/// # Arguments
///
/// * `expr` - Expression matrix
/// * `tf_genes` - Indices of transcription factor genes (regulators)
/// * `method` - Correlation method
/// * `threshold` - Minimum score to include an edge
pub fn infer_grn(
    expr: &ExpressionMatrix,
    tf_genes: &[usize],
    method: CorrelationMethod,
    threshold: f64,
) -> Result<GrnResult> {
    let n_genes = expr.num_genes();

    for &tf in tf_genes {
        if tf >= n_genes {
            return Err(CyaneaError::InvalidInput(
                format!("TF index {} out of range ({})", tf, n_genes),
            ));
        }
    }

    let mut graph = Graph::new(GraphType::Directed);
    let mut scores = HashMap::new();

    // Add all genes as nodes
    for (i, gene) in expr.genes.iter().enumerate() {
        let mut node = crate::graph::Node::new(gene, gene);
        if tf_genes.contains(&i) {
            node = node.with_attr("type", "TF");
        } else {
            node = node.with_attr("type", "target");
        }
        graph.add_node_struct(node)?;
    }

    // Compute pairwise scores between TFs and all genes
    for &tf_idx in tf_genes {
        let tf_expr = expr.gene_expression(tf_idx);

        for target_idx in 0..n_genes {
            if target_idx == tf_idx {
                continue;
            }

            let target_expr = expr.gene_expression(target_idx);

            let score = match method {
                CorrelationMethod::Pearson => pearson(tf_expr, target_expr).abs(),
                CorrelationMethod::Spearman => spearman(tf_expr, target_expr).abs(),
                CorrelationMethod::MutualInformation => {
                    mutual_information(tf_expr, target_expr)
                }
            };

            if score >= threshold {
                let edge = Edge::weighted(
                    &expr.genes[tf_idx],
                    &expr.genes[target_idx],
                    score,
                )
                .with_type("regulation");
                graph.add_edge_struct(edge)?;
                scores.insert(
                    (expr.genes[tf_idx].clone(), expr.genes[target_idx].clone()),
                    score,
                );
            }
        }
    }

    let method_name = match method {
        CorrelationMethod::Pearson => "pearson",
        CorrelationMethod::Spearman => "spearman",
        CorrelationMethod::MutualInformation => "mutual_information",
    };

    Ok(GrnResult {
        network: graph,
        scores,
        method: method_name.to_string(),
    })
}

/// Context Likelihood of Relatedness (CLR) algorithm.
///
/// Enhances mutual information scores by comparing each MI value to the
/// background distribution of MI values for each gene. This reduces
/// false positives from highly variable genes.
///
/// CLR score = sqrt(z_tf² + z_target²) where z = (MI - mean) / std.
pub fn clr(
    expr: &ExpressionMatrix,
    tf_genes: &[usize],
    threshold: f64,
) -> Result<GrnResult> {
    let n_genes = expr.num_genes();

    // Compute full MI matrix between TFs and all genes
    let mut mi_matrix: Vec<Vec<f64>> = vec![vec![0.0; n_genes]; tf_genes.len()];

    for (ti, &tf_idx) in tf_genes.iter().enumerate() {
        let tf_expr = expr.gene_expression(tf_idx);
        for target_idx in 0..n_genes {
            if target_idx == tf_idx {
                continue;
            }
            mi_matrix[ti][target_idx] = mutual_information(tf_expr, expr.gene_expression(target_idx));
        }
    }

    // Compute z-scores for each TF row
    let mut z_tf: Vec<Vec<f64>> = vec![vec![0.0; n_genes]; tf_genes.len()];
    for ti in 0..tf_genes.len() {
        let vals: Vec<f64> = mi_matrix[ti].iter().copied().filter(|&v| v > 0.0).collect();
        let mean = if vals.is_empty() {
            0.0
        } else {
            vals.iter().sum::<f64>() / vals.len() as f64
        };
        let std = if vals.len() < 2 {
            1.0
        } else {
            let var = vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (vals.len() - 1) as f64;
            var.sqrt().max(1e-15)
        };
        for j in 0..n_genes {
            z_tf[ti][j] = ((mi_matrix[ti][j] - mean) / std).max(0.0);
        }
    }

    // Compute z-scores for each target column
    let mut z_target: Vec<Vec<f64>> = vec![vec![0.0; n_genes]; tf_genes.len()];
    for j in 0..n_genes {
        let vals: Vec<f64> = (0..tf_genes.len())
            .map(|ti| mi_matrix[ti][j])
            .filter(|&v| v > 0.0)
            .collect();
        let mean = if vals.is_empty() {
            0.0
        } else {
            vals.iter().sum::<f64>() / vals.len() as f64
        };
        let std = if vals.len() < 2 {
            1.0
        } else {
            let var = vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (vals.len() - 1) as f64;
            var.sqrt().max(1e-15)
        };
        for ti in 0..tf_genes.len() {
            z_target[ti][j] = ((mi_matrix[ti][j] - mean) / std).max(0.0);
        }
    }

    // Build CLR network
    let mut graph = Graph::new(GraphType::Directed);
    let mut scores_map = HashMap::new();

    for gene in &expr.genes {
        graph.add_node(gene, gene)?;
    }

    for (ti, &tf_idx) in tf_genes.iter().enumerate() {
        for target_idx in 0..n_genes {
            if target_idx == tf_idx {
                continue;
            }
            let clr_score = (z_tf[ti][target_idx].powi(2) + z_target[ti][target_idx].powi(2)).sqrt();
            if clr_score >= threshold {
                let edge = Edge::weighted(
                    &expr.genes[tf_idx],
                    &expr.genes[target_idx],
                    clr_score,
                )
                .with_type("regulation");
                graph.add_edge_struct(edge)?;
                scores_map.insert(
                    (expr.genes[tf_idx].clone(), expr.genes[target_idx].clone()),
                    clr_score,
                );
            }
        }
    }

    Ok(GrnResult {
        network: graph,
        scores: scores_map,
        method: "CLR".to_string(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_expression() -> ExpressionMatrix {
        // TF1 activates G1 and G2, TF2 activates G3
        // 10 samples with correlated patterns
        ExpressionMatrix::new(
            vec![
                "TF1".into(),
                "TF2".into(),
                "G1".into(),
                "G2".into(),
                "G3".into(),
            ],
            (0..10).map(|i| format!("S{}", i)).collect(),
            vec![
                vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], // TF1
                vec![10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0], // TF2
                vec![1.1, 2.2, 2.8, 4.1, 5.2, 5.9, 7.1, 8.0, 8.8, 10.2], // G1 (corr w/ TF1)
                vec![0.9, 1.8, 3.2, 3.9, 4.8, 6.1, 7.2, 7.9, 9.1, 9.8],  // G2 (corr w/ TF1)
                vec![9.8, 9.1, 7.9, 7.2, 6.1, 4.8, 3.9, 3.2, 1.8, 0.9],  // G3 (corr w/ TF2)
            ],
        )
        .unwrap()
    }

    #[test]
    fn test_pearson_perfect() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![2.0, 4.0, 6.0, 8.0, 10.0];
        assert!((pearson(&x, &y) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_pearson_negative() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![10.0, 8.0, 6.0, 4.0, 2.0];
        assert!((pearson(&x, &y) - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_spearman_rank() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![1.0, 4.0, 9.0, 16.0, 25.0]; // monotonic
        assert!(spearman(&x, &y) > 0.99);
    }

    #[test]
    fn test_mutual_information_correlated() {
        let x: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let y: Vec<f64> = x.iter().map(|v| v * 2.0 + 1.0).collect();
        let mi = mutual_information(&x, &y);
        assert!(mi > 0.0);
    }

    #[test]
    fn test_expression_matrix() {
        let expr = sample_expression();
        assert_eq!(expr.num_genes(), 5);
        assert_eq!(expr.num_samples(), 10);
    }

    #[test]
    fn test_expression_matrix_invalid() {
        let result = ExpressionMatrix::new(
            vec!["A".into()],
            vec!["S1".into(), "S2".into()],
            vec![vec![1.0]], // wrong column count
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_infer_grn_pearson() {
        let expr = sample_expression();
        let result = infer_grn(&expr, &[0, 1], CorrelationMethod::Pearson, 0.9).unwrap();

        // TF1→G1 and TF1→G2 should be strong positive correlations
        assert!(result.scores.contains_key(&("TF1".into(), "G1".into())));
        assert!(result.scores.contains_key(&("TF1".into(), "G2".into())));

        // TF2→G3 should also be strong (negative becomes abs)
        assert!(result.scores.contains_key(&("TF2".into(), "G3".into())));

        assert_eq!(result.method, "pearson");
    }

    #[test]
    fn test_infer_grn_spearman() {
        let expr = sample_expression();
        let result = infer_grn(&expr, &[0, 1], CorrelationMethod::Spearman, 0.9).unwrap();
        assert!(!result.scores.is_empty());
    }

    #[test]
    fn test_infer_grn_mi() {
        let expr = sample_expression();
        let result = infer_grn(
            &expr,
            &[0, 1],
            CorrelationMethod::MutualInformation,
            0.0,
        )
        .unwrap();
        assert!(!result.scores.is_empty());
    }

    #[test]
    fn test_infer_grn_high_threshold() {
        let expr = sample_expression();
        let result = infer_grn(&expr, &[0, 1], CorrelationMethod::Pearson, 1.1).unwrap();
        assert!(result.scores.is_empty());
    }

    #[test]
    fn test_clr() {
        let expr = sample_expression();
        let result = clr(&expr, &[0, 1], 0.0).unwrap();
        assert!(!result.scores.is_empty());
        assert_eq!(result.method, "CLR");
    }

    #[test]
    fn test_grn_node_types() {
        let expr = sample_expression();
        let result = infer_grn(&expr, &[0, 1], CorrelationMethod::Pearson, 0.0).unwrap();
        let tf1 = result.network.get_node("TF1").unwrap();
        assert_eq!(tf1.attributes.get("type").unwrap(), "TF");
        let g1 = result.network.get_node("G1").unwrap();
        assert_eq!(g1.attributes.get("type").unwrap(), "target");
    }
}
