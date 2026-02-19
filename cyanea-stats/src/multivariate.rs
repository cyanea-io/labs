//! Multivariate statistical tests for community ecology.
//!
//! - **PERMANOVA** — permutational multivariate analysis of variance
//! - **ANOSIM** — analysis of similarities
//! - **Mantel test** — correlation between distance matrices
//! - **AMOVA** — analysis of molecular variance
//! - **BIOENV** — best subset of environmental variables explaining community structure

use cyanea_core::{CyaneaError, Result};

use crate::correlation;
use crate::rank::{rank, RankMethod};

// ── PERMANOVA ───────────────────────────────────────────────────────────────

/// Result of PERMANOVA.
#[derive(Debug, Clone)]
pub struct PermanovaResult {
    /// Pseudo-F statistic.
    pub f_statistic: f64,
    /// Permutation p-value.
    pub p_value: f64,
    /// Number of permutations performed.
    pub n_permutations: usize,
    /// R² (proportion of variance explained by grouping).
    pub r_squared: f64,
    /// Number of groups.
    pub n_groups: usize,
    /// Method name.
    pub method: String,
}

/// PERMANOVA: test whether groups differ in multivariate dispersion.
///
/// Uses squared distances to partition total sum of squares into among-group
/// and within-group components, then tests significance by permuting group labels.
///
/// # Arguments
///
/// * `distances` — symmetric pairwise distance matrix
/// * `groups` — group label for each sample (0-indexed integers)
/// * `n_permutations` — number of permutations for p-value estimation
/// * `seed` — random seed for reproducibility
///
/// # Errors
///
/// Returns an error if the distance matrix is not square, groups don't match,
/// or fewer than 2 groups exist.
pub fn permanova(
    distances: &[Vec<f64>],
    groups: &[usize],
    n_permutations: usize,
    seed: u64,
) -> Result<PermanovaResult> {
    let n = distances.len();
    validate_distance_matrix(distances, n)?;
    if groups.len() != n {
        return Err(CyaneaError::InvalidInput(format!(
            "permanova: groups length ({}) != distance matrix size ({})",
            groups.len(),
            n
        )));
    }
    if n_permutations == 0 {
        return Err(CyaneaError::InvalidInput(
            "permanova: n_permutations must be > 0".into(),
        ));
    }

    let unique_groups: std::collections::HashSet<usize> = groups.iter().copied().collect();
    let k = unique_groups.len();
    if k < 2 {
        return Err(CyaneaError::InvalidInput(
            "permanova: at least 2 groups required".into(),
        ));
    }

    let observed_f = compute_pseudo_f(distances, groups, n, k);

    // Permutation test
    let mut rng = Xorshift64::new(seed);
    let mut perm_groups: Vec<usize> = groups.to_vec();
    let mut n_extreme = 0usize;

    for _ in 0..n_permutations {
        fisher_yates_shuffle_usize(&mut perm_groups, &mut rng);
        let perm_f = compute_pseudo_f(distances, &perm_groups, n, k);
        if perm_f >= observed_f {
            n_extreme += 1;
        }
    }

    let p_value = (n_extreme as f64 + 1.0) / (n_permutations as f64 + 1.0);

    // R² = SS_among / SS_total
    let ss_total = compute_ss_total(distances, n);
    let ss_within = compute_ss_within(distances, groups, n);
    let ss_among = ss_total - ss_within;
    let r_squared = if ss_total > 0.0 {
        ss_among / ss_total
    } else {
        0.0
    };

    Ok(PermanovaResult {
        f_statistic: observed_f,
        p_value,
        n_permutations,
        r_squared,
        n_groups: k,
        method: "PERMANOVA".to_string(),
    })
}

fn compute_pseudo_f(distances: &[Vec<f64>], groups: &[usize], n: usize, k: usize) -> f64 {
    let ss_total = compute_ss_total(distances, n);
    let ss_within = compute_ss_within(distances, groups, n);
    let ss_among = ss_total - ss_within;

    let df_among = (k - 1) as f64;
    let df_within = (n - k) as f64;

    if df_within <= 0.0 || ss_within <= 0.0 {
        return 0.0;
    }

    (ss_among / df_among) / (ss_within / df_within)
}

fn compute_ss_total(distances: &[Vec<f64>], n: usize) -> f64 {
    let mut ss = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            ss += distances[i][j] * distances[i][j];
        }
    }
    ss / n as f64
}

fn compute_ss_within(distances: &[Vec<f64>], groups: &[usize], _n: usize) -> f64 {
    // Group samples by label
    let mut group_map: std::collections::HashMap<usize, Vec<usize>> =
        std::collections::HashMap::new();
    for (i, &g) in groups.iter().enumerate() {
        group_map.entry(g).or_default().push(i);
    }

    let mut ss_within = 0.0;
    for members in group_map.values() {
        let ng = members.len();
        if ng < 2 {
            continue;
        }
        let mut ss_g = 0.0;
        for ii in 0..ng {
            for jj in (ii + 1)..ng {
                let d = distances[members[ii]][members[jj]];
                ss_g += d * d;
            }
        }
        ss_within += ss_g / ng as f64;
    }
    ss_within
}

// ── ANOSIM ──────────────────────────────────────────────────────────────────

/// Result of ANOSIM.
#[derive(Debug, Clone)]
pub struct AnosimResult {
    /// ANOSIM R statistic (range: -1 to 1).
    pub r_statistic: f64,
    /// Permutation p-value.
    pub p_value: f64,
    /// Number of permutations performed.
    pub n_permutations: usize,
}

/// ANOSIM: analysis of similarities.
///
/// Compares mean ranks of between-group distances to mean ranks of within-group
/// distances. `R = (r_between - r_within) / (n(n-1)/4)`.
///
/// # Errors
///
/// Returns an error if the distance matrix is invalid or fewer than 2 groups.
pub fn anosim(
    distances: &[Vec<f64>],
    groups: &[usize],
    n_permutations: usize,
    seed: u64,
) -> Result<AnosimResult> {
    let n = distances.len();
    validate_distance_matrix(distances, n)?;
    if groups.len() != n {
        return Err(CyaneaError::InvalidInput(format!(
            "anosim: groups length ({}) != distance matrix size ({})",
            groups.len(),
            n
        )));
    }
    if n_permutations == 0 {
        return Err(CyaneaError::InvalidInput(
            "anosim: n_permutations must be > 0".into(),
        ));
    }

    let unique_groups: std::collections::HashSet<usize> = groups.iter().copied().collect();
    if unique_groups.len() < 2 {
        return Err(CyaneaError::InvalidInput(
            "anosim: at least 2 groups required".into(),
        ));
    }

    // Rank all pairwise distances
    let ranks = rank_distances(distances, n);

    let observed_r = compute_anosim_r(&ranks, groups, n);

    // Permutation test
    let mut rng = Xorshift64::new(seed);
    let mut perm_groups: Vec<usize> = groups.to_vec();
    let mut n_extreme = 0usize;

    for _ in 0..n_permutations {
        fisher_yates_shuffle_usize(&mut perm_groups, &mut rng);
        let perm_r = compute_anosim_r(&ranks, &perm_groups, n);
        if perm_r >= observed_r {
            n_extreme += 1;
        }
    }

    let p_value = (n_extreme as f64 + 1.0) / (n_permutations as f64 + 1.0);

    Ok(AnosimResult {
        r_statistic: observed_r,
        p_value,
        n_permutations,
    })
}

fn rank_distances(distances: &[Vec<f64>], n: usize) -> Vec<Vec<f64>> {
    // Collect upper triangle, rank them, put back
    let n_pairs = n * (n - 1) / 2;
    let mut flat = Vec::with_capacity(n_pairs);
    for i in 0..n {
        for j in (i + 1)..n {
            flat.push(distances[i][j]);
        }
    }
    let ranked = rank(&flat, RankMethod::Average);

    let mut result = vec![vec![0.0; n]; n];
    let mut idx = 0;
    for i in 0..n {
        for j in (i + 1)..n {
            result[i][j] = ranked[idx];
            result[j][i] = ranked[idx];
            idx += 1;
        }
    }
    result
}

fn compute_anosim_r(ranks: &[Vec<f64>], groups: &[usize], n: usize) -> f64 {
    let mut r_between_sum = 0.0;
    let mut r_within_sum = 0.0;
    let mut n_between = 0usize;
    let mut n_within = 0usize;

    for i in 0..n {
        for j in (i + 1)..n {
            if groups[i] == groups[j] {
                r_within_sum += ranks[i][j];
                n_within += 1;
            } else {
                r_between_sum += ranks[i][j];
                n_between += 1;
            }
        }
    }

    let r_between = if n_between > 0 {
        r_between_sum / n_between as f64
    } else {
        0.0
    };
    let r_within = if n_within > 0 {
        r_within_sum / n_within as f64
    } else {
        0.0
    };

    let m = (n * (n - 1)) as f64 / 4.0;
    if m == 0.0 {
        return 0.0;
    }

    (r_between - r_within) / m
}

// ── Mantel test ─────────────────────────────────────────────────────────────

/// Result of Mantel test.
#[derive(Debug, Clone)]
pub struct MantelResult {
    /// Correlation statistic (Pearson or Spearman).
    pub statistic: f64,
    /// Permutation p-value.
    pub p_value: f64,
    /// Number of permutations performed.
    pub n_permutations: usize,
    /// Method used ("pearson" or "spearman").
    pub method: String,
}

/// Mantel test: correlation between two distance matrices.
///
/// Flattens the upper triangles of both matrices, computes the correlation
/// (Pearson or Spearman), and tests significance by permuting rows/columns
/// of one matrix.
///
/// # Arguments
///
/// * `matrix_a`, `matrix_b` — symmetric distance matrices of the same size
/// * `n_permutations` — number of permutations
/// * `seed` — random seed
/// * `method` — "pearson" or "spearman"
///
/// # Errors
///
/// Returns an error if matrices have different sizes or are not square.
pub fn mantel_test(
    matrix_a: &[Vec<f64>],
    matrix_b: &[Vec<f64>],
    n_permutations: usize,
    seed: u64,
    method: &str,
) -> Result<MantelResult> {
    let n = matrix_a.len();
    validate_distance_matrix(matrix_a, n)?;
    validate_distance_matrix(matrix_b, n)?;
    if matrix_b.len() != n {
        return Err(CyaneaError::InvalidInput(format!(
            "mantel: matrices must have same size ({} vs {})",
            n,
            matrix_b.len()
        )));
    }
    if n_permutations == 0 {
        return Err(CyaneaError::InvalidInput(
            "mantel: n_permutations must be > 0".into(),
        ));
    }
    if method != "pearson" && method != "spearman" {
        return Err(CyaneaError::InvalidInput(format!(
            "mantel: method must be 'pearson' or 'spearman', got '{}'",
            method
        )));
    }

    // Flatten upper triangles
    let flat_a = upper_triangle_vec(matrix_a, n);
    let flat_b = upper_triangle_vec(matrix_b, n);

    let observed = if method == "pearson" {
        correlation::pearson(&flat_a, &flat_b)?
    } else {
        correlation::spearman(&flat_a, &flat_b)?
    };

    // Permute rows/columns of matrix_a
    let mut rng = Xorshift64::new(seed);
    let mut perm_indices: Vec<usize> = (0..n).collect();
    let mut n_extreme = 0usize;

    for _ in 0..n_permutations {
        fisher_yates_shuffle_usize(&mut perm_indices, &mut rng);
        let perm_flat: Vec<f64> = {
            let mut v = Vec::with_capacity(flat_a.len());
            for i in 0..n {
                for j in (i + 1)..n {
                    v.push(matrix_a[perm_indices[i]][perm_indices[j]]);
                }
            }
            v
        };
        let perm_stat = if method == "pearson" {
            correlation::pearson(&perm_flat, &flat_b)?
        } else {
            correlation::spearman(&perm_flat, &flat_b)?
        };
        if perm_stat >= observed {
            n_extreme += 1;
        }
    }

    let p_value = (n_extreme as f64 + 1.0) / (n_permutations as f64 + 1.0);

    Ok(MantelResult {
        statistic: observed,
        p_value,
        n_permutations,
        method: method.to_string(),
    })
}

// ── AMOVA ───────────────────────────────────────────────────────────────────

/// Result of AMOVA.
#[derive(Debug, Clone)]
pub struct AmovaResult {
    /// Sum of squared deviations among groups.
    pub ss_among: f64,
    /// Sum of squared deviations within groups.
    pub ss_within: f64,
    /// Total sum of squared deviations.
    pub ss_total: f64,
    /// Degrees of freedom (among, within, total).
    pub df: (usize, usize, usize),
    /// Mean squares (among, within).
    pub ms: (f64, f64),
    /// F-statistic.
    pub f_statistic: f64,
    /// Permutation p-value.
    pub p_value: f64,
    /// Phi-statistic (proportion of variance among groups).
    pub phi_statistic: f64,
    /// Number of permutations performed.
    pub n_permutations: usize,
}

/// AMOVA: analysis of molecular variance.
///
/// Partitions squared-distance variance into among- and within-group components.
/// Tests significance by permuting group labels.
///
/// # Errors
///
/// Returns an error if the distance matrix is invalid or fewer than 2 groups.
pub fn amova(
    distances: &[Vec<f64>],
    groups: &[usize],
    n_permutations: usize,
    seed: u64,
) -> Result<AmovaResult> {
    let n = distances.len();
    validate_distance_matrix(distances, n)?;
    if groups.len() != n {
        return Err(CyaneaError::InvalidInput(format!(
            "amova: groups length ({}) != distance matrix size ({})",
            groups.len(),
            n
        )));
    }
    if n_permutations == 0 {
        return Err(CyaneaError::InvalidInput(
            "amova: n_permutations must be > 0".into(),
        ));
    }

    let unique_groups: std::collections::HashSet<usize> = groups.iter().copied().collect();
    let k = unique_groups.len();
    if k < 2 {
        return Err(CyaneaError::InvalidInput(
            "amova: at least 2 groups required".into(),
        ));
    }

    let (ss_among, ss_within, ss_total) = compute_amova_ss(distances, groups, n);

    let df_among = k - 1;
    let df_within = n - k;
    let df_total = n - 1;

    let ms_among = if df_among > 0 {
        ss_among / df_among as f64
    } else {
        0.0
    };
    let ms_within = if df_within > 0 {
        ss_within / df_within as f64
    } else {
        0.0
    };

    let f_statistic = if ms_within > 0.0 {
        ms_among / ms_within
    } else {
        0.0
    };

    // Phi_ST = sigma²_among / sigma²_total
    // sigma²_among = (MS_among - MS_within) / n_0  (n_0 is avg group size)
    let group_sizes: Vec<usize> = {
        let mut map: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
        for &g in groups {
            *map.entry(g).or_insert(0) += 1;
        }
        map.values().copied().collect()
    };
    let n0 = {
        let sum_n: usize = group_sizes.iter().sum();
        let sum_n2: usize = group_sizes.iter().map(|&ng| ng * ng).sum();
        if k > 1 {
            (sum_n as f64 - sum_n2 as f64 / sum_n as f64) / (k - 1) as f64
        } else {
            sum_n as f64
        }
    };

    let sigma2_within = ms_within;
    let sigma2_among = if n0 > 0.0 {
        ((ms_among - ms_within) / n0).max(0.0)
    } else {
        0.0
    };
    let sigma2_total = sigma2_among + sigma2_within;
    let phi_statistic = if sigma2_total > 0.0 {
        sigma2_among / sigma2_total
    } else {
        0.0
    };

    // Permutation test
    let mut rng = Xorshift64::new(seed);
    let mut perm_groups: Vec<usize> = groups.to_vec();
    let mut n_extreme = 0usize;

    for _ in 0..n_permutations {
        fisher_yates_shuffle_usize(&mut perm_groups, &mut rng);
        let (perm_ss_among, perm_ss_within, _) =
            compute_amova_ss(distances, &perm_groups, n);
        let perm_ms_among = if df_among > 0 {
            perm_ss_among / df_among as f64
        } else {
            0.0
        };
        let perm_ms_within = if df_within > 0 {
            perm_ss_within / df_within as f64
        } else {
            0.0
        };
        let perm_f = if perm_ms_within > 0.0 {
            perm_ms_among / perm_ms_within
        } else {
            0.0
        };
        if perm_f >= f_statistic {
            n_extreme += 1;
        }
    }

    let p_value = (n_extreme as f64 + 1.0) / (n_permutations as f64 + 1.0);

    Ok(AmovaResult {
        ss_among,
        ss_within,
        ss_total,
        df: (df_among, df_within, df_total),
        ms: (ms_among, ms_within),
        f_statistic,
        p_value,
        phi_statistic,
        n_permutations,
    })
}

fn compute_amova_ss(
    distances: &[Vec<f64>],
    groups: &[usize],
    n: usize,
) -> (f64, f64, f64) {
    // SSD_total = (1/n) Σ_{i<j} d²_{ij}
    let mut ssd_total = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            ssd_total += distances[i][j] * distances[i][j];
        }
    }
    let ss_total = ssd_total / n as f64;

    // SSD_within = Σ_g (1/n_g) Σ_{i<j in g} d²_{ij}
    let mut group_map: std::collections::HashMap<usize, Vec<usize>> =
        std::collections::HashMap::new();
    for (i, &g) in groups.iter().enumerate() {
        group_map.entry(g).or_default().push(i);
    }

    let mut ss_within = 0.0;
    for members in group_map.values() {
        let ng = members.len();
        if ng < 2 {
            continue;
        }
        let mut ss_g = 0.0;
        for ii in 0..ng {
            for jj in (ii + 1)..ng {
                let d = distances[members[ii]][members[jj]];
                ss_g += d * d;
            }
        }
        ss_within += ss_g / ng as f64;
    }

    let ss_among = ss_total - ss_within;
    (ss_among, ss_within, ss_total)
}

// ── BIOENV ──────────────────────────────────────────────────────────────────

/// Result of BIOENV analysis.
#[derive(Debug, Clone)]
pub struct BioenvResult {
    /// Indices of the best subset of environmental variables (0-based).
    pub best_variables: Vec<usize>,
    /// Spearman correlation for the best subset.
    pub best_correlation: f64,
    /// All subsets tested with their correlations: `(variable indices, correlation)`.
    pub all_results: Vec<(Vec<usize>, f64)>,
}

/// BIOENV: find the best subset of environmental variables explaining
/// community distance structure.
///
/// For each subset of environmental variables (sizes 1 to `max_vars`), computes
/// a Euclidean distance matrix and correlates it with the community distance
/// matrix using Spearman correlation. Returns the best-correlated subset.
///
/// # Arguments
///
/// * `community_distances` — pairwise community distance matrix (n × n)
/// * `env_variables` — flat row-major environmental data (n_samples × n_vars)
/// * `n_samples` — number of samples
/// * `n_vars` — number of environmental variables
/// * `max_vars` — maximum subset size to test (capped at `n_vars`)
///
/// # Errors
///
/// Returns an error if dimensions are inconsistent.
pub fn bioenv(
    community_distances: &[Vec<f64>],
    env_variables: &[f64],
    n_samples: usize,
    n_vars: usize,
    max_vars: usize,
) -> Result<BioenvResult> {
    let n = community_distances.len();
    validate_distance_matrix(community_distances, n)?;
    if n != n_samples {
        return Err(CyaneaError::InvalidInput(format!(
            "bioenv: community_distances size ({}) != n_samples ({})",
            n, n_samples
        )));
    }
    if env_variables.len() != n_samples * n_vars {
        return Err(CyaneaError::InvalidInput(format!(
            "bioenv: env_variables length ({}) != n_samples ({}) * n_vars ({})",
            env_variables.len(),
            n_samples,
            n_vars
        )));
    }
    if n_vars == 0 {
        return Err(CyaneaError::InvalidInput(
            "bioenv: n_vars must be > 0".into(),
        ));
    }

    let max_k = max_vars.min(n_vars);

    // Flatten upper triangle of community distances for Spearman
    let comm_flat = upper_triangle_vec(community_distances, n);

    let mut best_corr = f64::NEG_INFINITY;
    let mut best_vars = Vec::new();
    let mut all_results = Vec::new();

    // Enumerate all subsets of size 1..=max_k
    for size in 1..=max_k {
        let subsets = combinations(n_vars, size);
        for subset in subsets {
            // Compute Euclidean distance matrix for this subset
            let env_dists = euclidean_distance_subset(env_variables, n_samples, n_vars, &subset);
            let env_flat = upper_triangle_vec_owned(&env_dists, n);

            // Spearman correlation
            let rho = correlation::spearman(&comm_flat, &env_flat).unwrap_or(0.0);

            all_results.push((subset.clone(), rho));
            if rho > best_corr {
                best_corr = rho;
                best_vars = subset;
            }
        }
    }

    Ok(BioenvResult {
        best_variables: best_vars,
        best_correlation: best_corr,
        all_results,
    })
}

fn euclidean_distance_subset(
    env: &[f64],
    n: usize,
    n_vars: usize,
    vars: &[usize],
) -> Vec<Vec<f64>> {
    let mut dists = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let mut d2 = 0.0;
            for &v in vars {
                let diff = env[i * n_vars + v] - env[j * n_vars + v];
                d2 += diff * diff;
            }
            let d = d2.sqrt();
            dists[i][j] = d;
            dists[j][i] = d;
        }
    }
    dists
}

fn combinations(n: usize, k: usize) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    let mut current = Vec::with_capacity(k);
    fn helper(
        start: usize,
        n: usize,
        k: usize,
        current: &mut Vec<usize>,
        result: &mut Vec<Vec<usize>>,
    ) {
        if current.len() == k {
            result.push(current.clone());
            return;
        }
        for i in start..n {
            current.push(i);
            helper(i + 1, n, k, current, result);
            current.pop();
        }
    }
    helper(0, n, k, &mut current, &mut result);
    result
}

// ── Common helpers ──────────────────────────────────────────────────────────

struct Xorshift64 {
    state: u64,
}

impl Xorshift64 {
    fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 { 1 } else { seed },
        }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    fn next_usize(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}

fn fisher_yates_shuffle_usize(slice: &mut [usize], rng: &mut Xorshift64) {
    let n = slice.len();
    for i in (1..n).rev() {
        let j = rng.next_usize(i + 1);
        slice.swap(i, j);
    }
}

fn validate_distance_matrix(distances: &[Vec<f64>], n: usize) -> Result<()> {
    if n < 2 {
        return Err(CyaneaError::InvalidInput(
            "distance matrix must have at least 2 samples".into(),
        ));
    }
    for row in distances {
        if row.len() != n {
            return Err(CyaneaError::InvalidInput(
                "distance matrix must be square".into(),
            ));
        }
    }
    Ok(())
}

fn upper_triangle_vec(matrix: &[Vec<f64>], n: usize) -> Vec<f64> {
    let mut v = Vec::with_capacity(n * (n - 1) / 2);
    for i in 0..n {
        for j in (i + 1)..n {
            v.push(matrix[i][j]);
        }
    }
    v
}

fn upper_triangle_vec_owned(matrix: &[Vec<f64>], n: usize) -> Vec<f64> {
    upper_triangle_vec(matrix, n)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_grouped_distances() -> (Vec<Vec<f64>>, Vec<usize>) {
        // Two well-separated groups: {0..5} and {6..11}
        // Within-group distances small, between-group distances large
        let n = 12;
        let mut dists = vec![vec![0.0; n]; n];
        // Within group 0: small distances
        for i in 0..6 {
            for j in (i + 1)..6 {
                dists[i][j] = 0.1;
                dists[j][i] = 0.1;
            }
        }
        // Within group 1: small distances
        for i in 6..12 {
            for j in (i + 1)..12 {
                dists[i][j] = 0.1;
                dists[j][i] = 0.1;
            }
        }
        // Between groups: large distances
        for i in 0..6 {
            for j in 6..12 {
                dists[i][j] = 5.0;
                dists[j][i] = 5.0;
            }
        }
        let groups = vec![0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1];
        (dists, groups)
    }

    fn make_uniform_distances() -> (Vec<Vec<f64>>, Vec<usize>) {
        // All pairwise distances equal → no group structure
        let n = 12;
        let mut dists = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in (i + 1)..n {
                dists[i][j] = 1.0;
                dists[j][i] = 1.0;
            }
        }
        let groups = vec![0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1];
        (dists, groups)
    }

    // ── PERMANOVA tests ──────────────────────────────────────────────

    #[test]
    fn permanova_separated_groups_significant() {
        let (dists, groups) = make_grouped_distances();
        let result = permanova(&dists, &groups, 999, 42).unwrap();
        assert!(result.f_statistic > 1.0, "F={}", result.f_statistic);
        assert!(result.p_value < 0.05, "p={}", result.p_value);
        assert!(result.r_squared > 0.5, "R²={}", result.r_squared);
    }

    #[test]
    fn permanova_uniform_not_significant() {
        let (dists, groups) = make_uniform_distances();
        let result = permanova(&dists, &groups, 999, 42).unwrap();
        assert!(result.p_value > 0.05, "p={}", result.p_value);
    }

    #[test]
    fn permanova_r_squared_in_range() {
        let (dists, groups) = make_grouped_distances();
        let result = permanova(&dists, &groups, 99, 42).unwrap();
        assert!(result.r_squared >= 0.0 && result.r_squared <= 1.0);
    }

    #[test]
    fn permanova_too_few_groups_error() {
        let dists = vec![vec![0.0, 1.0], vec![1.0, 0.0]];
        let groups = vec![0, 0]; // Only 1 group
        assert!(permanova(&dists, &groups, 99, 42).is_err());
    }

    #[test]
    fn permanova_mismatched_error() {
        let dists = vec![vec![0.0, 1.0], vec![1.0, 0.0]];
        let groups = vec![0, 1, 2]; // Wrong length
        assert!(permanova(&dists, &groups, 99, 42).is_err());
    }

    // ── ANOSIM tests ─────────────────────────────────────────────────

    #[test]
    fn anosim_separated_groups() {
        let (dists, groups) = make_grouped_distances();
        let result = anosim(&dists, &groups, 999, 42).unwrap();
        assert!(result.r_statistic > 0.0, "R={}", result.r_statistic);
        assert!(result.p_value < 0.05, "p={}", result.p_value);
    }

    #[test]
    fn anosim_r_in_valid_range() {
        let (dists, groups) = make_grouped_distances();
        let result = anosim(&dists, &groups, 99, 42).unwrap();
        assert!(
            result.r_statistic >= -1.0 && result.r_statistic <= 1.0,
            "R={}",
            result.r_statistic
        );
    }

    #[test]
    fn anosim_uniform_not_significant() {
        let (dists, groups) = make_uniform_distances();
        let result = anosim(&dists, &groups, 999, 42).unwrap();
        // R should be near 0 (no group structure in uniform distances)
        assert!(result.r_statistic.abs() < 0.5, "R={}", result.r_statistic);
    }

    #[test]
    fn anosim_too_few_groups_error() {
        let dists = vec![vec![0.0, 1.0], vec![1.0, 0.0]];
        assert!(anosim(&dists, &[0, 0], 99, 42).is_err());
    }

    // ── Mantel test ──────────────────────────────────────────────────

    #[test]
    fn mantel_identical_matrices_correlation_one() {
        let mat = vec![
            vec![0.0, 1.0, 2.0, 3.0],
            vec![1.0, 0.0, 1.5, 2.5],
            vec![2.0, 1.5, 0.0, 1.0],
            vec![3.0, 2.5, 1.0, 0.0],
        ];
        let result = mantel_test(&mat, &mat, 99, 42, "pearson").unwrap();
        assert!(
            (result.statistic - 1.0).abs() < 1e-10,
            "r={}",
            result.statistic
        );
    }

    #[test]
    fn mantel_pearson_vs_spearman() {
        let mat_a = vec![
            vec![0.0, 1.0, 3.0],
            vec![1.0, 0.0, 2.0],
            vec![3.0, 2.0, 0.0],
        ];
        let mat_b = vec![
            vec![0.0, 1.0, 3.0],
            vec![1.0, 0.0, 2.0],
            vec![3.0, 2.0, 0.0],
        ];
        let r_p = mantel_test(&mat_a, &mat_b, 99, 42, "pearson").unwrap();
        let r_s = mantel_test(&mat_a, &mat_b, 99, 42, "spearman").unwrap();
        // Both should be high for identical matrices
        assert!(r_p.statistic > 0.9);
        assert!(r_s.statistic > 0.9);
    }

    #[test]
    fn mantel_invalid_method_error() {
        let mat = vec![vec![0.0, 1.0], vec![1.0, 0.0]];
        assert!(mantel_test(&mat, &mat, 99, 42, "invalid").is_err());
    }

    #[test]
    fn mantel_size_mismatch_error() {
        let a = vec![vec![0.0, 1.0], vec![1.0, 0.0]];
        let b = vec![vec![0.0, 1.0, 2.0], vec![1.0, 0.0, 1.5], vec![2.0, 1.5, 0.0]];
        assert!(mantel_test(&a, &b, 99, 42, "pearson").is_err());
    }

    // ── AMOVA tests ──────────────────────────────────────────────────

    #[test]
    fn amova_ss_additivity() {
        let (dists, groups) = make_grouped_distances();
        let result = amova(&dists, &groups, 99, 42).unwrap();
        let ss_sum = result.ss_among + result.ss_within;
        assert!(
            (ss_sum - result.ss_total).abs() < 1e-10,
            "SS_among + SS_within = {} != SS_total = {}",
            ss_sum,
            result.ss_total
        );
    }

    #[test]
    fn amova_separated_groups_significant() {
        let (dists, groups) = make_grouped_distances();
        let result = amova(&dists, &groups, 999, 42).unwrap();
        assert!(result.f_statistic > 1.0, "F={}", result.f_statistic);
        assert!(result.p_value < 0.05, "p={}", result.p_value);
        assert!(result.phi_statistic > 0.0, "Phi={}", result.phi_statistic);
    }

    #[test]
    fn amova_phi_in_range() {
        let (dists, groups) = make_grouped_distances();
        let result = amova(&dists, &groups, 99, 42).unwrap();
        assert!(
            result.phi_statistic >= 0.0 && result.phi_statistic <= 1.0,
            "Phi={}",
            result.phi_statistic
        );
    }

    #[test]
    fn amova_too_few_groups_error() {
        let dists = vec![vec![0.0, 1.0], vec![1.0, 0.0]];
        assert!(amova(&dists, &[0, 0], 99, 42).is_err());
    }

    // ── BIOENV tests ─────────────────────────────────────────────────

    #[test]
    fn bioenv_finds_correlated_variable() {
        // Community distances correlated with env variable 0
        let comm_dists = vec![
            vec![0.0, 1.0, 2.0, 3.0],
            vec![1.0, 0.0, 1.0, 2.0],
            vec![2.0, 1.0, 0.0, 1.0],
            vec![3.0, 2.0, 1.0, 0.0],
        ];
        // Env var 0: gradient matching community; var 1: noise
        let env = vec![
            0.0, 5.0, // sample 0
            1.0, 3.0, // sample 1
            2.0, 7.0, // sample 2
            3.0, 1.0, // sample 3
        ];
        let result = bioenv(&comm_dists, &env, 4, 2, 2).unwrap();
        // Best subset should include variable 0
        assert!(
            result.best_variables.contains(&0),
            "best_vars={:?}",
            result.best_variables
        );
        assert!(result.best_correlation > 0.0);
    }

    #[test]
    fn bioenv_single_variable() {
        let comm_dists = vec![
            vec![0.0, 1.0, 2.0],
            vec![1.0, 0.0, 1.0],
            vec![2.0, 1.0, 0.0],
        ];
        let env = vec![0.0, 1.0, 2.0];
        let result = bioenv(&comm_dists, &env, 3, 1, 1).unwrap();
        assert_eq!(result.best_variables, vec![0]);
        assert_eq!(result.all_results.len(), 1);
    }

    #[test]
    fn bioenv_dimension_mismatch_error() {
        let comm_dists = vec![vec![0.0, 1.0], vec![1.0, 0.0]];
        let env = vec![1.0]; // Wrong size
        assert!(bioenv(&comm_dists, &env, 2, 1, 1).is_err());
    }

    #[test]
    fn bioenv_zero_vars_error() {
        let comm_dists = vec![vec![0.0, 1.0], vec![1.0, 0.0]];
        assert!(bioenv(&comm_dists, &[], 2, 0, 1).is_err());
    }
}
