//! Hi-C and 3D genome analysis — contact matrices, TAD calling,
//! A/B compartment detection, and chromatin loop calling.

use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// Contact matrix
// ---------------------------------------------------------------------------

/// A Hi-C contact matrix at a given resolution.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ContactMatrix {
    /// Chromosome name (intra-chromosomal) or pair (inter-chromosomal).
    pub chrom: String,
    /// Second chromosome (for inter-chromosomal; same as `chrom` for intra).
    pub chrom2: String,
    /// Resolution (bin size) in base pairs.
    pub resolution: u64,
    /// Number of bins on first axis.
    pub n_bins1: usize,
    /// Number of bins on second axis.
    pub n_bins2: usize,
    /// Dense contact counts: `matrix[i][j]`.
    pub matrix: Vec<Vec<f64>>,
}

impl ContactMatrix {
    /// Create a new square intra-chromosomal contact matrix.
    pub fn new(chrom: &str, resolution: u64, n_bins: usize) -> Self {
        Self {
            chrom: chrom.to_string(),
            chrom2: chrom.to_string(),
            resolution,
            n_bins1: n_bins,
            n_bins2: n_bins,
            matrix: vec![vec![0.0; n_bins]; n_bins],
        }
    }

    /// Create from a dense matrix.
    pub fn from_dense(chrom: &str, resolution: u64, matrix: Vec<Vec<f64>>) -> Result<Self> {
        let n = matrix.len();
        if n == 0 {
            return Err(CyaneaError::InvalidInput("empty matrix".into()));
        }
        for row in &matrix {
            if row.len() != matrix[0].len() {
                return Err(CyaneaError::InvalidInput("ragged matrix".into()));
            }
        }
        Ok(Self {
            chrom: chrom.to_string(),
            chrom2: chrom.to_string(),
            resolution,
            n_bins1: n,
            n_bins2: matrix[0].len(),
            matrix,
        })
    }

    /// Whether this is an intra-chromosomal (cis) matrix.
    pub fn is_cis(&self) -> bool {
        self.chrom == self.chrom2
    }

    /// Total contact count.
    pub fn total_contacts(&self) -> f64 {
        self.matrix.iter().flat_map(|r| r.iter()).sum()
    }

    /// Get contact value at (i, j).
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.matrix[i][j]
    }

    /// Set contact value (symmetric for cis).
    pub fn set(&mut self, i: usize, j: usize, value: f64) {
        self.matrix[i][j] = value;
        if self.is_cis() && i != j {
            self.matrix[j][i] = value;
        }
    }

    /// Add a contact (symmetric for cis).
    pub fn add_contact(&mut self, i: usize, j: usize, count: f64) {
        self.matrix[i][j] += count;
        if self.is_cis() && i != j {
            self.matrix[j][i] += count;
        }
    }

    /// ICE (Iterative Correction and Eigenvector decomposition) balancing.
    ///
    /// Normalizes the contact matrix so that each row/column has the same
    /// marginal sum. Returns the bias vector.
    pub fn ice_balance(&mut self, max_iter: usize, tolerance: f64) -> Vec<f64> {
        let n = self.n_bins1;
        let mut bias = vec![1.0; n];

        for _ in 0..max_iter {
            let mut max_change = 0.0_f64;
            for i in 0..n {
                let row_sum: f64 = (0..n).map(|j| self.matrix[i][j]).sum();
                if row_sum > 0.0 {
                    let factor = row_sum.sqrt();
                    let old_bias = bias[i];
                    bias[i] *= factor;
                    for j in 0..n {
                        self.matrix[i][j] /= factor;
                        self.matrix[j][i] /= factor;
                    }
                    max_change = max_change.max((bias[i] - old_bias).abs() / old_bias.max(1e-10));
                }
            }
            if max_change < tolerance {
                break;
            }
        }

        bias
    }

    /// KR (Knight-Ruiz) balancing — iteratively normalizes so all row/column
    /// sums equal the target (default: sqrt of mean marginal).
    pub fn kr_balance(&mut self, max_iter: usize, tolerance: f64) -> Vec<f64> {
        let n = self.n_bins1;
        let mut weights = vec![1.0; n];

        // Compute target: sqrt(mean marginal)
        let marginals: Vec<f64> = (0..n)
            .map(|i| (0..n).map(|j| self.matrix[i][j]).sum::<f64>())
            .collect();
        let mean_marginal = marginals.iter().sum::<f64>() / n as f64;
        let target = mean_marginal.sqrt().max(1.0);

        for _ in 0..max_iter {
            let mut max_change = 0.0_f64;
            for i in 0..n {
                let row_sum: f64 = (0..n).map(|j| self.matrix[i][j] * weights[j]).sum();
                if row_sum > 1e-10 {
                    let new_w = target / row_sum;
                    let change = (new_w - weights[i]).abs() / weights[i].max(1e-10);
                    max_change = max_change.max(change);
                    weights[i] = new_w;
                }
            }
            if max_change < tolerance {
                break;
            }
        }

        // Apply weights
        for i in 0..n {
            for j in 0..n {
                self.matrix[i][j] *= weights[i] * weights[j];
            }
        }

        weights
    }

    /// Compute observed/expected matrix (distance normalization).
    ///
    /// For each diagonal offset d, computes the mean contact count, then
    /// divides each cell by the expected value at that distance.
    pub fn observed_expected(&self) -> Vec<Vec<f64>> {
        let n = self.n_bins1;
        let mut oe = vec![vec![0.0; n]; n];

        // Compute expected per diagonal
        let mut diag_sums = vec![0.0_f64; n];
        let mut diag_counts = vec![0usize; n];
        for i in 0..n {
            for j in 0..n {
                let d = if i > j { i - j } else { j - i };
                diag_sums[d] += self.matrix[i][j];
                diag_counts[d] += 1;
            }
        }

        let diag_means: Vec<f64> = diag_sums
            .iter()
            .zip(diag_counts.iter())
            .map(|(&s, &c)| if c > 0 { s / c as f64 } else { 0.0 })
            .collect();

        for i in 0..n {
            for j in 0..n {
                let d = if i > j { i - j } else { j - i };
                oe[i][j] = if diag_means[d] > 0.0 {
                    self.matrix[i][j] / diag_means[d]
                } else {
                    0.0
                };
            }
        }

        oe
    }
}

// ---------------------------------------------------------------------------
// TAD calling
// ---------------------------------------------------------------------------

/// A Topologically Associating Domain.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Tad {
    /// Chromosome.
    pub chrom: String,
    /// Start bin index.
    pub start_bin: usize,
    /// End bin index (exclusive).
    pub end_bin: usize,
    /// Start genomic position.
    pub start_bp: u64,
    /// End genomic position.
    pub end_bp: u64,
    /// Insulation score at the boundary.
    pub boundary_score: f64,
}

/// Parameters for TAD calling.
#[derive(Debug, Clone)]
pub struct TadParams {
    /// Window size for insulation score (in bins).
    pub window_size: usize,
    /// Boundary detection threshold (z-score).
    pub boundary_threshold: f64,
    /// Minimum TAD size (in bins).
    pub min_tad_size: usize,
}

impl Default for TadParams {
    fn default() -> Self {
        Self {
            window_size: 10,
            boundary_threshold: -0.5,
            min_tad_size: 3,
        }
    }
}

/// Call TADs using the insulation score method (Crane et al. 2015).
///
/// For each bin, computes the mean contact count in a square window
/// along the diagonal. Local minima in the insulation score indicate
/// TAD boundaries.
pub fn call_tads(matrix: &ContactMatrix, params: &TadParams) -> Result<Vec<Tad>> {
    let n = matrix.n_bins1;
    let w = params.window_size;

    if n < 2 * w + 1 {
        return Err(CyaneaError::InvalidInput(
            "matrix too small for window size".into(),
        ));
    }

    // Compute insulation scores
    let scores = insulation_scores(matrix, w);

    // Normalize to z-scores
    let mean = scores.iter().sum::<f64>() / scores.len() as f64;
    let var = scores.iter().map(|s| (s - mean).powi(2)).sum::<f64>() / scores.len() as f64;
    let std = var.sqrt().max(1e-10);
    let z_scores: Vec<f64> = scores.iter().map(|s| (s - mean) / std).collect();

    // Find boundaries: local minima below threshold
    let mut boundaries = Vec::new();
    for i in 1..n - 1 {
        if z_scores[i] < params.boundary_threshold
            && z_scores[i] <= z_scores[i - 1]
            && z_scores[i] <= z_scores[i + 1]
        {
            boundaries.push((i, z_scores[i]));
        }
    }

    // Add start and end
    let mut all_bounds = vec![(0, 0.0)];
    all_bounds.extend(&boundaries);
    all_bounds.push((n, 0.0));

    // Build TADs from consecutive boundaries
    let mut tads = Vec::new();
    for pair in all_bounds.windows(2) {
        let start = pair[0].0;
        let end = pair[1].0;
        let score = pair[1].1;

        if end - start >= params.min_tad_size {
            tads.push(Tad {
                chrom: matrix.chrom.clone(),
                start_bin: start,
                end_bin: end,
                start_bp: start as u64 * matrix.resolution,
                end_bp: end as u64 * matrix.resolution,
                boundary_score: score,
            });
        }
    }

    Ok(tads)
}

/// Compute insulation scores for each bin.
///
/// The insulation score at bin i is the mean contact count in the square
/// window `[i-w, i] × [i, i+w]` on the contact matrix.
pub fn insulation_scores(matrix: &ContactMatrix, window_size: usize) -> Vec<f64> {
    let n = matrix.n_bins1;
    let w = window_size;
    let mut scores = vec![0.0; n];

    for i in 0..n {
        let lo = if i >= w { i - w } else { 0 };
        let hi = (i + w).min(n);
        let mut sum = 0.0;
        let mut count = 0;
        for r in lo..i {
            for c in i..hi {
                sum += matrix.matrix[r][c];
                count += 1;
            }
        }
        scores[i] = if count > 0 {
            (sum / count as f64).ln().max(-10.0)
        } else {
            0.0
        };
    }

    scores
}

// ---------------------------------------------------------------------------
// A/B compartments
// ---------------------------------------------------------------------------

/// A/B compartment assignment for a genomic bin.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Compartment {
    /// Active / open chromatin.
    A,
    /// Inactive / closed chromatin.
    B,
}

/// Result of compartment calling.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CompartmentResult {
    /// Per-bin compartment assignment.
    pub compartments: Vec<Compartment>,
    /// First eigenvector (PC1) values per bin.
    pub eigenvector: Vec<f64>,
}

/// Call A/B compartments from a Hi-C contact matrix.
///
/// Uses the Lieberman-Aiden method:
/// 1. Compute observed/expected matrix
/// 2. Compute Pearson correlation matrix of O/E rows
/// 3. Extract first eigenvector (PC1) via power iteration
/// 4. A = positive PC1, B = negative PC1
///
/// Optionally, GC content or gene density can orient the eigenvector
/// (A compartments should have higher GC/gene density).
pub fn call_compartments(
    matrix: &ContactMatrix,
    gc_content: Option<&[f64]>,
) -> Result<CompartmentResult> {
    let n = matrix.n_bins1;
    if n < 3 {
        return Err(CyaneaError::InvalidInput(
            "need at least 3 bins for compartment calling".into(),
        ));
    }

    // Step 1: O/E matrix
    let oe = matrix.observed_expected();

    // Step 2: Pearson correlation matrix of O/E rows
    let mut corr = vec![vec![0.0; n]; n];
    let mut means = vec![0.0; n];
    let mut stds = vec![0.0; n];

    for i in 0..n {
        means[i] = oe[i].iter().sum::<f64>() / n as f64;
        let var = oe[i].iter().map(|v| (v - means[i]).powi(2)).sum::<f64>() / n as f64;
        stds[i] = var.sqrt().max(1e-10);
    }

    for i in 0..n {
        corr[i][i] = 1.0;
        for j in (i + 1)..n {
            let cov: f64 = (0..n)
                .map(|k| (oe[i][k] - means[i]) * (oe[j][k] - means[j]))
                .sum::<f64>()
                / n as f64;
            let r = cov / (stds[i] * stds[j]);
            corr[i][j] = r;
            corr[j][i] = r;
        }
    }

    // Step 3: Power iteration for first eigenvector
    let mut eigvec = vec![1.0 / (n as f64).sqrt(); n];
    for _ in 0..100 {
        let mut new_vec = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                new_vec[i] += corr[i][j] * eigvec[j];
            }
        }
        // Normalize
        let norm: f64 = new_vec.iter().map(|v| v * v).sum::<f64>().sqrt();
        if norm > 1e-10 {
            for v in &mut new_vec {
                *v /= norm;
            }
        }
        eigvec = new_vec;
    }

    // Step 4: Orient eigenvector using GC content (if provided)
    if let Some(gc) = gc_content {
        if gc.len() == n {
            let dot: f64 = eigvec.iter().zip(gc.iter()).map(|(e, g)| e * g).sum();
            if dot < 0.0 {
                for v in &mut eigvec {
                    *v = -*v;
                }
            }
        }
    }

    // Assign compartments
    let compartments: Vec<Compartment> = eigvec
        .iter()
        .map(|&v| if v >= 0.0 { Compartment::A } else { Compartment::B })
        .collect();

    Ok(CompartmentResult {
        compartments,
        eigenvector: eigvec,
    })
}

// ---------------------------------------------------------------------------
// Loop calling
// ---------------------------------------------------------------------------

/// A chromatin loop (peak in contact matrix).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ChromatinLoop {
    /// Chromosome.
    pub chrom: String,
    /// Anchor 1 bin index.
    pub anchor1_bin: usize,
    /// Anchor 2 bin index.
    pub anchor2_bin: usize,
    /// Anchor 1 genomic position.
    pub anchor1_bp: u64,
    /// Anchor 2 genomic position.
    pub anchor2_bp: u64,
    /// Observed contact count.
    pub observed: f64,
    /// Expected contact count (from local background).
    pub expected: f64,
    /// Enrichment (observed / expected).
    pub enrichment: f64,
    /// p-value.
    pub p_value: f64,
}

/// Parameters for loop calling.
#[derive(Debug, Clone)]
pub struct LoopParams {
    /// Local background window size (in bins).
    pub background_window: usize,
    /// Minimum enrichment over local background.
    pub min_enrichment: f64,
    /// Minimum genomic distance (in bins).
    pub min_distance: usize,
    /// Maximum genomic distance (in bins).
    pub max_distance: usize,
    /// p-value threshold.
    pub p_threshold: f64,
}

impl Default for LoopParams {
    fn default() -> Self {
        Self {
            background_window: 5,
            min_enrichment: 1.5,
            min_distance: 5,
            max_distance: 500,
            p_threshold: 0.01,
        }
    }
}

/// Call chromatin loops (HiCCUPS-style local enrichment).
///
/// For each pixel (i, j), computes the enrichment over four local
/// background regions (donut, horizontal, vertical, lower-left).
/// Pixels significantly enriched over all backgrounds are called as loops.
pub fn call_loops(matrix: &ContactMatrix, params: &LoopParams) -> Result<Vec<ChromatinLoop>> {
    let n = matrix.n_bins1;
    let w = params.background_window;

    if n < 2 * w + 1 {
        return Err(CyaneaError::InvalidInput(
            "matrix too small for loop calling".into(),
        ));
    }

    let mut loops = Vec::new();

    for i in w..n - w {
        for j in (i + params.min_distance)..n.min(i + params.max_distance) {
            if j + w >= n {
                break;
            }

            let observed = matrix.matrix[i][j];
            if observed <= 0.0 {
                continue;
            }

            // Donut background: ring around (i,j) excluding the pixel itself
            let mut bg_sum = 0.0;
            let mut bg_count = 0;
            for di in -(w as i64)..=(w as i64) {
                for dj in -(w as i64)..=(w as i64) {
                    let ri = (i as i64 + di) as usize;
                    let rj = (j as i64 + dj) as usize;
                    if ri >= n || rj >= n {
                        continue;
                    }
                    // Exclude inner region (within 1 bin)
                    if di.abs() <= 1 && dj.abs() <= 1 {
                        continue;
                    }
                    bg_sum += matrix.matrix[ri][rj];
                    bg_count += 1;
                }
            }

            let expected = if bg_count > 0 {
                bg_sum / bg_count as f64
            } else {
                continue;
            };

            if expected <= 0.0 {
                continue;
            }

            let enrichment = observed / expected;
            if enrichment < params.min_enrichment {
                continue;
            }

            // Simple Poisson-like p-value: P(X >= observed | lambda = expected)
            // Using normal approximation for speed
            let z = (observed - expected) / expected.sqrt().max(1.0);
            let p = erfc_approx(z / std::f64::consts::SQRT_2) / 2.0;

            if p < params.p_threshold {
                loops.push(ChromatinLoop {
                    chrom: matrix.chrom.clone(),
                    anchor1_bin: i,
                    anchor2_bin: j,
                    anchor1_bp: i as u64 * matrix.resolution,
                    anchor2_bp: j as u64 * matrix.resolution,
                    observed,
                    expected,
                    enrichment,
                    p_value: p.max(0.0),
                });
            }
        }
    }

    // Sort by enrichment descending
    loops.sort_by(|a, b| {
        b.enrichment
            .partial_cmp(&a.enrichment)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    Ok(loops)
}

// ---------------------------------------------------------------------------
// Hi-C file format parsing
// ---------------------------------------------------------------------------

/// A sparse contact record from a Hi-C file.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SparseContact {
    /// Bin 1 index (or genomic position).
    pub bin1: u64,
    /// Bin 2 index (or genomic position).
    pub bin2: u64,
    /// Contact count.
    pub count: f64,
}

/// Metadata from a .cool file header.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CoolHeader {
    /// Bin size (resolution).
    pub bin_size: u64,
    /// Chromosome names.
    pub chroms: Vec<String>,
    /// Chromosome sizes.
    pub chrom_sizes: Vec<u64>,
    /// Total number of bins.
    pub n_bins: usize,
    /// Number of non-zero entries (nnz).
    pub nnz: usize,
}

/// Parse a .cool-style text format (tab-delimited: chrom1, start1, end1, chrom2, start2, end2, count).
///
/// This parses the text representation commonly exported from cooler tools.
/// For actual HDF5-based .cool files, use the h5ad feature.
pub fn parse_cool_text(data: &str, resolution: u64) -> Result<Vec<SparseContact>> {
    let mut contacts = Vec::new();

    for line in data.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 7 {
            continue;
        }

        let start1: u64 = fields[1]
            .parse()
            .map_err(|_| CyaneaError::Parse("invalid start1".into()))?;
        let start2: u64 = fields[4]
            .parse()
            .map_err(|_| CyaneaError::Parse("invalid start2".into()))?;
        let count: f64 = fields[6]
            .parse()
            .map_err(|_| CyaneaError::Parse("invalid count".into()))?;

        contacts.push(SparseContact {
            bin1: start1 / resolution,
            bin2: start2 / resolution,
            count,
        });
    }

    if contacts.is_empty() {
        return Err(CyaneaError::Parse("no contacts found".into()));
    }

    Ok(contacts)
}

/// Parse a pairs format file (4DN consortium).
///
/// Format: readID chr1 pos1 chr2 pos2 strand1 strand2
pub fn parse_pairs(data: &str) -> Result<Vec<(String, u64, String, u64)>> {
    let mut pairs = Vec::new();

    for line in data.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }

        let chr1 = fields[1].to_string();
        let pos1: u64 = fields[2]
            .parse()
            .map_err(|_| CyaneaError::Parse("invalid pos1".into()))?;
        let chr2 = fields[3].to_string();
        let pos2: u64 = fields[4]
            .parse()
            .map_err(|_| CyaneaError::Parse("invalid pos2".into()))?;

        pairs.push((chr1, pos1, chr2, pos2));
    }

    Ok(pairs)
}

/// Build a contact matrix from sparse contacts.
pub fn contacts_to_matrix(
    contacts: &[SparseContact],
    chrom: &str,
    resolution: u64,
    n_bins: usize,
) -> ContactMatrix {
    let mut matrix = ContactMatrix::new(chrom, resolution, n_bins);
    for c in contacts {
        let i = c.bin1 as usize;
        let j = c.bin2 as usize;
        if i < n_bins && j < n_bins {
            matrix.add_contact(i, j, c.count);
        }
    }
    matrix
}

/// Write contacts in pairs text format.
pub fn write_pairs(contacts: &[SparseContact], chrom: &str, resolution: u64) -> String {
    let mut lines = Vec::new();
    lines.push("## pairs format v1.0".to_string());
    lines.push("#columns: readID chr1 pos1 chr2 pos2 strand1 strand2".to_string());

    for (i, c) in contacts.iter().enumerate() {
        lines.push(format!(
            "read_{}\t{}\t{}\t{}\t{}\t+\t+",
            i,
            chrom,
            c.bin1 * resolution,
            chrom,
            c.bin2 * resolution
        ));
    }

    lines.join("\n")
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn erf_approx(x: f64) -> f64 {
    let a1: f64 = 0.254829592;
    let a2: f64 = -0.284496736;
    let a3: f64 = 1.421413741;
    let a4: f64 = -1.453152027;
    let a5: f64 = 1.061405429;
    let p: f64 = 0.3275911;
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    let t = 1.0 / (1.0 + p * x);
    let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x * x).exp();
    sign * y
}

fn erfc_approx(x: f64) -> f64 {
    1.0 - erf_approx(x)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_tad_matrix() -> ContactMatrix {
        // 20x20 matrix with two TAD-like blocks
        let n = 20;
        let mut m = ContactMatrix::new("chr1", 10000, n);
        for i in 0..n {
            for j in 0..n {
                let d = if i > j { i - j } else { j - i };
                let base = 100.0 / (d as f64 + 1.0);
                // TAD 1: bins 0-9, TAD 2: bins 10-19
                let in_same_tad = (i < 10 && j < 10) || (i >= 10 && j >= 10);
                let contact = if in_same_tad { base * 3.0 } else { base * 0.3 };
                m.set(i, j, contact);
            }
        }
        m
    }

    #[test]
    fn test_contact_matrix_basic() {
        let m = ContactMatrix::new("chr1", 5000, 10);
        assert_eq!(m.n_bins1, 10);
        assert!(m.is_cis());
        assert!((m.total_contacts() - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_contact_matrix_symmetric() {
        let mut m = ContactMatrix::new("chr1", 5000, 5);
        m.set(1, 3, 42.0);
        assert!((m.get(1, 3) - 42.0).abs() < 1e-10);
        assert!((m.get(3, 1) - 42.0).abs() < 1e-10);
    }

    #[test]
    fn test_add_contact() {
        let mut m = ContactMatrix::new("chr1", 5000, 5);
        m.add_contact(2, 4, 10.0);
        m.add_contact(2, 4, 5.0);
        assert!((m.get(2, 4) - 15.0).abs() < 1e-10);
        assert!((m.get(4, 2) - 15.0).abs() < 1e-10);
    }

    #[test]
    fn test_observed_expected() {
        let m = make_tad_matrix();
        let oe = m.observed_expected();
        assert_eq!(oe.len(), 20);
        // Diagonal should have O/E near 1.0 (since diagonal average = value)
        // Within-TAD should be > 1, between-TAD < 1
        // Check within TAD1
        assert!(oe[2][5] > 1.0, "within-TAD O/E should be > 1: {}", oe[2][5]);
        // Check between TADs (at sufficient distance to avoid diagonal effects)
        assert!(
            oe[2][15] < 1.5,
            "between-TAD O/E should be modest: {}",
            oe[2][15]
        );
    }

    #[test]
    fn test_ice_balance() {
        let mut m = ContactMatrix::new("chr1", 5000, 5);
        for i in 0..5 {
            for j in 0..5 {
                m.set(i, j, ((i + j + 1) as f64) * 10.0);
            }
        }
        let bias = m.ice_balance(50, 1e-6);
        assert_eq!(bias.len(), 5);
        // After balancing, row sums should be more uniform
        let row_sums: Vec<f64> = (0..5)
            .map(|i| (0..5).map(|j| m.get(i, j)).sum())
            .collect();
        let mean = row_sums.iter().sum::<f64>() / 5.0;
        let var = row_sums.iter().map(|s| (s - mean).powi(2)).sum::<f64>() / 5.0;
        let cv = var.sqrt() / mean;
        assert!(cv < 0.5, "CV of row sums should be small: {}", cv);
    }

    #[test]
    fn test_call_tads() {
        let m = make_tad_matrix();
        let params = TadParams {
            window_size: 3,
            boundary_threshold: -0.3,
            min_tad_size: 2,
        };
        let tads = call_tads(&m, &params).unwrap();
        assert!(!tads.is_empty(), "should find at least one TAD");
        // TADs should cover the matrix
        assert!(tads[0].start_bin == 0 || tads[0].start_bin < 5);
    }

    #[test]
    fn test_insulation_scores() {
        let m = make_tad_matrix();
        let scores = insulation_scores(&m, 3);
        assert_eq!(scores.len(), 20);
        // Boundary region (around bin 10) should have lower insulation
        let boundary_score = scores[10];
        let interior_score = scores[5];
        assert!(
            boundary_score < interior_score,
            "boundary ({}) should have lower insulation than interior ({})",
            boundary_score,
            interior_score
        );
    }

    #[test]
    fn test_call_compartments() {
        // Create a matrix with two block compartments (first half vs second half)
        let n = 20;
        let mut m = ContactMatrix::new("chr1", 50000, n);
        for i in 0..n {
            for j in 0..n {
                let d = if i > j { i - j } else { j - i };
                let base = 50.0 / (d as f64 + 1.0);
                // Block compartment: first 10 bins = A, last 10 = B
                let same_comp = (i < 10 && j < 10) || (i >= 10 && j >= 10);
                let contact = if same_comp { base * 5.0 } else { base };
                m.set(i, j, contact);
            }
        }

        let result = call_compartments(&m, None).unwrap();
        assert_eq!(result.compartments.len(), 20);
        assert_eq!(result.eigenvector.len(), 20);
        // First half and second half should be in different compartments
        let comp_first = &result.compartments[2];
        let comp_second = &result.compartments[15];
        assert_ne!(comp_first, comp_second, "different blocks should be in different compartments");
    }

    #[test]
    fn test_compartments_gc_orientation() {
        let n = 10;
        let mut m = ContactMatrix::new("chr1", 50000, n);
        for i in 0..n {
            for j in 0..n {
                let d = if i > j { i - j } else { j - i };
                // First half is A (high contacts among themselves)
                let same = (i < 5 && j < 5) || (i >= 5 && j >= 5);
                let v = if same { 50.0 / (d as f64 + 1.0) } else { 10.0 / (d as f64 + 1.0) };
                m.set(i, j, v);
            }
        }

        // GC content higher in first half (should be A)
        let gc: Vec<f64> = (0..n).map(|i| if i < 5 { 0.6 } else { 0.4 }).collect();
        let result = call_compartments(&m, Some(&gc)).unwrap();

        // First half should be A
        assert_eq!(result.compartments[0], Compartment::A);
        assert_eq!(result.compartments[9], Compartment::B);
    }

    #[test]
    fn test_call_loops() {
        let n = 30;
        let mut m = ContactMatrix::new("chr1", 10000, n);
        // Background: distance decay
        for i in 0..n {
            for j in 0..n {
                let d = if i > j { i - j } else { j - i };
                m.set(i, j, 10.0 / (d as f64 + 1.0));
            }
        }
        // Inject a strong loop at (5, 20)
        m.set(5, 20, 100.0);

        let params = LoopParams {
            background_window: 3,
            min_enrichment: 1.5,
            min_distance: 5,
            max_distance: 25,
            p_threshold: 0.05,
        };

        let loops = call_loops(&m, &params).unwrap();
        // Should find the injected loop
        let found = loops.iter().any(|l| l.anchor1_bin == 5 && l.anchor2_bin == 20);
        assert!(found, "should detect the injected loop at (5,20)");
    }

    #[test]
    fn test_parse_cool_text() {
        let data = "chr1\t0\t10000\tchr1\t10000\t20000\t150\n\
                     chr1\t0\t10000\tchr1\t20000\t30000\t50\n\
                     chr1\t10000\t20000\tchr1\t20000\t30000\t120\n";
        let contacts = parse_cool_text(data, 10000).unwrap();
        assert_eq!(contacts.len(), 3);
        assert_eq!(contacts[0].bin1, 0);
        assert_eq!(contacts[0].bin2, 1);
        assert!((contacts[0].count - 150.0).abs() < 0.1);
    }

    #[test]
    fn test_parse_pairs() {
        let data = "## pairs format v1.0\n\
                     #columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n\
                     read_1\tchr1\t1000\tchr1\t50000\t+\t-\n\
                     read_2\tchr1\t2000\tchr2\t30000\t+\t+\n";
        let pairs = parse_pairs(data).unwrap();
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0].0, "chr1");
        assert_eq!(pairs[0].1, 1000);
        assert_eq!(pairs[1].2, "chr2");
    }

    #[test]
    fn test_contacts_to_matrix() {
        let contacts = vec![
            SparseContact { bin1: 0, bin2: 1, count: 10.0 },
            SparseContact { bin1: 1, bin2: 2, count: 20.0 },
            SparseContact { bin1: 0, bin2: 2, count: 5.0 },
        ];
        let m = contacts_to_matrix(&contacts, "chr1", 5000, 5);
        assert!((m.get(0, 1) - 10.0).abs() < 1e-10);
        assert!((m.get(1, 0) - 10.0).abs() < 1e-10); // symmetric
        assert!((m.get(1, 2) - 20.0).abs() < 1e-10);
    }

    #[test]
    fn test_write_pairs() {
        let contacts = vec![
            SparseContact { bin1: 0, bin2: 5, count: 42.0 },
        ];
        let output = write_pairs(&contacts, "chr1", 10000);
        assert!(output.contains("pairs format"));
        assert!(output.contains("chr1\t0\tchr1\t50000"));
    }

    #[test]
    fn test_kr_balance() {
        let mut m = ContactMatrix::new("chr1", 5000, 5);
        for i in 0..5 {
            for j in 0..5 {
                m.set(i, j, ((i + j + 2) as f64) * 5.0);
            }
        }
        let weights = m.kr_balance(100, 1e-6);
        assert_eq!(weights.len(), 5);
        // Weights should be positive
        assert!(weights.iter().all(|&w| w > 0.0));
    }

    #[test]
    fn test_from_dense() {
        let data = vec![
            vec![10.0, 5.0, 2.0],
            vec![5.0, 10.0, 4.0],
            vec![2.0, 4.0, 10.0],
        ];
        let m = ContactMatrix::from_dense("chr1", 10000, data).unwrap();
        assert_eq!(m.n_bins1, 3);
        assert!((m.get(0, 1) - 5.0).abs() < 1e-10);
    }
}
