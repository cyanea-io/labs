//! Multiple sequence alignment via progressive alignment.
//!
//! Implements a simplified ClustalW-style progressive alignment:
//! 1. Compute all pairwise alignment scores
//! 2. Build a UPGMA guide tree from pairwise distances
//! 3. Progressively align sequences/profiles following the guide tree
//!
//! # Example
//!
//! ```
//! use cyanea_align::msa::progressive_msa;
//! use cyanea_align::scoring::{ScoringMatrix, ScoringScheme};
//!
//! let seqs: Vec<&[u8]> = vec![b"ACGT", b"ACGT", b"ACGA"];
//! let scoring = ScoringScheme::Simple(ScoringMatrix::dna_default());
//! let msa = progressive_msa(&seqs, &scoring).unwrap();
//! assert_eq!(msa.n_sequences(), 3);
//! ```

use cyanea_core::{CyaneaError, Result};

use crate::scoring::ScoringScheme;

/// A multiple sequence alignment result.
#[derive(Debug, Clone)]
pub struct MsaResult {
    /// Aligned sequences (with `-` gap characters), all the same length.
    pub aligned: Vec<Vec<u8>>,
    /// Number of columns in the alignment.
    pub n_columns: usize,
}

impl MsaResult {
    /// Number of sequences in the alignment.
    pub fn n_sequences(&self) -> usize {
        self.aligned.len()
    }

    /// Get a column of the alignment (one byte per sequence).
    pub fn column(&self, col: usize) -> Option<Vec<u8>> {
        if col >= self.n_columns {
            return None;
        }
        Some(self.aligned.iter().map(|s| s[col]).collect())
    }

    /// Fraction of columns where all sequences agree (no gaps, identical base).
    pub fn conservation(&self) -> f64 {
        if self.n_columns == 0 || self.aligned.is_empty() {
            return 0.0;
        }
        let conserved = (0..self.n_columns)
            .filter(|&c| {
                let first = self.aligned[0][c];
                first != b'-'
                    && self
                        .aligned
                        .iter()
                        .all(|s| s[c].to_ascii_uppercase() == first.to_ascii_uppercase())
            })
            .count();
        conserved as f64 / self.n_columns as f64
    }
}

/// Perform progressive multiple sequence alignment.
///
/// # Errors
///
/// Returns an error if fewer than 2 sequences are provided or any sequence is empty.
pub fn progressive_msa(sequences: &[&[u8]], scoring: &ScoringScheme) -> Result<MsaResult> {
    let n = sequences.len();
    if n < 2 {
        return Err(CyaneaError::InvalidInput(
            "need at least 2 sequences for MSA".into(),
        ));
    }
    for (i, seq) in sequences.iter().enumerate() {
        if seq.is_empty() {
            return Err(CyaneaError::InvalidInput(format!(
                "sequence {} is empty",
                i
            )));
        }
    }

    // Step 1: pairwise distances
    let mut dist = vec![vec![0.0f64; n]; n];
    let max_pair_score = scoring.score_pair(b'A', b'A') as f64;
    for i in 0..n {
        for j in (i + 1)..n {
            let score = pairwise_nw_score(sequences[i], sequences[j], scoring);
            let max_len = sequences[i].len().max(sequences[j].len()) as f64;
            let d = 1.0 - (score as f64 / (max_len * max_pair_score)).clamp(0.0, 1.0);
            dist[i][j] = d;
            dist[j][i] = d;
        }
    }

    // Step 2: UPGMA guide tree
    let merge_order = upgma_guide_tree(&dist);

    // Step 3: progressive alignment
    let mut profiles: Vec<Option<Vec<Vec<u8>>>> = sequences
        .iter()
        .map(|s| Some(vec![s.to_vec()]))
        .collect();

    for (i, j) in merge_order {
        let prof_i = profiles[i].take().unwrap();
        let prof_j = profiles[j].take().unwrap();
        let merged = align_profiles(&prof_i, &prof_j, scoring);
        profiles[i] = Some(merged);
    }

    let final_profile = profiles.into_iter().find_map(|p| p).unwrap();
    let n_columns = final_profile[0].len();

    Ok(MsaResult {
        aligned: final_profile,
        n_columns,
    })
}

/// Simple Needleman-Wunsch score (no traceback).
fn pairwise_nw_score(a: &[u8], b: &[u8], scoring: &ScoringScheme) -> i32 {
    let m = a.len();
    let n = b.len();
    let mut prev = vec![0i32; n + 1];
    let mut curr = vec![0i32; n + 1];

    for j in 1..=n {
        prev[j] = scoring.gap_open() + (j as i32 - 1) * scoring.gap_extend();
    }

    for i in 1..=m {
        curr[0] = scoring.gap_open() + (i as i32 - 1) * scoring.gap_extend();
        for j in 1..=n {
            let sub = prev[j - 1] + scoring.score_pair(a[i - 1], b[j - 1]);
            let del = prev[j] + scoring.gap_open();
            let ins = curr[j - 1] + scoring.gap_open();
            curr[j] = sub.max(del).max(ins);
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    prev[n]
}

/// Build UPGMA guide tree, returning merge order as `(keep, remove)` pairs.
fn upgma_guide_tree(dist: &[Vec<f64>]) -> Vec<(usize, usize)> {
    let n = dist.len();
    let mut active: Vec<bool> = vec![true; n];
    let mut sizes: Vec<usize> = vec![1; n];
    let mut d = dist.to_vec();
    let mut merges = Vec::new();

    for _ in 0..(n - 1) {
        let mut best_d = f64::MAX;
        let mut bi = 0;
        let mut bj = 0;
        for i in 0..n {
            if !active[i] {
                continue;
            }
            for j in (i + 1)..n {
                if !active[j] && d[i][j] < best_d {
                    // intentionally skip inactive
                }
                if active[j] && d[i][j] < best_d {
                    best_d = d[i][j];
                    bi = i;
                    bj = j;
                }
            }
        }

        merges.push((bi, bj));
        let si = sizes[bi] as f64;
        let sj = sizes[bj] as f64;

        for k in 0..n {
            if !active[k] || k == bi || k == bj {
                continue;
            }
            let new_d = (d[bi][k] * si + d[bj][k] * sj) / (si + sj);
            d[bi][k] = new_d;
            d[k][bi] = new_d;
        }

        sizes[bi] += sizes[bj];
        active[bj] = false;
    }

    merges
}

/// Align two profiles using profile-profile Needleman-Wunsch.
fn align_profiles(
    prof_a: &[Vec<u8>],
    prof_b: &[Vec<u8>],
    scoring: &ScoringScheme,
) -> Vec<Vec<u8>> {
    let m = prof_a[0].len();
    let n = prof_b[0].len();

    let mut h = vec![vec![0i32; n + 1]; m + 1];
    let mut trace = vec![vec![0u8; n + 1]; m + 1];

    for i in 1..=m {
        h[i][0] = scoring.gap_open() + (i as i32 - 1) * scoring.gap_extend();
        trace[i][0] = 1;
    }
    for j in 1..=n {
        h[0][j] = scoring.gap_open() + (j as i32 - 1) * scoring.gap_extend();
        trace[0][j] = 2;
    }

    for i in 1..=m {
        for j in 1..=n {
            let sub = profile_column_score(prof_a, i - 1, prof_b, j - 1, scoring);
            let diag = h[i - 1][j - 1] + sub;
            let up = h[i - 1][j] + scoring.gap_open();
            let left = h[i][j - 1] + scoring.gap_open();

            let best = diag.max(up).max(left);
            h[i][j] = best;
            trace[i][j] = if best == diag {
                0
            } else if best == up {
                1
            } else {
                2
            };
        }
    }

    // Traceback
    let mut cols_a: Vec<Option<usize>> = Vec::new();
    let mut cols_b: Vec<Option<usize>> = Vec::new();
    let mut i = m;
    let mut j = n;

    while i > 0 || j > 0 {
        if i > 0 && j > 0 && trace[i][j] == 0 {
            cols_a.push(Some(i - 1));
            cols_b.push(Some(j - 1));
            i -= 1;
            j -= 1;
        } else if i > 0 && (j == 0 || trace[i][j] == 1) {
            cols_a.push(Some(i - 1));
            cols_b.push(None);
            i -= 1;
        } else {
            cols_a.push(None);
            cols_b.push(Some(j - 1));
            j -= 1;
        }
    }

    cols_a.reverse();
    cols_b.reverse();
    let aln_len = cols_a.len();

    let mut result: Vec<Vec<u8>> = Vec::new();

    for seq in prof_a {
        let mut aligned = Vec::with_capacity(aln_len);
        for ca in &cols_a {
            match ca {
                Some(col) => aligned.push(seq[*col]),
                None => aligned.push(b'-'),
            }
        }
        result.push(aligned);
    }

    for seq in prof_b {
        let mut aligned = Vec::with_capacity(aln_len);
        for cb in &cols_b {
            match cb {
                Some(col) => aligned.push(seq[*col]),
                None => aligned.push(b'-'),
            }
        }
        result.push(aligned);
    }

    result
}

/// Average pairwise score between column `col_a` of profile A and `col_b` of profile B.
fn profile_column_score(
    prof_a: &[Vec<u8>],
    col_a: usize,
    prof_b: &[Vec<u8>],
    col_b: usize,
    scoring: &ScoringScheme,
) -> i32 {
    let mut total = 0i64;
    let mut count = 0i64;
    for sa in prof_a {
        let a = sa[col_a];
        if a == b'-' {
            continue;
        }
        for sb in prof_b {
            let b = sb[col_b];
            if b == b'-' {
                continue;
            }
            total += scoring.score_pair(a, b) as i64;
            count += 1;
        }
    }
    if count == 0 {
        return scoring.gap_open();
    }
    (total / count) as i32
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::{ScoringMatrix, ScoringScheme};

    fn dna_scheme() -> ScoringScheme {
        ScoringScheme::Simple(ScoringMatrix::dna_default())
    }

    #[test]
    fn msa_identical_sequences() {
        let seqs: Vec<&[u8]> = vec![b"ACGT", b"ACGT", b"ACGT"];
        let result = progressive_msa(&seqs, &dna_scheme()).unwrap();
        assert_eq!(result.n_sequences(), 3);
        assert_eq!(result.n_columns, 4);
        assert!((result.conservation() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn msa_two_sequences() {
        let seqs: Vec<&[u8]> = vec![b"ACGT", b"ACGA"];
        let result = progressive_msa(&seqs, &dna_scheme()).unwrap();
        assert_eq!(result.n_sequences(), 2);
        assert!(result.n_columns >= 4);
    }

    #[test]
    fn msa_different_lengths() {
        let seqs: Vec<&[u8]> = vec![b"ACGTACGT", b"ACGT", b"ACGTAC"];
        let result = progressive_msa(&seqs, &dna_scheme()).unwrap();
        assert_eq!(result.n_sequences(), 3);
        let len = result.aligned[0].len();
        assert!(result.aligned.iter().all(|s| s.len() == len));
    }

    #[test]
    fn msa_too_few_sequences() {
        let seqs: Vec<&[u8]> = vec![b"ACGT"];
        assert!(progressive_msa(&seqs, &dna_scheme()).is_err());
    }

    #[test]
    fn msa_empty_sequence_error() {
        let seqs: Vec<&[u8]> = vec![b"ACGT", b""];
        assert!(progressive_msa(&seqs, &dna_scheme()).is_err());
    }

    #[test]
    fn msa_conservation_partial() {
        let seqs: Vec<&[u8]> = vec![b"AAAA", b"AATT"];
        let result = progressive_msa(&seqs, &dna_scheme()).unwrap();
        assert!(result.conservation() > 0.0);
        assert!(result.conservation() <= 1.0);
    }

    #[test]
    fn msa_column_access() {
        let seqs: Vec<&[u8]> = vec![b"ACGT", b"ACGT"];
        let result = progressive_msa(&seqs, &dna_scheme()).unwrap();
        let col = result.column(0).unwrap();
        assert_eq!(col.len(), 2);
        assert!(result.column(100).is_none());
    }
}
