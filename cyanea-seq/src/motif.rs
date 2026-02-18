//! DNA motif discovery — PWM construction, scanning, and EM-based de novo discovery.
//!
//! Complements the generic [`Pssm`](crate::pssm::Pssm) with DNA-specific features:
//! reverse complement scanning, consensus sequences, and MEME-style EM discovery.

use cyanea_core::{CyaneaError, Result};

/// A Position Weight Matrix for DNA motifs (A=0, C=1, G=2, T=3).
#[derive(Debug, Clone)]
pub struct Pwm {
    /// Frequency matrix: `matrix[pos] = [p_A, p_C, p_G, p_T]`.
    pub matrix: Vec<[f64; 4]>,
    /// Motif length.
    pub length: usize,
}

/// A motif match found by scanning.
#[derive(Debug, Clone)]
pub struct MotifMatch {
    /// Position in the sequence where the match starts.
    pub position: usize,
    /// Log-odds score of the match.
    pub score: f64,
    /// Which strand the match is on.
    pub strand: Strand,
}

/// Strand orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

/// A motif discovered by EM.
#[derive(Debug, Clone)]
pub struct DiscoveredMotif {
    /// The discovered PWM.
    pub pwm: Pwm,
    /// Sites where the motif was found: `(sequence_index, position)`.
    pub sites: Vec<(usize, usize)>,
    /// Log-likelihood ratio score of the motif.
    pub score: f64,
}

fn base_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

impl Pwm {
    /// Build a PWM from a set of aligned sequences of equal length.
    ///
    /// Adds a pseudocount of 0.25 per base to avoid zero probabilities.
    ///
    /// # Errors
    ///
    /// Returns an error if `sequences` is empty or sequences differ in length.
    pub fn from_aligned(sequences: &[&[u8]]) -> Result<Self> {
        if sequences.is_empty() {
            return Err(CyaneaError::InvalidInput(
                "at least one sequence is required".into(),
            ));
        }
        let len = sequences[0].len();
        if len == 0 {
            return Err(CyaneaError::InvalidInput(
                "sequences must be non-empty".into(),
            ));
        }
        for s in sequences {
            if s.len() != len {
                return Err(CyaneaError::InvalidInput(
                    "all sequences must have the same length".into(),
                ));
            }
        }

        let n = sequences.len() as f64;
        let pseudocount = 0.25;
        let total = n + 4.0 * pseudocount;

        let mut matrix = vec![[0.0f64; 4]; len];
        for pos in 0..len {
            let mut counts = [pseudocount; 4];
            for seq in sequences {
                if let Some(idx) = base_index(seq[pos]) {
                    counts[idx] += 1.0;
                }
            }
            for j in 0..4 {
                matrix[pos][j] = counts[j] / total;
            }
        }

        Ok(Self {
            matrix,
            length: len,
        })
    }

    /// Build a PWM from raw base counts (no pseudocount added).
    pub fn from_counts(counts: &[[usize; 4]]) -> Self {
        let mut matrix = Vec::with_capacity(counts.len());
        for row in counts {
            let total: usize = row.iter().sum();
            let t = if total > 0 { total as f64 } else { 1.0 };
            matrix.push([
                row[0] as f64 / t,
                row[1] as f64 / t,
                row[2] as f64 / t,
                row[3] as f64 / t,
            ]);
        }
        let length = matrix.len();
        Self { matrix, length }
    }

    /// Score a sequence window against this PWM using log-odds.
    ///
    /// `background` is `[p_A, p_C, p_G, p_T]`. The window must be
    /// exactly `self.length` bases long.
    pub fn score_sequence(&self, seq: &[u8], background: &[f64; 4]) -> f64 {
        let mut score = 0.0;
        for (pos, &base) in seq.iter().enumerate().take(self.length) {
            if let Some(idx) = base_index(base) {
                let p = self.matrix[pos][idx];
                let bg = background[idx];
                if p > 0.0 && bg > 0.0 {
                    score += (p / bg).log2();
                }
            }
        }
        score
    }

    /// Scan a sequence for motif matches above a score threshold.
    ///
    /// Checks both forward and reverse complement strands.
    pub fn scan(
        &self,
        seq: &[u8],
        background: &[f64; 4],
        threshold: f64,
    ) -> Vec<MotifMatch> {
        let mut matches = Vec::new();
        if seq.len() < self.length {
            return matches;
        }

        let rc_pwm = self.reverse_complement();

        for i in 0..=seq.len() - self.length {
            let window = &seq[i..i + self.length];

            // Forward strand.
            let fwd_score = self.score_sequence(window, background);
            if fwd_score >= threshold {
                matches.push(MotifMatch {
                    position: i,
                    score: fwd_score,
                    strand: Strand::Forward,
                });
            }

            // Reverse strand.
            let rev_score = rc_pwm.score_sequence(window, background);
            if rev_score >= threshold {
                matches.push(MotifMatch {
                    position: i,
                    score: rev_score,
                    strand: Strand::Reverse,
                });
            }
        }
        matches
    }

    /// Information content at each position (in bits).
    ///
    /// IC = 2 - H, where H = -Σ p log2(p).
    pub fn information_content(&self) -> Vec<f64> {
        self.matrix
            .iter()
            .map(|row| {
                let entropy: f64 = row
                    .iter()
                    .filter(|&&p| p > 0.0)
                    .map(|&p| -p * p.log2())
                    .sum();
                2.0 - entropy
            })
            .collect()
    }

    /// Total information content of the motif (sum across positions).
    pub fn total_information(&self) -> f64 {
        self.information_content().iter().sum()
    }

    /// Consensus sequence (most frequent base at each position).
    pub fn consensus(&self) -> Vec<u8> {
        let bases = [b'A', b'C', b'G', b'T'];
        self.matrix
            .iter()
            .map(|row| {
                let max_idx = row
                    .iter()
                    .enumerate()
                    .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                    .unwrap()
                    .0;
                bases[max_idx]
            })
            .collect()
    }

    /// Reverse complement of the PWM.
    ///
    /// Reverses position order and swaps A↔T, C↔G columns.
    pub fn reverse_complement(&self) -> Self {
        let matrix: Vec<[f64; 4]> = self
            .matrix
            .iter()
            .rev()
            .map(|row| {
                // Swap: A(0)↔T(3), C(1)↔G(2)
                [row[3], row[2], row[1], row[0]]
            })
            .collect();
        Self {
            length: self.length,
            matrix,
        }
    }
}

/// Discover motifs in a set of sequences using EM (MEME-style).
///
/// Searches for `n_motifs` motifs of `motif_length` bases. Each motif is
/// found by expectation-maximization, then its sites are masked before
/// searching for the next.
///
/// # Errors
///
/// Returns an error if sequences are empty, any sequence is shorter than
/// `motif_length`, or `motif_length` is zero.
pub fn discover_motifs(
    sequences: &[&[u8]],
    motif_length: usize,
    n_motifs: usize,
    max_iter: usize,
) -> Result<Vec<DiscoveredMotif>> {
    if sequences.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "at least one sequence is required".into(),
        ));
    }
    if motif_length == 0 {
        return Err(CyaneaError::InvalidInput(
            "motif_length must be at least 1".into(),
        ));
    }
    for (i, seq) in sequences.iter().enumerate() {
        if seq.len() < motif_length {
            return Err(CyaneaError::InvalidInput(format!(
                "sequence {} (length {}) is shorter than motif_length {}",
                i,
                seq.len(),
                motif_length
            )));
        }
    }

    let background = [0.25f64; 4];

    // Work on owned copies so we can mask discovered sites.
    let mut working_seqs: Vec<Vec<u8>> = sequences.iter().map(|s| s.to_vec()).collect();
    let mut motifs = Vec::new();

    for _ in 0..n_motifs {
        let refs: Vec<&[u8]> = working_seqs.iter().map(|s| s.as_slice()).collect();
        if let Some(motif) = em_one_motif(&refs, motif_length, max_iter, &background) {
            // Mask discovered sites with N.
            for &(seq_idx, pos) in &motif.sites {
                for j in pos..pos + motif_length {
                    if j < working_seqs[seq_idx].len() {
                        working_seqs[seq_idx][j] = b'N';
                    }
                }
            }
            motifs.push(motif);
        } else {
            break;
        }
    }

    Ok(motifs)
}

/// Run one round of EM motif discovery.
fn em_one_motif(
    sequences: &[&[u8]],
    motif_length: usize,
    max_iter: usize,
    background: &[f64; 4],
) -> Option<DiscoveredMotif> {
    // Find the best seed: subsequence that yields highest initial score.
    let mut best_pwm: Option<Pwm> = None;
    let mut best_ll = f64::NEG_INFINITY;

    // Try seeds from the first few sequences.
    let n_seed_seqs = sequences.len().min(3);
    for si in 0..n_seed_seqs {
        if sequences[si].len() < motif_length {
            continue;
        }
        // Sample a few starting positions.
        let step = (sequences[si].len() - motif_length + 1).max(1);
        let n_seeds = step.min(10);
        let stride = step / n_seeds;
        for seed_start_idx in 0..n_seeds {
            let pos = seed_start_idx * stride;
            let seed = &sequences[si][pos..pos + motif_length];

            // Skip seeds with N.
            if seed.iter().any(|&b| base_index(b).is_none()) {
                continue;
            }

            // Initialize PWM from this single seed.
            let mut pwm = Pwm::from_aligned(&[seed]).ok()?;

            // Run EM iterations.
            for _ in 0..max_iter {
                // E-step: compute Z (posterior probability of motif start at each position).
                let mut weighted_counts = vec![[0.25f64; 4]; motif_length]; // pseudocounts
                let mut total_weight = 4.0 * 0.25 * motif_length as f64;

                for seq in sequences {
                    if seq.len() < motif_length {
                        continue;
                    }
                    let n_pos = seq.len() - motif_length + 1;

                    // Compute scores for all positions.
                    let mut scores: Vec<f64> = Vec::with_capacity(n_pos);
                    for j in 0..n_pos {
                        let window = &seq[j..j + motif_length];
                        if window.iter().any(|&b| base_index(b).is_none()) {
                            scores.push(f64::NEG_INFINITY);
                        } else {
                            scores.push(pwm.score_sequence(window, background));
                        }
                    }

                    // Convert to probabilities via softmax.
                    let max_score = scores
                        .iter()
                        .copied()
                        .filter(|s| s.is_finite())
                        .fold(f64::NEG_INFINITY, f64::max);
                    if !max_score.is_finite() {
                        continue;
                    }

                    let exp_scores: Vec<f64> = scores
                        .iter()
                        .map(|&s| if s.is_finite() { (s - max_score).exp() } else { 0.0 })
                        .collect();
                    let sum_exp: f64 = exp_scores.iter().sum();
                    if sum_exp <= 0.0 {
                        continue;
                    }

                    // M-step: accumulate weighted counts.
                    for j in 0..n_pos {
                        let z = exp_scores[j] / sum_exp;
                        if z < 1e-10 {
                            continue;
                        }
                        for p in 0..motif_length {
                            if let Some(idx) = base_index(seq[j + p]) {
                                weighted_counts[p][idx] += z;
                                total_weight += z;
                            }
                        }
                    }
                }

                // Update PWM from weighted counts.
                let _ = total_weight; // used implicitly in per-position normalization
                let mut new_matrix = vec![[0.0f64; 4]; motif_length];
                for p in 0..motif_length {
                    let row_total: f64 = weighted_counts[p].iter().sum();
                    if row_total > 0.0 {
                        for j in 0..4 {
                            new_matrix[p][j] = weighted_counts[p][j] / row_total;
                        }
                    } else {
                        new_matrix[p] = [0.25; 4];
                    }
                }
                pwm.matrix = new_matrix;
            }

            // Compute log-likelihood of this PWM.
            let ll = compute_ll(sequences, &pwm, background, motif_length);
            if ll > best_ll {
                best_ll = ll;
                best_pwm = Some(pwm);
            }
        }
    }

    let pwm = best_pwm?;

    // Find best site in each sequence.
    let mut sites = Vec::new();
    let mut total_score = 0.0;
    for (si, seq) in sequences.iter().enumerate() {
        if seq.len() < motif_length {
            continue;
        }
        let mut best_pos = 0;
        let mut best_score = f64::NEG_INFINITY;
        for j in 0..=seq.len() - motif_length {
            let window = &seq[j..j + motif_length];
            if window.iter().any(|&b| base_index(b).is_none()) {
                continue;
            }
            let s = pwm.score_sequence(window, background);
            if s > best_score {
                best_score = s;
                best_pos = j;
            }
        }
        if best_score.is_finite() && best_score > 0.0 {
            sites.push((si, best_pos));
            total_score += best_score;
        }
    }

    if sites.is_empty() {
        return None;
    }

    Some(DiscoveredMotif {
        pwm,
        sites,
        score: total_score,
    })
}

fn compute_ll(sequences: &[&[u8]], pwm: &Pwm, background: &[f64; 4], motif_length: usize) -> f64 {
    let mut ll = 0.0;
    for seq in sequences {
        if seq.len() < motif_length {
            continue;
        }
        let mut best = f64::NEG_INFINITY;
        for j in 0..=seq.len() - motif_length {
            let window = &seq[j..j + motif_length];
            if window.iter().any(|&b| base_index(b).is_none()) {
                continue;
            }
            let s = pwm.score_sequence(window, background);
            if s > best {
                best = s;
            }
        }
        if best.is_finite() {
            ll += best;
        }
    }
    ll
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pwm_from_aligned_sequences() {
        let seqs: Vec<&[u8]> = vec![b"ACGT", b"ACGT", b"ACGT"];
        let pwm = Pwm::from_aligned(&seqs).unwrap();
        assert_eq!(pwm.length, 4);
        // Position 0 should be heavily A.
        assert!(pwm.matrix[0][0] > pwm.matrix[0][1]);
        assert!(pwm.matrix[0][0] > pwm.matrix[0][2]);
        assert!(pwm.matrix[0][0] > pwm.matrix[0][3]);
    }

    #[test]
    fn pwm_score_perfect_match() {
        let seqs: Vec<&[u8]> = vec![b"ACGT", b"ACGT", b"ACGT"];
        let pwm = Pwm::from_aligned(&seqs).unwrap();
        let bg = [0.25; 4];
        let score = pwm.score_sequence(b"ACGT", &bg);
        // Perfect match should have positive score.
        assert!(score > 0.0);
    }

    #[test]
    fn pwm_scan_finds_motif() {
        let seqs: Vec<&[u8]> = vec![b"GATTACA", b"GATTACA"];
        let pwm = Pwm::from_aligned(&seqs).unwrap();
        let bg = [0.25; 4];
        let target = b"AAAGATTACAAAA";
        let matches = pwm.scan(target, &bg, 0.0);
        // Should find the motif on the forward strand.
        let fwd_matches: Vec<_> = matches
            .iter()
            .filter(|m| m.strand == Strand::Forward)
            .collect();
        assert!(!fwd_matches.is_empty());
        // Best forward match should be at position 3.
        let best = fwd_matches.iter().max_by(|a, b| a.score.partial_cmp(&b.score).unwrap()).unwrap();
        assert_eq!(best.position, 3);
    }

    #[test]
    fn pwm_information_content() {
        // Uniform distribution: IC = 0 at each position.
        let pwm = Pwm {
            matrix: vec![[0.25, 0.25, 0.25, 0.25]; 3],
            length: 3,
        };
        let ic = pwm.information_content();
        for &v in &ic {
            assert!(v.abs() < 1e-10);
        }

        // Perfect conservation: IC = 2 at each position.
        let pwm2 = Pwm {
            matrix: vec![[1.0, 0.0, 0.0, 0.0]; 3],
            length: 3,
        };
        let ic2 = pwm2.information_content();
        for &v in &ic2 {
            assert!((v - 2.0).abs() < 1e-10);
        }
    }

    #[test]
    fn pwm_reverse_complement() {
        // A motif "ACG" → reverse complement should be "CGT".
        let pwm = Pwm {
            matrix: vec![
                [1.0, 0.0, 0.0, 0.0], // A
                [0.0, 1.0, 0.0, 0.0], // C
                [0.0, 0.0, 1.0, 0.0], // G
            ],
            length: 3,
        };
        let rc = pwm.reverse_complement();
        // Position 0 of rc should represent T (complement of A at position 2 reversed).
        // Original pos 2: G → complement: C → rc pos 0 should be [0, 1, 0, 0]
        assert!((rc.matrix[0][1] - 1.0).abs() < 1e-10); // C
        assert!((rc.matrix[1][2] - 1.0).abs() < 1e-10); // G
        assert!((rc.matrix[2][3] - 1.0).abs() < 1e-10); // T
    }

    #[test]
    fn em_discovers_planted_motif() {
        // Plant a strong motif "ACGTAC" in random-ish backgrounds.
        let seqs: Vec<&[u8]> = vec![
            b"TTTTACGTACTTTT",
            b"GGGGACGTACGGGG",
            b"AAAACGTACAAAA",
            b"CCCCACGTACCCCC",
        ];
        let motifs = discover_motifs(&seqs, 6, 1, 20).unwrap();
        assert!(!motifs.is_empty());
        let m = &motifs[0];
        // The discovered motif should have consensus close to ACGTAC.
        let consensus = m.pwm.consensus();
        assert_eq!(&consensus, b"ACGTAC");
    }
}
