//! Shared host-side logic for GPU alignment: sequence encoding, traceback
//! reconstruction, and partitioning.

use crate::scoring::ScoringScheme;
use crate::types::{AlignmentMode, AlignmentResult, CigarOp};

/// Index entry for one sequence pair on the GPU.
#[derive(Debug, Clone, Copy)]
pub struct SeqIndex {
    pub query_offset: u32,
    pub query_len: u32,
    pub target_offset: u32,
    pub target_len: u32,
}

/// Configuration for GPU batch alignment.
#[derive(Debug, Clone)]
pub struct GpuAlignConfig {
    /// Maximum DP band width (default: 128). Pairs where `|m - n| > max_bandwidth`
    /// fall back to CPU.
    pub max_bandwidth: usize,
    /// Minimum number of pairs before GPU dispatch is worthwhile (default: 64).
    pub min_pairs_for_gpu: usize,
}

impl Default for GpuAlignConfig {
    fn default() -> Self {
        Self {
            max_bandwidth: 128,
            min_pairs_for_gpu: 64,
        }
    }
}

/// Traceback direction encoding (4 bits per cell):
/// - Bits [1:0]: H direction (diag=0, up=1, left=2, stop=3)
/// - Bit [2]: E from H+gap_open (1) vs E+gap_extend (0)
/// - Bit [3]: F from H+gap_open (1) vs F+gap_extend (0)
pub const TB_DIAG: u8 = 0;
pub const TB_UP: u8 = 1;
pub const TB_LEFT: u8 = 2;
pub const TB_STOP: u8 = 3;
pub const TB_H_MASK: u8 = 0x03;
pub const TB_E_OPEN: u8 = 0x04;
pub const TB_F_OPEN: u8 = 0x08;

/// Encode a nucleotide as a 2-bit value: A=0, C=1, G=2, T/U=3, other=0.
pub fn encode_base(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        _ => 0,
    }
}

/// Pack sequences into a contiguous buffer for GPU transfer.
/// Returns (packed bytes, index entries).
pub fn encode_pairs(pairs: &[(&[u8], &[u8])]) -> (Vec<u8>, Vec<SeqIndex>) {
    let total_len: usize = pairs.iter().map(|(q, t)| q.len() + t.len()).sum();
    let mut packed = Vec::with_capacity(total_len);
    let mut index = Vec::with_capacity(pairs.len());

    for (query, target) in pairs {
        let q_off = packed.len() as u32;
        packed.extend_from_slice(query);
        let t_off = packed.len() as u32;
        packed.extend_from_slice(target);

        index.push(SeqIndex {
            query_offset: q_off,
            query_len: query.len() as u32,
            target_offset: t_off,
            target_len: target.len() as u32,
        });
    }

    (packed, index)
}

/// Partition pairs into those suitable for GPU (within bandwidth) and those
/// that need CPU fallback.
pub fn partition_pairs<'a>(
    pairs: &[(&'a [u8], &'a [u8])],
    max_bandwidth: usize,
) -> (Vec<(usize, &'a [u8], &'a [u8])>, Vec<(usize, &'a [u8], &'a [u8])>) {
    let mut gpu_pairs = Vec::new();
    let mut cpu_pairs = Vec::new();

    for (i, (q, t)) in pairs.iter().enumerate() {
        let len_diff = if q.len() > t.len() {
            q.len() - t.len()
        } else {
            t.len() - q.len()
        };
        if len_diff > max_bandwidth || q.is_empty() || t.is_empty() {
            cpu_pairs.push((i, *q, *t));
        } else {
            gpu_pairs.push((i, *q, *t));
        }
    }

    (gpu_pairs, cpu_pairs)
}

/// Reconstruct an AlignmentResult from GPU traceback data for a single pair.
///
/// `traceback` is a flat array of traceback cells for the banded DP matrix.
/// `score` and `end_pos` come from the GPU output.
pub fn reconstruct_alignment(
    query: &[u8],
    target: &[u8],
    score: i32,
    end_i: usize,
    end_j: usize,
    traceback: &[u8],
    bandwidth: usize,
    mode: AlignmentMode,
    _scoring: &ScoringScheme,
) -> AlignmentResult {
    let m = query.len();
    let _n = target.len();

    if score == 0 && mode == AlignmentMode::Local {
        return AlignmentResult {
            score: 0,
            aligned_query: Vec::new(),
            aligned_target: Vec::new(),
            query_start: 0,
            query_end: 0,
            target_start: 0,
            target_end: 0,
            cigar: Vec::new(),
        };
    }

    let mut aligned_query = Vec::new();
    let mut aligned_target = Vec::new();
    let mut cigar_ops: Vec<CigarOp> = Vec::new();

    let mut i = end_i;
    let mut j = end_j;

    // Walk traceback
    while i > 0 || j > 0 {
        if i == 0 && j == 0 {
            break;
        }

        // Get band-relative position
        let band_center = (i as isize * _n as isize) / (m.max(1) as isize);
        let band_offset = j as isize - band_center + (bandwidth as isize / 2);
        if band_offset < 0 || band_offset >= bandwidth as isize {
            break;
        }

        let tb_idx = i * bandwidth + band_offset as usize;
        if tb_idx >= traceback.len() {
            break;
        }

        let tb = traceback[tb_idx];
        let h_dir = tb & TB_H_MASK;

        match h_dir {
            TB_STOP => break,
            TB_DIAG if i > 0 && j > 0 => {
                let op = if query[i - 1].to_ascii_uppercase()
                    == target[j - 1].to_ascii_uppercase()
                {
                    CigarOp::Match(1)
                } else {
                    CigarOp::Mismatch(1)
                };
                aligned_query.push(query[i - 1]);
                aligned_target.push(target[j - 1]);
                push_cigar(&mut cigar_ops, op);
                i -= 1;
                j -= 1;
            }
            TB_UP if i > 0 => {
                aligned_query.push(query[i - 1]);
                aligned_target.push(b'-');
                push_cigar(&mut cigar_ops, CigarOp::Deletion(1));
                i -= 1;
            }
            TB_LEFT if j > 0 => {
                aligned_query.push(b'-');
                aligned_target.push(target[j - 1]);
                push_cigar(&mut cigar_ops, CigarOp::Insertion(1));
                j -= 1;
            }
            _ => break,
        }

        if mode == AlignmentMode::Local {
            // Check if we've reached a zero-score cell
            let at_zero = h_dir == TB_STOP;
            if at_zero {
                break;
            }
        }
    }

    aligned_query.reverse();
    aligned_target.reverse();
    cigar_ops.reverse();

    AlignmentResult {
        score,
        aligned_query,
        aligned_target,
        query_start: i,
        query_end: end_i,
        target_start: j,
        target_end: end_j,
        cigar: cigar_ops,
    }
}

/// Merge a 1-length CIGAR op with the last op if same variant.
fn push_cigar(ops: &mut Vec<CigarOp>, op: CigarOp) {
    if let Some(last) = ops.last_mut() {
        match (last, &op) {
            (CigarOp::Match(ref mut n), CigarOp::Match(1)) => {
                *n += 1;
                return;
            }
            (CigarOp::Mismatch(ref mut n), CigarOp::Mismatch(1)) => {
                *n += 1;
                return;
            }
            (CigarOp::Insertion(ref mut n), CigarOp::Insertion(1)) => {
                *n += 1;
                return;
            }
            (CigarOp::Deletion(ref mut n), CigarOp::Deletion(1)) => {
                *n += 1;
                return;
            }
            _ => {}
        }
    }
    ops.push(op);
}

/// Extract scoring parameters for GPU transfer (simple nucleotide mode).
pub fn extract_scoring_params(scoring: &ScoringScheme) -> (i32, i32, i32, i32) {
    match scoring {
        ScoringScheme::Simple(m) => (m.match_score, m.mismatch_score, m.gap_open, m.gap_extend),
        ScoringScheme::Substitution(m) => {
            // For protein scoring, use gap params; match/mismatch are in the matrix
            (0, 0, m.gap_open, m.gap_extend)
        }
    }
}

/// Extract the 24x24 substitution matrix as a flat i32 array, or None for simple scoring.
pub fn extract_substitution_matrix(scoring: &ScoringScheme) -> Option<Vec<i32>> {
    match scoring {
        ScoringScheme::Simple(_) => None,
        ScoringScheme::Substitution(m) => {
            // Score all pairs to build the matrix
            let aa_chars = b"ARNDCQEGHILKMFPSTWYVBZX*";
            let mut matrix = vec![0i32; 24 * 24];
            for (i, &a) in aa_chars.iter().enumerate() {
                for (j, &b) in aa_chars.iter().enumerate() {
                    matrix[i * 24 + j] = m.score_pair(a, b);
                }
            }
            Some(matrix)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::ScoringMatrix;

    #[test]
    fn encode_base_values() {
        assert_eq!(encode_base(b'A'), 0);
        assert_eq!(encode_base(b'C'), 1);
        assert_eq!(encode_base(b'G'), 2);
        assert_eq!(encode_base(b'T'), 3);
        assert_eq!(encode_base(b'a'), 0);
    }

    #[test]
    fn encode_pairs_round_trip() {
        let pairs: Vec<(&[u8], &[u8])> = vec![
            (b"ACGT", b"TGCA"),
            (b"AA", b"CC"),
        ];
        let (packed, index) = encode_pairs(&pairs);
        assert_eq!(packed.len(), 4 + 4 + 2 + 2);
        assert_eq!(index.len(), 2);
        assert_eq!(index[0].query_len, 4);
        assert_eq!(index[0].target_len, 4);
        assert_eq!(index[1].query_len, 2);
        assert_eq!(index[1].target_len, 2);
    }

    #[test]
    fn partition_within_bandwidth() {
        let pairs: Vec<(&[u8], &[u8])> = vec![
            (b"ACGT", b"ACGT"),      // same length, within band
            (b"A", b"ACGTACGTACGT"),  // diff > 8
        ];
        let (gpu, cpu) = partition_pairs(&pairs, 8);
        assert_eq!(gpu.len(), 1);
        assert_eq!(cpu.len(), 1);
        assert_eq!(gpu[0].0, 0);
        assert_eq!(cpu[0].0, 1);
    }

    #[test]
    fn extract_dna_scoring() {
        let scoring = ScoringScheme::Simple(ScoringMatrix::dna_default());
        let (m, mm, go, ge) = extract_scoring_params(&scoring);
        assert_eq!(m, 2);
        assert_eq!(mm, -1);
        assert_eq!(go, -5);
        assert_eq!(ge, -2);
    }

    #[test]
    fn extract_no_substitution_for_simple() {
        let scoring = ScoringScheme::Simple(ScoringMatrix::dna_default());
        assert!(extract_substitution_matrix(&scoring).is_none());
    }
}
