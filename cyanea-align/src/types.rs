//! Core types for sequence alignment results.

use core::fmt;

/// The alignment strategy to use.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum AlignmentMode {
    /// Local alignment (Smith-Waterman) — finds the best-scoring local region.
    Local,
    /// Global alignment (Needleman-Wunsch) — aligns sequences end-to-end.
    Global,
    /// Semi-global alignment — one sequence aligned end-to-end, the other locally.
    SemiGlobal,
}

/// A single CIGAR operation describing how aligned sequences relate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum CigarOp {
    /// Matching bases (identical in query and target).
    Match(usize),
    /// Mismatching bases (different in query and target).
    Mismatch(usize),
    /// Insertion in the query (gap in target).
    Insertion(usize),
    /// Deletion from the query (gap in query).
    Deletion(usize),
}

impl CigarOp {
    /// Single-character CIGAR code: `=` (match), `X` (mismatch), `I`, `D`.
    pub fn code(&self) -> char {
        match self {
            CigarOp::Match(_) => '=',
            CigarOp::Mismatch(_) => 'X',
            CigarOp::Insertion(_) => 'I',
            CigarOp::Deletion(_) => 'D',
        }
    }

    /// Number of positions consumed by this operation.
    pub fn len(&self) -> usize {
        match self {
            CigarOp::Match(n)
            | CigarOp::Mismatch(n)
            | CigarOp::Insertion(n)
            | CigarOp::Deletion(n) => *n,
        }
    }

    /// Whether this operation has zero length.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl fmt::Display for CigarOp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", self.len(), self.code())
    }
}

/// The result of a pairwise sequence alignment.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct AlignmentResult {
    /// Alignment score.
    pub score: i32,
    /// Aligned query sequence (with `-` for gaps).
    pub aligned_query: Vec<u8>,
    /// Aligned target sequence (with `-` for gaps).
    pub aligned_target: Vec<u8>,
    /// Start position in the original query (0-based, inclusive).
    pub query_start: usize,
    /// End position in the original query (0-based, exclusive).
    pub query_end: usize,
    /// Start position in the original target (0-based, inclusive).
    pub target_start: usize,
    /// End position in the original target (0-based, exclusive).
    pub target_end: usize,
    /// CIGAR operations describing the alignment.
    pub cigar: Vec<CigarOp>,
}

impl AlignmentResult {
    /// Format the CIGAR vector as a compact string, e.g. `"4=1I3="`.
    pub fn cigar_string(&self) -> String {
        let mut s = String::new();
        for op in &self.cigar {
            s.push_str(&format!("{}{}", op.len(), op.code()));
        }
        s
    }

    /// Fraction of aligned columns that are exact matches, in `[0.0, 1.0]`.
    ///
    /// Returns 0.0 if the alignment is empty.
    pub fn identity(&self) -> f64 {
        let total = self.length();
        if total == 0 {
            return 0.0;
        }
        self.matches() as f64 / total as f64
    }

    /// Number of matching positions.
    pub fn matches(&self) -> usize {
        self.cigar
            .iter()
            .filter_map(|op| match op {
                CigarOp::Match(n) => Some(n),
                _ => None,
            })
            .sum()
    }

    /// Number of mismatching positions.
    pub fn mismatches(&self) -> usize {
        self.cigar
            .iter()
            .filter_map(|op| match op {
                CigarOp::Mismatch(n) => Some(n),
                _ => None,
            })
            .sum()
    }

    /// Number of gap positions (insertions + deletions).
    pub fn gaps(&self) -> usize {
        self.cigar
            .iter()
            .filter_map(|op| match op {
                CigarOp::Insertion(n) | CigarOp::Deletion(n) => Some(n),
                _ => None,
            })
            .sum()
    }

    /// Total length of the alignment (all CIGAR operation lengths).
    pub fn length(&self) -> usize {
        self.cigar.iter().map(|op| op.len()).sum()
    }
}

impl cyanea_core::Scored for AlignmentResult {
    fn score(&self) -> f64 {
        self.score as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cigar_string_formatting() {
        let result = AlignmentResult {
            score: 10,
            aligned_query: b"ACGT".to_vec(),
            aligned_target: b"ACGT".to_vec(),
            query_start: 0,
            query_end: 4,
            target_start: 0,
            target_end: 4,
            cigar: vec![CigarOp::Match(4)],
        };
        assert_eq!(result.cigar_string(), "4=");
    }

    #[test]
    fn cigar_string_mixed_ops() {
        let result = AlignmentResult {
            score: 5,
            aligned_query: Vec::new(),
            aligned_target: Vec::new(),
            query_start: 0,
            query_end: 0,
            target_start: 0,
            target_end: 0,
            cigar: vec![
                CigarOp::Match(3),
                CigarOp::Insertion(1),
                CigarOp::Match(2),
                CigarOp::Deletion(1),
                CigarOp::Mismatch(1),
            ],
        };
        assert_eq!(result.cigar_string(), "3=1I2=1D1X");
    }

    #[test]
    fn identity_perfect_match() {
        let result = AlignmentResult {
            score: 8,
            aligned_query: b"ACGT".to_vec(),
            aligned_target: b"ACGT".to_vec(),
            query_start: 0,
            query_end: 4,
            target_start: 0,
            target_end: 4,
            cigar: vec![CigarOp::Match(4)],
        };
        assert!((result.identity() - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn identity_with_mismatches_and_gaps() {
        let result = AlignmentResult {
            score: 0,
            aligned_query: Vec::new(),
            aligned_target: Vec::new(),
            query_start: 0,
            query_end: 0,
            target_start: 0,
            target_end: 0,
            cigar: vec![
                CigarOp::Match(3),
                CigarOp::Mismatch(1),
                CigarOp::Insertion(1),
            ],
        };
        assert!((result.identity() - 0.6).abs() < f64::EPSILON);
    }

    #[test]
    fn identity_empty_alignment() {
        let result = AlignmentResult {
            score: 0,
            aligned_query: Vec::new(),
            aligned_target: Vec::new(),
            query_start: 0,
            query_end: 0,
            target_start: 0,
            target_end: 0,
            cigar: vec![],
        };
        assert!((result.identity() - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn gap_count() {
        let result = AlignmentResult {
            score: 0,
            aligned_query: Vec::new(),
            aligned_target: Vec::new(),
            query_start: 0,
            query_end: 0,
            target_start: 0,
            target_end: 0,
            cigar: vec![
                CigarOp::Match(2),
                CigarOp::Insertion(3),
                CigarOp::Match(1),
                CigarOp::Deletion(2),
            ],
        };
        assert_eq!(result.gaps(), 5);
        assert_eq!(result.matches(), 3);
        assert_eq!(result.length(), 8);
    }

    #[test]
    fn scored_trait() {
        use cyanea_core::Scored;
        let result = AlignmentResult {
            score: 42,
            aligned_query: Vec::new(),
            aligned_target: Vec::new(),
            query_start: 0,
            query_end: 0,
            target_start: 0,
            target_end: 0,
            cigar: vec![],
        };
        assert!((result.score() - 42.0).abs() < f64::EPSILON);
    }

    #[test]
    fn cigar_op_display() {
        assert_eq!(format!("{}", CigarOp::Match(5)), "5=");
        assert_eq!(format!("{}", CigarOp::Mismatch(2)), "2X");
        assert_eq!(format!("{}", CigarOp::Insertion(1)), "1I");
        assert_eq!(format!("{}", CigarOp::Deletion(3)), "3D");
    }
}
