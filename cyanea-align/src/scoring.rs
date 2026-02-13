//! Scoring schemes for pairwise sequence alignment.
//!
//! Provides simple match/mismatch scoring for nucleotides ([`ScoringMatrix`]),
//! amino acid substitution matrices ([`SubstitutionMatrix`]) with BLOSUM and PAM
//! variants, and a unified [`ScoringScheme`] enum that the alignment algorithms accept.

use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// Simple scoring (nucleotides)
// ---------------------------------------------------------------------------

/// A simple match/mismatch scoring matrix with affine gap penalties.
///
/// Suitable for nucleotide alignments where all matches score the same
/// and all mismatches score the same.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ScoringMatrix {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

impl ScoringMatrix {
    /// Create a new scoring matrix.
    ///
    /// # Errors
    ///
    /// Returns an error if `match_score` is not positive, `mismatch_score` is not
    /// negative, or `gap_open`/`gap_extend` are not negative.
    pub fn new(
        match_score: i32,
        mismatch_score: i32,
        gap_open: i32,
        gap_extend: i32,
    ) -> Result<Self> {
        if match_score <= 0 {
            return Err(CyaneaError::InvalidInput(
                "match_score must be positive".into(),
            ));
        }
        if mismatch_score >= 0 {
            return Err(CyaneaError::InvalidInput(
                "mismatch_score must be negative".into(),
            ));
        }
        if gap_open >= 0 {
            return Err(CyaneaError::InvalidInput(
                "gap_open must be negative".into(),
            ));
        }
        if gap_extend >= 0 {
            return Err(CyaneaError::InvalidInput(
                "gap_extend must be negative".into(),
            ));
        }
        Ok(Self {
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
        })
    }

    /// Default scoring for DNA alignment: +2 match, -1 mismatch, -5 gap open, -2 gap extend.
    pub fn dna_default() -> Self {
        Self {
            match_score: 2,
            mismatch_score: -1,
            gap_open: -5,
            gap_extend: -2,
        }
    }

    /// Score a pair of bases. Case-insensitive.
    pub fn score_pair(&self, a: u8, b: u8) -> i32 {
        if a.to_ascii_uppercase() == b.to_ascii_uppercase() {
            self.match_score
        } else {
            self.mismatch_score
        }
    }
}

// ---------------------------------------------------------------------------
// Amino acid index mapping
// ---------------------------------------------------------------------------

/// Maps an amino acid letter to a 0-based index in substitution matrices.
///
/// Standard 20 amino acids + B (Asx), Z (Glx), X (unknown), * (stop).
/// Returns `None` for unrecognized characters.
fn aa_to_index(aa: u8) -> Option<usize> {
    match aa.to_ascii_uppercase() {
        b'A' => Some(0),
        b'R' => Some(1),
        b'N' => Some(2),
        b'D' => Some(3),
        b'C' => Some(4),
        b'Q' => Some(5),
        b'E' => Some(6),
        b'G' => Some(7),
        b'H' => Some(8),
        b'I' => Some(9),
        b'L' => Some(10),
        b'K' => Some(11),
        b'M' => Some(12),
        b'F' => Some(13),
        b'P' => Some(14),
        b'S' => Some(15),
        b'T' => Some(16),
        b'W' => Some(17),
        b'Y' => Some(18),
        b'V' => Some(19),
        b'B' => Some(20),
        b'Z' => Some(21),
        b'X' => Some(22),
        b'*' => Some(23),
        _ => None,
    }
}

/// Matrix dimension: 24 amino acid symbols.
const AA_DIM: usize = 24;

// ---------------------------------------------------------------------------
// Substitution matrices (protein)
// ---------------------------------------------------------------------------

/// An amino acid substitution matrix with affine gap penalties.
///
/// Stores a 24x24 lookup table covering the 20 standard amino acids plus
/// B (Asx), Z (Glx), X (unknown), and * (stop codon).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SubstitutionMatrix {
    /// 24x24 flattened score table (row-major).
    scores: Vec<i32>,
    pub gap_open: i32,
    pub gap_extend: i32,
    name: &'static str,
}

impl SubstitutionMatrix {
    /// Score a pair of amino acids. Case-insensitive.
    ///
    /// Returns the worst score in the matrix for unrecognised residues.
    pub fn score_pair(&self, a: u8, b: u8) -> i32 {
        match (aa_to_index(a), aa_to_index(b)) {
            (Some(i), Some(j)) => self.scores[i * AA_DIM + j],
            _ => self.worst_score(),
        }
    }

    fn worst_score(&self) -> i32 {
        self.scores.iter().copied().min().unwrap_or(-4)
    }

    /// Matrix name (e.g. "BLOSUM62").
    pub fn name(&self) -> &str {
        self.name
    }

    // -----------------------------------------------------------------------
    // Standard matrices — values from NCBI
    // Order: A R N D C Q E G H I L K M F P S T W Y V B Z X *
    // -----------------------------------------------------------------------

    /// BLOSUM62 substitution matrix. Gap penalties: -11 open, -1 extend.
    pub fn blosum62() -> Self {
        Self {
            scores: BLOSUM62.to_vec(),
            gap_open: -11,
            gap_extend: -1,
            name: "BLOSUM62",
        }
    }

    /// BLOSUM45 substitution matrix. Gap penalties: -13 open, -3 extend.
    pub fn blosum45() -> Self {
        Self {
            scores: BLOSUM45.to_vec(),
            gap_open: -13,
            gap_extend: -3,
            name: "BLOSUM45",
        }
    }

    /// BLOSUM80 substitution matrix. Gap penalties: -10 open, -1 extend.
    pub fn blosum80() -> Self {
        Self {
            scores: BLOSUM80.to_vec(),
            gap_open: -10,
            gap_extend: -1,
            name: "BLOSUM80",
        }
    }

    /// PAM250 substitution matrix. Gap penalties: -11 open, -1 extend.
    pub fn pam250() -> Self {
        Self {
            scores: PAM250.to_vec(),
            gap_open: -11,
            gap_extend: -1,
            name: "PAM250",
        }
    }

    /// PAM40 substitution matrix. Gap penalties: -10 open, -2 extend.
    pub fn pam40() -> Self {
        Self {
            scores: PAM40.to_vec(),
            gap_open: -10,
            gap_extend: -2,
            name: "PAM40",
        }
    }

    /// PAM120 substitution matrix. Gap penalties: -11 open, -1 extend.
    pub fn pam120() -> Self {
        Self {
            scores: PAM120.to_vec(),
            gap_open: -11,
            gap_extend: -1,
            name: "PAM120",
        }
    }

    /// PAM200 substitution matrix. Gap penalties: -11 open, -1 extend.
    pub fn pam200() -> Self {
        Self {
            scores: PAM200.to_vec(),
            gap_open: -11,
            gap_extend: -1,
            name: "PAM200",
        }
    }

    /// BLOSUM30 substitution matrix. Gap penalties: -14 open, -4 extend.
    pub fn blosum30() -> Self {
        Self {
            scores: BLOSUM30.to_vec(),
            gap_open: -14,
            gap_extend: -4,
            name: "BLOSUM30",
        }
    }
}

// ---------------------------------------------------------------------------
// Unified scoring scheme
// ---------------------------------------------------------------------------

/// A unified scoring scheme accepted by alignment algorithms.
#[derive(Debug, Clone)]
pub enum ScoringScheme {
    /// Simple match/mismatch scoring (typically for nucleotides).
    Simple(ScoringMatrix),
    /// Amino acid substitution matrix (BLOSUM, PAM, etc.).
    Substitution(SubstitutionMatrix),
}

impl ScoringScheme {
    /// Score a pair of residues under this scheme.
    pub fn score_pair(&self, a: u8, b: u8) -> i32 {
        match self {
            ScoringScheme::Simple(m) => m.score_pair(a, b),
            ScoringScheme::Substitution(m) => m.score_pair(a, b),
        }
    }

    /// Gap opening penalty (negative).
    pub fn gap_open(&self) -> i32 {
        match self {
            ScoringScheme::Simple(m) => m.gap_open,
            ScoringScheme::Substitution(m) => m.gap_open,
        }
    }

    /// Gap extension penalty (negative).
    pub fn gap_extend(&self) -> i32 {
        match self {
            ScoringScheme::Simple(m) => m.gap_extend,
            ScoringScheme::Substitution(m) => m.gap_extend,
        }
    }
}

impl From<ScoringMatrix> for ScoringScheme {
    fn from(m: ScoringMatrix) -> Self {
        ScoringScheme::Simple(m)
    }
}

impl From<SubstitutionMatrix> for ScoringScheme {
    fn from(m: SubstitutionMatrix) -> Self {
        ScoringScheme::Substitution(m)
    }
}

// ===========================================================================
// NCBI substitution matrix data
// Row/column order: A R N D C Q E G H I L K M F P S T W Y V B Z X *
// ===========================================================================

/// BLOSUM62 — 24x24 flattened, NCBI reference.
#[rustfmt::skip]
const BLOSUM62: [i32; AA_DIM * AA_DIM] = [
//   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
     4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4, // A
    -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4, // R
    -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4, // N
    -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4, // D
     0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4, // C
    -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4, // Q
    -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4, // E
     0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4, // G
    -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4, // H
    -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4, // I
    -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4, // L
    -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4, // K
    -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4, // M
    -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4, // F
    -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4, // P
     1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4, // S
     0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4, // T
    -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4, // W
    -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4, // Y
     0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4, // V
    -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4, // B
    -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4, // Z
     0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4, // X
    -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1, // *
];

/// BLOSUM45 — 24x24 flattened, NCBI reference.
#[rustfmt::skip]
const BLOSUM45: [i32; AA_DIM * AA_DIM] = [
//   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
     5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -2, -2,  0, -1, -1,  0, -5, // A
    -2,  7,  0, -1, -3,  1,  0, -2,  0, -3, -2,  3, -1, -2, -2, -1, -1, -2, -1, -2, -1,  0, -1, -5, // R
    -1,  0,  6,  2, -2,  0,  0,  0,  1, -2, -3,  0, -2, -2, -2,  1,  0, -4, -2, -3,  4,  0, -1, -5, // N
    -2, -1,  2,  7, -3,  0,  2, -1,  0, -4, -3,  0, -3, -4, -1,  0, -1, -4, -2, -3,  5,  1, -1, -5, // D
    -1, -3, -2, -3, 12, -3, -3, -3, -3, -3, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -2, -3, -2, -5, // C
    -1,  1,  0,  0, -3,  6,  2, -2,  1, -2, -2,  1,  0, -4, -1,  0, -1, -2, -1, -3,  0,  4, -1, -5, // Q
    -1,  0,  0,  2, -3,  2,  6, -2,  0, -3, -2,  1, -2, -3,  0,  0, -1, -3, -2, -3,  1,  4, -1, -5, // E
     0, -2,  0, -1, -3, -2, -2,  7, -2, -4, -3, -2, -2, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -5, // G
    -2,  0,  1,  0, -3,  1,  0, -2, 10, -3, -2, -1,  0, -2, -2, -1, -2, -3,  2, -3,  0,  0, -1, -5, // H
    -1, -3, -2, -4, -3, -2, -3, -4, -3,  5,  2, -3,  2,  0, -2, -2, -1, -2,  0,  3, -3, -3, -1, -5, // I
    -1, -2, -3, -3, -2, -2, -2, -3, -2,  2,  5, -3,  2,  1, -3, -3, -1, -2,  0,  1, -3, -2, -1, -5, // L
    -1,  3,  0,  0, -3,  1,  1, -2, -1, -3, -3,  5, -1, -3, -1, -1, -1, -2, -1, -2,  0,  1, -1, -5, // K
    -1, -1, -2, -3, -2,  0, -2, -2,  0,  2,  2, -1,  6,  0, -2, -2, -1, -2,  0,  1, -2, -1, -1, -5, // M
    -2, -2, -2, -4, -2, -4, -3, -3, -2,  0,  1, -3,  0,  8, -3, -2, -1,  1,  3,  0, -3, -3, -1, -5, // F
    -1, -2, -2, -1, -4, -1,  0, -2, -2, -2, -3, -1, -2, -3,  9, -1, -1, -3, -3, -3, -2, -1, -1, -5, // P
     1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -3, -1, -2, -2, -1,  4,  2, -4, -2, -1,  0,  0,  0, -5, // S
     0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -1, -1,  2,  5, -3, -1,  0,  0, -1,  0, -5, // T
    -2, -2, -4, -4, -5, -2, -3, -2, -3, -2, -2, -2, -2,  1, -3, -4, -3, 15,  3, -3, -4, -2, -2, -5, // W
    -2, -1, -2, -2, -3, -1, -2, -3,  2,  0,  0, -1,  0,  3, -3, -2, -1,  3,  8, -1, -2, -2, -1, -5, // Y
     0, -2, -3, -3, -1, -3, -3, -3, -3,  3,  1, -2,  1,  0, -3, -1,  0, -3, -1,  5, -3, -3, -1, -5, // V
    -1, -1,  4,  5, -2,  0,  1, -1,  0, -3, -3,  0, -2, -3, -2,  0,  0, -4, -2, -3,  4,  2, -1, -5, // B
    -1,  0,  0,  1, -3,  4,  4, -2,  0, -3, -2,  1, -1, -3, -1,  0, -1, -2, -2, -3,  2,  4, -1, -5, // Z
     0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0, -2, -1, -1, -1, -1, -1, -5, // X
    -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  1, // *
];

/// BLOSUM80 — 24x24 flattened, NCBI reference.
#[rustfmt::skip]
const BLOSUM80: [i32; AA_DIM * AA_DIM] = [
//   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
     7, -3, -3, -3, -1, -2, -2,  0, -3, -3, -3, -1, -2, -4, -1,  2,  0, -5, -4, -1, -3, -2, -1, -8, // A
    -3,  9, -1, -3, -6,  1, -1, -4,  0, -5, -4,  3, -3, -5, -3, -2, -2, -5, -4, -4, -2,  0, -2, -8, // R
    -3, -1,  9,  2, -5,  0, -1, -1,  1, -6, -6,  0, -4, -6, -4,  1,  0, -7, -4, -5,  5,  0, -2, -8, // N
    -3, -3,  2, 10, -7, -1,  2, -3, -2, -7, -7, -2, -6, -6, -3, -1, -2, -8, -6, -6,  6,  1, -3, -8, // D
    -1, -6, -5, -7, 13, -5, -7, -6, -7, -2, -3, -6, -3, -4, -6, -2, -2, -5, -5, -2, -6, -7, -4, -8, // C
    -2,  1,  0, -1, -5,  9,  3, -4,  1, -5, -4,  2, -1, -5, -3, -1, -1, -4, -3, -4, -1,  5, -2, -8, // Q
    -2, -1, -1,  2, -7,  3,  8, -4,  0, -6, -6,  1, -4, -6, -2,  0, -2, -6, -5, -4,  1,  6, -2, -8, // E
     0, -4, -1, -3, -6, -4, -4,  9, -4, -7, -7, -3, -5, -6, -5, -1, -3, -6, -6, -6, -2, -4, -3, -8, // G
    -3,  0,  1, -2, -7,  1,  0, -4, 12, -6, -5, -1, -4, -2, -4, -2, -3, -4,  3, -5, -1,  0, -2, -8, // H
    -3, -5, -6, -7, -2, -5, -6, -7, -6,  7,  2, -5,  2, -1, -5, -4, -2, -5, -3,  4, -6, -6, -2, -8, // I
    -3, -4, -6, -7, -3, -4, -6, -7, -5,  2,  6, -4,  3,  0, -5, -4, -3, -4, -2,  1, -7, -5, -2, -8, // L
    -1,  3,  0, -2, -6,  2,  1, -3, -1, -5, -4,  8, -3, -5, -2, -1, -1, -6, -4, -4, -1,  1, -2, -8, // K
    -2, -3, -4, -6, -3, -1, -4, -5, -4,  2,  3, -3,  9, -1, -4, -3, -1, -3, -3,  1, -5, -3, -2, -8, // M
    -4, -5, -6, -6, -4, -5, -6, -6, -2, -1,  0, -5, -1, 10, -6, -4, -4,  0,  4, -2, -6, -6, -3, -8, // F
    -1, -3, -4, -3, -6, -3, -2, -5, -4, -5, -5, -2, -4, -6, 12, -2, -3, -7, -6, -4, -4, -2, -3, -8, // P
     2, -2,  1, -1, -2, -1,  0, -1, -2, -4, -4, -1, -3, -4, -2,  7,  2, -6, -3, -3,  0, -1, -1, -8, // S
     0, -2,  0, -2, -2, -1, -2, -3, -3, -2, -3, -1, -1, -4, -3,  2,  8, -5, -3,  0, -1, -2, -1, -8, // T
    -5, -5, -7, -8, -5, -4, -6, -6, -4, -5, -4, -6, -3,  0, -7, -6, -5, 16,  3, -5, -8, -5, -5, -8, // W
    -4, -4, -4, -6, -5, -3, -5, -6,  3, -3, -2, -4, -3,  4, -6, -3, -3,  3, 11, -3, -5, -4, -3, -8, // Y
    -1, -4, -5, -6, -2, -4, -4, -6, -5,  4,  1, -4,  1, -2, -4, -3,  0, -5, -3,  7, -6, -4, -2, -8, // V
    -3, -2,  5,  6, -6, -1,  1, -2, -1, -6, -7, -1, -5, -6, -4,  0, -1, -8, -5, -6,  6,  0, -3, -8, // B
    -2,  0,  0,  1, -7,  5,  6, -4,  0, -6, -5,  1, -3, -6, -2, -1, -2, -5, -4, -4,  0,  6, -1, -8, // Z
    -1, -2, -2, -3, -4, -2, -2, -3, -2, -2, -2, -2, -2, -3, -3, -1, -1, -5, -3, -2, -3, -1, -2, -8, // X
    -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8,  1, // *
];

/// PAM250 — 24x24 flattened, NCBI/Dayhoff reference.
#[rustfmt::skip]
const PAM250: [i32; AA_DIM * AA_DIM] = [
//   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
     2, -2,  0,  0, -2,  0,  0,  1, -1, -1, -2, -1, -1, -3,  1,  1,  1, -6, -3,  0,  0,  0,  0, -8, // A
    -2,  6,  0, -1, -4,  1, -1, -3,  2, -2, -3,  3,  0, -4,  0,  0, -1,  2, -4, -2, -1,  0, -1, -8, // R
     0,  0,  2,  2, -4,  1,  1,  0,  2, -2, -3,  1, -2, -3,  0,  1,  0, -4, -2, -2,  2,  1,  0, -8, // N
     0, -1,  2,  4, -5,  2,  3,  1,  1, -2, -4,  0, -3, -6, -1,  0,  0, -7, -4, -2,  3,  3, -1, -8, // D
    -2, -4, -4, -5, 12, -5, -5, -3, -3, -2, -6, -5, -5, -4, -3,  0, -2, -8,  0, -2, -4, -5, -3, -8, // C
     0,  1,  1,  2, -5,  4,  2, -1,  3, -2, -2,  1, -1, -5,  0, -1, -1, -5, -4, -2,  1,  3, -1, -8, // Q
     0, -1,  1,  3, -5,  2,  4,  0,  1, -2, -3,  0, -2, -5, -1,  0,  0, -7, -4, -2,  3,  3, -1, -8, // E
     1, -3,  0,  1, -3, -1,  0,  5, -2, -3, -4, -2, -3, -5,  0,  1,  0, -7, -5, -1,  0,  0, -1, -8, // G
    -1,  2,  2,  1, -3,  3,  1, -2,  6, -2, -2,  0, -2, -2,  0, -1, -1, -3,  0, -2,  1,  2, -1, -8, // H
    -1, -2, -2, -2, -2, -2, -2, -3, -2,  5,  2, -2,  2,  1, -2, -1,  0, -5, -1,  4, -2, -2, -1, -8, // I
    -2, -3, -3, -4, -6, -2, -3, -4, -2,  2,  6, -3,  4,  2, -3, -3, -2, -2, -1,  2, -3, -3, -1, -8, // L
    -1,  3,  1,  0, -5,  1,  0, -2,  0, -2, -3,  5,  0, -5, -1,  0,  0, -3, -4, -2,  1,  0, -1, -8, // K
    -1,  0, -2, -3, -5, -1, -2, -3, -2,  2,  4,  0,  6,  0, -2, -2, -1, -4, -2,  2, -2, -2, -1, -8, // M
    -3, -4, -3, -6, -4, -5, -5, -5, -2,  1,  2, -5,  0,  9, -5, -3, -3,  0,  7, -1, -4, -5, -2, -8, // F
     1,  0,  0, -1, -3,  0, -1,  0,  0, -2, -3, -1, -2, -5,  6,  1,  0, -6, -5, -1, -1,  0, -1, -8, // P
     1,  0,  1,  0,  0, -1,  0,  1, -1, -1, -3,  0, -2, -3,  1,  2,  1, -2, -3, -1,  0,  0,  0, -8, // S
     1, -1,  0,  0, -2, -1,  0,  0, -1,  0, -2,  0, -1, -3,  0,  1,  3, -5, -3,  0,  0, -1,  0, -8, // T
    -6,  2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4,  0, -6, -2, -5, 17,  0, -6, -5, -6, -4, -8, // W
    -3, -4, -2, -4,  0, -4, -4, -5,  0, -1, -1, -4, -2,  7, -5, -3, -3,  0, 10, -2, -3, -4, -2, -8, // Y
     0, -2, -2, -2, -2, -2, -2, -1, -2,  4,  2, -2,  2, -1, -1, -1,  0, -6, -2,  4, -2, -2, -1, -8, // V
     0, -1,  2,  3, -4,  1,  3,  0,  1, -2, -3,  1, -2, -4, -1,  0,  0, -5, -3, -2,  3,  2, -1, -8, // B
     0,  0,  1,  3, -5,  3,  3,  0,  2, -2, -3,  0, -2, -5,  0,  0, -1, -6, -4, -2,  2,  3, -1, -8, // Z
     0, -1,  0, -1, -3, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1,  0,  0, -4, -2, -1, -1, -1, -1, -8, // X
    -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8,  1, // *
];

/// PAM40 — 24x24 flattened, NCBI reference.
#[rustfmt::skip]
const PAM40: [i32; AA_DIM * AA_DIM] = [
//   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
     6, -4, -2, -1, -4, -2, -1,  0, -4, -3, -4, -4, -3, -6,  0,  1,  1, -9, -6, -1, -1, -1, -2,-11, // A
    -4,  8, -2, -4, -6,  0, -3, -6,  1, -3, -5,  2, -2, -7, -2, -1, -4,  0, -7, -5, -3, -2, -3,-11, // R
    -2, -2,  6,  3, -7, -1,  0, -1,  2, -4, -5,  0, -5, -6, -3,  1, -1, -6, -3, -5,  5, -1, -2,-11, // N
    -1, -4,  3,  7, -9,  0,  4, -2, -1, -5, -8, -2, -6, -9, -4, -1, -2,-10, -7, -5,  5,  3, -3,-11, // D
    -4, -6, -7, -9, 15, -9, -9, -6, -5, -4, -9, -9, -9, -8, -5, -1, -5,-11, -2, -4, -8, -9, -6,-11, // C
    -2,  0, -1,  0, -9,  7,  2, -4,  2, -5, -3,  0, -2, -8, -1, -3, -3, -8, -7, -4,  0,  5, -2,-11, // Q
    -1, -3,  0,  4, -9,  2,  7, -2, -2, -4, -6, -2, -4, -9, -3, -2, -3,-11, -6, -4,  3,  5, -3,-11, // E
     0, -6, -1, -2, -6, -4, -2,  7, -5, -6, -7, -4, -5, -7, -3,  0, -3,-10, -8, -4, -1, -3, -3,-11, // G
    -4,  1,  2, -1, -5,  2, -2, -5,  8, -5, -4, -3, -5, -4, -2, -3, -4, -5, -1, -4,  1,  0, -3,-11, // H
    -3, -3, -4, -5, -4, -5, -4, -6, -5,  8,  1, -4,  1,  0, -5, -4, -1, -9, -3,  3, -4, -4, -3,-11, // I
    -4, -5, -5, -8, -9, -3, -6, -7, -4,  1,  7, -5,  3,  0, -5, -5, -4, -4, -3,  0, -6, -4, -4,-11, // L
    -4,  2,  0, -2, -9,  0, -2, -4, -3, -4, -5,  6, -1, -9, -4, -2, -2, -7, -7, -5, -1, -1, -3,-11, // K
    -3, -2, -5, -6, -9, -2, -4, -5, -5,  1,  3, -1, 10, -2, -5, -3, -2, -8, -6,  0, -5, -3, -3,-11, // M
    -6, -7, -6, -9, -8, -8, -9, -7, -4,  0,  0, -9, -2, 11, -7, -4, -6, -3,  4, -4, -7, -8, -5,-11, // F
     0, -2, -3, -4, -5, -1, -3, -3, -2, -5, -5, -4, -5, -7,  8,  0, -2, -9, -8, -3, -4, -2, -3,-11, // P
     1, -1,  1, -1, -1, -3, -2,  0, -3, -4, -5, -2, -3, -4,  0,  5,  2, -4, -4, -3,  0, -2, -2,-11, // S
     1, -4, -1, -2, -5, -3, -3, -3, -4, -1, -4, -2, -2, -6, -2,  2,  6, -8, -4, -1, -2, -3, -2,-11, // T
    -9,  0, -6,-10,-11, -8,-11,-10, -5, -9, -4, -7, -8, -3, -9, -4, -8, 17, -3, -9, -7,-10, -7,-11, // W
    -6, -7, -3, -7, -2, -7, -6, -8, -1, -3, -3, -7, -6,  4, -8, -4, -4, -3, 12, -5, -4, -7, -5,-11, // Y
    -1, -5, -5, -5, -4, -4, -4, -4, -4,  3,  0, -5,  0, -4, -3, -3, -1, -9, -5,  6, -5, -4, -3,-11, // V
    -1, -3,  5,  5, -8,  0,  3, -1,  1, -4, -6, -1, -5, -7, -4,  0, -2, -7, -4, -5,  5,  2, -2,-11, // B
    -1, -2, -1,  3, -9,  5,  5, -3,  0, -4, -4, -1, -3, -8, -2, -2, -3,-10, -7, -4,  2,  5, -2,-11, // Z
    -2, -3, -2, -3, -6, -2, -3, -3, -3, -3, -4, -3, -3, -5, -3, -2, -2, -7, -5, -3, -2, -2, -3,-11, // X
   -11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,  1, // *
];

/// PAM120 — 24x24 flattened, NCBI reference.
#[rustfmt::skip]
const PAM120: [i32; AA_DIM * AA_DIM] = [
//   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
     3, -3,  0,  0, -3, -1,  0,  1, -3, -1, -3, -2, -2, -4,  1,  1,  1, -7, -4,  0,  0,  0, -1, -8, // A
    -3,  6, -1, -3, -4,  1, -3, -4,  1, -2, -4,  2, -1, -5, -1, -1, -2,  1, -5, -3, -2, -1, -2, -8, // R
     0, -1,  4,  2, -5,  0,  1,  0,  2, -2, -4,  1, -3, -4, -2,  1,  0, -4, -2, -3,  3,  0, -1, -8, // N
     0, -3,  2,  5, -7,  1,  3,  0,  0, -3, -5, -1, -4, -7, -3,  0, -1, -8, -5, -3,  4,  3, -2, -8, // D
    -3, -4, -5, -7, 12, -7, -7, -4, -4, -3, -7, -7, -6, -6, -4,  0, -3, -8, -1, -3, -6, -7, -4, -8, // C
    -1,  1,  0,  1, -7,  6,  2, -3,  3, -3, -2,  0, -1, -6,  0, -2, -2, -6, -5, -3,  0,  4, -1, -8, // Q
     0, -3,  1,  3, -7,  2,  5, -1, -1, -3, -4, -1, -3, -7, -2,  0, -2, -8, -5, -3,  3,  4, -1, -8, // E
     1, -4,  0,  0, -4, -3, -1,  5, -4, -4, -5, -3, -4, -5, -2,  1, -1, -8, -6, -2,  0, -2, -2, -8, // G
    -3,  1,  2,  0, -4,  3, -1, -4,  7, -4, -3, -2, -4, -3, -1, -2, -3, -3,  0, -3,  1,  1, -2, -8, // H
    -1, -2, -2, -3, -3, -3, -3, -4, -4,  6,  1, -3,  1,  0, -3, -2,  0, -6, -2,  3, -3, -3, -1, -8, // I
    -3, -4, -4, -5, -7, -2, -4, -5, -3,  1,  5, -4,  3,  0, -4, -4, -3, -3, -2,  1, -5, -3, -2, -8, // L
    -2,  2,  1, -1, -7,  0, -1, -3, -2, -3, -4,  5,  0, -7, -2, -1, -1, -5, -5, -4,  0, -1, -2, -8, // K
    -2, -1, -3, -4, -6, -1, -3, -4, -4,  1,  3,  0,  8, -1, -3, -2, -1, -6, -4,  1, -4, -2, -2, -8, // M
    -4, -5, -4, -7, -6, -6, -7, -5, -3,  0,  1, -7, -1,  8, -5, -3, -4,  0,  5, -3, -5, -6, -3, -8, // F
     1, -1, -2, -3, -4,  0, -2, -2, -1, -3, -4, -2, -3, -5,  6,  1,  0, -7, -6, -2, -2, -1, -2, -8, // P
     1, -1,  1,  0,  0, -2,  0,  1, -2, -2, -4, -1, -2, -3,  1,  3,  2, -2, -3, -2,  0, -1, -1, -8, // S
     1, -2,  0, -1, -3, -2, -2, -1, -3,  0, -3, -1, -1, -4,  0,  2,  4, -6, -3,  0, -1, -2, -1, -8, // T
    -7,  1, -4, -8, -8, -6, -8, -8, -3, -6, -3, -5, -6,  0, -7, -2, -6, 17, -1, -8, -6, -7, -5, -8, // W
    -4, -5, -2, -5, -1, -5, -5, -6,  0, -2, -2, -5, -4,  5, -6, -3, -3, -1, 10, -3, -3, -5, -3, -8, // Y
     0, -3, -3, -3, -3, -3, -3, -2, -3,  3,  1, -4,  1, -3, -2, -2,  0, -8, -3,  5, -3, -3, -1, -8, // V
     0, -2,  3,  4, -6,  0,  3,  0,  1, -3, -5,  0, -4, -5, -2,  0, -1, -6, -3, -3,  4,  2, -1, -8, // B
     0, -1,  0,  3, -7,  4,  4, -2,  1, -3, -3, -1, -2, -6, -1, -1, -2, -7, -5, -3,  2,  4, -1, -8, // Z
    -1, -2, -1, -2, -4, -1, -1, -2, -2, -1, -2, -2, -2, -3, -2, -1, -1, -5, -3, -1, -1, -1, -2, -8, // X
    -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8,  1, // *
];

/// PAM200 — 24x24 flattened, NCBI reference.
#[rustfmt::skip]
const PAM200: [i32; AA_DIM * AA_DIM] = [
//   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
     2, -2,  0,  0, -2,  0,  0,  1, -1, -1, -2, -1, -1, -3,  1,  1,  1, -6, -3,  0,  0,  0,  0, -8, // A
    -2,  6,  0, -1, -4,  1, -1, -3,  2, -2, -3,  3,  0, -4,  0,  0, -1,  2, -4, -2, -1,  0, -1, -8, // R
     0,  0,  2,  2, -4,  1,  1,  0,  2, -2, -3,  1, -2, -3,  0,  1,  0, -4, -2, -2,  2,  1,  0, -8, // N
     0, -1,  2,  4, -5,  2,  3,  1,  1, -2, -4,  0, -3, -6, -1,  0,  0, -7, -4, -2,  3,  3, -1, -8, // D
    -2, -4, -4, -5, 12, -5, -5, -3, -3, -2, -6, -5, -5, -4, -3,  0, -2, -8,  0, -2, -4, -5, -3, -8, // C
     0,  1,  1,  2, -5,  4,  2, -1,  3, -2, -2,  1, -1, -5,  0, -1, -1, -5, -4, -2,  1,  3, -1, -8, // Q
     0, -1,  1,  3, -5,  2,  4,  0,  1, -2, -3,  0, -2, -5, -1,  0,  0, -7, -4, -2,  3,  3, -1, -8, // E
     1, -3,  0,  1, -3, -1,  0,  5, -2, -3, -4, -2, -3, -5,  0,  1,  0, -7, -5, -1,  0,  0, -1, -8, // G
    -1,  2,  2,  1, -3,  3,  1, -2,  6, -2, -2,  0, -2, -2,  0, -1, -1, -3,  0, -2,  1,  2, -1, -8, // H
    -1, -2, -2, -2, -2, -2, -2, -3, -2,  5,  2, -2,  2,  1, -2, -1,  0, -5, -1,  4, -2, -2, -1, -8, // I
    -2, -3, -3, -4, -6, -2, -3, -4, -2,  2,  6, -3,  4,  2, -3, -3, -2, -2, -1,  2, -3, -3, -1, -8, // L
    -1,  3,  1,  0, -5,  1,  0, -2,  0, -2, -3,  5,  0, -5, -1,  0,  0, -3, -4, -2,  1,  0, -1, -8, // K
    -1,  0, -2, -3, -5, -1, -2, -3, -2,  2,  4,  0,  6,  0, -2, -2, -1, -4, -2,  2, -2, -2, -1, -8, // M
    -3, -4, -3, -6, -4, -5, -5, -5, -2,  1,  2, -5,  0,  9, -5, -3, -3,  0,  7, -1, -4, -5, -2, -8, // F
     1,  0,  0, -1, -3,  0, -1,  0,  0, -2, -3, -1, -2, -5,  6,  1,  0, -6, -5, -1, -1,  0, -1, -8, // P
     1,  0,  1,  0,  0, -1,  0,  1, -1, -1, -3,  0, -2, -3,  1,  2,  1, -2, -3, -1,  0,  0,  0, -8, // S
     1, -1,  0,  0, -2, -1,  0,  0, -1,  0, -2,  0, -1, -3,  0,  1,  3, -5, -3,  0,  0, -1,  0, -8, // T
    -6,  2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4,  0, -6, -2, -5, 17,  0, -6, -5, -6, -4, -8, // W
    -3, -4, -2, -4,  0, -4, -4, -5,  0, -1, -1, -4, -2,  7, -5, -3, -3,  0, 10, -2, -3, -4, -2, -8, // Y
     0, -2, -2, -2, -2, -2, -2, -1, -2,  4,  2, -2,  2, -1, -1, -1,  0, -6, -2,  4, -2, -2, -1, -8, // V
     0, -1,  2,  3, -4,  1,  3,  0,  1, -2, -3,  1, -2, -4, -1,  0,  0, -5, -3, -2,  3,  2, -1, -8, // B
     0,  0,  1,  3, -5,  3,  3,  0,  2, -2, -3,  0, -2, -5,  0,  0, -1, -6, -4, -2,  2,  3, -1, -8, // Z
     0, -1,  0, -1, -3, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1,  0,  0, -4, -2, -1, -1, -1, -1, -8, // X
    -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8,  1, // *
];

/// BLOSUM30 — 24x24 flattened, NCBI reference.
#[rustfmt::skip]
const BLOSUM30: [i32; AA_DIM * AA_DIM] = [
//   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
     4, -1,  0,  0, -3,  1,  0,  0, -2,  0, -1,  0,  1, -2, -1,  1,  1, -5, -4,  1,  0,  0,  0, -7, // A
    -1,  8, -2, -1, -2,  3, -1, -2, -1, -3, -2,  1,  0, -1, -1, -1, -3,  0,  0, -1, -2,  0, -1, -7, // R
     0, -2,  8,  1, -1, -1, -1,  0,  0, -2, -2,  0, -1, -1, -3,  0,  1, -7, -4, -2,  5, -1, -1, -7, // N
     0, -1,  1,  9, -3, -1,  1, -1, -2, -4, -1,  0, -3, -5, -1,  0, -1, -4, -1, -2,  5,  0, -1, -7, // D
    -3, -2, -1, -3, 17, -2,  1, -4, -5,  0, -2,  0, -2, -3, -3, -2, -2, -2, -6, -2, -2,  0, -2, -7, // C
     1,  3, -1, -1, -2,  8,  2, -2,  0, -2,  0,  0, -1, -3,  0, -1,  0, -1, -1, -3, -1,  4, -1, -7, // Q
     0, -1, -1,  1,  1,  2,  6, -2,  0, -3, -1,  1, -1, -4,  1,  0, -2, -1, -2, -3,  0,  5, -1, -7, // E
     0, -2,  0, -1, -4, -2, -2,  8, -3, -1, -2, -1, -2, -3, -1,  0, -2,  1, -3, -3,  0, -2, -1, -7, // G
    -2, -1,  0, -2, -5,  0,  0, -3, 14, -2, -1, -2,  2, -3,  1, -1, -2, -5,  0, -3, -1,  0, -1, -7, // H
     0, -3, -2, -4,  0, -2, -3, -1, -2,  6,  2, -2,  1,  0, -3, -1,  0, -3, -1,  4, -3, -3,  0, -7, // I
    -1, -2, -2, -1, -2,  0, -1, -2, -1,  2,  4, -2,  2,  2, -3, -2,  0, -2,  3,  1, -1, -1, -1, -7, // L
     0,  1,  0,  0,  0,  0,  1, -1, -2, -2, -2,  4,  2, -1,  0,  0, -1, -2, -1, -2,  0,  0,  0, -7, // K
     1,  0, -1, -3, -2, -1, -1, -2,  2,  1,  2,  2,  6,  0, -4, -2,  0, -3, -1,  0, -2, -1,  0, -7, // M
    -2, -1, -1, -5, -3, -3, -4, -3, -3,  0,  2, -1,  0, 10, -4, -1, -2,  0,  3, -1, -3, -4, -1, -7, // F
    -1, -1, -3, -1, -3,  0,  1, -1,  1, -3, -3,  0, -4, -4, 11, -1,  0, -3, -2, -4, -2,  0, -1, -7, // P
     1, -1,  0,  0, -2, -1,  0,  0, -1, -1, -2,  0, -2, -1, -1,  4,  2, -3, -2, -1,  0,  0,  0, -7, // S
     1, -3,  1, -1, -2,  0, -2, -2, -2,  0,  0, -1,  0, -2,  0,  2,  5, -5, -1,  1,  0, -1,  0, -7, // T
    -5,  0, -7, -4, -2, -1, -1,  1, -5, -3, -2, -2, -3,  0, -3, -3, -5, 20,  5, -3, -5, -1, -2, -7, // W
    -4,  0, -4, -1, -6, -1, -2, -3,  0, -1,  3, -1, -1,  3, -2, -2, -1,  5,  9,  1, -3, -2, -1, -7, // Y
     1, -1, -2, -2, -2, -3, -3, -3, -3,  4,  1, -2,  0, -1, -4, -1,  1, -3,  1,  5, -2, -3,  0, -7, // V
     0, -2,  5,  5, -2, -1,  0,  0, -1, -3, -1,  0, -2, -3, -2,  0,  0, -5, -3, -2,  5,  0, -1, -7, // B
     0,  0, -1,  0,  0,  4,  5, -2,  0, -3, -1,  0, -1, -4,  0,  0, -1, -1, -2, -3,  0,  4, -1, -7, // Z
     0, -1, -1, -1, -2, -1, -1, -1, -1,  0, -1,  0,  0, -1, -1,  0,  0, -2, -1,  0, -1, -1, -1, -7, // X
    -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7,  1, // *
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dna_default_values() {
        let m = ScoringMatrix::dna_default();
        assert_eq!(m.match_score, 2);
        assert_eq!(m.mismatch_score, -1);
        assert_eq!(m.gap_open, -5);
        assert_eq!(m.gap_extend, -2);
    }

    #[test]
    fn dna_score_pair_case_insensitive() {
        let m = ScoringMatrix::dna_default();
        assert_eq!(m.score_pair(b'A', b'A'), 2);
        assert_eq!(m.score_pair(b'a', b'A'), 2);
        assert_eq!(m.score_pair(b'A', b'T'), -1);
    }

    #[test]
    fn scoring_matrix_validation() {
        assert!(ScoringMatrix::new(0, -1, -5, -2).is_err());
        assert!(ScoringMatrix::new(2, 0, -5, -2).is_err());
        assert!(ScoringMatrix::new(2, -1, 0, -2).is_err());
        assert!(ScoringMatrix::new(2, -1, -5, 0).is_err());
        assert!(ScoringMatrix::new(2, -1, -5, -2).is_ok());
    }

    #[test]
    fn blosum62_diagonal_spot_checks() {
        let m = SubstitutionMatrix::blosum62();
        // A-A = 4
        assert_eq!(m.score_pair(b'A', b'A'), 4);
        // W-W = 11
        assert_eq!(m.score_pair(b'W', b'W'), 11);
        // R-R = 5
        assert_eq!(m.score_pair(b'R', b'R'), 5);
        // Case insensitive
        assert_eq!(m.score_pair(b'a', b'a'), 4);
    }

    #[test]
    fn blosum62_off_diagonal() {
        let m = SubstitutionMatrix::blosum62();
        // A-R = -1
        assert_eq!(m.score_pair(b'A', b'R'), -1);
        // Symmetric
        assert_eq!(m.score_pair(b'R', b'A'), -1);
    }

    #[test]
    fn blosum62_gap_penalties() {
        let m = SubstitutionMatrix::blosum62();
        assert_eq!(m.gap_open, -11);
        assert_eq!(m.gap_extend, -1);
    }

    #[test]
    fn blosum45_diagonal() {
        let m = SubstitutionMatrix::blosum45();
        assert_eq!(m.score_pair(b'A', b'A'), 5);
        assert_eq!(m.score_pair(b'W', b'W'), 15);
    }

    #[test]
    fn blosum80_diagonal() {
        let m = SubstitutionMatrix::blosum80();
        assert_eq!(m.score_pair(b'A', b'A'), 7);
        assert_eq!(m.score_pair(b'W', b'W'), 16);
    }

    #[test]
    fn pam250_diagonal() {
        let m = SubstitutionMatrix::pam250();
        assert_eq!(m.score_pair(b'A', b'A'), 2);
        assert_eq!(m.score_pair(b'W', b'W'), 17);
    }

    #[test]
    fn unrecognised_residue_returns_worst() {
        let m = SubstitutionMatrix::blosum62();
        let worst = m.worst_score();
        assert_eq!(m.score_pair(b'?', b'A'), worst);
    }

    #[test]
    fn scoring_scheme_delegates() {
        let dna = ScoringScheme::Simple(ScoringMatrix::dna_default());
        assert_eq!(dna.score_pair(b'A', b'A'), 2);
        assert_eq!(dna.gap_open(), -5);
        assert_eq!(dna.gap_extend(), -2);

        let protein = ScoringScheme::Substitution(SubstitutionMatrix::blosum62());
        assert_eq!(protein.score_pair(b'W', b'W'), 11);
        assert_eq!(protein.gap_open(), -11);
        assert_eq!(protein.gap_extend(), -1);
    }

    #[test]
    fn from_conversions() {
        let _scheme: ScoringScheme = ScoringMatrix::dna_default().into();
        let _scheme: ScoringScheme = SubstitutionMatrix::blosum62().into();
    }

    #[test]
    fn pam40_diagonal_spot_checks() {
        let m = SubstitutionMatrix::pam40();
        assert_eq!(m.score_pair(b'A', b'A'), 6);
        assert_eq!(m.score_pair(b'W', b'W'), 17);
        assert_eq!(m.score_pair(b'C', b'C'), 15);
    }

    #[test]
    fn pam40_symmetry() {
        let m = SubstitutionMatrix::pam40();
        assert_eq!(m.score_pair(b'A', b'R'), m.score_pair(b'R', b'A'));
        assert_eq!(m.score_pair(b'D', b'E'), m.score_pair(b'E', b'D'));
    }

    #[test]
    fn pam40_gap_penalties() {
        let m = SubstitutionMatrix::pam40();
        assert_eq!(m.gap_open, -10);
        assert_eq!(m.gap_extend, -2);
        assert_eq!(m.name(), "PAM40");
    }

    #[test]
    fn pam120_diagonal_spot_checks() {
        let m = SubstitutionMatrix::pam120();
        assert_eq!(m.score_pair(b'A', b'A'), 3);
        assert_eq!(m.score_pair(b'W', b'W'), 17);
        assert_eq!(m.score_pair(b'C', b'C'), 12);
    }

    #[test]
    fn pam120_gap_penalties() {
        let m = SubstitutionMatrix::pam120();
        assert_eq!(m.gap_open, -11);
        assert_eq!(m.gap_extend, -1);
        assert_eq!(m.name(), "PAM120");
    }

    #[test]
    fn pam200_diagonal_spot_checks() {
        let m = SubstitutionMatrix::pam200();
        assert_eq!(m.score_pair(b'A', b'A'), 2);
        assert_eq!(m.score_pair(b'W', b'W'), 17);
        assert_eq!(m.score_pair(b'Y', b'Y'), 10);
    }

    #[test]
    fn pam200_gap_penalties() {
        let m = SubstitutionMatrix::pam200();
        assert_eq!(m.gap_open, -11);
        assert_eq!(m.gap_extend, -1);
        assert_eq!(m.name(), "PAM200");
    }

    #[test]
    fn blosum30_diagonal_spot_checks() {
        let m = SubstitutionMatrix::blosum30();
        assert_eq!(m.score_pair(b'A', b'A'), 4);
        assert_eq!(m.score_pair(b'W', b'W'), 20);
        assert_eq!(m.score_pair(b'H', b'H'), 14);
        assert_eq!(m.score_pair(b'C', b'C'), 17);
    }

    #[test]
    fn blosum30_gap_penalties() {
        let m = SubstitutionMatrix::blosum30();
        assert_eq!(m.gap_open, -14);
        assert_eq!(m.gap_extend, -4);
        assert_eq!(m.name(), "BLOSUM30");
    }
}
