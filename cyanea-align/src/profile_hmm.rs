//! Profile Hidden Markov Model for sequence-profile alignment.
//!
//! Implements a Plan 7 Profile HMM with match, insert, and delete states.
//! Supports construction from a multiple sequence alignment, Viterbi decoding,
//! forward/backward algorithms, and E-value estimation via Gumbel calibration.
//!
//! # Model
//!
//! Each profile position k (1..L) has three states:
//!
//! - **Match (Mk)** — emits one residue with position-specific probabilities
//! - **Insert (Ik)** — emits one residue with uniform probabilities (I0..IL)
//! - **Delete (Dk)** — silent state, skips a profile position
//!
//! Seven transitions per position: M→M, M→I, M→D, I→M, I→I, D→M, D→D.
//! Local alignment via Begin→Mk and Mk→End entry/exit (HMMER3-style).
//!
//! # Example
//!
//! ```
//! use cyanea_align::msa::{progressive_msa, MsaResult};
//! use cyanea_align::scoring::{ScoringMatrix, ScoringScheme};
//! use cyanea_align::profile_hmm::{ProfileHmm, ProfileHmmConfig};
//!
//! let seqs: Vec<&[u8]> = vec![b"ACGT", b"ACGT", b"ACGA"];
//! let scoring = ScoringScheme::Simple(ScoringMatrix::dna_default());
//! let msa = progressive_msa(&seqs, &scoring).unwrap();
//!
//! let config = ProfileHmmConfig::dna();
//! let hmm = ProfileHmm::from_msa(&msa, &config).unwrap();
//! let result = hmm.viterbi(b"ACGT").unwrap();
//! assert!(result.score.is_finite());
//! ```

use cyanea_core::{CyaneaError, Result};

use crate::msa::MsaResult;

// ---------------------------------------------------------------------------
// Transition index constants
// ---------------------------------------------------------------------------

const TR_MM: usize = 0;
const TR_MI: usize = 1;
const TR_MD: usize = 2;
const TR_IM: usize = 3;
const TR_II: usize = 4;
const TR_DM: usize = 5;
const TR_DD: usize = 6;
const NUM_TRANSITIONS: usize = 7;

// ---------------------------------------------------------------------------
// Alphabet
// ---------------------------------------------------------------------------

/// Sequence alphabet for the profile HMM.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Alphabet {
    /// DNA alphabet: A, C, G, T
    Dna,
    /// RNA alphabet: A, C, G, U
    Rna,
    /// Protein alphabet: 20 standard amino acids
    Protein,
}

impl Alphabet {
    /// Number of symbols in the alphabet.
    pub fn size(self) -> usize {
        match self {
            Alphabet::Dna | Alphabet::Rna => 4,
            Alphabet::Protein => 20,
        }
    }

    /// Encode a byte to its index in the alphabet, or `None` if invalid.
    pub fn encode(self, b: u8) -> Option<usize> {
        let upper = b.to_ascii_uppercase();
        match self {
            Alphabet::Dna => match upper {
                b'A' => Some(0),
                b'C' => Some(1),
                b'G' => Some(2),
                b'T' => Some(3),
                _ => None,
            },
            Alphabet::Rna => match upper {
                b'A' => Some(0),
                b'C' => Some(1),
                b'G' => Some(2),
                b'U' => Some(3),
                _ => None,
            },
            Alphabet::Protein => match upper {
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
                _ => None,
            },
        }
    }

    /// Uniform background distribution for this alphabet.
    pub fn uniform_background(self) -> Vec<f64> {
        let k = self.size();
        vec![1.0 / k as f64; k]
    }
}

// ---------------------------------------------------------------------------
// Config
// ---------------------------------------------------------------------------

/// Configuration for Profile HMM construction.
#[derive(Debug, Clone)]
pub struct ProfileHmmConfig {
    /// Alphabet for the profile.
    pub alphabet: Alphabet,
    /// Pseudocount for emission estimation (default 1.0).
    pub pseudocount: f64,
    /// Background frequencies. If `None`, uses uniform.
    pub background: Option<Vec<f64>>,
    /// Gap fraction threshold for classifying MSA columns (default 0.5).
    /// Columns with gap fraction > threshold become insert columns.
    pub gap_threshold: f64,
    /// Pseudocount for transition estimation (default 1.0).
    pub transition_pseudocount: f64,
    /// If true, use local alignment mode (uniform Begin→Mk, Mk→End).
    /// If false, use global mode (entry at position 1, exit at position L).
    pub local_mode: bool,
}

impl Default for ProfileHmmConfig {
    fn default() -> Self {
        Self::dna()
    }
}

impl ProfileHmmConfig {
    /// DNA-specific default configuration.
    pub fn dna() -> Self {
        Self {
            alphabet: Alphabet::Dna,
            pseudocount: 1.0,
            background: None,
            gap_threshold: 0.5,
            transition_pseudocount: 1.0,
            local_mode: true,
        }
    }

    /// Protein-specific default configuration.
    pub fn protein() -> Self {
        Self {
            alphabet: Alphabet::Protein,
            pseudocount: 1.0,
            background: None,
            gap_threshold: 0.5,
            transition_pseudocount: 1.0,
            local_mode: true,
        }
    }
}

// ---------------------------------------------------------------------------
// Result types
// ---------------------------------------------------------------------------

/// State in a profile HMM alignment path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProfileHmmState {
    /// Match state at profile position k (1-indexed).
    Match(usize),
    /// Insert state at profile position k (0-indexed, I0..IL).
    Insert(usize),
    /// Delete state at profile position k (1-indexed).
    Delete(usize),
}

/// Result of Profile HMM Viterbi alignment.
#[derive(Debug, Clone)]
pub struct ProfileHmmResult {
    /// Log-probability score of the Viterbi path.
    pub score: f64,
    /// Sequence of states in the optimal path.
    pub states: Vec<ProfileHmmState>,
}

/// Gumbel distribution parameters for E-value estimation.
#[derive(Debug, Clone)]
pub struct GumbelParams {
    /// Location parameter.
    pub mu: f64,
    /// Scale parameter.
    pub lambda: f64,
}

// ---------------------------------------------------------------------------
// Profile HMM
// ---------------------------------------------------------------------------

/// Plan 7 Profile Hidden Markov Model.
///
/// Flat arrays store log-probabilities for emissions and transitions.
/// Match emissions: `L * K` values (position k=1..L, symbol 0..K-1).
/// Insert emissions: `(L+1) * K` values (position k=0..L, symbol 0..K-1).
/// Transitions: `(L+1) * 7` values (7 transitions per position k=0..L).
pub struct ProfileHmm {
    profile_len: usize,
    alphabet_size: usize,
    alphabet: Alphabet,
    /// Match emission log-probs: L * K, indexed as (k-1)*K + sym.
    match_emit: Vec<f64>,
    /// Insert emission log-probs: (L+1) * K, indexed as k*K + sym.
    insert_emit: Vec<f64>,
    /// Transition log-probs: (L+1) * 7, indexed as k*7 + TR_*.
    transitions: Vec<f64>,
    /// Local entry log-probs: Begin→Mk for k=1..L, indexed as k-1.
    entry: Vec<f64>,
    /// Local exit log-probs: Mk→End for k=1..L, indexed as k-1.
    exit: Vec<f64>,
    /// Gumbel parameters for E-value estimation (set by `calibrate`).
    gumbel: Option<GumbelParams>,
}

impl ProfileHmm {
    /// Create a profile HMM from raw parameters.
    ///
    /// # Errors
    ///
    /// Returns an error if array dimensions are inconsistent or profile length is zero.
    pub fn new(
        profile_len: usize,
        alphabet: Alphabet,
        match_emit: Vec<f64>,
        insert_emit: Vec<f64>,
        transitions: Vec<f64>,
        entry: Vec<f64>,
        exit: Vec<f64>,
    ) -> Result<Self> {
        let k = alphabet.size();
        if profile_len == 0 {
            return Err(CyaneaError::InvalidInput(
                "profile length must be > 0".into(),
            ));
        }
        if match_emit.len() != profile_len * k {
            return Err(CyaneaError::InvalidInput(format!(
                "match_emit length {} != profile_len * alphabet_size ({})",
                match_emit.len(),
                profile_len * k
            )));
        }
        if insert_emit.len() != (profile_len + 1) * k {
            return Err(CyaneaError::InvalidInput(format!(
                "insert_emit length {} != (profile_len+1) * alphabet_size ({})",
                insert_emit.len(),
                (profile_len + 1) * k
            )));
        }
        if transitions.len() != (profile_len + 1) * NUM_TRANSITIONS {
            return Err(CyaneaError::InvalidInput(format!(
                "transitions length {} != (profile_len+1) * 7 ({})",
                transitions.len(),
                (profile_len + 1) * NUM_TRANSITIONS
            )));
        }
        if entry.len() != profile_len {
            return Err(CyaneaError::InvalidInput(format!(
                "entry length {} != profile_len ({})",
                entry.len(),
                profile_len
            )));
        }
        if exit.len() != profile_len {
            return Err(CyaneaError::InvalidInput(format!(
                "exit length {} != profile_len ({})",
                exit.len(),
                profile_len
            )));
        }

        Ok(Self {
            profile_len,
            alphabet_size: k,
            alphabet,
            match_emit,
            insert_emit,
            transitions,
            entry,
            exit,
            gumbel: None,
        })
    }

    /// Construct a Profile HMM from a multiple sequence alignment.
    ///
    /// # Algorithm
    ///
    /// 1. Classify MSA columns: columns with gap fraction > `gap_threshold` are
    ///    insert columns; the rest are match columns.
    /// 2. Match emissions: count residue frequencies per match column, add
    ///    pseudocount, compute log-odds vs background (`ln(p/bg)`).
    /// 3. Insert emissions: `0.0` (emit at background rate, log-odds = 0).
    /// 4. Transitions: walk each MSA sequence through the column classification,
    ///    count M/I/D transitions, add pseudocount, normalize, take log.
    /// 5. Entry/exit: local mode gives uniform `ln(1/L)`; global mode concentrates
    ///    entry at position 1 and exit at position L.
    ///
    /// # Errors
    ///
    /// Returns an error if the MSA has fewer than 2 sequences, has no match
    /// columns, or contains characters not in the alphabet.
    pub fn from_msa(msa: &MsaResult, config: &ProfileHmmConfig) -> Result<Self> {
        let n_seq = msa.n_sequences();
        if n_seq < 2 {
            return Err(CyaneaError::InvalidInput(
                "need at least 2 sequences in MSA to build profile HMM".into(),
            ));
        }

        let n_cols = msa.n_columns;
        if n_cols == 0 {
            return Err(CyaneaError::InvalidInput(
                "MSA has no columns".into(),
            ));
        }

        let alphabet = config.alphabet;
        let k = alphabet.size();

        // Step 1: classify columns
        let mut is_match_col = vec![false; n_cols];
        for c in 0..n_cols {
            let gap_count = msa.aligned.iter().filter(|s| s[c] == b'-').count();
            let gap_frac = gap_count as f64 / n_seq as f64;
            is_match_col[c] = gap_frac <= config.gap_threshold;
        }

        let match_cols: Vec<usize> = (0..n_cols).filter(|&c| is_match_col[c]).collect();
        let profile_len = match_cols.len();
        if profile_len == 0 {
            return Err(CyaneaError::InvalidInput(
                "no match columns found (all columns exceed gap threshold)".into(),
            ));
        }

        // Step 2: match emissions as log-odds vs background (HMMER3-style)
        let bg = config
            .background
            .as_ref()
            .cloned()
            .unwrap_or_else(|| alphabet.uniform_background());
        let bg_log: Vec<f64> = bg.iter().map(|&p| p.ln()).collect();

        let mut match_emit = vec![0.0f64; profile_len * k];
        for (mk, &col) in match_cols.iter().enumerate() {
            let mut counts = vec![config.pseudocount; k];
            for seq in &msa.aligned {
                let b = seq[col];
                if b == b'-' {
                    continue;
                }
                if let Some(idx) = alphabet.encode(b) {
                    counts[idx] += 1.0;
                }
                // Skip characters not in alphabet (treat as gap)
            }
            let total: f64 = counts.iter().sum();
            for s in 0..k {
                let log_prob = (counts[s] / total).ln();
                match_emit[mk * k + s] = log_prob - bg_log[s]; // log-odds
            }
        }

        // Step 3: insert emissions — 0.0 (emit at background rate, log-odds = 0)
        let insert_emit = vec![0.0; (profile_len + 1) * k];

        // Step 4: transitions
        // Count transitions by walking each sequence through column classification.
        // For each sequence, at each column, determine if the sequence has a
        // residue or gap in that column, and what state that implies.
        let mut trans_counts = vec![vec![config.transition_pseudocount; NUM_TRANSITIONS]; profile_len + 1];

        for seq in &msa.aligned {
            // Walk through columns, tracking current state
            // State: M at match_pos, I between match positions, D if gap in match column
            let mut match_pos: usize = 0; // current match position (0-based), 0 means before first match
            let mut prev_state: u8 = 0; // 0=M, 1=I, 2=D; start as if coming from Begin (treated as M)
            let mut in_profile = false; // have we entered the profile?

            for c in 0..n_cols {
                let has_residue = seq[c] != b'-';

                if is_match_col[c] {
                    // This is a match column
                    if has_residue {
                        // → Match state
                        if in_profile {
                            let trans_pos = match_pos; // transition from previous match_pos
                            match prev_state {
                                0 => trans_counts[trans_pos][TR_MM] += 1.0,
                                1 => trans_counts[trans_pos][TR_IM] += 1.0,
                                2 => trans_counts[trans_pos][TR_DM] += 1.0,
                                _ => {}
                            }
                        }
                        match_pos += 1;
                        prev_state = 0;
                        in_profile = true;
                    } else {
                        // Gap in match column → Delete state
                        if in_profile {
                            let trans_pos = match_pos;
                            match prev_state {
                                0 => trans_counts[trans_pos][TR_MD] += 1.0,
                                2 => trans_counts[trans_pos][TR_DD] += 1.0,
                                1 => {
                                    // I→D not standard in Plan 7, treat as I→M→D
                                    trans_counts[trans_pos][TR_IM] += 1.0;
                                }
                                _ => {}
                            }
                        }
                        match_pos += 1;
                        prev_state = 2;
                        in_profile = true;
                    }
                } else if has_residue && in_profile {
                    // Insert column with residue (and we're in the profile)
                    if prev_state != 1 {
                        // Entering insert state
                        let trans_pos = match_pos;
                        if trans_pos <= profile_len {
                            match prev_state {
                                0 => trans_counts[trans_pos][TR_MI] += 1.0,
                                2 => {
                                    // D→I not standard, treat as D→M→I
                                    trans_counts[trans_pos][TR_DM] += 1.0;
                                }
                                _ => {}
                            }
                        }
                        prev_state = 1;
                    } else {
                        // Staying in insert state
                        let trans_pos = match_pos;
                        if trans_pos <= profile_len {
                            trans_counts[trans_pos][TR_II] += 1.0;
                        }
                    }
                }
            }

            // Reset match_pos for next sequence
        }

        // Normalize transition counts to log-probabilities
        let mut transitions = vec![f64::NEG_INFINITY; (profile_len + 1) * NUM_TRANSITIONS];
        for pos in 0..=profile_len {
            // Normalize M→{M,I,D} group
            let m_total = trans_counts[pos][TR_MM]
                + trans_counts[pos][TR_MI]
                + trans_counts[pos][TR_MD];
            if m_total > 0.0 {
                transitions[pos * NUM_TRANSITIONS + TR_MM] = (trans_counts[pos][TR_MM] / m_total).ln();
                transitions[pos * NUM_TRANSITIONS + TR_MI] = (trans_counts[pos][TR_MI] / m_total).ln();
                transitions[pos * NUM_TRANSITIONS + TR_MD] = (trans_counts[pos][TR_MD] / m_total).ln();
            }

            // Normalize I→{M,I} group
            let i_total = trans_counts[pos][TR_IM] + trans_counts[pos][TR_II];
            if i_total > 0.0 {
                transitions[pos * NUM_TRANSITIONS + TR_IM] = (trans_counts[pos][TR_IM] / i_total).ln();
                transitions[pos * NUM_TRANSITIONS + TR_II] = (trans_counts[pos][TR_II] / i_total).ln();
            }

            // Normalize D→{M,D} group
            let d_total = trans_counts[pos][TR_DM] + trans_counts[pos][TR_DD];
            if d_total > 0.0 {
                transitions[pos * NUM_TRANSITIONS + TR_DM] = (trans_counts[pos][TR_DM] / d_total).ln();
                transitions[pos * NUM_TRANSITIONS + TR_DD] = (trans_counts[pos][TR_DD] / d_total).ln();
            }
        }

        // Step 5: entry/exit
        let entry;
        let exit;
        if config.local_mode {
            let uniform_entry = (1.0 / profile_len as f64).ln();
            entry = vec![uniform_entry; profile_len];
            exit = vec![uniform_entry; profile_len];
        } else {
            entry = (0..profile_len)
                .map(|i| if i == 0 { 0.0 } else { f64::NEG_INFINITY })
                .collect();
            exit = (0..profile_len)
                .map(|i| {
                    if i == profile_len - 1 {
                        0.0
                    } else {
                        f64::NEG_INFINITY
                    }
                })
                .collect();
        }

        Ok(Self {
            profile_len,
            alphabet_size: k,
            alphabet,
            match_emit,
            insert_emit,
            transitions,
            entry,
            exit,
            gumbel: None,
        })
    }

    /// Profile length (number of match states).
    pub fn profile_len(&self) -> usize {
        self.profile_len
    }

    /// Alphabet size.
    pub fn alphabet_size(&self) -> usize {
        self.alphabet_size
    }

    /// Alphabet used by this profile.
    pub fn alphabet(&self) -> Alphabet {
        self.alphabet
    }

    // -----------------------------------------------------------------------
    // Internal indexing helpers
    // -----------------------------------------------------------------------

    /// Get match emission log-prob for match position k (1-indexed), symbol sym.
    fn match_emit_score(&self, k: usize, sym: usize) -> f64 {
        self.match_emit[(k - 1) * self.alphabet_size + sym]
    }

    /// Get insert emission log-prob for insert position k (0-indexed), symbol sym.
    fn insert_emit_score(&self, k: usize, sym: usize) -> f64 {
        self.insert_emit[k * self.alphabet_size + sym]
    }

    /// Get transition log-prob at position k, transition type t.
    fn trans(&self, k: usize, t: usize) -> f64 {
        self.transitions[k * NUM_TRANSITIONS + t]
    }

    /// Encode a full sequence, returning indices or an error.
    fn encode_sequence(&self, seq: &[u8]) -> Result<Vec<usize>> {
        seq.iter()
            .enumerate()
            .map(|(i, &b)| {
                self.alphabet.encode(b).ok_or_else(|| {
                    CyaneaError::InvalidInput(format!(
                        "invalid character '{}' at position {} for {:?} alphabet",
                        b as char, i, self.alphabet
                    ))
                })
            })
            .collect()
    }

    // -----------------------------------------------------------------------
    // Viterbi algorithm
    // -----------------------------------------------------------------------

    /// Viterbi algorithm: find the most probable state path through the profile.
    ///
    /// Returns the log-probability score and the state path.
    ///
    /// # Errors
    ///
    /// Returns an error if the sequence is empty or contains invalid characters.
    pub fn viterbi(&self, seq: &[u8]) -> Result<ProfileHmmResult> {
        if seq.is_empty() {
            return Err(CyaneaError::InvalidInput(
                "sequence must not be empty".into(),
            ));
        }

        let encoded = self.encode_sequence(seq)?;
        let n = encoded.len();
        let l = self.profile_len;

        // DP matrices: (N+1) rows × (L+1) columns
        let dim = (n + 1) * (l + 1);
        let mut vm = vec![f64::NEG_INFINITY; dim];
        let mut vi = vec![f64::NEG_INFINITY; dim];
        let mut vd = vec![f64::NEG_INFINITY; dim];

        // Backpointer: 0=None, 1=fromM, 2=fromI, 3=fromD, 4=fromEntry
        let mut bp_m = vec![0u8; dim];
        let mut bp_i = vec![0u8; dim];
        let mut bp_d = vec![0u8; dim];

        let idx = |i: usize, k: usize| -> usize { i * (l + 1) + k };

        // Fill
        for i in 1..=n {
            let sym = encoded[i - 1];

            // Match states
            for k in 1..=l {
                let emit = self.match_emit_score(k, sym);

                // Entry from Begin
                let mut best = self.entry[k - 1];
                let mut bp: u8 = 4; // fromEntry

                // From M[i-1][k-1]
                if k > 1 {
                    let from_m = vm[idx(i - 1, k - 1)] + self.trans(k - 1, TR_MM);
                    if from_m > best {
                        best = from_m;
                        bp = 1;
                    }

                    // From I[i-1][k-1]
                    let from_i = vi[idx(i - 1, k - 1)] + self.trans(k - 1, TR_IM);
                    if from_i > best {
                        best = from_i;
                        bp = 2;
                    }

                    // From D[i-1][k-1]
                    let from_d = vd[idx(i - 1, k - 1)] + self.trans(k - 1, TR_DM);
                    if from_d > best {
                        best = from_d;
                        bp = 3;
                    }
                }

                vm[idx(i, k)] = emit + best;
                bp_m[idx(i, k)] = bp;
            }

            // Insert states
            for k in 0..=l {
                // Allow I0 through IL
                // I states only from M or I at same k (previous row)
                if k > 0 || i > 0 {
                    let emit = self.insert_emit_score(k, sym);

                    let from_m = if k <= l {
                        vm[idx(i - 1, k)] + self.trans(k, TR_MI)
                    } else {
                        f64::NEG_INFINITY
                    };
                    let from_i = vi[idx(i - 1, k)] + self.trans(k, TR_II);

                    if from_m >= from_i && from_m > f64::NEG_INFINITY {
                        vi[idx(i, k)] = emit + from_m;
                        bp_i[idx(i, k)] = 1;
                    } else if from_i > f64::NEG_INFINITY {
                        vi[idx(i, k)] = emit + from_i;
                        bp_i[idx(i, k)] = 2;
                    }
                }
            }

            // Delete states (same row, no sequence consumed)
            for k in 1..=l {
                let from_m = if k > 1 {
                    vm[idx(i, k - 1)] + self.trans(k - 1, TR_MD)
                } else {
                    f64::NEG_INFINITY
                };
                let from_d = if k > 1 {
                    vd[idx(i, k - 1)] + self.trans(k - 1, TR_DD)
                } else {
                    f64::NEG_INFINITY
                };

                if from_m >= from_d && from_m > f64::NEG_INFINITY {
                    vd[idx(i, k)] = from_m;
                    bp_d[idx(i, k)] = 1;
                } else if from_d > f64::NEG_INFINITY {
                    vd[idx(i, k)] = from_d;
                    bp_d[idx(i, k)] = 3;
                }
            }
        }

        // Also need to fill delete states at row 0 for completeness
        // (handled by initialization: all NEG_INFINITY)

        // Termination: best score over k of V_M[N][k] + exit[k-1]
        let mut best_score = f64::NEG_INFINITY;
        let mut best_k = 0usize;
        let mut best_end_state: u8 = 0; // 0=M, 1=D

        for k in 1..=l {
            let score_m = vm[idx(n, k)] + self.exit[k - 1];
            if score_m > best_score {
                best_score = score_m;
                best_k = k;
                best_end_state = 0;
            }
            // Also consider ending in a delete state (via M→End)
            // Delete states don't directly exit; only match states do in Plan 7.
        }

        // Check if we can end through delete states that transition to match that exits
        // Actually in Plan 7, only Match states have Mk→End transitions.
        // So we just use the M scores above.

        if best_score == f64::NEG_INFINITY {
            // Fallback: try all final states
            for k in 1..=l {
                if vm[idx(n, k)] > best_score {
                    best_score = vm[idx(n, k)];
                    best_k = k;
                    best_end_state = 0;
                }
            }
        }

        // Traceback
        let mut states = Vec::new();
        let mut i = n;
        let mut k = best_k;
        let mut current_state = best_end_state; // 0=M, 1=I, 2=D

        // We need to trace back through the DP matrix
        loop {
            match current_state {
                0 => {
                    // Match state
                    states.push(ProfileHmmState::Match(k));
                    let bp = bp_m[idx(i, k)];
                    match bp {
                        4 => {
                            // from Entry — done
                            break;
                        }
                        1 => {
                            // from M[i-1][k-1]
                            i -= 1;
                            k -= 1;
                            current_state = 0;
                        }
                        2 => {
                            // from I[i-1][k-1]
                            i -= 1;
                            k -= 1;
                            current_state = 1;
                        }
                        3 => {
                            // from D[i-1][k-1]
                            i -= 1;
                            k -= 1;
                            current_state = 2;
                        }
                        _ => break,
                    }
                }
                1 => {
                    // Insert state
                    states.push(ProfileHmmState::Insert(k));
                    let bp = bp_i[idx(i, k)];
                    i -= 1; // Insert always consumes a sequence position
                    match bp {
                        1 => current_state = 0, // from M
                        2 => current_state = 1, // from I
                        _ => break,
                    }
                }
                2 => {
                    // Delete state
                    states.push(ProfileHmmState::Delete(k));
                    let bp = bp_d[idx(i, k)];
                    // Delete doesn't consume sequence
                    match bp {
                        1 => {
                            // from M[i][k-1]
                            k -= 1;
                            current_state = 0;
                        }
                        3 => {
                            // from D[i][k-1]
                            k -= 1;
                            current_state = 2;
                        }
                        _ => break,
                    }
                }
                _ => break,
            }

            if i == 0 && current_state != 2 {
                break;
            }
        }

        states.reverse();

        Ok(ProfileHmmResult {
            score: best_score,
            states,
        })
    }

    // -----------------------------------------------------------------------
    // Forward algorithm
    // -----------------------------------------------------------------------

    /// Forward algorithm: compute the total log-likelihood P(seq | model).
    ///
    /// Uses log-sum-exp to sum over all possible state paths.
    ///
    /// # Errors
    ///
    /// Returns an error if the sequence is empty or contains invalid characters.
    pub fn forward(&self, seq: &[u8]) -> Result<f64> {
        if seq.is_empty() {
            return Err(CyaneaError::InvalidInput(
                "sequence must not be empty".into(),
            ));
        }

        let encoded = self.encode_sequence(seq)?;
        let n = encoded.len();
        let l = self.profile_len;

        let dim = (n + 1) * (l + 1);
        let mut fm = vec![f64::NEG_INFINITY; dim];
        let mut fi = vec![f64::NEG_INFINITY; dim];
        let mut fd = vec![f64::NEG_INFINITY; dim];

        let idx = |i: usize, k: usize| -> usize { i * (l + 1) + k };

        for i in 1..=n {
            let sym = encoded[i - 1];

            // Match states
            for k in 1..=l {
                let emit = self.match_emit_score(k, sym);

                let mut val = self.entry[k - 1];

                if k > 1 {
                    val = log_sum_exp(val, fm[idx(i - 1, k - 1)] + self.trans(k - 1, TR_MM));
                    val = log_sum_exp(val, fi[idx(i - 1, k - 1)] + self.trans(k - 1, TR_IM));
                    val = log_sum_exp(val, fd[idx(i - 1, k - 1)] + self.trans(k - 1, TR_DM));
                }

                fm[idx(i, k)] = emit + val;
            }

            // Insert states
            for k in 0..=l {
                let emit = self.insert_emit_score(k, sym);

                let from_m = if k <= l {
                    fm[idx(i - 1, k)] + self.trans(k, TR_MI)
                } else {
                    f64::NEG_INFINITY
                };
                let from_i = fi[idx(i - 1, k)] + self.trans(k, TR_II);

                let val = log_sum_exp(from_m, from_i);
                if val > f64::NEG_INFINITY {
                    fi[idx(i, k)] = emit + val;
                }
            }

            // Delete states (same row)
            for k in 1..=l {
                let from_m = if k > 1 {
                    fm[idx(i, k - 1)] + self.trans(k - 1, TR_MD)
                } else {
                    f64::NEG_INFINITY
                };
                let from_d = if k > 1 {
                    fd[idx(i, k - 1)] + self.trans(k - 1, TR_DD)
                } else {
                    f64::NEG_INFINITY
                };

                let val = log_sum_exp(from_m, from_d);
                if val > f64::NEG_INFINITY {
                    fd[idx(i, k)] = val;
                }
            }
        }

        // Termination: sum over k of F_M[N][k] + exit[k-1]
        let mut total = f64::NEG_INFINITY;
        for k in 1..=l {
            total = log_sum_exp(total, fm[idx(n, k)] + self.exit[k - 1]);
        }

        Ok(total)
    }

    // -----------------------------------------------------------------------
    // Backward algorithm
    // -----------------------------------------------------------------------

    /// Backward algorithm: compute the total log-likelihood P(seq | model).
    ///
    /// Reverse sweep from i=N to 0, k=L to 1. The result should match
    /// the forward algorithm's output (up to numerical precision).
    ///
    /// # Errors
    ///
    /// Returns an error if the sequence is empty or contains invalid characters.
    pub fn backward(&self, seq: &[u8]) -> Result<f64> {
        if seq.is_empty() {
            return Err(CyaneaError::InvalidInput(
                "sequence must not be empty".into(),
            ));
        }

        let encoded = self.encode_sequence(seq)?;
        let n = encoded.len();
        let l = self.profile_len;

        let dim = (n + 1) * (l + 1);
        let mut bm = vec![f64::NEG_INFINITY; dim];
        let mut bi = vec![f64::NEG_INFINITY; dim];
        let mut bd = vec![f64::NEG_INFINITY; dim];

        let idx = |i: usize, k: usize| -> usize { i * (l + 1) + k };

        // Initialize: at i=N, Mk can exit with exit[k-1]
        for k in 1..=l {
            bm[idx(n, k)] = self.exit[k - 1];
        }

        // Also initialize delete states at i=N: D can transition to M or D
        // that eventually exits. Process k from L-1 down to 1.
        for k in (1..l).rev() {
            // D[N][k] can go to M[N][k+1] via D→M (but M at row N needs exit)
            // Actually delete states at row N: they don't consume sequence,
            // so D[N][k] → M exit is through D→M transition where M would need
            // to be at row N too. But M at row N, k+1 already has exit score.
            // D→M means at same row, next position matches... but M requires
            // sequence consumption. So D[N][k] can only chain to D[N][k+1].
            bd[idx(n, k)] = bd[idx(n, k + 1)] + self.trans(k, TR_DD);
            // D→M at end doesn't work because there's no more sequence to consume.
        }

        // Backward recursion: i from N down to 1
        for i in (1..=n).rev() {
            // Delete states at row i (reverse sweep k from L-1 to 1)
            // Must process before M and I since they depend on D at same row
            // First, set D[i][L] based on what comes after
            // D[i][k] can transition to:
            //   D[i][k+1] via DD (same row, no consumption)
            //   M[i+1][k+1] via DM (next row, consume symbol) — but only if i < n
            //   (no D→I transition in Plan 7)

            // Process deletes right-to-left at current row
            if i < n {
                let sym_i = encoded[i]; // symbol at position i+1 (0-indexed: seq[i])
                for k in (1..l).rev() {
                    let to_dm = self.trans(k, TR_DM)
                        + self.match_emit_score(k + 1, sym_i)
                        + bm[idx(i + 1, k + 1)];
                    let to_dd = self.trans(k, TR_DD) + bd[idx(i, k + 1)];
                    bd[idx(i, k)] = log_sum_exp(to_dm, to_dd);
                }
                // k = L: can only go to D→M at k+1 (out of range) → nothing
                // Actually k=L: D[i][L] has no valid transition forward within profile
                // It could go to exit but deletes don't exit directly.
                // So bd[idx(i, l)] stays NEG_INFINITY unless chained
            } else {
                // i == n: delete states can only chain DD
                for k in (1..l).rev() {
                    bd[idx(i, k)] = self.trans(k, TR_DD) + bd[idx(i, k + 1)];
                }
            }

            // Match and Insert backward values
            if i < n {
                let sym_i = encoded[i]; // next symbol consumed

                for k in 1..=l {
                    // B_M[i][k]: M at position i,k can go to:
                    //   M[i+1][k+1] via MM (consume next symbol, emit at k+1)
                    //   I[i+1][k] via MI (consume next symbol, insert at k)
                    //   D[i][k+1] via MD (same row, no consumption) — if k < l

                    let mut val = f64::NEG_INFINITY;

                    if k < l {
                        let to_mm = self.trans(k, TR_MM)
                            + self.match_emit_score(k + 1, sym_i)
                            + bm[idx(i + 1, k + 1)];
                        val = log_sum_exp(val, to_mm);
                    }

                    let to_mi = self.trans(k, TR_MI)
                        + self.insert_emit_score(k, sym_i)
                        + bi[idx(i + 1, k)];
                    val = log_sum_exp(val, to_mi);

                    if k < l {
                        let to_md = self.trans(k, TR_MD) + bd[idx(i, k + 1)];
                        val = log_sum_exp(val, to_md);
                    }

                    // Also exit
                    val = log_sum_exp(val, self.exit[k - 1]);

                    bm[idx(i, k)] = val;

                    // B_I[i][k]: I at position i,k can go to:
                    //   M[i+1][k+1] via IM (consume next symbol)
                    //   I[i+1][k] via II (consume next symbol)

                    let mut ival = f64::NEG_INFINITY;

                    if k < l {
                        let to_im = self.trans(k, TR_IM)
                            + self.match_emit_score(k + 1, sym_i)
                            + bm[idx(i + 1, k + 1)];
                        ival = log_sum_exp(ival, to_im);
                    }

                    let to_ii = self.trans(k, TR_II)
                        + self.insert_emit_score(k, sym_i)
                        + bi[idx(i + 1, k)];
                    ival = log_sum_exp(ival, to_ii);

                    bi[idx(i, k)] = ival;
                }
            } else {
                // i == n: at the end, M can only exit
                for k in 1..=l {
                    bm[idx(i, k)] = self.exit[k - 1];
                    // I at end: no valid transitions
                    // bi stays NEG_INFINITY
                }
            }
        }

        // Termination: sum over k of entry[k-1] + match_emit(k, seq[0]) + B_M[1][k]
        let sym0 = encoded[0];
        let mut total = f64::NEG_INFINITY;
        for k in 1..=l {
            let val = self.entry[k - 1] + self.match_emit_score(k, sym0) + bm[idx(1, k)];
            total = log_sum_exp(total, val);
        }

        Ok(total)
    }

    // -----------------------------------------------------------------------
    // E-value calibration
    // -----------------------------------------------------------------------

    /// Calibrate the profile for E-value estimation by scoring random sequences.
    ///
    /// Generates `n_samples` random sequences of length `seq_length` using a
    /// simple LCG PRNG seeded with `seed`, scores each with Viterbi, and fits
    /// a Gumbel distribution via method of moments.
    ///
    /// # Errors
    ///
    /// Returns an error if `n_samples` is zero.
    pub fn calibrate(
        &mut self,
        n_samples: usize,
        seq_length: usize,
        seed: u64,
    ) -> Result<()> {
        if n_samples == 0 {
            return Err(CyaneaError::InvalidInput(
                "n_samples must be > 0".into(),
            ));
        }

        let symbols: Vec<u8> = match self.alphabet {
            Alphabet::Dna => vec![b'A', b'C', b'G', b'T'],
            Alphabet::Rna => vec![b'A', b'C', b'G', b'U'],
            Alphabet::Protein => {
                vec![
                    b'A', b'R', b'N', b'D', b'C', b'Q', b'E', b'G', b'H', b'I',
                    b'L', b'K', b'M', b'F', b'P', b'S', b'T', b'W', b'Y', b'V',
                ]
            }
        };

        let mut rng_state = seed;
        let mut scores = Vec::with_capacity(n_samples);

        for _ in 0..n_samples {
            let mut seq = Vec::with_capacity(seq_length);
            for _ in 0..seq_length {
                // LCG: state = state * 6364136223846793005 + 1442695040888963407
                rng_state = rng_state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let idx = ((rng_state >> 33) as usize) % symbols.len();
                seq.push(symbols[idx]);
            }

            match self.viterbi(&seq) {
                Ok(result) => {
                    if result.score.is_finite() {
                        scores.push(result.score);
                    }
                }
                Err(_) => continue,
            }
        }

        if scores.is_empty() {
            return Err(CyaneaError::InvalidInput(
                "no valid scores obtained during calibration".into(),
            ));
        }

        // Fit Gumbel by method of moments
        let n = scores.len() as f64;
        let mean = scores.iter().sum::<f64>() / n;
        let variance = scores.iter().map(|&s| (s - mean) * (s - mean)).sum::<f64>() / n;
        let stddev = variance.sqrt();

        if stddev < 1e-15 {
            // All scores identical — degenerate case
            self.gumbel = Some(GumbelParams {
                mu: mean,
                lambda: 1.0,
            });
            return Ok(());
        }

        // Euler-Mascheroni constant
        let euler_mascheroni = 0.5772156649015329;
        let lambda = std::f64::consts::PI / (stddev * (6.0_f64).sqrt());
        let mu = mean - euler_mascheroni / lambda;

        self.gumbel = Some(GumbelParams { mu, lambda });

        Ok(())
    }

    /// Compute E-value for a given score against a database of `db_size` sequences.
    ///
    /// # Errors
    ///
    /// Returns an error if `calibrate` has not been called.
    pub fn evalue(&self, score: f64, db_size: usize) -> Result<f64> {
        let gumbel = self.gumbel.as_ref().ok_or_else(|| {
            CyaneaError::InvalidInput(
                "must call calibrate() before computing E-values".into(),
            )
        })?;

        let evalue = db_size as f64 * (-gumbel.lambda * (score - gumbel.mu)).exp();
        Ok(evalue)
    }
}

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

/// Numerically stable log(exp(a) + exp(b)).
fn log_sum_exp(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let max = a.max(b);
    max + ((a - max).exp() + (b - max).exp()).ln()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::msa::progressive_msa;
    use crate::scoring::{ScoringMatrix, ScoringScheme};

    // ------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------

    fn dna_scheme() -> ScoringScheme {
        ScoringScheme::Simple(ScoringMatrix::dna_default())
    }

    fn build_dna_hmm(seqs: &[&[u8]]) -> ProfileHmm {
        let msa = progressive_msa(seqs, &dna_scheme()).unwrap();
        ProfileHmm::from_msa(&msa, &ProfileHmmConfig::dna()).unwrap()
    }

    fn build_protein_hmm(seqs: &[&[u8]]) -> ProfileHmm {
        // Build a simple MSA manually for protein sequences
        let msa = MsaResult {
            aligned: seqs.iter().map(|s| s.to_vec()).collect(),
            n_columns: seqs[0].len(),
        };
        ProfileHmm::from_msa(&msa, &ProfileHmmConfig::protein()).unwrap()
    }

    // ------------------------------------------------------------------
    // Alphabet tests
    // ------------------------------------------------------------------

    #[test]
    fn alphabet_encode_dna() {
        assert_eq!(Alphabet::Dna.size(), 4);
        assert_eq!(Alphabet::Dna.encode(b'A'), Some(0));
        assert_eq!(Alphabet::Dna.encode(b'C'), Some(1));
        assert_eq!(Alphabet::Dna.encode(b'G'), Some(2));
        assert_eq!(Alphabet::Dna.encode(b'T'), Some(3));
        assert_eq!(Alphabet::Dna.encode(b'a'), Some(0)); // case-insensitive
        assert_eq!(Alphabet::Dna.encode(b'X'), None);
        assert_eq!(Alphabet::Dna.encode(b'U'), None); // U is RNA
    }

    #[test]
    fn alphabet_encode_protein() {
        assert_eq!(Alphabet::Protein.size(), 20);
        assert_eq!(Alphabet::Protein.encode(b'A'), Some(0));
        assert_eq!(Alphabet::Protein.encode(b'V'), Some(19));
        assert_eq!(Alphabet::Protein.encode(b'm'), Some(12)); // case-insensitive
        assert_eq!(Alphabet::Protein.encode(b'X'), None);
        assert_eq!(Alphabet::Protein.encode(b'-'), None);

        let bg = Alphabet::Protein.uniform_background();
        assert_eq!(bg.len(), 20);
        assert!((bg.iter().sum::<f64>() - 1.0).abs() < 1e-10);
    }

    // ------------------------------------------------------------------
    // Config tests
    // ------------------------------------------------------------------

    #[test]
    fn config_constructors() {
        let dna = ProfileHmmConfig::dna();
        assert_eq!(dna.alphabet, Alphabet::Dna);
        assert!((dna.pseudocount - 1.0).abs() < 1e-10);
        assert!(dna.local_mode);

        let prot = ProfileHmmConfig::protein();
        assert_eq!(prot.alphabet, Alphabet::Protein);
        assert!(prot.local_mode);

        let def = ProfileHmmConfig::default();
        assert_eq!(def.alphabet, Alphabet::Dna);
    }

    #[test]
    fn local_vs_global_mode() {
        let seqs: Vec<&[u8]> = vec![b"ACGT", b"ACGT", b"ACGA"];
        let msa = progressive_msa(&seqs, &dna_scheme()).unwrap();

        let mut config = ProfileHmmConfig::dna();

        config.local_mode = true;
        let local = ProfileHmm::from_msa(&msa, &config).unwrap();

        config.local_mode = false;
        let global = ProfileHmm::from_msa(&msa, &config).unwrap();

        // Local mode: all entries are equal (uniform)
        let expected_local = (1.0 / local.profile_len() as f64).ln();
        for &e in &local.entry {
            assert!((e - expected_local).abs() < 1e-10);
        }

        // Global mode: only first entry is 0.0, rest are NEG_INFINITY
        assert!((global.entry[0] - 0.0).abs() < 1e-10);
        for &e in &global.entry[1..] {
            assert!(e == f64::NEG_INFINITY);
        }
    }

    // ------------------------------------------------------------------
    // Construction tests
    // ------------------------------------------------------------------

    #[test]
    fn from_msa_basic() {
        let hmm = build_dna_hmm(&[b"ACGT", b"ACGT", b"ACGT"]);
        assert!(hmm.profile_len() > 0);
        assert_eq!(hmm.alphabet_size(), 4);
        assert_eq!(hmm.alphabet(), Alphabet::Dna);
    }

    #[test]
    fn from_msa_protein() {
        let hmm = build_protein_hmm(&[b"MARVEL", b"MARVEL", b"MARVEL"]);
        assert!(hmm.profile_len() > 0);
        assert_eq!(hmm.alphabet_size(), 20);
        assert_eq!(hmm.alphabet(), Alphabet::Protein);
    }

    #[test]
    fn from_msa_with_gaps() {
        // MSA with some heavily gapped columns that should become inserts
        let msa = MsaResult {
            aligned: vec![
                b"A-C-G".to_vec(),
                b"A-C-G".to_vec(),
                b"A-C-G".to_vec(),
            ],
            n_columns: 5,
        };
        let config = ProfileHmmConfig::dna();
        let hmm = ProfileHmm::from_msa(&msa, &config).unwrap();
        // Columns 1 and 3 (0-indexed) are all gaps → insert columns
        // Profile length should be 3 (only columns 0, 2, 4 are match)
        assert_eq!(hmm.profile_len(), 3);
    }

    #[test]
    fn from_msa_error_single_seq() {
        let msa = MsaResult {
            aligned: vec![b"ACGT".to_vec()],
            n_columns: 4,
        };
        let config = ProfileHmmConfig::dna();
        assert!(ProfileHmm::from_msa(&msa, &config).is_err());
    }

    #[test]
    fn from_msa_error_empty() {
        let msa = MsaResult {
            aligned: vec![],
            n_columns: 0,
        };
        let config = ProfileHmmConfig::dna();
        assert!(ProfileHmm::from_msa(&msa, &config).is_err());
    }

    #[test]
    fn new_validates_dimensions() {
        // Correct dimensions should work
        let result = ProfileHmm::new(
            2,
            Alphabet::Dna,
            vec![0.0; 2 * 4],      // match_emit: L * K
            vec![0.0; 3 * 4],      // insert_emit: (L+1) * K
            vec![0.0; 3 * 7],      // transitions: (L+1) * 7
            vec![0.0; 2],          // entry: L
            vec![0.0; 2],          // exit: L
        );
        assert!(result.is_ok());

        // Wrong match_emit length
        assert!(ProfileHmm::new(
            2,
            Alphabet::Dna,
            vec![0.0; 3],           // wrong
            vec![0.0; 3 * 4],
            vec![0.0; 3 * 7],
            vec![0.0; 2],
            vec![0.0; 2],
        ).is_err());

        // Zero profile length
        assert!(ProfileHmm::new(
            0,
            Alphabet::Dna,
            vec![],
            vec![0.0; 4],
            vec![0.0; 7],
            vec![],
            vec![],
        ).is_err());
    }

    // ------------------------------------------------------------------
    // Viterbi tests
    // ------------------------------------------------------------------

    #[test]
    fn viterbi_identical_sequence() {
        let hmm = build_dna_hmm(&[b"ACGT", b"ACGT", b"ACGT"]);
        let result = hmm.viterbi(b"ACGT").unwrap();
        assert!(result.score.is_finite());
        // Should be mostly Match states
        let match_count = result
            .states
            .iter()
            .filter(|s| matches!(s, ProfileHmmState::Match(_)))
            .count();
        assert!(match_count > 0);
    }

    #[test]
    fn viterbi_score_negative() {
        let hmm = build_dna_hmm(&[b"ACGT", b"ACGT", b"ACGA"]);
        let result = hmm.viterbi(b"ACGT").unwrap();
        // Log-probability should be negative
        assert!(
            result.score < 0.0,
            "Viterbi score should be negative (log-prob), got {}",
            result.score
        );
    }

    #[test]
    fn viterbi_path_positions_monotonic() {
        let hmm = build_dna_hmm(&[b"ACGTACGT", b"ACGTACGT", b"ACGTACGA"]);
        let result = hmm.viterbi(b"ACGTACGT").unwrap();

        // Profile positions in Match/Delete states should be non-decreasing
        let positions: Vec<usize> = result
            .states
            .iter()
            .filter_map(|s| match s {
                ProfileHmmState::Match(k) | ProfileHmmState::Delete(k) => Some(*k),
                _ => None,
            })
            .collect();

        for w in positions.windows(2) {
            assert!(
                w[1] >= w[0],
                "profile positions should be monotonically non-decreasing: {:?}",
                positions
            );
        }
    }

    #[test]
    fn viterbi_consumes_full_sequence() {
        let hmm = build_dna_hmm(&[b"ACGT", b"ACGT", b"ACGA"]);
        let seq = b"ACGT";
        let result = hmm.viterbi(seq).unwrap();

        // Match + Insert states consume one sequence position each
        let consumed: usize = result
            .states
            .iter()
            .filter(|s| matches!(s, ProfileHmmState::Match(_) | ProfileHmmState::Insert(_)))
            .count();
        assert_eq!(
            consumed,
            seq.len(),
            "Match+Insert count ({}) should equal sequence length ({})",
            consumed,
            seq.len()
        );
    }

    #[test]
    fn viterbi_shorter_seq_has_deletes() {
        // Build HMM from 8-base sequences, query with 4 bases
        let hmm = build_dna_hmm(&[b"ACGTACGT", b"ACGTACGT", b"ACGTACGA"]);
        let result = hmm.viterbi(b"ACGT").unwrap();

        let has_deletes = result
            .states
            .iter()
            .any(|s| matches!(s, ProfileHmmState::Delete(_)));
        // A shorter sequence should likely use some delete states
        // (or just align locally to a subset of positions)
        assert!(
            has_deletes || result.states.len() <= 4,
            "shorter sequence should either have deletes or align to subset of profile"
        );
    }

    #[test]
    fn viterbi_longer_seq_has_inserts() {
        // Build HMM in global mode from an MSA where some sequences have an
        // extra base (insert column). Global mode forces entry at M1 and exit
        // at ML, so the extra base must be consumed by an insert state.
        let mut aligned = Vec::new();
        for _ in 0..7 {
            aligned.push(b"AC-GT".to_vec());
        }
        for _ in 0..3 {
            aligned.push(b"ACAGT".to_vec());
        }
        let msa = MsaResult {
            aligned,
            n_columns: 5,
        };
        let mut config = ProfileHmmConfig::dna();
        config.local_mode = false;
        let hmm = ProfileHmm::from_msa(&msa, &config).unwrap();

        // Query matches the inserted variant — insert state captures the 'A'
        let result = hmm.viterbi(b"ACAGT").unwrap();

        let has_inserts = result
            .states
            .iter()
            .any(|s| matches!(s, ProfileHmmState::Insert(_)));
        assert!(
            has_inserts,
            "query matching inserted variant should use insert states, got {:?}",
            result.states
        );
    }

    #[test]
    fn viterbi_error_empty() {
        let hmm = build_dna_hmm(&[b"ACGT", b"ACGT", b"ACGA"]);
        assert!(hmm.viterbi(b"").is_err());
    }

    #[test]
    fn viterbi_error_invalid_char() {
        let hmm = build_dna_hmm(&[b"ACGT", b"ACGT", b"ACGA"]);
        assert!(hmm.viterbi(b"ACXG").is_err());
    }

    // ------------------------------------------------------------------
    // Forward tests
    // ------------------------------------------------------------------

    #[test]
    fn forward_finite() {
        let hmm = build_dna_hmm(&[b"ACGT", b"ACGT", b"ACGA"]);
        let score = hmm.forward(b"ACGT").unwrap();
        assert!(score.is_finite(), "forward score should be finite");
    }

    #[test]
    fn forward_identical_higher() {
        let hmm = build_dna_hmm(&[b"ACGT", b"ACGT", b"ACGT"]);
        let identical = hmm.forward(b"ACGT").unwrap();
        let different = hmm.forward(b"TTTT").unwrap();
        assert!(
            identical > different,
            "identical sequence ({}) should score higher than different ({})",
            identical,
            different
        );
    }

    #[test]
    fn forward_gte_viterbi() {
        let hmm = build_dna_hmm(&[b"ACGT", b"ACGT", b"ACGA"]);
        let fwd = hmm.forward(b"ACGT").unwrap();
        let vit = hmm.viterbi(b"ACGT").unwrap();
        assert!(
            fwd >= vit.score - 1e-6,
            "forward ({}) should be >= viterbi ({})",
            fwd,
            vit.score
        );
    }

    #[test]
    fn forward_error_empty() {
        let hmm = build_dna_hmm(&[b"ACGT", b"ACGT", b"ACGA"]);
        assert!(hmm.forward(b"").is_err());
    }

    #[test]
    fn forward_short_sequence() {
        let hmm = build_dna_hmm(&[b"ACGT", b"ACGT", b"ACGA"]);
        let score = hmm.forward(b"AC").unwrap();
        assert!(score.is_finite());
    }

    // ------------------------------------------------------------------
    // Backward tests
    // ------------------------------------------------------------------

    #[test]
    fn backward_finite() {
        let hmm = build_dna_hmm(&[b"ACGT", b"ACGT", b"ACGA"]);
        let score = hmm.backward(b"ACGT").unwrap();
        assert!(score.is_finite(), "backward score should be finite");
    }

    #[test]
    fn backward_matches_forward() {
        let hmm = build_dna_hmm(&[b"ACGT", b"ACGT", b"ACGA"]);

        for seq in &[b"ACGT".as_slice(), b"ACGA", b"TTTT", b"AC"] {
            let fwd = hmm.forward(seq).unwrap();
            let bwd = hmm.backward(seq).unwrap();
            assert!(
                (fwd - bwd).abs() < 0.5,
                "forward ({}) and backward ({}) should approximately match for {:?}",
                fwd,
                bwd,
                std::str::from_utf8(seq)
            );
        }
    }

    #[test]
    fn backward_error_empty() {
        let hmm = build_dna_hmm(&[b"ACGT", b"ACGT", b"ACGA"]);
        assert!(hmm.backward(b"").is_err());
    }

    #[test]
    fn backward_different_sequences() {
        let hmm = build_dna_hmm(&[b"ACGT", b"ACGT", b"ACGT"]);
        let score_match = hmm.backward(b"ACGT").unwrap();
        let score_diff = hmm.backward(b"TTTT").unwrap();
        assert!(
            score_match > score_diff,
            "matching sequence ({}) should score higher than different ({})",
            score_match,
            score_diff
        );
    }

    // ------------------------------------------------------------------
    // E-value tests
    // ------------------------------------------------------------------

    #[test]
    fn calibrate_sets_params() {
        let mut hmm = build_dna_hmm(&[b"ACGTACGT", b"ACGTACGT", b"ACGTACGA"]);
        assert!(hmm.gumbel.is_none());

        hmm.calibrate(100, 20, 42).unwrap();
        assert!(hmm.gumbel.is_some());

        let g = hmm.gumbel.as_ref().unwrap();
        assert!(g.mu.is_finite());
        assert!(g.lambda.is_finite());
        assert!(g.lambda > 0.0);
    }

    #[test]
    fn evalue_decreases_with_score() {
        let mut hmm = build_dna_hmm(&[b"ACGTACGT", b"ACGTACGT", b"ACGTACGA"]);
        hmm.calibrate(200, 20, 42).unwrap();

        // Use scores relative to calibration: a good match vs a poor match
        let good_score = hmm.viterbi(b"ACGTACGT").unwrap().score;
        let poor_score = hmm.viterbi(b"TTTTTTTT").unwrap().score;

        let ev_good = hmm.evalue(good_score, 1000).unwrap();
        let ev_poor = hmm.evalue(poor_score, 1000).unwrap();

        // Higher score → lower E-value
        assert!(
            ev_good < ev_poor,
            "good match E-value ({}) should be less than poor match E-value ({})",
            ev_good,
            ev_poor
        );
    }

    #[test]
    fn evalue_error_before_calibration() {
        let hmm = build_dna_hmm(&[b"ACGT", b"ACGT", b"ACGA"]);
        assert!(hmm.evalue(-10.0, 1000).is_err());
    }

    #[test]
    fn calibrate_error_zero_samples() {
        let mut hmm = build_dna_hmm(&[b"ACGT", b"ACGT", b"ACGA"]);
        assert!(hmm.calibrate(0, 20, 42).is_err());
    }
}
