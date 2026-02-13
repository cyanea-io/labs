//! X-drop and Z-drop seed extension alignment.
//!
//! Implements banded extension alignment with early termination for
//! efficient seed-and-extend workflows. X-drop terminates when all
//! active cells fall below `best - x_drop`. Z-drop tracks the best
//! score per anti-diagonal and is more permissive.

use cyanea_core::{CyaneaError, Result};

use crate::scoring::ScoringScheme;
use crate::types::{AlignmentResult, CigarOp};

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

/// Parameters for X-drop extension alignment.
#[derive(Debug, Clone)]
pub struct XDropParams {
    /// Drop threshold: terminate when all active cells score < best - x_drop (default: 50).
    pub x_drop: i32,
    /// Maximum bandwidth around the diagonal (default: 100).
    pub bandwidth: usize,
}

impl Default for XDropParams {
    fn default() -> Self {
        Self {
            x_drop: 50,
            bandwidth: 100,
        }
    }
}

/// Parameters for Z-drop extension alignment.
#[derive(Debug, Clone)]
pub struct ZDropParams {
    /// Drop threshold per anti-diagonal (default: 200).
    pub z_drop: i32,
    /// Maximum bandwidth (default: 500).
    pub bandwidth: usize,
}

impl Default for ZDropParams {
    fn default() -> Self {
        Self {
            z_drop: 200,
            bandwidth: 500,
        }
    }
}

/// Result of an extension alignment.
#[derive(Debug, Clone)]
pub struct ExtensionResult {
    /// Best alignment score achieved.
    pub score: i32,
    /// Number of query bases consumed.
    pub query_len: usize,
    /// Number of target bases consumed.
    pub target_len: usize,
    /// CIGAR operations describing the alignment.
    pub cigar: Vec<CigarOp>,
    /// Whether the extension was terminated by the drop condition.
    pub dropped: bool,
}

// ---------------------------------------------------------------------------
// CIGAR helper
// ---------------------------------------------------------------------------

fn push_cigar(cigar: &mut Vec<CigarOp>, op: CigarOp) {
    if let Some(last) = cigar.last_mut() {
        match (last, &op) {
            (CigarOp::Match(n), CigarOp::Match(m)) => {
                *n += m;
                return;
            }
            (CigarOp::Mismatch(n), CigarOp::Mismatch(m)) => {
                *n += m;
                return;
            }
            (CigarOp::Insertion(n), CigarOp::Insertion(m)) => {
                *n += m;
                return;
            }
            (CigarOp::Deletion(n), CigarOp::Deletion(m)) => {
                *n += m;
                return;
            }
            _ => {}
        }
    }
    cigar.push(op);
}

// ---------------------------------------------------------------------------
// Core DP engine
// ---------------------------------------------------------------------------

/// Drop mode: X-drop vs Z-drop
enum DropMode {
    XDrop(i32),
    ZDrop(i32),
}

const NEG_INF: i32 = i32::MIN / 2;

/// Core banded extension DP.
///
/// Extends from position (0,0) in the given query/target slices.
/// Returns the extension result with traceback.
fn extend_dp(
    query: &[u8],
    target: &[u8],
    scoring: &ScoringScheme,
    bandwidth: usize,
    drop: DropMode,
) -> ExtensionResult {
    let m = query.len();
    let n = target.len();

    if m == 0 || n == 0 {
        return ExtensionResult {
            score: 0,
            query_len: 0,
            target_len: 0,
            cigar: Vec::new(),
            dropped: false,
        };
    }

    let gap_open = scoring.gap_open();
    let gap_extend = scoring.gap_extend();
    let bw = bandwidth;

    // DP with affine gaps: H (match), E (gap in target/insertion), F (gap in query/deletion)
    // Use row-by-row storage, banded
    // Traceback: 0=diag(match/mismatch), 1=left(insertion from E), 2=up(deletion from F)
    // Store full matrix for traceback (bounded by bw)

    // For simplicity with bandwidth, we store full rows but only compute within band
    let cols = n + 1;
    let rows = m + 1;

    // Full DP (for correctness with traceback) — bounded by bandwidth
    let mut h = vec![vec![NEG_INF; cols]; rows];
    let mut e = vec![vec![NEG_INF; cols]; rows];
    let mut f = vec![vec![NEG_INF; cols]; rows];
    // Traceback: 0=diag, 1=ins(E), 2=del(F)
    let mut tb = vec![vec![0u8; cols]; rows];

    h[0][0] = 0;
    // First row: gaps in target (deletions)
    for j in 1..=n.min(bw) {
        f[0][j] = gap_open + gap_extend * j as i32;
        h[0][j] = f[0][j];
        tb[0][j] = 2;
    }
    // First column: gaps in query (insertions)
    for i in 1..=m.min(bw) {
        e[i][0] = gap_open + gap_extend * i as i32;
        h[i][0] = e[i][0];
        tb[i][0] = 1;
    }

    let mut best_score = 0i32;
    let mut best_i = 0usize;
    let mut best_j = 0usize;
    let mut dropped = false;

    // For z-drop: track best score on each anti-diagonal
    let max_ad = m + n;
    let mut best_per_ad = vec![NEG_INF; max_ad + 1];
    best_per_ad[0] = 0;

    for i in 1..=m {
        let j_center = i; // diagonal
        let j_lo = if j_center > bw { j_center - bw } else { 1 };
        let j_hi = (j_center + bw).min(n);

        let mut row_max = NEG_INF;

        for j in j_lo..=j_hi {
            // E: insertion (gap in target, consumes query)
            let e_open = h[i - 1][j].saturating_add(gap_open + gap_extend);
            let e_ext = e[i - 1][j].saturating_add(gap_extend);
            e[i][j] = e_open.max(e_ext);

            // F: deletion (gap in query, consumes target)
            let f_open = h[i][j - 1].saturating_add(gap_open + gap_extend);
            let f_ext = f[i][j - 1].saturating_add(gap_extend);
            f[i][j] = f_open.max(f_ext);

            // H: match/mismatch
            let sub = scoring.score_pair(query[i - 1], target[j - 1]);
            let diag = h[i - 1][j - 1].saturating_add(sub);

            let h_val = diag.max(e[i][j]).max(f[i][j]);
            h[i][j] = h_val;

            if h_val == diag {
                tb[i][j] = 0;
            } else if h_val == e[i][j] {
                tb[i][j] = 1;
            } else {
                tb[i][j] = 2;
            }

            row_max = row_max.max(h_val);

            if h_val > best_score {
                best_score = h_val;
                best_i = i;
                best_j = j;
            }

            let ad = i + j;
            if ad <= max_ad {
                best_per_ad[ad] = best_per_ad[ad].max(h_val);
            }
        }

        // Check drop condition
        match drop {
            DropMode::XDrop(x) => {
                if row_max < best_score - x {
                    dropped = true;
                    break;
                }
            }
            DropMode::ZDrop(z) => {
                let ad = i + (i.min(n)); // current anti-diagonal
                if ad <= max_ad && best_score - best_per_ad[ad] > z {
                    dropped = true;
                    break;
                }
            }
        }
    }

    // Traceback from (best_i, best_j)
    let mut cigar = Vec::new();
    let mut ci = best_i;
    let mut cj = best_j;

    while ci > 0 && cj > 0 {
        match tb[ci][cj] {
            0 => {
                // diagonal
                let op = if query[ci - 1].to_ascii_uppercase()
                    == target[cj - 1].to_ascii_uppercase()
                {
                    CigarOp::Match(1)
                } else {
                    CigarOp::Mismatch(1)
                };
                push_cigar(&mut cigar, op);
                ci -= 1;
                cj -= 1;
            }
            1 => {
                // insertion (gap in target)
                push_cigar(&mut cigar, CigarOp::Insertion(1));
                ci -= 1;
            }
            2 => {
                // deletion (gap in query)
                push_cigar(&mut cigar, CigarOp::Deletion(1));
                cj -= 1;
            }
            _ => break,
        }
    }

    // Handle remaining bases at start
    while ci > 0 {
        push_cigar(&mut cigar, CigarOp::Insertion(1));
        ci -= 1;
    }
    while cj > 0 {
        push_cigar(&mut cigar, CigarOp::Deletion(1));
        cj -= 1;
    }

    cigar.reverse();

    ExtensionResult {
        score: best_score,
        query_len: best_i,
        target_len: best_j,
        cigar,
        dropped,
    }
}

// ---------------------------------------------------------------------------
// X-drop public API
// ---------------------------------------------------------------------------

/// Extend alignment to the right using X-drop.
///
/// Aligns `query` against `target` starting from position 0, terminating
/// when all active cells fall below `best_score - x_drop`.
///
/// # Errors
///
/// Returns an error if either sequence is empty.
pub fn xdrop_extend_right(
    query: &[u8],
    target: &[u8],
    scoring: &ScoringScheme,
    params: &XDropParams,
) -> Result<ExtensionResult> {
    if query.is_empty() || target.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "query and target must be non-empty".into(),
        ));
    }
    Ok(extend_dp(
        query,
        target,
        scoring,
        params.bandwidth,
        DropMode::XDrop(params.x_drop),
    ))
}

/// Extend alignment to the left using X-drop.
///
/// Reverses both sequences, runs forward extension, then reverses the CIGAR.
///
/// # Errors
///
/// Returns an error if either sequence is empty.
pub fn xdrop_extend_left(
    query: &[u8],
    target: &[u8],
    scoring: &ScoringScheme,
    params: &XDropParams,
) -> Result<ExtensionResult> {
    if query.is_empty() || target.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "query and target must be non-empty".into(),
        ));
    }
    let rev_q: Vec<u8> = query.iter().rev().copied().collect();
    let rev_t: Vec<u8> = target.iter().rev().copied().collect();
    let mut result = extend_dp(
        &rev_q,
        &rev_t,
        scoring,
        params.bandwidth,
        DropMode::XDrop(params.x_drop),
    );
    result.cigar.reverse();
    Ok(result)
}

/// Seed-and-extend alignment using X-drop.
///
/// Given a seed match at `(query_pos, target_pos)` of length `seed_len`,
/// extends left and right, combines into a full alignment result.
///
/// # Errors
///
/// Returns an error if positions are out of bounds.
pub fn xdrop_seed_extend(
    query: &[u8],
    target: &[u8],
    query_pos: usize,
    target_pos: usize,
    seed_len: usize,
    scoring: &ScoringScheme,
    params: &XDropParams,
) -> Result<AlignmentResult> {
    seed_extend_impl(
        query,
        target,
        query_pos,
        target_pos,
        seed_len,
        scoring,
        params.bandwidth,
        DropMode::XDrop(params.x_drop),
    )
}

// ---------------------------------------------------------------------------
// Z-drop public API
// ---------------------------------------------------------------------------

/// Extend alignment to the right using Z-drop.
///
/// Z-drop tracks the best score per anti-diagonal and is more permissive
/// than X-drop, allowing brief score dips.
///
/// # Errors
///
/// Returns an error if either sequence is empty.
pub fn zdrop_extend_right(
    query: &[u8],
    target: &[u8],
    scoring: &ScoringScheme,
    params: &ZDropParams,
) -> Result<ExtensionResult> {
    if query.is_empty() || target.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "query and target must be non-empty".into(),
        ));
    }
    Ok(extend_dp(
        query,
        target,
        scoring,
        params.bandwidth,
        DropMode::ZDrop(params.z_drop),
    ))
}

/// Extend alignment to the left using Z-drop.
///
/// # Errors
///
/// Returns an error if either sequence is empty.
pub fn zdrop_extend_left(
    query: &[u8],
    target: &[u8],
    scoring: &ScoringScheme,
    params: &ZDropParams,
) -> Result<ExtensionResult> {
    if query.is_empty() || target.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "query and target must be non-empty".into(),
        ));
    }
    let rev_q: Vec<u8> = query.iter().rev().copied().collect();
    let rev_t: Vec<u8> = target.iter().rev().copied().collect();
    let mut result = extend_dp(
        &rev_q,
        &rev_t,
        scoring,
        params.bandwidth,
        DropMode::ZDrop(params.z_drop),
    );
    result.cigar.reverse();
    Ok(result)
}

/// Seed-and-extend alignment using Z-drop.
///
/// # Errors
///
/// Returns an error if positions are out of bounds.
pub fn zdrop_seed_extend(
    query: &[u8],
    target: &[u8],
    query_pos: usize,
    target_pos: usize,
    seed_len: usize,
    scoring: &ScoringScheme,
    params: &ZDropParams,
) -> Result<AlignmentResult> {
    seed_extend_impl(
        query,
        target,
        query_pos,
        target_pos,
        seed_len,
        scoring,
        params.bandwidth,
        DropMode::ZDrop(params.z_drop),
    )
}

// ---------------------------------------------------------------------------
// Shared seed-extend
// ---------------------------------------------------------------------------

fn seed_extend_impl(
    query: &[u8],
    target: &[u8],
    query_pos: usize,
    target_pos: usize,
    seed_len: usize,
    scoring: &ScoringScheme,
    bandwidth: usize,
    drop: DropMode,
) -> Result<AlignmentResult> {
    if query_pos + seed_len > query.len() || target_pos + seed_len > target.len() {
        return Err(CyaneaError::InvalidInput(
            "seed position + length exceeds sequence bounds".into(),
        ));
    }

    // Score the seed
    let mut seed_score = 0i32;
    let mut seed_cigar = Vec::new();
    for k in 0..seed_len {
        let s = scoring.score_pair(query[query_pos + k], target[target_pos + k]);
        seed_score += s;
        let op = if query[query_pos + k].to_ascii_uppercase()
            == target[target_pos + k].to_ascii_uppercase()
        {
            CigarOp::Match(1)
        } else {
            CigarOp::Mismatch(1)
        };
        push_cigar(&mut seed_cigar, op);
    }

    // Extend left
    let left = if query_pos > 0 && target_pos > 0 {
        let lq: Vec<u8> = query[..query_pos].iter().rev().copied().collect();
        let lt: Vec<u8> = target[..target_pos].iter().rev().copied().collect();
        let drop_left = match &drop {
            DropMode::XDrop(x) => DropMode::XDrop(*x),
            DropMode::ZDrop(z) => DropMode::ZDrop(*z),
        };
        let mut r = extend_dp(&lq, &lt, scoring, bandwidth, drop_left);
        r.cigar.reverse();
        r
    } else {
        ExtensionResult {
            score: 0,
            query_len: 0,
            target_len: 0,
            cigar: Vec::new(),
            dropped: false,
        }
    };

    // Extend right
    let right_q_start = query_pos + seed_len;
    let right_t_start = target_pos + seed_len;
    let right = if right_q_start < query.len() && right_t_start < target.len() {
        let drop_right = match &drop {
            DropMode::XDrop(x) => DropMode::XDrop(*x),
            DropMode::ZDrop(z) => DropMode::ZDrop(*z),
        };
        extend_dp(
            &query[right_q_start..],
            &target[right_t_start..],
            scoring,
            bandwidth,
            drop_right,
        )
    } else {
        ExtensionResult {
            score: 0,
            query_len: 0,
            target_len: 0,
            cigar: Vec::new(),
            dropped: false,
        }
    };

    // Combine
    let total_score = left.score + seed_score + right.score;
    let q_start = query_pos - left.query_len;
    let q_end = query_pos + seed_len + right.query_len;
    let t_start = target_pos - left.target_len;
    let t_end = target_pos + seed_len + right.target_len;

    let mut combined_cigar = left.cigar;
    for op in seed_cigar {
        push_cigar(&mut combined_cigar, op);
    }
    for op in right.cigar {
        push_cigar(&mut combined_cigar, op);
    }

    // Build aligned sequences
    let mut aligned_query = Vec::new();
    let mut aligned_target = Vec::new();
    let mut qi = q_start;
    let mut ti = t_start;

    for op in &combined_cigar {
        match op {
            CigarOp::Match(n) | CigarOp::Mismatch(n) => {
                for _ in 0..*n {
                    aligned_query.push(query[qi]);
                    aligned_target.push(target[ti]);
                    qi += 1;
                    ti += 1;
                }
            }
            CigarOp::Insertion(n) => {
                for _ in 0..*n {
                    aligned_query.push(query[qi]);
                    aligned_target.push(b'-');
                    qi += 1;
                }
            }
            CigarOp::Deletion(n) => {
                for _ in 0..*n {
                    aligned_query.push(b'-');
                    aligned_target.push(target[ti]);
                    ti += 1;
                }
            }
            _ => {}
        }
    }

    Ok(AlignmentResult {
        score: total_score,
        aligned_query,
        aligned_target,
        query_start: q_start,
        query_end: q_end,
        target_start: t_start,
        target_end: t_end,
        cigar: combined_cigar,
    })
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::ScoringMatrix;

    fn dna_scheme() -> ScoringScheme {
        ScoringScheme::Simple(ScoringMatrix::dna_default())
    }

    #[test]
    fn identical_sequences() {
        let result =
            xdrop_extend_right(b"ACGTACGT", b"ACGTACGT", &dna_scheme(), &XDropParams::default())
                .unwrap();
        assert_eq!(result.score, 16); // 8 matches * 2
        assert_eq!(result.query_len, 8);
        assert_eq!(result.target_len, 8);
        assert!(!result.dropped);
    }

    #[test]
    fn early_termination() {
        // Good match followed by many mismatches
        let query = b"ACGTACGTTTTTTTTTTTTTTTTTTT";
        let target = b"ACGTACGTAAAAAAAAAAAAAAAAAAA";
        let params = XDropParams {
            x_drop: 10,
            bandwidth: 100,
        };
        let result = xdrop_extend_right(query, target, &dna_scheme(), &params).unwrap();
        // Should terminate before reaching the end
        assert!(result.dropped || result.query_len < query.len());
    }

    #[test]
    fn matches_sw_for_short() {
        // For short identical sequences, should match SW
        let result =
            xdrop_extend_right(b"ACGT", b"ACGT", &dna_scheme(), &XDropParams::default()).unwrap();
        assert_eq!(result.score, 8); // 4 * 2
    }

    #[test]
    fn empty_error() {
        let result = xdrop_extend_right(b"", b"ACGT", &dna_scheme(), &XDropParams::default());
        assert!(result.is_err());
        let result = xdrop_extend_right(b"ACGT", b"", &dna_scheme(), &XDropParams::default());
        assert!(result.is_err());
    }

    #[test]
    fn single_base() {
        let result =
            xdrop_extend_right(b"A", b"A", &dna_scheme(), &XDropParams::default()).unwrap();
        assert_eq!(result.score, 2);
        assert_eq!(result.query_len, 1);
    }

    #[test]
    fn left_extension_reversal() {
        let result =
            xdrop_extend_left(b"ACGTACGT", b"ACGTACGT", &dna_scheme(), &XDropParams::default())
                .unwrap();
        assert_eq!(result.score, 16);
        assert_eq!(result.query_len, 8);
    }

    #[test]
    fn seed_extend_combine() {
        let query = b"AAAACGTACGTAAAA";
        let target = b"AAAACGTACGTAAAA";
        let result = xdrop_seed_extend(
            query,
            target,
            4,
            4,
            7,
            &dna_scheme(),
            &XDropParams::default(),
        )
        .unwrap();
        assert_eq!(result.score, 30); // 15 * 2
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 15);
    }

    #[test]
    fn seed_extend_out_of_bounds() {
        let result =
            xdrop_seed_extend(b"ACGT", b"ACGT", 3, 3, 5, &dna_scheme(), &XDropParams::default());
        assert!(result.is_err());
    }

    #[test]
    fn zdrop_extend() {
        let result =
            zdrop_extend_right(b"ACGTACGT", b"ACGTACGT", &dna_scheme(), &ZDropParams::default())
                .unwrap();
        assert_eq!(result.score, 16);
        assert!(!result.dropped);
    }

    #[test]
    fn zdrop_more_permissive() {
        // Z-drop should extend further than X-drop for same threshold
        let query = b"ACGTACGTTTACGTACGT";
        let target = b"ACGTACGTAAACGTACGT";
        let x_params = XDropParams {
            x_drop: 5,
            bandwidth: 100,
        };
        let z_params = ZDropParams {
            z_drop: 5,
            bandwidth: 100,
        };
        let x_result = xdrop_extend_right(query, target, &dna_scheme(), &x_params).unwrap();
        let z_result = zdrop_extend_right(query, target, &dna_scheme(), &z_params).unwrap();
        // Z-drop should generally score at least as well
        assert!(
            z_result.score >= x_result.score,
            "z_drop score {} should be >= x_drop score {}",
            z_result.score,
            x_result.score
        );
    }

    #[test]
    fn zdrop_seed_extend_combine() {
        let query = b"AAAACGTACGTAAAA";
        let target = b"AAAACGTACGTAAAA";
        let result = super::zdrop_seed_extend(
            query,
            target,
            4,
            4,
            7,
            &dna_scheme(),
            &ZDropParams::default(),
        )
        .unwrap();
        assert_eq!(result.score, 30);
    }

    #[test]
    fn dropped_flag_set() {
        let query = b"ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        let target = b"ACGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let params = XDropParams {
            x_drop: 5,
            bandwidth: 100,
        };
        let result = xdrop_extend_right(query, target, &dna_scheme(), &params).unwrap();
        assert!(result.dropped, "should have dropped");
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use crate::scoring::ScoringMatrix;
    use proptest::prelude::*;

    fn dna_scheme() -> ScoringScheme {
        ScoringScheme::Simple(ScoringMatrix::dna_default())
    }

    fn dna_seq(max_len: usize) -> impl Strategy<Value = Vec<u8>> {
        proptest::collection::vec(
            prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')],
            1..=max_len,
        )
    }

    proptest! {
        #[test]
        fn xdrop_deterministic(q in dna_seq(30), t in dna_seq(30)) {
            let r1 = xdrop_extend_right(&q, &t, &dna_scheme(), &XDropParams::default()).unwrap();
            let r2 = xdrop_extend_right(&q, &t, &dna_scheme(), &XDropParams::default()).unwrap();
            prop_assert_eq!(r1.score, r2.score);
        }
    }
}
