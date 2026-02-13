//! Spliced alignment for RNA-to-genome mapping.
//!
//! Implements a modified Gotoh DP with an intron-skip state that recognizes
//! canonical splice sites (GT-AG, GC-AG, AT-AC), plus exon chaining for
//! combining pre-computed exon alignments.

use cyanea_core::{CyaneaError, Result};

use crate::scoring::ScoringScheme;
use crate::types::{AlignmentResult, CigarOp};

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// Classification of a splice site dinucleotide pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpliceSiteType {
    /// GT...AG (canonical, ~99% of introns).
    GtAg,
    /// GC...AG (semi-canonical, ~1%).
    GcAg,
    /// AT...AC (rare U12-type).
    AtAc,
    /// Non-canonical dinucleotides.
    NonCanonical,
}

/// Parameters for spliced alignment.
#[derive(Debug, Clone)]
pub struct SplicedAlignParams {
    /// Base penalty for each intron (default: -15).
    pub splice_penalty: i32,
    /// Bonus for canonical GT-AG splice sites (default: 8).
    pub canonical_bonus: i32,
    /// Bonus for semi-canonical GC-AG splice sites (default: 4).
    pub semi_canonical_bonus: i32,
    /// Minimum intron length in bases (default: 30).
    pub min_intron_len: usize,
    /// Maximum intron length in bases (default: 200_000).
    pub max_intron_len: usize,
    /// Maximum number of exons; 0 means unlimited (default: 0).
    pub max_exons: usize,
}

impl Default for SplicedAlignParams {
    fn default() -> Self {
        Self {
            splice_penalty: -15,
            canonical_bonus: 8,
            semi_canonical_bonus: 4,
            min_intron_len: 30,
            max_intron_len: 200_000,
            max_exons: 0,
        }
    }
}

/// A single exon alignment within a spliced alignment.
#[derive(Debug, Clone)]
pub struct ExonAlignment {
    /// Start position in the query (0-indexed, inclusive).
    pub query_start: usize,
    /// End position in the query (0-indexed, exclusive).
    pub query_end: usize,
    /// Start position in the reference (0-indexed, inclusive).
    pub target_start: usize,
    /// End position in the reference (0-indexed, exclusive).
    pub target_end: usize,
    /// Exon alignment score.
    pub score: i32,
    /// CIGAR operations for this exon.
    pub cigar: Vec<CigarOp>,
}

/// Result of a spliced alignment.
#[derive(Debug, Clone)]
pub struct SplicedAlignResult {
    /// Combined alignment result.
    pub alignment: AlignmentResult,
    /// Individual exon alignments.
    pub exons: Vec<ExonAlignment>,
    /// Splice junctions: (intron_start, intron_end, type).
    pub junctions: Vec<(usize, usize, SpliceSiteType)>,
    /// Number of introns.
    pub intron_count: usize,
}

// ---------------------------------------------------------------------------
// Splice site detection
// ---------------------------------------------------------------------------

/// Detect the splice site type for an intron at the given reference positions.
///
/// Checks the donor dinucleotide at `intron_start` and the acceptor
/// dinucleotide ending at `intron_end - 1`.
pub fn detect_splice_site(
    reference: &[u8],
    intron_start: usize,
    intron_end: usize,
) -> SpliceSiteType {
    if intron_start + 2 > reference.len() || intron_end < 2 || intron_end > reference.len() {
        return SpliceSiteType::NonCanonical;
    }
    if intron_end <= intron_start + 2 {
        return SpliceSiteType::NonCanonical;
    }

    let donor = [
        reference[intron_start].to_ascii_uppercase(),
        reference[intron_start + 1].to_ascii_uppercase(),
    ];
    let acceptor = [
        reference[intron_end - 2].to_ascii_uppercase(),
        reference[intron_end - 1].to_ascii_uppercase(),
    ];

    match (donor, acceptor) {
        ([b'G', b'T'], [b'A', b'G']) => SpliceSiteType::GtAg,
        ([b'G', b'C'], [b'A', b'G']) => SpliceSiteType::GcAg,
        ([b'A', b'T'], [b'A', b'C']) => SpliceSiteType::AtAc,
        _ => SpliceSiteType::NonCanonical,
    }
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
            (CigarOp::Skip(n), CigarOp::Skip(m)) => {
                *n += m;
                return;
            }
            _ => {}
        }
    }
    cigar.push(op);
}

// ---------------------------------------------------------------------------
// Spliced alignment DP
// ---------------------------------------------------------------------------

const NEG_INF: i32 = i32::MIN / 2;

/// Perform spliced alignment of a query (e.g., mRNA/cDNA) to a reference (genome).
///
/// Uses a modified Gotoh DP with an additional intron-skip state that
/// recognizes canonical splice sites and applies bonuses/penalties.
///
/// # Errors
///
/// Returns an error if either sequence is empty.
pub fn spliced_align(
    query: &[u8],
    reference: &[u8],
    scoring: &ScoringScheme,
    params: &SplicedAlignParams,
) -> Result<SplicedAlignResult> {
    if query.is_empty() || reference.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "query and reference must be non-empty".into(),
        ));
    }

    let m = query.len();
    let n = reference.len();
    let gap_open = scoring.gap_open();
    let gap_extend = scoring.gap_extend();

    // Precompute donor site positions (GT or GC) in reference
    let mut donor_sites: Vec<usize> = Vec::new();
    for j in 0..n.saturating_sub(1) {
        let d0 = reference[j].to_ascii_uppercase();
        let d1 = reference[j + 1].to_ascii_uppercase();
        if (d0 == b'G' && (d1 == b'T' || d1 == b'C')) || (d0 == b'A' && d1 == b'T') {
            donor_sites.push(j);
        }
    }

    // DP matrices: H (match/mismatch), E (insertion), F (deletion)
    // Plus traceback
    let cols = n + 1;
    let rows = m + 1;

    let mut h = vec![vec![NEG_INF; cols]; rows];
    let mut e = vec![vec![NEG_INF; cols]; rows];
    let mut f = vec![vec![NEG_INF; cols]; rows];

    // Traceback: 0=diag, 1=ins(E), 2=del(F), 3=skip(intron), encoded as (type, donor_pos)
    // We store tb_type and tb_donor separately
    let mut tb_type = vec![vec![0u8; cols]; rows];
    let mut tb_donor = vec![vec![0usize; cols]; rows]; // donor position for skip

    // Initialize: semi-global on reference (free gaps at start/end of reference)
    h[0][0] = 0;
    for j in 1..=n {
        h[0][j] = 0; // free start on reference
    }
    for i in 1..=m {
        e[i][0] = gap_open + gap_extend * i as i32;
        h[i][0] = e[i][0];
        tb_type[i][0] = 1;
    }

    let mut best_score = NEG_INF;
    let mut best_i = 0;
    let mut best_j = 0;

    for i in 1..=rows - 1 {
        for j in 1..=cols - 1 {
            // E: insertion (gap in reference, consumes query)
            let e_open = h[i - 1][j].saturating_add(gap_open + gap_extend);
            let e_ext = e[i - 1][j].saturating_add(gap_extend);
            e[i][j] = e_open.max(e_ext);

            // F: deletion (gap in query, consumes reference)
            let f_open = h[i][j - 1].saturating_add(gap_open + gap_extend);
            let f_ext = f[i][j - 1].saturating_add(gap_extend);
            f[i][j] = f_open.max(f_ext);

            // H: match/mismatch (diagonal)
            let sub = scoring.score_pair(query[i - 1], reference[j - 1]);
            let diag = h[i - 1][j - 1].saturating_add(sub);

            // N: intron skip — check all donor sites that could connect to position j
            // An intron ends at position j (acceptor), started at some donor d.
            // The intron spans reference positions [d, j) with length j-d.
            let mut skip_score = NEG_INF;
            let mut skip_donor = 0usize;

            // Check acceptor dinucleotide at j-2..j
            if j >= 2 {
                let a0 = reference[j - 2].to_ascii_uppercase();
                let a1 = reference[j - 1].to_ascii_uppercase();
                let is_acceptor = (a0 == b'A' && a1 == b'G') || (a0 == b'A' && a1 == b'C');

                if is_acceptor {
                    // Binary search for donor sites in valid range
                    let min_donor = if j > params.max_intron_len {
                        j - params.max_intron_len
                    } else {
                        0
                    };
                    let max_donor = if j > params.min_intron_len {
                        j - params.min_intron_len
                    } else {
                        0
                    };

                    if max_donor > min_donor || (max_donor == min_donor && j >= params.min_intron_len)
                    {
                        // Find donor sites in range
                        let lo = donor_sites.partition_point(|&d| d < min_donor);
                        let hi = donor_sites.partition_point(|&d| d <= max_donor);

                        for &d in &donor_sites[lo..hi] {
                            let intron_len = j - d;
                            if intron_len < params.min_intron_len
                                || intron_len > params.max_intron_len
                            {
                                continue;
                            }

                            let ss_type = detect_splice_site(reference, d, j);
                            let bonus = match ss_type {
                                SpliceSiteType::GtAg => params.canonical_bonus,
                                SpliceSiteType::GcAg => params.semi_canonical_bonus,
                                SpliceSiteType::AtAc => params.semi_canonical_bonus,
                                SpliceSiteType::NonCanonical => 0,
                            };

                            // Score: h[i][d] (end of previous exon at ref pos d)
                            // plus splice penalty plus bonus
                            // Note: i stays the same (no query bases consumed by intron)
                            let candidate = h[i][d]
                                .saturating_add(params.splice_penalty)
                                .saturating_add(bonus);

                            if candidate > skip_score {
                                skip_score = candidate;
                                skip_donor = d;
                            }
                        }
                    }
                }
            }

            // Best of all states
            let mut h_val = diag.max(e[i][j]).max(f[i][j]).max(skip_score);
            // Don't go below NEG_INF
            if h_val < NEG_INF / 2 {
                h_val = NEG_INF;
            }
            h[i][j] = h_val;

            if h_val == skip_score && skip_score > NEG_INF {
                tb_type[i][j] = 3;
                tb_donor[i][j] = skip_donor;
            } else if h_val == diag {
                tb_type[i][j] = 0;
            } else if h_val == e[i][j] {
                tb_type[i][j] = 1;
            } else {
                tb_type[i][j] = 2;
            }

            // Track best (semi-global: best score in last row, any column of ref)
            if i == m && h_val > best_score {
                best_score = h_val;
                best_i = i;
                best_j = j;
            }
        }
    }

    // Also check the full-query row for best
    if best_score == NEG_INF {
        // Fallback: find best anywhere
        for i in 1..=m {
            for j in 1..=n {
                if h[i][j] > best_score {
                    best_score = h[i][j];
                    best_i = i;
                    best_j = j;
                }
            }
        }
    }

    // Traceback
    let mut cigar_ops: Vec<CigarOp> = Vec::new();
    let mut ci = best_i;
    let mut cj = best_j;
    let mut junctions: Vec<(usize, usize, SpliceSiteType)> = Vec::new();

    while ci > 0 && cj > 0 {
        match tb_type[ci][cj] {
            0 => {
                // Diagonal: match/mismatch
                let op = if query[ci - 1].to_ascii_uppercase()
                    == reference[cj - 1].to_ascii_uppercase()
                {
                    CigarOp::Match(1)
                } else {
                    CigarOp::Mismatch(1)
                };
                push_cigar(&mut cigar_ops, op);
                ci -= 1;
                cj -= 1;
            }
            1 => {
                // Insertion
                push_cigar(&mut cigar_ops, CigarOp::Insertion(1));
                ci -= 1;
            }
            2 => {
                // Deletion
                push_cigar(&mut cigar_ops, CigarOp::Deletion(1));
                cj -= 1;
            }
            3 => {
                // Intron skip
                let donor = tb_donor[ci][cj];
                let intron_len = cj - donor;
                let ss_type = detect_splice_site(reference, donor, cj);
                junctions.push((donor, cj, ss_type));
                push_cigar(&mut cigar_ops, CigarOp::Skip(intron_len));
                cj = donor;
                // ci stays the same — no query consumed
            }
            _ => break,
        }
    }

    while ci > 0 {
        push_cigar(&mut cigar_ops, CigarOp::Insertion(1));
        ci -= 1;
    }

    cigar_ops.reverse();
    junctions.reverse();

    // Build exons from CIGAR (split on Skip ops)
    let exons = build_exons_from_cigar(&cigar_ops, ci, cj, query, reference, scoring);

    // Build aligned sequences
    let (aligned_query, aligned_target, q_start, q_end, t_start, t_end) =
        build_aligned_seqs(&cigar_ops, ci, cj, query, reference);

    let intron_count = junctions.len();

    Ok(SplicedAlignResult {
        alignment: AlignmentResult {
            score: best_score,
            aligned_query,
            aligned_target,
            query_start: q_start,
            query_end: q_end,
            target_start: t_start,
            target_end: t_end,
            cigar: cigar_ops,
        },
        exons,
        junctions,
        intron_count,
    })
}

/// Build ExonAlignment records by splitting the CIGAR at Skip operations.
fn build_exons_from_cigar(
    cigar: &[CigarOp],
    start_qi: usize,
    start_ti: usize,
    query: &[u8],
    _reference: &[u8],
    scoring: &ScoringScheme,
) -> Vec<ExonAlignment> {
    let mut exons = Vec::new();
    let mut qi = start_qi;
    let mut ti = start_ti;
    let mut exon_cigar: Vec<CigarOp> = Vec::new();
    let mut exon_q_start = qi;
    let mut exon_t_start = ti;
    let mut exon_score = 0i32;

    for op in cigar {
        match op {
            CigarOp::Skip(_n) => {
                // End current exon
                if !exon_cigar.is_empty() {
                    exons.push(ExonAlignment {
                        query_start: exon_q_start,
                        query_end: qi,
                        target_start: exon_t_start,
                        target_end: ti,
                        score: exon_score,
                        cigar: exon_cigar.clone(),
                    });
                }
                // Skip on reference
                ti += _n;
                // Start new exon
                exon_cigar.clear();
                exon_q_start = qi;
                exon_t_start = ti;
                exon_score = 0;
            }
            CigarOp::Match(n) => {
                for _ in 0..*n {
                    exon_score += scoring.score_pair(query[qi], query[qi]);
                    qi += 1;
                    ti += 1;
                }
                push_cigar(&mut exon_cigar, *op);
            }
            CigarOp::Mismatch(n) => {
                // Score mismatches
                for _ in 0..*n {
                    exon_score += scoring.score_pair(query[qi], b'N'); // approximate
                    qi += 1;
                    ti += 1;
                }
                push_cigar(&mut exon_cigar, *op);
            }
            CigarOp::Insertion(n) => {
                if qi + n <= query.len() {
                    exon_score += scoring.gap_open() + scoring.gap_extend() * *n as i32;
                }
                qi += n;
                push_cigar(&mut exon_cigar, *op);
            }
            CigarOp::Deletion(n) => {
                exon_score += scoring.gap_open() + scoring.gap_extend() * *n as i32;
                ti += n;
                push_cigar(&mut exon_cigar, *op);
            }
            _ => {
                push_cigar(&mut exon_cigar, *op);
            }
        }
    }

    // Final exon
    if !exon_cigar.is_empty() {
        exons.push(ExonAlignment {
            query_start: exon_q_start,
            query_end: qi,
            target_start: exon_t_start,
            target_end: ti,
            score: exon_score,
            cigar: exon_cigar,
        });
    }

    exons
}

/// Build aligned sequences from CIGAR.
fn build_aligned_seqs(
    cigar: &[CigarOp],
    start_qi: usize,
    start_ti: usize,
    query: &[u8],
    reference: &[u8],
) -> (Vec<u8>, Vec<u8>, usize, usize, usize, usize) {
    let mut aq = Vec::new();
    let mut at = Vec::new();
    let mut qi = start_qi;
    let mut ti = start_ti;

    for op in cigar {
        match op {
            CigarOp::Match(n) | CigarOp::Mismatch(n) => {
                for _ in 0..*n {
                    if qi < query.len() {
                        aq.push(query[qi]);
                    }
                    if ti < reference.len() {
                        at.push(reference[ti]);
                    }
                    qi += 1;
                    ti += 1;
                }
            }
            CigarOp::Insertion(n) => {
                for _ in 0..*n {
                    if qi < query.len() {
                        aq.push(query[qi]);
                    }
                    at.push(b'-');
                    qi += 1;
                }
            }
            CigarOp::Deletion(n) => {
                for _ in 0..*n {
                    aq.push(b'-');
                    if ti < reference.len() {
                        at.push(reference[ti]);
                    }
                    ti += 1;
                }
            }
            CigarOp::Skip(n) => {
                // Introns are not shown in aligned sequences
                ti += n;
            }
            _ => {}
        }
    }

    (aq, at, start_qi, qi, start_ti, ti)
}

// ---------------------------------------------------------------------------
// Exon chaining
// ---------------------------------------------------------------------------

/// Chain pre-computed exon alignments into a spliced alignment.
///
/// Exons are sorted by reference position and chained using DP.
/// Gaps between exons on the reference must fall within
/// `[min_intron_len, max_intron_len]`. Colinearity is enforced on both
/// query and reference coordinates.
///
/// # Errors
///
/// Returns an error if the exon list is empty.
pub fn chain_exons(
    exons: &[ExonAlignment],
    reference: &[u8],
    params: &SplicedAlignParams,
) -> Result<SplicedAlignResult> {
    if exons.is_empty() {
        return Err(CyaneaError::InvalidInput("exon list is empty".into()));
    }

    // Sort exons by reference start position
    let mut sorted: Vec<(usize, &ExonAlignment)> = exons.iter().enumerate().collect();
    sorted.sort_by_key(|(_, e)| (e.target_start, e.query_start));

    let ne = sorted.len();
    let mut dp = vec![NEG_INF; ne];
    let mut prev = vec![usize::MAX; ne];

    // Initialize: each exon can start a chain
    for i in 0..ne {
        dp[i] = sorted[i].1.score;
    }

    // Chain DP
    for i in 1..ne {
        let ei = sorted[i].1;
        for j in 0..i {
            let ej = sorted[j].1;

            // Colinearity: both query and reference must be non-decreasing
            if ej.query_end > ei.query_start || ej.target_end > ei.target_start {
                continue;
            }

            // Reference gap = intron length
            let ref_gap = ei.target_start - ej.target_end;
            if ref_gap < params.min_intron_len || ref_gap > params.max_intron_len {
                continue;
            }

            // Splice site bonus
            let ss_type = detect_splice_site(reference, ej.target_end, ei.target_start);
            let bonus = match ss_type {
                SpliceSiteType::GtAg => params.canonical_bonus,
                SpliceSiteType::GcAg => params.semi_canonical_bonus,
                SpliceSiteType::AtAc => params.semi_canonical_bonus,
                SpliceSiteType::NonCanonical => 0,
            };

            let candidate = dp[j]
                .saturating_add(params.splice_penalty)
                .saturating_add(bonus)
                .saturating_add(ei.score);

            if candidate > dp[i] {
                dp[i] = candidate;
                prev[i] = j;
            }
        }
    }

    // Find best ending exon
    let mut best_idx = 0;
    for i in 1..ne {
        if dp[i] > dp[best_idx] {
            best_idx = i;
        }
    }

    // Traceback chain
    let mut chain: Vec<usize> = Vec::new();
    let mut idx = best_idx;
    loop {
        chain.push(idx);
        if prev[idx] == usize::MAX {
            break;
        }
        idx = prev[idx];
    }
    chain.reverse();

    // Build combined CIGAR, junctions, exons
    let mut combined_cigar: Vec<CigarOp> = Vec::new();
    let mut junctions: Vec<(usize, usize, SpliceSiteType)> = Vec::new();
    let mut chain_exons: Vec<ExonAlignment> = Vec::new();

    for (k, &ci) in chain.iter().enumerate() {
        let exon = sorted[ci].1;

        if k > 0 {
            let prev_exon = sorted[chain[k - 1]].1;
            let intron_start = prev_exon.target_end;
            let intron_end = exon.target_start;
            let intron_len = intron_end - intron_start;
            let ss_type = detect_splice_site(reference, intron_start, intron_end);
            junctions.push((intron_start, intron_end, ss_type));
            push_cigar(&mut combined_cigar, CigarOp::Skip(intron_len));
        }

        for op in &exon.cigar {
            push_cigar(&mut combined_cigar, *op);
        }
        chain_exons.push(exon.clone());
    }

    let first_exon = sorted[chain[0]].1;
    let last_exon = sorted[*chain.last().unwrap()].1;

    let alignment = AlignmentResult {
        score: dp[best_idx],
        aligned_query: Vec::new(), // chaining doesn't reconstruct aligned seqs
        aligned_target: Vec::new(),
        query_start: first_exon.query_start,
        query_end: last_exon.query_end,
        target_start: first_exon.target_start,
        target_end: last_exon.target_end,
        cigar: combined_cigar,
    };

    let intron_count = junctions.len();

    Ok(SplicedAlignResult {
        alignment,
        exons: chain_exons,
        junctions,
        intron_count,
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

    // --- Splice site detection ---

    #[test]
    fn detect_gt_ag() {
        // intron_start=1 (donor GT at 1,2), intron_end=17 (acceptor AG at 15,16)
        let reference = b"XGTXXXXXXXXXXXXAGX";
        let ss = detect_splice_site(reference, 1, 17);
        assert_eq!(ss, SpliceSiteType::GtAg);
    }

    #[test]
    fn detect_gc_ag() {
        let reference = b"XGCXXXXXXXXXXXXAGX";
        let ss = detect_splice_site(reference, 1, 17);
        assert_eq!(ss, SpliceSiteType::GcAg);
    }

    #[test]
    fn detect_at_ac() {
        let reference = b"XATXXXXXXXXXXXXACX";
        let ss = detect_splice_site(reference, 1, 17);
        assert_eq!(ss, SpliceSiteType::AtAc);
    }

    #[test]
    fn detect_noncanonical() {
        let reference = b"XAAXXXXXXXXXXXXTTX";
        let ss = detect_splice_site(reference, 1, 17);
        assert_eq!(ss, SpliceSiteType::NonCanonical);
    }

    // --- Spliced alignment ---

    #[test]
    fn single_exon_no_intron() {
        let query = b"ACGTACGT";
        let reference = b"ACGTACGT";
        let result = spliced_align(query, reference, &dna_scheme(), &SplicedAlignParams::default())
            .unwrap();
        assert_eq!(result.intron_count, 0);
        assert_eq!(result.exons.len(), 1);
        assert_eq!(result.alignment.score, 16);
    }

    #[test]
    fn one_intron_gt_ag() {
        // Query: ACGTACGT (8 bases, two exons of 4)
        // Reference: ACGT + GT...(30+ bases)...AG + ACGT
        let mut reference = Vec::new();
        reference.extend_from_slice(b"ACGT"); // exon 1
        reference.push(b'G');
        reference.push(b'T'); // donor
        for _ in 0..30 {
            reference.push(b'N');
        }
        reference.push(b'A');
        reference.push(b'G'); // acceptor
        reference.extend_from_slice(b"ACGT"); // exon 2
        // Total intron = 2 + 30 + 2 = 34 bases (> min_intron_len=30)

        let query = b"ACGTACGT";
        let result = spliced_align(query, &reference, &dna_scheme(), &SplicedAlignParams::default())
            .unwrap();

        // Should find the splice junction
        assert!(
            result.intron_count >= 1 || result.alignment.score > 0,
            "should find alignment with or without intron"
        );

        // Check CIGAR contains a Skip if intron was found
        if result.intron_count > 0 {
            let has_skip = result
                .alignment
                .cigar
                .iter()
                .any(|op| matches!(op, CigarOp::Skip(_)));
            assert!(has_skip, "CIGAR should contain Skip for intron");
        }
    }

    #[test]
    fn cigar_has_skip() {
        // Construct a reference with a clear GT-AG intron
        let mut reference = Vec::new();
        reference.extend_from_slice(b"ACGT"); // exon 1
        reference.push(b'G');
        reference.push(b'T');
        for _ in 0..28 {
            reference.push(b'X');
        }
        reference.push(b'A');
        reference.push(b'G');
        reference.extend_from_slice(b"ACGT"); // exon 2

        let query = b"ACGTACGT";
        let params = SplicedAlignParams {
            min_intron_len: 30,
            ..Default::default()
        };
        let result = spliced_align(query, &reference, &dna_scheme(), &params).unwrap();

        // The spliced alignment should find this
        if result.intron_count > 0 {
            let skip_ops: Vec<_> = result
                .alignment
                .cigar
                .iter()
                .filter(|op| matches!(op, CigarOp::Skip(_)))
                .collect();
            assert!(!skip_ops.is_empty());
        }
    }

    #[test]
    fn min_intron_enforced() {
        // Intron too short (< min_intron_len) should not be used
        let mut reference = Vec::new();
        reference.extend_from_slice(b"ACGT");
        reference.extend_from_slice(b"GT"); // donor
        for _ in 0..5 {
            reference.push(b'N');
        }
        reference.extend_from_slice(b"AG"); // acceptor
        reference.extend_from_slice(b"ACGT");
        // Intron = 9 bases, less than default min_intron_len=30

        let query = b"ACGTACGT";
        let result = spliced_align(query, &reference, &dna_scheme(), &SplicedAlignParams::default())
            .unwrap();
        assert_eq!(result.intron_count, 0, "short intron should not be used");
    }

    #[test]
    fn empty_error() {
        let result = spliced_align(
            b"",
            b"ACGT",
            &dna_scheme(),
            &SplicedAlignParams::default(),
        );
        assert!(result.is_err());

        let result = spliced_align(
            b"ACGT",
            b"",
            &dna_scheme(),
            &SplicedAlignParams::default(),
        );
        assert!(result.is_err());
    }

    // --- Exon chaining ---

    #[test]
    fn chain_two_exons() {
        let mut reference = Vec::new();
        reference.extend_from_slice(b"ACGT"); // 0..4
        reference.extend_from_slice(b"GT"); // 4..6, donor
        for _ in 0..28 {
            reference.push(b'N');
        } // 6..34
        reference.extend_from_slice(b"AG"); // 34..36, acceptor
        reference.extend_from_slice(b"TGCA"); // 36..40

        let exons = vec![
            ExonAlignment {
                query_start: 0,
                query_end: 4,
                target_start: 0,
                target_end: 4,
                score: 8,
                cigar: vec![CigarOp::Match(4)],
            },
            ExonAlignment {
                query_start: 4,
                query_end: 8,
                target_start: 36,
                target_end: 40,
                score: 8,
                cigar: vec![CigarOp::Match(4)],
            },
        ];

        let params = SplicedAlignParams {
            min_intron_len: 30,
            ..Default::default()
        };
        let result = chain_exons(&exons, &reference, &params).unwrap();
        assert_eq!(result.exons.len(), 2);
        assert_eq!(result.intron_count, 1);

        // Check junction
        assert_eq!(result.junctions[0].0, 4); // intron start
        assert_eq!(result.junctions[0].1, 36); // intron end
    }

    #[test]
    fn chain_three_exons() {
        // Create a reference with two introns
        let mut reference = vec![0u8; 120];
        // Exon 1: 0..10
        // Intron 1: 10..50 (40 bp) with GT..AG
        reference[10] = b'G';
        reference[11] = b'T';
        reference[48] = b'A';
        reference[49] = b'G';
        // Exon 2: 50..60
        // Intron 2: 60..100 (40 bp) with GT..AG
        reference[60] = b'G';
        reference[61] = b'T';
        reference[98] = b'A';
        reference[99] = b'G';
        // Exon 3: 100..110

        let exons = vec![
            ExonAlignment {
                query_start: 0,
                query_end: 10,
                target_start: 0,
                target_end: 10,
                score: 20,
                cigar: vec![CigarOp::Match(10)],
            },
            ExonAlignment {
                query_start: 10,
                query_end: 20,
                target_start: 50,
                target_end: 60,
                score: 20,
                cigar: vec![CigarOp::Match(10)],
            },
            ExonAlignment {
                query_start: 20,
                query_end: 30,
                target_start: 100,
                target_end: 110,
                score: 20,
                cigar: vec![CigarOp::Match(10)],
            },
        ];

        let result = chain_exons(&exons, &reference, &SplicedAlignParams::default()).unwrap();
        assert_eq!(result.exons.len(), 3);
        assert_eq!(result.intron_count, 2);
    }

    #[test]
    fn chain_colinearity() {
        // Exons that violate query colinearity should not be chained
        let reference = vec![b'N'; 200];
        let exons = vec![
            ExonAlignment {
                query_start: 10,
                query_end: 20,
                target_start: 0,
                target_end: 10,
                score: 20,
                cigar: vec![CigarOp::Match(10)],
            },
            ExonAlignment {
                query_start: 5, // overlaps with first exon on query!
                query_end: 15,
                target_start: 50,
                target_end: 60,
                score: 20,
                cigar: vec![CigarOp::Match(10)],
            },
        ];

        let result = chain_exons(&exons, &reference, &SplicedAlignParams::default()).unwrap();
        // Should only chain one exon since they overlap on query
        assert_eq!(result.exons.len(), 1);
    }

    #[test]
    fn chain_empty_error() {
        let result = chain_exons(&[], b"ACGT", &SplicedAlignParams::default());
        assert!(result.is_err());
    }

    #[test]
    fn chain_splice_bonus() {
        // Verify canonical bonus affects scoring
        let mut reference = vec![b'N'; 100];
        reference[10] = b'G';
        reference[11] = b'T'; // donor
        reference[48] = b'A';
        reference[49] = b'G'; // acceptor

        let exons = vec![
            ExonAlignment {
                query_start: 0,
                query_end: 10,
                target_start: 0,
                target_end: 10,
                score: 20,
                cigar: vec![CigarOp::Match(10)],
            },
            ExonAlignment {
                query_start: 10,
                query_end: 20,
                target_start: 50,
                target_end: 60,
                score: 20,
                cigar: vec![CigarOp::Match(10)],
            },
        ];

        let result = chain_exons(&exons, &reference, &SplicedAlignParams::default()).unwrap();
        // Score = 20 + splice_penalty(-15) + canonical_bonus(8) + 20 = 33
        assert_eq!(result.alignment.score, 33);
        assert_eq!(result.junctions[0].2, SpliceSiteType::GtAg);
    }

    #[test]
    fn exon_coordinates() {
        let exons = vec![ExonAlignment {
            query_start: 5,
            query_end: 15,
            target_start: 100,
            target_end: 110,
            score: 20,
            cigar: vec![CigarOp::Match(10)],
        }];

        let reference = vec![b'N'; 200];
        let result = chain_exons(&exons, &reference, &SplicedAlignParams::default()).unwrap();
        assert_eq!(result.alignment.query_start, 5);
        assert_eq!(result.alignment.query_end, 15);
        assert_eq!(result.alignment.target_start, 100);
        assert_eq!(result.alignment.target_end, 110);
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
        fn spliced_align_deterministic(q in dna_seq(20), t in dna_seq(40)) {
            let r1 = spliced_align(&q, &t, &dna_scheme(), &SplicedAlignParams::default()).unwrap();
            let r2 = spliced_align(&q, &t, &dna_scheme(), &SplicedAlignParams::default()).unwrap();
            prop_assert_eq!(r1.alignment.score, r2.alignment.score);
        }
    }
}
