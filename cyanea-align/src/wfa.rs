//! Wavefront Alignment (WFA) for gap-affine scoring.
//!
//! WFA computes optimal global alignment in O(ns) time and space, where n is
//! sequence length and s is the optimal alignment score (penalty).  This is
//! dramatically faster than O(nm) Needleman-Wunsch for similar sequences.
//!
//! Reference: Santiago Marco-Sola et al., "Fast gap-affine pairwise alignment
//! using the wavefront algorithm", *Bioinformatics*, 2021.
//!
//! # Example
//!
//! ```
//! use cyanea_align::{wfa_align, ScoringMatrix, ScoringScheme};
//!
//! let scoring = ScoringScheme::Simple(ScoringMatrix::dna_default());
//! let result = wfa_align(b"ACGTACGT", b"ACGTACGT", &scoring).unwrap();
//! assert_eq!(result.score, 16); // 8 matches * 2
//! ```

use crate::scoring::ScoringScheme;
use crate::types::{AlignmentResult, CigarOp};
use cyanea_core::{CyaneaError, Result};

// ---------------------------------------------------------------------------
// Wavefront data structures
// ---------------------------------------------------------------------------

/// A single wavefront: offsets along diagonals `lo..=hi`.
///
/// Convention: diagonal `k = i - j` where `i` indexes the query and `j` the
/// target.  The offset at diagonal `k` is the furthest query position `i`
/// reached.  Target position is `j = i - k`.
#[derive(Debug, Clone)]
struct Wavefront {
    lo: i32,
    hi: i32,
    /// Offset for each diagonal in `[lo, hi]`.  Index with `(k - lo) as usize`.
    offsets: Vec<i32>,
}

impl Wavefront {
    fn new(lo: i32, hi: i32) -> Self {
        let len = (hi - lo + 1).max(0) as usize;
        Self {
            lo,
            hi,
            offsets: vec![-1; len],
        }
    }

    #[inline]
    fn get(&self, k: i32) -> i32 {
        if k < self.lo || k > self.hi {
            return -1;
        }
        self.offsets[(k - self.lo) as usize]
    }

    #[inline]
    fn set(&mut self, k: i32, val: i32) {
        if k >= self.lo && k <= self.hi {
            self.offsets[(k - self.lo) as usize] = val;
        }
    }
}

/// A set of three wavefronts (M, I, D) at a given score.
#[derive(Debug, Clone)]
struct WavefrontSet {
    m: Option<Wavefront>,
    i: Option<Wavefront>,
    d: Option<Wavefront>,
}

// ---------------------------------------------------------------------------
// Backtrace encoding
// ---------------------------------------------------------------------------

/// Source of an M-wavefront entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MSource {
    /// Mismatch from M[s-x][k].
    Mismatch,
    /// Came from I[s][k].
    FromI,
    /// Came from D[s][k].
    FromD,
}

/// Source of an I-wavefront entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ISource {
    /// Gap open from M[s-o-e][k+1].
    Open,
    /// Gap extend from I[s-e][k+1].
    Extend,
}

/// Source of a D-wavefront entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DSource {
    /// Gap open from M[s-o-e][k-1].
    Open,
    /// Gap extend from D[s-e][k-1].
    Extend,
}

/// Which wavefront type we're tracing through.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum WfType {
    M,
    I,
    D,
}

/// Traceback info for one score level.
#[derive(Debug, Clone)]
struct TraceLevel {
    lo: i32,
    m: Vec<Option<MSource>>,
    i: Vec<Option<ISource>>,
    d: Vec<Option<DSource>>,
}

impl TraceLevel {
    fn new(lo: i32, len: usize) -> Self {
        Self {
            lo,
            m: vec![None; len],
            i: vec![None; len],
            d: vec![None; len],
        }
    }

    fn idx(&self, k: i32) -> usize {
        (k - self.lo) as usize
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Perform global alignment using the Wavefront Alignment algorithm with
/// gap-affine scoring.
///
/// WFA is optimal for aligning similar sequences where the edit distance is
/// small relative to sequence length.  For divergent sequences, classic
/// Needleman-Wunsch may be more memory-efficient.
///
/// # Errors
///
/// - Returns an error if either sequence is empty.
/// - Returns an error if a `ScoringScheme::Substitution` matrix is provided
///   (WFA currently only supports simple match/mismatch scoring).
pub fn wfa_align(
    query: &[u8],
    target: &[u8],
    scoring: &ScoringScheme,
) -> Result<AlignmentResult> {
    if query.is_empty() || target.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "sequences must not be empty".into(),
        ));
    }

    let sm = match scoring {
        ScoringScheme::Simple(m) => m,
        ScoringScheme::Substitution(_) => {
            return Err(CyaneaError::InvalidInput(
                "WFA only supports simple scoring".into(),
            ));
        }
    };

    // WFA adjusted penalties (positive values).
    //
    // Traditional score = matches*m + mismatches*mm + gap_cost.
    // We rewrite as: score = n*m - penalty, where penalty accounts for
    // mismatches (lost match vs mismatch), deletions (lost match + gap cost),
    // and insertions (gap cost only).  This requires asymmetric gap costs:
    //
    //   x  = match_score - mismatch_score  (mismatch penalty)
    //   o  = |gap_open| - |gap_extend|     (gap-open component, shared)
    //   e_d = |gap_extend| + match_score   (deletion extend: gap + lost match)
    //   e_i = |gap_extend|                 (insertion extend: gap only)
    //
    // After WFA: traditional_score = n * match_score - wfa_penalty.
    let x = (sm.match_score - sm.mismatch_score) as u32;
    let abs_ge = (-sm.gap_extend) as u32;
    let abs_go = (-sm.gap_open) as u32;
    let o = abs_go.saturating_sub(abs_ge);
    let e_d = abs_ge + sm.match_score as u32; // deletion extend
    let e_i = abs_ge; // insertion extend
    let match_score = sm.match_score;

    let n = query.len();
    let m = target.len();
    let target_k = n as i32 - m as i32;
    let target_h = n as i32;

    // Upper bound on penalty score
    let max_score = ((n + m) as u32) * (x + o + e_d + e_i);
    let cap = max_score as usize + 1;

    // Storage indexed by penalty score
    let mut wf: Vec<Option<WavefrontSet>> = vec![None; cap];
    let mut trace: Vec<Option<TraceLevel>> = vec![None; cap];

    // Initialise score 0: M[0] on diagonal 0 with offset 0
    let mut m0 = Wavefront::new(0, 0);
    m0.set(0, 0);
    extend_wf(&mut m0, query, target);

    wf[0] = Some(WavefrontSet {
        m: Some(m0),
        i: None,
        d: None,
    });
    trace[0] = Some(TraceLevel::new(0, 1));

    if done(&wf[0], target_k, target_h) {
        return build_result(query, target, &wf, &trace, 0, target_k, match_score, sm);
    }

    // Main WFA loop over increasing penalty scores
    for s in 1..cap {
        let (lo, hi) = diag_range_asym(s, &wf, x, o, e_d, e_i);
        if lo > hi {
            continue;
        }
        let len = (hi - lo + 1) as usize;
        let mut tl = TraceLevel::new(lo, len);

        // ---- I wavefront (insertion = gap in query, target consumed) ----
        // Uses e_i (no lost match for insertions).
        let mut i_wf = Wavefront::new(lo, hi);
        let oe_i = o + e_i;
        if s as u32 >= oe_i {
            let ps = (s as u32 - oe_i) as usize;
            if let Some(ref ws) = wf.get(ps).and_then(|w| w.as_ref()) {
                if let Some(ref mw) = ws.m {
                    for k in lo..=hi {
                        let v = mw.get(k + 1);
                        if v >= 0 && v > i_wf.get(k) {
                            i_wf.set(k, v);
                            tl.i[(k - tl.lo) as usize] = Some(ISource::Open);
                        }
                    }
                }
            }
        }
        if s as u32 >= e_i {
            let ps = (s as u32 - e_i) as usize;
            if let Some(ref ws) = wf.get(ps).and_then(|w| w.as_ref()) {
                if let Some(ref iw) = ws.i {
                    for k in lo..=hi {
                        let v = iw.get(k + 1);
                        if v >= 0 && v > i_wf.get(k) {
                            i_wf.set(k, v);
                            tl.i[(k - tl.lo) as usize] = Some(ISource::Extend);
                        }
                    }
                }
            }
        }

        // ---- D wavefront (deletion = gap in target, query consumed) ----
        // Uses e_d (includes lost match reward per deletion base).
        let mut d_wf = Wavefront::new(lo, hi);
        let oe_d = o + e_d;
        if s as u32 >= oe_d {
            let ps = (s as u32 - oe_d) as usize;
            if let Some(ref ws) = wf.get(ps).and_then(|w| w.as_ref()) {
                if let Some(ref mw) = ws.m {
                    for k in lo..=hi {
                        let v = mw.get(k - 1);
                        if v >= 0 {
                            let nv = v + 1;
                            if nv > d_wf.get(k) {
                                d_wf.set(k, nv);
                                tl.d[(k - tl.lo) as usize] = Some(DSource::Open);
                            }
                        }
                    }
                }
            }
        }
        if s as u32 >= e_d {
            let ps = (s as u32 - e_d) as usize;
            if let Some(ref ws) = wf.get(ps).and_then(|w| w.as_ref()) {
                if let Some(ref dw) = ws.d {
                    for k in lo..=hi {
                        let v = dw.get(k - 1);
                        if v >= 0 {
                            let nv = v + 1;
                            if nv > d_wf.get(k) {
                                d_wf.set(k, nv);
                                tl.d[(k - tl.lo) as usize] = Some(DSource::Extend);
                            }
                        }
                    }
                }
            }
        }

        // ---- M wavefront ----
        // M[s][k] = max(M[s-x][k]+1, I[s][k], D[s][k])
        let mut m_wf = Wavefront::new(lo, hi);

        // From mismatch
        if s as u32 >= x {
            let ps = (s as u32 - x) as usize;
            if let Some(ref ws) = wf.get(ps).and_then(|w| w.as_ref()) {
                if let Some(ref mw) = ws.m {
                    for k in lo..=hi {
                        let v = mw.get(k);
                        if v >= 0 {
                            let nv = v + 1;
                            if nv > m_wf.get(k) {
                                m_wf.set(k, nv);
                                tl.m[(k - tl.lo) as usize] = Some(MSource::Mismatch);
                            }
                        }
                    }
                }
            }
        }

        // From I
        for k in lo..=hi {
            let v = i_wf.get(k);
            if v >= 0 && v > m_wf.get(k) {
                m_wf.set(k, v);
                tl.m[(k - tl.lo) as usize] = Some(MSource::FromI);
            }
        }

        // From D
        for k in lo..=hi {
            let v = d_wf.get(k);
            if v >= 0 && v > m_wf.get(k) {
                m_wf.set(k, v);
                tl.m[(k - tl.lo) as usize] = Some(MSource::FromD);
            }
        }

        // Extend M along matching characters
        extend_wf(&mut m_wf, query, target);

        wf[s] = Some(WavefrontSet {
            m: Some(m_wf),
            i: Some(i_wf),
            d: Some(d_wf),
        });
        trace[s] = Some(tl);

        if done(&wf[s], target_k, target_h) {
            return build_result(
                query,
                target,
                &wf,
                &trace,
                s as u32,
                target_k,
                match_score,
                sm,
            );
        }
    }

    Err(CyaneaError::InvalidInput(
        "WFA failed to converge".into(),
    ))
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Extend all diagonals in a wavefront as far as characters match.
fn extend_wf(wf: &mut Wavefront, query: &[u8], target: &[u8]) {
    let n = query.len() as i32;
    let m = target.len() as i32;
    for k in wf.lo..=wf.hi {
        let idx = (k - wf.lo) as usize;
        let mut h = wf.offsets[idx];
        if h < 0 {
            continue;
        }
        loop {
            let j = h - k;
            if h >= n || j < 0 || j >= m {
                break;
            }
            if query[h as usize].to_ascii_uppercase()
                != target[j as usize].to_ascii_uppercase()
            {
                break;
            }
            h += 1;
        }
        wf.offsets[idx] = h;
    }
}

/// Check if the M wavefront reaches the target.
fn done(wfs: &Option<WavefrontSet>, target_k: i32, target_h: i32) -> bool {
    if let Some(ref ws) = wfs {
        if let Some(ref m) = ws.m {
            return m.get(target_k) >= target_h;
        }
    }
    false
}

/// Compute diagonal range needed at score `s` with asymmetric gap costs.
fn diag_range_asym(
    s: usize,
    wf: &[Option<WavefrontSet>],
    x: u32,
    o: u32,
    e_d: u32,
    e_i: u32,
) -> (i32, i32) {
    let mut lo = i32::MAX;
    let mut hi = i32::MIN;

    let expand = |score: usize, lo: &mut i32, hi: &mut i32, delta: i32| {
        if let Some(ref ws) = wf.get(score).and_then(|w| w.as_ref()) {
            for wfopt in [&ws.m, &ws.i, &ws.d] {
                if let Some(ref w) = wfopt {
                    *lo = (*lo).min(w.lo + delta);
                    *hi = (*hi).max(w.hi + delta);
                }
            }
        }
    };

    // From mismatch: M[s-x][k]
    if s as u32 >= x {
        expand((s as u32 - x) as usize, &mut lo, &mut hi, 0);
    }
    // From I gap open: M[s-o-e_i][k+1] → delta -1 for lo
    let oe_i = o + e_i;
    if s as u32 >= oe_i {
        expand((s as u32 - oe_i) as usize, &mut lo, &mut hi, -1);
    }
    // From I gap extend: I[s-e_i][k+1] → delta -1
    if s as u32 >= e_i {
        expand((s as u32 - e_i) as usize, &mut lo, &mut hi, -1);
    }
    // From D gap open: M[s-o-e_d][k-1] → delta +1 for hi
    let oe_d = o + e_d;
    if s as u32 >= oe_d {
        expand((s as u32 - oe_d) as usize, &mut lo, &mut hi, 1);
    }
    // From D gap extend: D[s-e_d][k-1] → delta +1
    if s as u32 >= e_d {
        expand((s as u32 - e_d) as usize, &mut lo, &mut hi, 1);
    }

    (lo, hi)
}

/// Reconstruct alignment from the traceback.
fn build_result(
    query: &[u8],
    target: &[u8],
    wf: &[Option<WavefrontSet>],
    trace: &[Option<TraceLevel>],
    final_score: u32,
    target_k: i32,
    match_score: i32,
    sm: &crate::scoring::ScoringMatrix,
) -> Result<AlignmentResult> {
    let n = query.len();
    let m = target.len();

    // We'll collect operations in reverse order, then reverse at the end.
    let mut ops: Vec<CigarOp> = Vec::new();

    // WFA penalties (must match forward pass)
    let x = (sm.match_score - sm.mismatch_score) as u32;
    let abs_ge = (-sm.gap_extend) as u32;
    let abs_go = (-sm.gap_open) as u32;
    let o = abs_go.saturating_sub(abs_ge);
    let e_d = abs_ge + sm.match_score as u32;
    let e_i = abs_ge;

    // Current position in the traceback
    let mut s = final_score;
    let mut k = target_k;
    let mut wf_type = WfType::M;

    loop {
        if s == 0 && wf_type == WfType::M {
            // Remaining offset at score 0 is all initial matches.
            let h = get_offset(wf, 0, k, WfType::M);
            for _ in 0..h {
                ops.push(CigarOp::Match(1));
            }
            break;
        }

        match wf_type {
            WfType::M => {
                let tl = trace[s as usize].as_ref().unwrap();
                let idx = tl.idx(k);
                let src = tl.m.get(idx).copied().flatten();
                let src = match src {
                    Some(src) => src,
                    None => break,
                };

                // Post-extend offset (stored in wf)
                let h_post = get_offset(wf, s as usize, k, WfType::M);
                // Pre-extend offset depends on source
                let h_pre = match src {
                    MSource::Mismatch => {
                        let ps = (s - x) as usize;
                        get_offset(wf, ps, k, WfType::M) + 1
                    }
                    MSource::FromI => get_offset(wf, s as usize, k, WfType::I),
                    MSource::FromD => get_offset(wf, s as usize, k, WfType::D),
                };

                // Record matches from extension
                let num_matches = h_post - h_pre;
                for _ in 0..num_matches {
                    ops.push(CigarOp::Match(1));
                }

                // Follow the source
                match src {
                    MSource::Mismatch => {
                        ops.push(CigarOp::Mismatch(1));
                        s -= x;
                        // k stays, wf_type stays M
                    }
                    MSource::FromI => {
                        wf_type = WfType::I;
                        // s and k stay — we descend into I[s][k]
                    }
                    MSource::FromD => {
                        wf_type = WfType::D;
                    }
                }
            }
            WfType::I => {
                let tl = trace[s as usize].as_ref().unwrap();
                let idx = tl.idx(k);
                let src = tl.i.get(idx).copied().flatten();
                let src = match src {
                    Some(src) => src,
                    None => break,
                };

                // Record an insertion (target consumed, gap in query)
                ops.push(CigarOp::Insertion(1));

                match src {
                    ISource::Open => {
                        s -= o + e_i;
                        k += 1;
                        wf_type = WfType::M;
                    }
                    ISource::Extend => {
                        s -= e_i;
                        k += 1;
                    }
                }
            }
            WfType::D => {
                let tl = trace[s as usize].as_ref().unwrap();
                let idx = tl.idx(k);
                let src = tl.d.get(idx).copied().flatten();
                let src = match src {
                    Some(src) => src,
                    None => break,
                };

                ops.push(CigarOp::Deletion(1));

                match src {
                    DSource::Open => {
                        s -= o + e_d;
                        k -= 1;
                        wf_type = WfType::M;
                    }
                    DSource::Extend => {
                        s -= e_d;
                        k -= 1;
                    }
                }
            }
        }
    }

    // Reverse to forward order
    ops.reverse();

    // Merge consecutive same-type CIGAR ops
    let mut cigar: Vec<CigarOp> = Vec::new();
    for op in &ops {
        push_cigar(&mut cigar, *op);
    }

    // Build aligned sequences from CIGAR
    let mut aligned_query = Vec::with_capacity(ops.len());
    let mut aligned_target = Vec::with_capacity(ops.len());
    let mut qi = 0usize;
    let mut ti = 0usize;
    for op in &ops {
        match op {
            CigarOp::Match(1) => {
                aligned_query.push(query[qi]);
                aligned_target.push(target[ti]);
                qi += 1;
                ti += 1;
            }
            CigarOp::Mismatch(1) => {
                aligned_query.push(query[qi]);
                aligned_target.push(target[ti]);
                qi += 1;
                ti += 1;
            }
            CigarOp::Insertion(1) => {
                aligned_query.push(b'-');
                aligned_target.push(target[ti]);
                ti += 1;
            }
            CigarOp::Deletion(1) => {
                aligned_query.push(query[qi]);
                aligned_target.push(b'-');
                qi += 1;
            }
            _ => {}
        }
    }

    // Compute traditional alignment score
    let num_matches = cigar
        .iter()
        .map(|op| match op {
            CigarOp::Match(n) => *n as i32,
            _ => 0,
        })
        .sum::<i32>();
    let num_mismatches = cigar
        .iter()
        .map(|op| match op {
            CigarOp::Mismatch(n) => *n as i32,
            _ => 0,
        })
        .sum::<i32>();
    let mut gap_opens = 0i32;
    let mut total_gap_bases = 0i32;
    for op in &cigar {
        match op {
            CigarOp::Insertion(n) | CigarOp::Deletion(n) => {
                gap_opens += 1;
                total_gap_bases += *n as i32;
            }
            _ => {}
        }
    }

    // NW gap model: gap(L) = gap_open + (L-1) * gap_extend
    // = gap_opens * gap_open + (total_gap_bases - gap_opens) * gap_extend
    let score = num_matches * match_score
        + num_mismatches * sm.mismatch_score
        + gap_opens * sm.gap_open
        + (total_gap_bases - gap_opens) * sm.gap_extend;

    Ok(AlignmentResult {
        score,
        aligned_query,
        aligned_target,
        query_start: 0,
        query_end: n,
        target_start: 0,
        target_end: m,
        cigar,
    })
}

/// Get the offset stored in a wavefront.
fn get_offset(wf: &[Option<WavefrontSet>], s: usize, k: i32, wf_type: WfType) -> i32 {
    if let Some(ref ws) = wf.get(s).and_then(|w| w.as_ref()) {
        match wf_type {
            WfType::M => ws.m.as_ref().map_or(-1, |w| w.get(k)),
            WfType::I => ws.i.as_ref().map_or(-1, |w| w.get(k)),
            WfType::D => ws.d.as_ref().map_or(-1, |w| w.get(k)),
        }
    } else {
        -1
    }
}

/// Merge a new 1-length CIGAR op with the last op if they are the same variant.
fn push_cigar(ops: &mut Vec<CigarOp>, op: CigarOp) {
    if let Some(last) = ops.last_mut() {
        let merged = match (last, &op) {
            (CigarOp::Match(n), CigarOp::Match(_)) => {
                *n += 1;
                true
            }
            (CigarOp::Mismatch(n), CigarOp::Mismatch(_)) => {
                *n += 1;
                true
            }
            (CigarOp::Insertion(n), CigarOp::Insertion(_)) => {
                *n += 1;
                true
            }
            (CigarOp::Deletion(n), CigarOp::Deletion(_)) => {
                *n += 1;
                true
            }
            _ => false,
        };
        if merged {
            return;
        }
    }
    ops.push(op);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::needleman_wunsch::needleman_wunsch;
    use crate::scoring::{ScoringMatrix, SubstitutionMatrix};

    fn dna_scheme() -> ScoringScheme {
        ScoringScheme::Simple(ScoringMatrix::dna_default())
    }

    #[test]
    fn identical_sequences() {
        let result = wfa_align(b"ACGT", b"ACGT", &dna_scheme()).unwrap();
        assert_eq!(result.score, 8); // 4 matches * 2
        assert_eq!(result.cigar_string(), "4=");
        assert!((result.identity() - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn identical_longer() {
        let seq = b"ACGTACGTACGT";
        let result = wfa_align(seq, seq, &dna_scheme()).unwrap();
        assert_eq!(result.score, 24); // 12 matches * 2
        assert_eq!(result.cigar_string(), "12=");
    }

    #[test]
    fn single_mismatch() {
        let result = wfa_align(b"ACGT", b"ACAT", &dna_scheme()).unwrap();
        let nw = needleman_wunsch(b"ACGT", b"ACAT", &dna_scheme()).unwrap();
        assert_eq!(result.score, nw.score);
        assert_eq!(result.matches(), 3);
        assert_eq!(result.mismatches(), 1);
    }

    #[test]
    fn single_insertion() {
        // query = ACGT, target = ACGGT (extra G in target = insertion)
        let result = wfa_align(b"ACGT", b"ACGGT", &dna_scheme()).unwrap();
        let nw = needleman_wunsch(b"ACGT", b"ACGGT", &dna_scheme()).unwrap();
        assert_eq!(result.score, nw.score);
        assert!(result.gaps() > 0);
    }

    #[test]
    fn single_deletion() {
        // query = ACGGT, target = ACGT (extra G in query = deletion)
        let result = wfa_align(b"ACGGT", b"ACGT", &dna_scheme()).unwrap();
        let nw = needleman_wunsch(b"ACGGT", b"ACGT", &dna_scheme()).unwrap();
        assert_eq!(result.score, nw.score);
        assert!(result.gaps() > 0);
    }

    #[test]
    fn empty_sequence_errors() {
        assert!(wfa_align(b"", b"ACGT", &dna_scheme()).is_err());
        assert!(wfa_align(b"ACGT", b"", &dna_scheme()).is_err());
    }

    #[test]
    fn substitution_matrix_rejected() {
        let scheme = ScoringScheme::Substitution(SubstitutionMatrix::blosum62());
        let err = wfa_align(b"HEAG", b"HEAG", &scheme).unwrap_err();
        assert!(err.to_string().contains("WFA only supports simple scoring"));
    }

    #[test]
    fn score_matches_nw() {
        let cases: Vec<(&[u8], &[u8])> = vec![
            (b"ACGT", b"ACGT"),
            (b"AAAA", b"TTTT"),
            (b"ACGTACGT", b"ACGACGT"),
            (b"GATTACA", b"GCATGCU"),
            (b"A", b"A"),
            (b"AC", b"AG"),
            (b"A", b"AA"),
            (b"AA", b"A"),
        ];
        let scheme = dna_scheme();
        for (q, t) in cases {
            let wfa = wfa_align(q, t, &scheme).unwrap();
            let nw = needleman_wunsch(q, t, &scheme).unwrap();
            assert_eq!(
                wfa.score, nw.score,
                "WFA score {} != NW score {} for {:?} vs {:?}",
                wfa.score,
                nw.score,
                core::str::from_utf8(q).unwrap_or("?"),
                core::str::from_utf8(t).unwrap_or("?")
            );
        }
    }

    #[test]
    fn long_similar_sequences() {
        let query = vec![b'A'; 200];
        let mut target = query.clone();
        target[50] = b'T';
        target[100] = b'G';
        target[150] = b'C';

        let result = wfa_align(&query, &target, &dna_scheme()).unwrap();
        let nw = needleman_wunsch(&query, &target, &dna_scheme()).unwrap();
        assert_eq!(result.score, nw.score);
        assert_eq!(result.matches(), 197);
        assert_eq!(result.mismatches(), 3);
    }

    #[test]
    fn single_base_match() {
        let result = wfa_align(b"A", b"A", &dna_scheme()).unwrap();
        assert_eq!(result.score, 2);
        assert_eq!(result.cigar_string(), "1=");
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use crate::needleman_wunsch::needleman_wunsch;
    use crate::scoring::ScoringMatrix;
    use proptest::prelude::*;

    fn dna_seq(max_len: usize) -> impl Strategy<Value = Vec<u8>> {
        proptest::collection::vec(
            prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')],
            1..=max_len,
        )
    }

    proptest! {
        #[test]
        fn wfa_score_equals_nw_score(
            q in dna_seq(30),
            t in dna_seq(30),
        ) {
            let scheme = ScoringScheme::Simple(ScoringMatrix::dna_default());
            let wfa = wfa_align(&q, &t, &scheme).unwrap();
            let nw = needleman_wunsch(&q, &t, &scheme).unwrap();
            prop_assert_eq!(wfa.score, nw.score,
                "WFA={} NW={} for q={:?} t={:?}",
                wfa.score, nw.score,
                String::from_utf8_lossy(&q),
                String::from_utf8_lossy(&t));
        }
    }
}
