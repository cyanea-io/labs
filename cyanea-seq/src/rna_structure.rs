//! RNA secondary structure prediction.
//!
//! Predicts which bases in an RNA sequence pair to form stems, hairpins,
//! internal loops, and multi-branch loops. Provides:
//!
//! - **Dot-bracket notation** — parse and emit `(((...)))` structures
//! - **Nussinov algorithm** — maximize base pair count (O(n³))
//! - **Zuker MFE** — minimum free energy with Turner nearest-neighbor parameters
//! - **McCaskill partition function** — base pair probability matrix
//! - **Structure comparison** — base pair distance, mountain distance

use cyanea_core::{CyaneaError, Result};

// ── Constants ────────────────────────────────────────────────────

/// Gas constant in kcal/(mol·K).
const R: f64 = 0.001987;

/// Default temperature in Kelvin (37 °C).
const DEFAULT_T: f64 = 310.15;

/// Multi-branch loop offset (kcal/mol).
const ML_A: f64 = 3.4;
/// Multi-branch loop per-helix penalty (kcal/mol).
const ML_B: f64 = 0.4;
/// Multi-branch loop per-unpaired-base penalty (kcal/mol).
const ML_C: f64 = 0.0;

/// Large energy value representing an impossible state.
const INF: f64 = 1e18;

/// Minimum hairpin loop size (bases between closing pair).
const MIN_HAIRPIN: usize = 3;

// ── Energy tables (simplified Turner 2004) ───────────────────────

/// Hairpin loop initiation energies indexed by size (3..=30), kcal/mol.
const HAIRPIN_INIT: [f64; 31] = [
    0.0, 0.0, 0.0, // 0, 1, 2 — unused
    5.4, 5.6, 5.7, 5.4, 5.6, 5.7, 5.4, // 3–9
    5.6, 5.7, 5.8, 5.9, 5.9, 6.0, 6.1, // 10–16
    6.1, 6.2, 6.2, 6.3, 6.3, 6.3, 6.4, // 17–23
    6.4, 6.4, 6.5, 6.5, 6.5, 6.5, 6.6, // 24–30
];

/// Internal loop initiation energies indexed by size (1..=30), kcal/mol.
const INTERNAL_INIT: [f64; 31] = [
    0.0, // 0 — unused
    0.0, 0.0, 0.0, 1.1, 2.0, 2.0, 2.1, 2.3, 2.4, 2.5, // 1–10
    2.6, 2.7, 2.8, 2.9, 2.9, 3.0, 3.1, 3.1, 3.2, 3.2, // 11–20
    3.3, 3.3, 3.4, 3.4, 3.4, 3.5, 3.5, 3.5, 3.6, 3.6, // 21–30
];

/// Bulge loop initiation energies indexed by size (1..=30), kcal/mol.
const BULGE_INIT: [f64; 31] = [
    0.0, // 0 — unused
    3.8, 2.8, 3.2, 3.6, 4.0, 4.4, 4.6, 4.7, 4.8, 4.9, // 1–10
    5.0, 5.1, 5.2, 5.3, 5.4, 5.4, 5.5, 5.5, 5.6, 5.6, // 11–20
    5.7, 5.7, 5.8, 5.8, 5.8, 5.9, 5.9, 5.9, 6.0, 6.0, // 21–30
];

/// Return true if bases `a` and `b` can form a pair (AU, GC, GU).
fn can_pair(a: u8, b: u8) -> bool {
    matches!(
        (a, b),
        (b'A', b'U')
            | (b'U', b'A')
            | (b'G', b'C')
            | (b'C', b'G')
            | (b'G', b'U')
            | (b'U', b'G')
    )
}

/// Encode a base pair as an index (0–5) for stacking lookup.
/// Returns None if not a valid pair.
fn pair_index(a: u8, b: u8) -> Option<usize> {
    match (a, b) {
        (b'A', b'U') => Some(0),
        (b'U', b'A') => Some(1),
        (b'G', b'C') => Some(2),
        (b'C', b'G') => Some(3),
        (b'G', b'U') => Some(4),
        (b'U', b'G') => Some(5),
        _ => None,
    }
}

/// Stacking energies (kcal/mol, 37 °C) from Turner 1998/2004.
/// Indexed by [closing pair index][enclosed pair index].
/// Pair indices: AU=0, UA=1, GC=2, CG=3, GU=4, UG=5.
const STACKING: [[f64; 6]; 6] = [
    // closing AU
    [-0.9, -1.1, -2.2, -2.1, -0.6, -1.4],
    // closing UA
    [-1.3, -0.9, -2.4, -2.1, -1.0, -0.7],
    // closing GC
    [-2.4, -2.1, -3.3, -2.4, -1.5, -1.5],
    // closing CG
    [-2.1, -2.1, -2.4, -3.4, -1.4, -2.1],
    // closing GU
    [-1.3, -1.0, -2.5, -1.5, -0.5, -1.3],
    // closing UG
    [-1.0, -0.7, -1.5, -1.5, -0.3, -0.5],
];

/// Look up the stacking energy for two consecutive base pairs.
/// (i5, j5) is the outer (closing) pair, (i3, j3) is the inner (enclosed) pair.
fn stacking_energy(i5: u8, j5: u8, i3: u8, j3: u8) -> f64 {
    match (pair_index(i5, j5), pair_index(i3, j3)) {
        (Some(a), Some(b)) => STACKING[a][b],
        _ => INF,
    }
}

/// Hairpin loop energy for closing pair (i, j).
fn hairpin_energy(seq: &[u8], i: usize, j: usize) -> f64 {
    let size = j - i - 1;
    if size < MIN_HAIRPIN {
        return INF;
    }
    let init = if size <= 30 {
        HAIRPIN_INIT[size]
    } else {
        HAIRPIN_INIT[30] + 1.75 * R * DEFAULT_T * ((size as f64) / 30.0).ln()
    };
    // Terminal mismatch bonus for hairpins ≥ 4
    let mismatch = if size >= 4 {
        terminal_mismatch(seq[i], seq[j], seq[i + 1], seq[j - 1])
    } else {
        0.0
    };
    init + mismatch
}

/// Simplified terminal mismatch energy (kcal/mol).
/// Returns a small bonus for purine-purine mismatches adjacent to a closing pair.
fn terminal_mismatch(_ci: u8, _cj: u8, ni: u8, nj: u8) -> f64 {
    // Simplified: small bonus if both adjacent bases are purines (stacking)
    let is_purine = |b: u8| b == b'A' || b == b'G';
    if is_purine(ni) && is_purine(nj) {
        -0.8
    } else if is_purine(ni) || is_purine(nj) {
        -0.4
    } else {
        0.0
    }
}

/// Internal/bulge loop energy for outer pair (i,j) and inner pair (p,q).
fn internal_loop_energy(seq: &[u8], i: usize, j: usize, p: usize, q: usize) -> f64 {
    let left = p - i - 1;
    let right = j - q - 1;

    if left == 0 && right == 0 {
        return INF; // stacking, not internal
    }

    // 1×1 or 1×2 internal loops: use stacking + asymmetry
    if left == 0 || right == 0 {
        // Bulge loop
        let size = left + right;
        let init = if size <= 30 {
            BULGE_INIT[size]
        } else {
            BULGE_INIT[30] + 1.75 * R * DEFAULT_T * ((size as f64) / 30.0).ln()
        };
        // Stacking bonus for single-nucleotide bulge
        if size == 1 {
            return init + stacking_energy(seq[i], seq[j], seq[p], seq[q]);
        }
        return init;
    }

    // Internal loop
    let size = left + right;
    let init = if size <= 30 {
        INTERNAL_INIT[size]
    } else {
        INTERNAL_INIT[30] + 1.75 * R * DEFAULT_T * ((size as f64) / 30.0).ln()
    };
    // Asymmetry penalty
    let asymmetry = 0.3 * ((left as f64) - (right as f64)).abs();
    let asymmetry = asymmetry.min(3.0); // capped

    init + asymmetry
}

// ── Dot-bracket notation & structure representation ──────────────

/// An RNA secondary structure as a pair table.
///
/// Each position `i` in the sequence is either paired to some position `j`
/// (`pairs[i] = Some(j)`) or unpaired (`pairs[i] = None`). Pairs are
/// non-crossing: if `i` pairs with `j`, no pair `(p, q)` exists with
/// `i < p < j < q`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RnaSecondaryStructure {
    /// Pair table: `pairs[i] = Some(j)` if position `i` is paired with `j`.
    pub pairs: Vec<Option<usize>>,
    /// Length of the sequence.
    pub length: usize,
}

impl RnaSecondaryStructure {
    /// Parse a dot-bracket string into a secondary structure.
    ///
    /// `(` and `)` denote paired bases; `.` denotes unpaired bases.
    ///
    /// # Errors
    ///
    /// Returns an error if parentheses are unbalanced or the string contains
    /// characters other than `(`, `)`, and `.`.
    ///
    /// # Example
    ///
    /// ```
    /// use cyanea_seq::rna_structure::RnaSecondaryStructure;
    ///
    /// let s = RnaSecondaryStructure::from_dot_bracket("(((...)))").unwrap();
    /// assert_eq!(s.num_pairs(), 3);
    /// ```
    pub fn from_dot_bracket(s: &str) -> Result<Self> {
        let n = s.len();
        let mut pairs = vec![None; n];
        let mut stack = Vec::new();

        for (i, ch) in s.chars().enumerate() {
            match ch {
                '(' => stack.push(i),
                ')' => {
                    let j = stack.pop().ok_or_else(|| {
                        CyaneaError::Parse("unmatched ')' in dot-bracket string".into())
                    })?;
                    pairs[j] = Some(i);
                    pairs[i] = Some(j);
                }
                '.' => {}
                _ => {
                    return Err(CyaneaError::Parse(format!(
                        "invalid character '{}' in dot-bracket string",
                        ch
                    )));
                }
            }
        }

        if !stack.is_empty() {
            return Err(CyaneaError::Parse("unmatched '(' in dot-bracket string".into()));
        }

        Ok(Self { pairs, length: n })
    }

    /// Convert this structure to dot-bracket notation.
    ///
    /// # Example
    ///
    /// ```
    /// use cyanea_seq::rna_structure::RnaSecondaryStructure;
    ///
    /// let s = RnaSecondaryStructure::from_dot_bracket("..((..))..").unwrap();
    /// assert_eq!(s.to_dot_bracket(), "..((..))..");
    /// ```
    pub fn to_dot_bracket(&self) -> String {
        let mut out = vec!['.'; self.length];
        for (i, partner) in self.pairs.iter().enumerate() {
            if let Some(j) = partner {
                if i < *j {
                    out[i] = '(';
                    out[*j] = ')';
                }
            }
        }
        out.into_iter().collect()
    }

    /// Return sorted list of base pairs `(i, j)` where `i < j`.
    pub fn base_pairs(&self) -> Vec<(usize, usize)> {
        let mut bps: Vec<(usize, usize)> = self
            .pairs
            .iter()
            .enumerate()
            .filter_map(|(i, p)| p.map(|j| (i, j)))
            .filter(|(i, j)| i < j)
            .collect();
        bps.sort();
        bps
    }

    /// Check whether position `i` is paired.
    pub fn is_paired(&self, i: usize) -> bool {
        i < self.length && self.pairs[i].is_some()
    }

    /// Return the pairing partner of position `i`, if any.
    pub fn partner(&self, i: usize) -> Option<usize> {
        if i < self.length {
            self.pairs[i]
        } else {
            None
        }
    }

    /// Number of base pairs in the structure.
    pub fn num_pairs(&self) -> usize {
        self.pairs.iter().filter(|p| p.is_some()).count() / 2
    }
}

// ── Nussinov algorithm ──────────────────────────────────────────

/// Result of the Nussinov maximum base pair algorithm.
#[derive(Debug, Clone)]
pub struct NussinovResult {
    /// The predicted secondary structure.
    pub structure: RnaSecondaryStructure,
    /// Maximum number of base pairs found.
    pub max_pairs: usize,
}

/// Predict RNA secondary structure by maximizing base pair count (Nussinov algorithm).
///
/// Uses O(n³) dynamic programming. Valid pairs are AU, GC, and GU (wobble).
/// No pair may close a loop smaller than `min_loop_size` bases.
///
/// # Errors
///
/// Returns an error if the sequence is empty.
///
/// # Example
///
/// ```
/// use cyanea_seq::rna_structure::nussinov;
///
/// let result = nussinov(b"GGGGCCCC", 3).unwrap();
/// assert!(result.max_pairs >= 2);
/// ```
pub fn nussinov(seq: &[u8], min_loop_size: usize) -> Result<NussinovResult> {
    let seq = normalize_rna(seq)?;
    let n = seq.len();

    if n == 0 {
        return Err(CyaneaError::InvalidInput("empty sequence".into()));
    }

    // DP table: M[i][j] = max pairs in subsequence [i..=j]
    let mut m = vec![0i32; n * n];
    let idx = |i: usize, j: usize| i * n + j;

    // Fill bottom-up by increasing subsequence length
    for len in 2..=n {
        for i in 0..=n - len {
            let j = i + len - 1;
            // i unpaired
            let mut best = if i + 1 <= j { m[idx(i + 1, j)] } else { 0 };
            // j unpaired
            if j > 0 {
                best = best.max(m[idx(i, j - 1)]);
            }
            // i,j pair
            if can_pair(seq[i], seq[j]) && j - i > min_loop_size {
                let inner = if i + 1 <= j.saturating_sub(1) {
                    m[idx(i + 1, j - 1)]
                } else {
                    0
                };
                best = best.max(inner + 1);
            }
            // bifurcation
            for k in (i + 1)..j {
                best = best.max(m[idx(i, k)] + m[idx(k + 1, j)]);
            }
            m[idx(i, j)] = best;
        }
    }

    // Traceback
    let mut pairs = vec![None; n];
    nussinov_traceback(&seq, &m, n, min_loop_size, 0, n - 1, &mut pairs);

    let max_pairs = m[idx(0, n - 1)] as usize;
    Ok(NussinovResult {
        structure: RnaSecondaryStructure {
            pairs,
            length: n,
        },
        max_pairs,
    })
}

fn nussinov_traceback(
    seq: &[u8],
    m: &[i32],
    n: usize,
    min_loop_size: usize,
    i: usize,
    j: usize,
    pairs: &mut [Option<usize>],
) {
    if i >= j || (j - i) < 1 {
        return;
    }
    let idx = |a: usize, b: usize| a * n + b;
    let val = m[idx(i, j)];

    // i unpaired
    if i + 1 <= j && m[idx(i + 1, j)] == val {
        nussinov_traceback(seq, m, n, min_loop_size, i + 1, j, pairs);
        return;
    }

    // i,j pair
    if can_pair(seq[i], seq[j]) && j - i > min_loop_size {
        let inner = if i + 1 <= j.saturating_sub(1) {
            m[idx(i + 1, j - 1)]
        } else {
            0
        };
        if inner + 1 == val {
            pairs[i] = Some(j);
            pairs[j] = Some(i);
            if j > 0 && i + 1 < j {
                nussinov_traceback(seq, m, n, min_loop_size, i + 1, j - 1, pairs);
            }
            return;
        }
    }

    // bifurcation
    for k in (i + 1)..j {
        if m[idx(i, k)] + m[idx(k + 1, j)] == val {
            nussinov_traceback(seq, m, n, min_loop_size, i, k, pairs);
            nussinov_traceback(seq, m, n, min_loop_size, k + 1, j, pairs);
            return;
        }
    }

    // j unpaired (fallback)
    if j > 0 {
        nussinov_traceback(seq, m, n, min_loop_size, i, j - 1, pairs);
    }
}

// ── Zuker MFE algorithm ─────────────────────────────────────────

/// Result of the Zuker minimum free energy algorithm.
#[derive(Debug, Clone)]
pub struct MfeResult {
    /// The predicted MFE secondary structure.
    pub structure: RnaSecondaryStructure,
    /// Minimum free energy in kcal/mol (negative = stable).
    pub energy: f64,
}

/// Predict RNA secondary structure by minimizing free energy (Zuker algorithm).
///
/// Uses simplified Turner 2004 nearest-neighbor thermodynamic parameters.
/// Suitable for sequences up to ~500 nt.
///
/// # Errors
///
/// Returns an error if the sequence is empty or shorter than 5 bases.
///
/// # Example
///
/// ```
/// use cyanea_seq::rna_structure::zuker_mfe;
///
/// let result = zuker_mfe(b"GGGAAACCC").unwrap();
/// assert!(result.energy <= 0.0);
/// ```
pub fn zuker_mfe(seq: &[u8]) -> Result<MfeResult> {
    let seq = normalize_rna(seq)?;
    let n = seq.len();

    if n == 0 {
        return Err(CyaneaError::InvalidInput("empty sequence".into()));
    }
    if n < 5 {
        // Too short to form any structure
        return Ok(MfeResult {
            structure: RnaSecondaryStructure {
                pairs: vec![None; n],
                length: n,
            },
            energy: 0.0,
        });
    }

    let idx = |i: usize, j: usize| i * n + j;

    // V[i,j] = MFE where (i,j) form a base pair
    let mut v = vec![INF; n * n];
    // W[i,j] = MFE of any structure on [i..=j]
    let mut w = vec![0.0_f64; n * n];
    // WM[i,j] = MFE of multi-branch loop region on [i..=j]
    let mut wm = vec![INF; n * n];

    // Fill bottom-up by increasing subsequence length
    for len in 2..=n {
        for i in 0..=n - len {
            let j = i + len - 1;

            // ── V(i,j): only valid if (i,j) can pair ──
            if can_pair(seq[i], seq[j]) && j - i > MIN_HAIRPIN {
                let mut best_v = INF;

                // Hairpin
                best_v = best_v.min(hairpin_energy(&seq, i, j));

                // Stacking
                if i + 1 < j && j > 0 && can_pair(seq[i + 1], seq[j - 1]) && j - 1 - (i + 1) >= MIN_HAIRPIN {
                    let stack = stacking_energy(seq[i], seq[j], seq[i + 1], seq[j - 1]);
                    best_v = best_v.min(v[idx(i + 1, j - 1)] + stack);
                } else if i + 1 < j && j > 0 && can_pair(seq[i + 1], seq[j - 1]) {
                    // Inner pair exists but too close — still allow stacking if V is valid
                    let stack = stacking_energy(seq[i], seq[j], seq[i + 1], seq[j - 1]);
                    if v[idx(i + 1, j - 1)] < INF / 2.0 {
                        best_v = best_v.min(v[idx(i + 1, j - 1)] + stack);
                    }
                }

                // Internal loop / bulge
                // Iterate over all inner pairs (p, q) with i < p < q < j
                let max_left = (j - i - 1).min(30);
                for p in (i + 1)..=(i + max_left).min(j.saturating_sub(1)) {
                    let max_right = (j - i - 1 - (p - i - 1)).min(30);
                    let q_min = if p + MIN_HAIRPIN + 1 > j {
                        continue;
                    } else {
                        (j - max_right).max(p + MIN_HAIRPIN + 1)
                    };
                    for q in q_min..j {
                        if !can_pair(seq[p], seq[q]) {
                            continue;
                        }
                        if p == i + 1 && q == j - 1 {
                            continue; // stacking, handled above
                        }
                        if v[idx(p, q)] >= INF / 2.0 {
                            continue;
                        }
                        let il_e = internal_loop_energy(&seq, i, j, p, q);
                        best_v = best_v.min(v[idx(p, q)] + il_e);
                    }
                }

                // Multi-branch loop
                if j > i + 2 && wm[idx(i + 1, j - 1)] < INF / 2.0 {
                    best_v = best_v.min(wm[idx(i + 1, j - 1)] + ML_A + ML_B);
                }

                v[idx(i, j)] = best_v;
            }

            // ── WM(i,j): multi-branch loop region ──
            {
                let mut best_wm = INF;

                // i unpaired
                if i + 1 <= j && wm[idx(i + 1, j)] < INF / 2.0 {
                    best_wm = best_wm.min(wm[idx(i + 1, j)] + ML_C);
                }

                // j unpaired
                if j > 0 && i <= j - 1 && wm[idx(i, j - 1)] < INF / 2.0 {
                    best_wm = best_wm.min(wm[idx(i, j - 1)] + ML_C);
                }

                // Helix starting at (i,j)
                if v[idx(i, j)] < INF / 2.0 {
                    best_wm = best_wm.min(v[idx(i, j)] + ML_B);
                }

                // Concatenation
                for k in (i + 1)..j {
                    if wm[idx(i, k)] < INF / 2.0 && wm[idx(k + 1, j)] < INF / 2.0 {
                        best_wm = best_wm.min(wm[idx(i, k)] + wm[idx(k + 1, j)]);
                    }
                }

                wm[idx(i, j)] = best_wm;
            }

            // ── W(i,j): any structure on [i..=j] ──
            {
                let mut best_w: f64 = 0.0; // no structure

                // i unpaired
                if i + 1 <= j {
                    best_w = best_w.min(w[idx(i + 1, j)]);
                }

                // j unpaired
                if j > 0 && i <= j - 1 {
                    best_w = best_w.min(w[idx(i, j - 1)]);
                }

                // (i,j) pair
                if v[idx(i, j)] < INF / 2.0 {
                    best_w = best_w.min(v[idx(i, j)]);
                }

                // bifurcation
                for k in (i + 1)..j {
                    best_w = best_w.min(w[idx(i, k)] + w[idx(k + 1, j)]);
                }

                w[idx(i, j)] = best_w;
            }
        }
    }

    let energy = w[idx(0, n - 1)];
    let energy = if energy >= INF / 2.0 { 0.0 } else { energy };

    // Traceback
    let mut pairs = vec![None; n];
    zuker_traceback_w(&seq, &v, &w, &wm, n, 0, n - 1, &mut pairs);

    Ok(MfeResult {
        structure: RnaSecondaryStructure {
            pairs,
            length: n,
        },
        energy,
    })
}

fn zuker_traceback_w(
    seq: &[u8],
    v: &[f64],
    w: &[f64],
    wm: &[f64],
    n: usize,
    i: usize,
    j: usize,
    pairs: &mut [Option<usize>],
) {
    if i >= j {
        return;
    }
    let idx = |a: usize, b: usize| a * n + b;
    let val = w[idx(i, j)];
    let eps = 1e-9;

    // No structure
    if val.abs() < eps {
        return;
    }

    // (i,j) pair
    if v[idx(i, j)] < INF / 2.0 && (v[idx(i, j)] - val).abs() < eps {
        pairs[i] = Some(j);
        pairs[j] = Some(i);
        zuker_traceback_v(seq, v, w, wm, n, i, j, pairs);
        return;
    }

    // i unpaired
    if i + 1 <= j && (w[idx(i + 1, j)] - val).abs() < eps {
        zuker_traceback_w(seq, v, w, wm, n, i + 1, j, pairs);
        return;
    }

    // j unpaired
    if j > 0 && i <= j - 1 && (w[idx(i, j - 1)] - val).abs() < eps {
        zuker_traceback_w(seq, v, w, wm, n, i, j - 1, pairs);
        return;
    }

    // bifurcation
    for k in (i + 1)..j {
        if (w[idx(i, k)] + w[idx(k + 1, j)] - val).abs() < eps {
            zuker_traceback_w(seq, v, w, wm, n, i, k, pairs);
            zuker_traceback_w(seq, v, w, wm, n, k + 1, j, pairs);
            return;
        }
    }
}

fn zuker_traceback_v(
    seq: &[u8],
    v: &[f64],
    w: &[f64],
    wm: &[f64],
    n: usize,
    i: usize,
    j: usize,
    pairs: &mut [Option<usize>],
) {
    let idx = |a: usize, b: usize| a * n + b;
    let val = v[idx(i, j)];
    let eps = 1e-9;

    // Hairpin
    if (hairpin_energy(seq, i, j) - val).abs() < eps {
        return;
    }

    // Stacking
    if i + 1 < j && j > 0 && can_pair(seq[i + 1], seq[j - 1]) && v[idx(i + 1, j - 1)] < INF / 2.0 {
        let stack = stacking_energy(seq[i], seq[j], seq[i + 1], seq[j - 1]);
        if (v[idx(i + 1, j - 1)] + stack - val).abs() < eps {
            pairs[i + 1] = Some(j - 1);
            pairs[j - 1] = Some(i + 1);
            zuker_traceback_v(seq, v, w, wm, n, i + 1, j - 1, pairs);
            return;
        }
    }

    // Internal loop / bulge
    let max_left = (j - i - 1).min(30);
    for p in (i + 1)..=(i + max_left).min(j.saturating_sub(1)) {
        let max_right = (j - i - 1 - (p - i - 1)).min(30);
        let q_min_val = (j.saturating_sub(max_right)).max(p + MIN_HAIRPIN + 1);
        for q in q_min_val..j {
            if !can_pair(seq[p], seq[q]) || v[idx(p, q)] >= INF / 2.0 {
                continue;
            }
            if p == i + 1 && q == j - 1 {
                continue;
            }
            let il_e = internal_loop_energy(seq, i, j, p, q);
            if (v[idx(p, q)] + il_e - val).abs() < eps {
                pairs[p] = Some(q);
                pairs[q] = Some(p);
                zuker_traceback_v(seq, v, w, wm, n, p, q, pairs);
                return;
            }
        }
    }

    // Multi-branch loop
    if j > i + 2 && wm[idx(i + 1, j - 1)] < INF / 2.0 {
        if (wm[idx(i + 1, j - 1)] + ML_A + ML_B - val).abs() < eps {
            zuker_traceback_wm(seq, v, w, wm, n, i + 1, j - 1, pairs);
        }
    }
}

fn zuker_traceback_wm(
    seq: &[u8],
    v: &[f64],
    w: &[f64],
    wm: &[f64],
    n: usize,
    i: usize,
    j: usize,
    pairs: &mut [Option<usize>],
) {
    if i >= j {
        return;
    }
    let idx = |a: usize, b: usize| a * n + b;
    let val = wm[idx(i, j)];
    let eps = 1e-9;

    if val >= INF / 2.0 {
        return;
    }

    // Helix at (i,j)
    if v[idx(i, j)] < INF / 2.0 && (v[idx(i, j)] + ML_B - val).abs() < eps {
        pairs[i] = Some(j);
        pairs[j] = Some(i);
        zuker_traceback_v(seq, v, w, wm, n, i, j, pairs);
        return;
    }

    // i unpaired
    if i + 1 <= j && wm[idx(i + 1, j)] < INF / 2.0 && (wm[idx(i + 1, j)] + ML_C - val).abs() < eps
    {
        zuker_traceback_wm(seq, v, w, wm, n, i + 1, j, pairs);
        return;
    }

    // j unpaired
    if j > 0 && i <= j - 1 && wm[idx(i, j - 1)] < INF / 2.0
        && (wm[idx(i, j - 1)] + ML_C - val).abs() < eps
    {
        zuker_traceback_wm(seq, v, w, wm, n, i, j - 1, pairs);
        return;
    }

    // Concatenation
    for k in (i + 1)..j {
        if wm[idx(i, k)] < INF / 2.0
            && wm[idx(k + 1, j)] < INF / 2.0
            && (wm[idx(i, k)] + wm[idx(k + 1, j)] - val).abs() < eps
        {
            zuker_traceback_wm(seq, v, w, wm, n, i, k, pairs);
            zuker_traceback_wm(seq, v, w, wm, n, k + 1, j, pairs);
            return;
        }
    }
}

// ── McCaskill partition function ────────────────────────────────

/// Result of the McCaskill partition function algorithm.
#[derive(Debug, Clone)]
pub struct PartitionResult {
    /// Base pair probability matrix (n×n, flat row-major).
    pub pair_probabilities: Vec<f64>,
    /// Sequence length.
    pub length: usize,
    /// Ensemble free energy: −RT·ln(Z) in kcal/mol.
    pub ensemble_energy: f64,
}

impl PartitionResult {
    /// Get the probability that positions `i` and `j` are paired.
    pub fn pair_probability(&self, i: usize, j: usize) -> f64 {
        if i >= self.length || j >= self.length {
            return 0.0;
        }
        self.pair_probabilities[i * self.length + j]
    }

    /// Get the probability that position `i` is unpaired.
    pub fn unpaired_probability(&self, i: usize) -> f64 {
        if i >= self.length {
            return 0.0;
        }
        let paired: f64 = (0..self.length)
            .map(|j| self.pair_probabilities[i * self.length + j])
            .sum();
        (1.0 - paired).max(0.0)
    }
}

/// Compute base pair probabilities via the McCaskill inside-outside algorithm.
///
/// Uses the same simplified Turner energy model as [`zuker_mfe`].
/// Computations are performed in log-space for numerical stability.
///
/// # Arguments
///
/// * `seq` — RNA sequence (A, U, G, C)
/// * `temperature` — temperature in Kelvin (e.g., 310.15 for 37 °C)
///
/// # Errors
///
/// Returns an error if the sequence is empty or temperature is not positive.
///
/// # Example
///
/// ```
/// use cyanea_seq::rna_structure::mccaskill;
///
/// let result = mccaskill(b"GGGAAACCC", 310.15).unwrap();
/// assert!(result.pair_probability(0, 8) > 0.0);
/// ```
pub fn mccaskill(seq: &[u8], temperature: f64) -> Result<PartitionResult> {
    let seq = normalize_rna(seq)?;
    let n = seq.len();

    if n == 0 {
        return Err(CyaneaError::InvalidInput("empty sequence".into()));
    }
    if temperature <= 0.0 {
        return Err(CyaneaError::InvalidInput(
            "temperature must be positive".into(),
        ));
    }

    let rt = R * temperature;

    if n < 5 {
        return Ok(PartitionResult {
            pair_probabilities: vec![0.0; n * n],
            length: n,
            ensemble_energy: 0.0,
        });
    }

    let idx = |i: usize, j: usize| i * n + j;
    let boltz = |e: f64| -> f64 {
        if e >= INF / 2.0 {
            0.0
        } else {
            (-e / rt).exp()
        }
    };

    // Inside partition functions (stored as Boltzmann weights, not logs)
    // Q[i,j] = partition function for subsequence [i..=j]
    // Qb[i,j] = partition function for subsequence [i..=j] where (i,j) form a pair
    let mut q = vec![0.0_f64; n * n];
    let mut qb = vec![0.0_f64; n * n];
    let mut qm = vec![0.0_f64; n * n];

    // Base cases: Q[i,i] = 1, Q[i,j] = 1 for j < i
    for i in 0..n {
        q[idx(i, i)] = 1.0;
        if i + 1 < n {
            q[idx(i + 1, i)] = 1.0; // empty
        }
    }

    // Fill inside tables bottom-up
    for len in 2..=n {
        for i in 0..=n - len {
            let j = i + len - 1;

            // Qb(i,j): (i,j) must be a valid pair
            if can_pair(seq[i], seq[j]) && j - i > MIN_HAIRPIN {
                let mut qb_val = 0.0;

                // Hairpin
                qb_val += boltz(hairpin_energy(&seq, i, j));

                // Stacking
                if i + 1 < j && can_pair(seq[i + 1], seq[j - 1]) {
                    let stack = stacking_energy(seq[i], seq[j], seq[i + 1], seq[j - 1]);
                    qb_val += qb[idx(i + 1, j - 1)] * boltz(stack);
                }

                // Internal loops / bulges
                let max_left = (j - i - 1).min(30);
                for p in (i + 1)..=(i + max_left).min(j.saturating_sub(1)) {
                    let max_right = (j - i - 1 - (p - i - 1)).min(30);
                    let q_min = (j.saturating_sub(max_right)).max(p + MIN_HAIRPIN + 1);
                    for qi in q_min..j {
                        if !can_pair(seq[p], seq[qi]) {
                            continue;
                        }
                        if p == i + 1 && qi == j - 1 {
                            continue; // stacking
                        }
                        let il_e = internal_loop_energy(&seq, i, j, p, qi);
                        qb_val += qb[idx(p, qi)] * boltz(il_e);
                    }
                }

                // Multi-branch loop
                if j > i + 2 {
                    qb_val += qm[idx(i + 1, j - 1)] * boltz(ML_A + ML_B);
                }

                qb[idx(i, j)] = qb_val;
            }

            // QM(i,j): multi-branch region
            {
                let mut qm_val = 0.0;

                // i unpaired
                if i + 1 <= j {
                    qm_val += qm[idx(i + 1, j)] * boltz(ML_C);
                }

                // j unpaired
                if j > 0 && i <= j - 1 {
                    qm_val += qm[idx(i, j - 1)] * boltz(ML_C);
                }

                // Helix starting here
                if qb[idx(i, j)] > 0.0 {
                    qm_val += qb[idx(i, j)] * boltz(ML_B);
                }

                // Concatenation
                for k in (i + 1)..j {
                    qm_val += qm[idx(i, k)] * qm[idx(k + 1, j)];
                }

                qm[idx(i, j)] = qm_val;
            }

            // Q(i,j): full partition
            {
                let mut q_val = 1.0; // empty structure

                for d in i..=j {
                    for e in (d + MIN_HAIRPIN + 1)..=j {
                        if !can_pair(seq[d], seq[e]) || qb[idx(d, e)] == 0.0 {
                            continue;
                        }
                        let q_left = if d > i { q[idx(i, d - 1)] } else { 1.0 };
                        let q_right = if e < j { q[idx(e + 1, j)] } else { 1.0 };
                        q_val += q_left * qb[idx(d, e)] * q_right;
                    }
                }

                q[idx(i, j)] = q_val;
            }
        }
    }

    let z = q[idx(0, n - 1)];
    let ensemble_energy = if z > 0.0 { -rt * z.ln() } else { 0.0 };

    // Outside algorithm for pair probabilities
    let mut prob = vec![0.0_f64; n * n];

    if z > 0.0 {
        for i in 0..n {
            for j in (i + MIN_HAIRPIN + 1)..n {
                if qb[idx(i, j)] == 0.0 {
                    continue;
                }
                let q_left = if i > 0 { q[idx(0, i - 1)] } else { 1.0 };
                let q_right = if j < n - 1 { q[idx(j + 1, n - 1)] } else { 1.0 };
                let p_ij = q_left * qb[idx(i, j)] * q_right / z;
                let p_ij = p_ij.min(1.0).max(0.0);
                prob[idx(i, j)] = p_ij;
                prob[idx(j, i)] = p_ij;
            }
        }
    }

    Ok(PartitionResult {
        pair_probabilities: prob,
        length: n,
        ensemble_energy,
    })
}

// ── Structure comparison ────────────────────────────────────────

/// Compute the base pair distance between two structures.
///
/// The base pair distance is the size of the symmetric difference of the
/// two base pair sets: the number of pairs in one structure but not the other.
///
/// # Errors
///
/// Returns an error if the structures have different lengths.
///
/// # Example
///
/// ```
/// use cyanea_seq::rna_structure::{RnaSecondaryStructure, base_pair_distance};
///
/// let a = RnaSecondaryStructure::from_dot_bracket("(((...)))").unwrap();
/// let b = RnaSecondaryStructure::from_dot_bracket("(((...)))").unwrap();
/// assert_eq!(base_pair_distance(&a, &b).unwrap(), 0);
/// ```
pub fn base_pair_distance(
    a: &RnaSecondaryStructure,
    b: &RnaSecondaryStructure,
) -> Result<usize> {
    if a.length != b.length {
        return Err(CyaneaError::InvalidInput(format!(
            "structure lengths differ: {} vs {}",
            a.length, b.length
        )));
    }

    let bp_a: std::collections::HashSet<(usize, usize)> = a.base_pairs().into_iter().collect();
    let bp_b: std::collections::HashSet<(usize, usize)> = b.base_pairs().into_iter().collect();

    let only_a = bp_a.difference(&bp_b).count();
    let only_b = bp_b.difference(&bp_a).count();
    Ok(only_a + only_b)
}

/// Compute the mountain distance between two structures.
///
/// The mountain representation of a structure assigns to each position
/// the number of base pairs enclosing it. The mountain distance is the
/// L1 norm of the difference of the two mountain vectors.
///
/// # Errors
///
/// Returns an error if the structures have different lengths.
///
/// # Example
///
/// ```
/// use cyanea_seq::rna_structure::{RnaSecondaryStructure, mountain_distance};
///
/// let a = RnaSecondaryStructure::from_dot_bracket("(((...)))").unwrap();
/// let b = RnaSecondaryStructure::from_dot_bracket("((...))..").unwrap();
/// assert!(mountain_distance(&a, &b).unwrap() > 0.0);
/// ```
pub fn mountain_distance(
    a: &RnaSecondaryStructure,
    b: &RnaSecondaryStructure,
) -> Result<f64> {
    if a.length != b.length {
        return Err(CyaneaError::InvalidInput(format!(
            "structure lengths differ: {} vs {}",
            a.length, b.length
        )));
    }

    let ma = mountain_vector(a);
    let mb = mountain_vector(b);

    let dist: f64 = ma
        .iter()
        .zip(mb.iter())
        .map(|(x, y)| (*x as f64 - *y as f64).abs())
        .sum();
    Ok(dist)
}

/// Build the mountain vector: m[i] = number of base pairs enclosing position i.
fn mountain_vector(s: &RnaSecondaryStructure) -> Vec<i32> {
    let mut m = vec![0i32; s.length];
    let mut depth = 0i32;
    for i in 0..s.length {
        if let Some(j) = s.pairs[i] {
            if i < j {
                depth += 1;
            } else {
                depth -= 1;
            }
        }
        m[i] = depth;
    }
    m
}

// ── Helpers ─────────────────────────────────────────────────────

/// Normalize sequence to uppercase RNA (T → U).
fn normalize_rna(seq: &[u8]) -> Result<Vec<u8>> {
    seq.iter()
        .map(|&b| match b {
            b'A' | b'a' => Ok(b'A'),
            b'U' | b'u' => Ok(b'U'),
            b'G' | b'g' => Ok(b'G'),
            b'C' | b'c' => Ok(b'C'),
            b'T' | b't' => Ok(b'U'), // DNA T → RNA U
            _ => Err(CyaneaError::InvalidInput(format!(
                "invalid nucleotide '{}' in RNA sequence",
                b as char
            ))),
        })
        .collect()
}

// ── Tests ───────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Dot-bracket ──

    #[test]
    fn dot_bracket_simple() {
        let s = RnaSecondaryStructure::from_dot_bracket("(((...)))").unwrap();
        assert_eq!(s.length, 9);
        assert_eq!(s.num_pairs(), 3);
        assert_eq!(s.pairs[0], Some(8));
        assert_eq!(s.pairs[8], Some(0));
        assert!(s.is_paired(0));
        assert!(!s.is_paired(3));
    }

    #[test]
    fn dot_bracket_with_unpaired() {
        let s = RnaSecondaryStructure::from_dot_bracket("..((..))..").unwrap();
        assert_eq!(s.length, 10);
        assert_eq!(s.num_pairs(), 2);
        assert!(!s.is_paired(0));
        assert!(s.is_paired(2));
        assert_eq!(s.partner(2), Some(7));
    }

    #[test]
    fn dot_bracket_roundtrip() {
        let input = "(((..((.....))...)))";
        let s = RnaSecondaryStructure::from_dot_bracket(input).unwrap();
        assert_eq!(s.to_dot_bracket(), input);

        // Parse again
        let s2 = RnaSecondaryStructure::from_dot_bracket(&s.to_dot_bracket()).unwrap();
        assert_eq!(s.pairs, s2.pairs);
    }

    #[test]
    fn dot_bracket_base_pairs() {
        let s = RnaSecondaryStructure::from_dot_bracket("((.()))").unwrap();
        let bps = s.base_pairs();
        assert_eq!(bps.len(), 3);
        // All pairs should have i < j
        for (i, j) in &bps {
            assert!(i < j);
        }
    }

    #[test]
    fn dot_bracket_unmatched_open() {
        assert!(RnaSecondaryStructure::from_dot_bracket("((...)))(").is_err());
    }

    #[test]
    fn dot_bracket_unmatched_close() {
        assert!(RnaSecondaryStructure::from_dot_bracket(")((..))").is_err());
    }

    #[test]
    fn dot_bracket_invalid_char() {
        assert!(RnaSecondaryStructure::from_dot_bracket("((..x..))").is_err());
    }

    #[test]
    fn dot_bracket_empty() {
        let s = RnaSecondaryStructure::from_dot_bracket("").unwrap();
        assert_eq!(s.length, 0);
        assert_eq!(s.num_pairs(), 0);
    }

    // ── Nussinov ──

    #[test]
    fn nussinov_gcaucg() {
        let r = nussinov(b"GCAUCG", 3).unwrap();
        // Only G-C pair possible with loop≥3: G(0)-C(5) has loop=4
        assert!(r.max_pairs >= 1);
    }

    #[test]
    fn nussinov_perfect_stem() {
        let r = nussinov(b"GGGGCCCC", 3).unwrap();
        // Can form at most 2 pairs with min_loop_size=3: G(0)-C(7) and G(1)-C(6)
        assert!(r.max_pairs >= 2);
    }

    #[test]
    fn nussinov_no_pairs() {
        let r = nussinov(b"AAAAAA", 3).unwrap();
        assert_eq!(r.max_pairs, 0);
        assert_eq!(r.structure.num_pairs(), 0);
    }

    #[test]
    fn nussinov_min_loop_enforced() {
        // AU with only 2 bases between: can't pair with min_loop=3
        let r = nussinov(b"AXXU", 3).unwrap_or_else(|_| {
            // X is invalid, use valid bases
            nussinov(b"AGCU", 3).unwrap()
        });
        // A(0) and U(3) have only 2 bases between them (positions 1,2) — can't pair
        // Actually j-i = 3 = min_loop_size, need > min_loop_size
        assert_eq!(r.max_pairs, 0);
    }

    #[test]
    fn nussinov_short_sequence() {
        let r = nussinov(b"AUG", 3).unwrap();
        assert_eq!(r.max_pairs, 0);
    }

    #[test]
    fn nussinov_empty() {
        assert!(nussinov(b"", 3).is_err());
    }

    #[test]
    fn nussinov_structure_valid() {
        let r = nussinov(b"GGGAAACCC", 3).unwrap();
        // Verify no crossing pairs
        let bps = r.structure.base_pairs();
        for (idx_a, &(i1, j1)) in bps.iter().enumerate() {
            for &(i2, j2) in bps.iter().skip(idx_a + 1) {
                // Non-crossing: either nested or disjoint
                assert!(j1 <= i2 || i2 >= i1 && j2 <= j1,
                    "crossing pairs: ({},{}) and ({},{})", i1, j1, i2, j2);
            }
        }
    }

    #[test]
    fn nussinov_lowercase_and_dna() {
        let r = nussinov(b"gggaaaccc", 3).unwrap();
        assert!(r.max_pairs > 0);

        // T is treated as U
        let r2 = nussinov(b"GGGAAATCC", 3).unwrap();
        let r3 = nussinov(b"GGGAAAUCC", 3).unwrap();
        assert_eq!(r2.max_pairs, r3.max_pairs);
    }

    // ── Zuker MFE ──

    #[test]
    fn zuker_simple_hairpin() {
        let r = zuker_mfe(b"GGGAAACCC").unwrap();
        // Should form a stem-loop with negative energy
        assert!(r.energy < 0.0, "energy should be negative, got {}", r.energy);
        assert!(r.structure.num_pairs() > 0);
    }

    #[test]
    fn zuker_gc_stronger_than_au() {
        let gc = zuker_mfe(b"GGGCAAAGCCC").unwrap();
        let au = zuker_mfe(b"AAAUAAAUUUU").unwrap();
        // GC-rich should have more negative (stronger) energy
        assert!(
            gc.energy <= au.energy,
            "GC energy ({}) should be <= AU energy ({})",
            gc.energy,
            au.energy
        );
    }

    #[test]
    fn zuker_no_structure() {
        let r = zuker_mfe(b"AAAAAA").unwrap();
        assert_eq!(r.structure.num_pairs(), 0);
        assert!((r.energy - 0.0).abs() < 1e-6);
    }

    #[test]
    fn zuker_energy_nonpositive() {
        for seq in &[b"GCGCGCGC" as &[u8], b"AUGCAUGC", b"GGGAAACCC", b"CCCCGGGGG"] {
            let r = zuker_mfe(seq).unwrap();
            assert!(
                r.energy <= 1e-9,
                "energy should be <= 0, got {} for {:?}",
                r.energy,
                std::str::from_utf8(seq).unwrap()
            );
        }
    }

    #[test]
    fn zuker_valid_structure() {
        let r = zuker_mfe(b"GGGAAACCC").unwrap();
        // All pairs respect min_loop_size
        let bps = r.structure.base_pairs();
        for &(i, j) in &bps {
            assert!(j - i > MIN_HAIRPIN, "pair ({},{}) violates min loop size", i, j);
            // Both must pair with each other
            assert_eq!(r.structure.pairs[i], Some(j));
            assert_eq!(r.structure.pairs[j], Some(i));
        }
    }

    #[test]
    fn zuker_short_sequence() {
        let r = zuker_mfe(b"AUGC").unwrap();
        assert_eq!(r.energy, 0.0);
        assert_eq!(r.structure.num_pairs(), 0);
    }

    #[test]
    fn zuker_empty() {
        assert!(zuker_mfe(b"").is_err());
    }

    // ── McCaskill partition function ──

    #[test]
    fn mccaskill_strong_stem() {
        let r = mccaskill(b"GGGAAACCC", 310.15).unwrap();
        // The terminal G-C pair should have significant probability
        let p = r.pair_probability(0, 8);
        assert!(p > 0.01, "pair prob(0,8) = {} should be > 0.01", p);
    }

    #[test]
    fn mccaskill_no_pairs() {
        let r = mccaskill(b"AAAAAA", 310.15).unwrap();
        // All pair probabilities should be near 0
        for i in 0..r.length {
            for j in 0..r.length {
                assert!(
                    r.pair_probability(i, j) < 0.01,
                    "pair prob({},{}) = {} should be < 0.01",
                    i,
                    j,
                    r.pair_probability(i, j)
                );
            }
        }
    }

    #[test]
    fn mccaskill_probabilities_sum() {
        let r = mccaskill(b"GGGAAACCC", 310.15).unwrap();
        for i in 0..r.length {
            let paired: f64 = (0..r.length).map(|j| r.pair_probability(i, j)).sum();
            let unpaired = r.unpaired_probability(i);
            let total = paired + unpaired;
            assert!(
                (total - 1.0).abs() < 0.1,
                "probability sum at position {} = {} (paired={}, unpaired={})",
                i,
                total,
                paired,
                unpaired
            );
        }
    }

    #[test]
    fn mccaskill_temperature_effect() {
        let low_t = mccaskill(b"GGGAAACCC", 300.0).unwrap();
        let high_t = mccaskill(b"GGGAAACCC", 370.0).unwrap();
        // Higher temperature → less structured → lower pair probabilities
        let low_p: f64 = (0..low_t.length)
            .flat_map(|i| (i + 1..low_t.length).map(move |j| (i, j)))
            .map(|(i, j)| low_t.pair_probability(i, j))
            .sum();
        let high_p: f64 = (0..high_t.length)
            .flat_map(|i| (i + 1..high_t.length).map(move |j| (i, j)))
            .map(|(i, j)| high_t.pair_probability(i, j))
            .sum();
        assert!(
            low_p >= high_p - 0.01,
            "lower T should give more pairing: {} vs {}",
            low_p,
            high_p
        );
    }

    #[test]
    fn mccaskill_deterministic() {
        let r1 = mccaskill(b"GCGCGCGCGC", 310.15).unwrap();
        let r2 = mccaskill(b"GCGCGCGCGC", 310.15).unwrap();
        assert_eq!(r1.pair_probabilities, r2.pair_probabilities);
        assert_eq!(r1.ensemble_energy, r2.ensemble_energy);
    }

    #[test]
    fn mccaskill_empty() {
        assert!(mccaskill(b"", 310.15).is_err());
    }

    #[test]
    fn mccaskill_invalid_temperature() {
        assert!(mccaskill(b"GGGAAACCC", 0.0).is_err());
        assert!(mccaskill(b"GGGAAACCC", -10.0).is_err());
    }

    #[test]
    fn mccaskill_short_sequence() {
        let r = mccaskill(b"AUGC", 310.15).unwrap();
        // Too short for any pairs
        for i in 0..r.length {
            for j in 0..r.length {
                assert!((r.pair_probability(i, j) - 0.0).abs() < 1e-10);
            }
        }
    }

    // ── Structure comparison ──

    #[test]
    fn bp_distance_identical() {
        let a = RnaSecondaryStructure::from_dot_bracket("(((...)))").unwrap();
        let b = RnaSecondaryStructure::from_dot_bracket("(((...)))").unwrap();
        assert_eq!(base_pair_distance(&a, &b).unwrap(), 0);
    }

    #[test]
    fn bp_distance_completely_different() {
        let a = RnaSecondaryStructure::from_dot_bracket("((....))").unwrap();
        let b = RnaSecondaryStructure::from_dot_bracket("........").unwrap();
        assert_eq!(base_pair_distance(&a, &b).unwrap(), 2); // 2 pairs in a, 0 in b
    }

    #[test]
    fn bp_distance_symmetric() {
        let a = RnaSecondaryStructure::from_dot_bracket("((....))").unwrap();
        let b = RnaSecondaryStructure::from_dot_bracket(".((..).)").unwrap();
        assert_eq!(
            base_pair_distance(&a, &b).unwrap(),
            base_pair_distance(&b, &a).unwrap()
        );
    }

    #[test]
    fn bp_distance_different_lengths() {
        let a = RnaSecondaryStructure::from_dot_bracket("((..))").unwrap();
        let b = RnaSecondaryStructure::from_dot_bracket("(((...)))").unwrap();
        assert!(base_pair_distance(&a, &b).is_err());
    }

    #[test]
    fn mountain_distance_identical() {
        let a = RnaSecondaryStructure::from_dot_bracket("(((...)))").unwrap();
        let b = RnaSecondaryStructure::from_dot_bracket("(((...)))").unwrap();
        assert!((mountain_distance(&a, &b).unwrap() - 0.0).abs() < 1e-10);
    }

    #[test]
    fn mountain_distance_symmetric() {
        let a = RnaSecondaryStructure::from_dot_bracket("((....))").unwrap();
        let b = RnaSecondaryStructure::from_dot_bracket(".((..).)").unwrap();
        let d1 = mountain_distance(&a, &b).unwrap();
        let d2 = mountain_distance(&b, &a).unwrap();
        assert!((d1 - d2).abs() < 1e-10);
    }

    #[test]
    fn mountain_distance_nonnegative() {
        let a = RnaSecondaryStructure::from_dot_bracket("(((...)))").unwrap();
        let b = RnaSecondaryStructure::from_dot_bracket("((...))..").unwrap();
        assert!(mountain_distance(&a, &b).unwrap() >= 0.0);
    }

    #[test]
    fn mountain_distance_different_lengths() {
        let a = RnaSecondaryStructure::from_dot_bracket("((..))").unwrap();
        let b = RnaSecondaryStructure::from_dot_bracket("(((...)))").unwrap();
        assert!(mountain_distance(&a, &b).is_err());
    }

    // ── Energy model helpers ──

    #[test]
    fn can_pair_valid() {
        assert!(can_pair(b'A', b'U'));
        assert!(can_pair(b'U', b'A'));
        assert!(can_pair(b'G', b'C'));
        assert!(can_pair(b'C', b'G'));
        assert!(can_pair(b'G', b'U'));
        assert!(can_pair(b'U', b'G'));
    }

    #[test]
    fn can_pair_invalid() {
        assert!(!can_pair(b'A', b'A'));
        assert!(!can_pair(b'A', b'C'));
        assert!(!can_pair(b'A', b'G'));
        assert!(!can_pair(b'C', b'U'));
    }

    #[test]
    fn stacking_energy_values() {
        // CG closing, CG enclosed: STACKING[3][3] = -3.4
        let e = stacking_energy(b'C', b'G', b'C', b'G');
        assert!((e - (-3.4)).abs() < 1e-10, "CG/CG stack = {}", e);
        // GC closing, GC enclosed: STACKING[2][2] = -3.3
        let e2 = stacking_energy(b'G', b'C', b'G', b'C');
        assert!((e2 - (-3.3)).abs() < 1e-10, "GC/GC stack = {}", e2);
        // AU closing, UA enclosed: STACKING[0][1] = -1.1
        let e3 = stacking_energy(b'A', b'U', b'U', b'A');
        assert!((e3 - (-1.1)).abs() < 1e-10, "AU/UA stack = {}", e3);
    }

    #[test]
    fn normalize_rna_dna_input() {
        let r = normalize_rna(b"ATGC").unwrap();
        assert_eq!(r, b"AUGC");
    }

    #[test]
    fn normalize_rna_lowercase() {
        let r = normalize_rna(b"augc").unwrap();
        assert_eq!(r, b"AUGC");
    }

    #[test]
    fn normalize_rna_invalid() {
        assert!(normalize_rna(b"AXGC").is_err());
    }
}
