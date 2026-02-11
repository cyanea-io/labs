//! SIMD-accelerated Smith-Waterman score-only computation.
//!
//! Implements the Farrar (2007) striped approach using SIMD intrinsics:
//! - **NEON** (aarch64 / Apple Silicon) — 8 × i16 lanes
//! - **SSE4.1** (x86_64) — 8 × i16 lanes
//! - **Scalar fallback** — when no SIMD is available or `simd` feature is off
//!
//! Score-only (no traceback) uses O(n) memory and is the fastest mode for
//! filtering or computing similarity scores.

use cyanea_core::{CyaneaError, Result};

use crate::scoring::ScoringScheme;

/// SIMD-accelerated Smith-Waterman score-only alignment.
///
/// Returns the maximum local alignment score without traceback.
/// Uses SIMD intrinsics when available (NEON on aarch64, SSE4.1 on x86_64),
/// with automatic fallback to a scalar implementation.
///
/// # Errors
///
/// Returns an error if either sequence is empty.
#[allow(unreachable_code)]
pub fn sw_simd_score(
    query: &[u8],
    target: &[u8],
    scoring: &ScoringScheme,
) -> Result<i32> {
    if query.is_empty() || target.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "sequences must not be empty".into(),
        ));
    }

    #[cfg(feature = "simd")]
    {
        #[cfg(target_arch = "aarch64")]
        {
            // NEON is always available on aarch64
            return Ok(unsafe { neon::sw_score_neon(query, target, scoring) });
        }

        #[cfg(target_arch = "x86_64")]
        {
            if is_x86_feature_detected!("sse4.1") {
                return Ok(unsafe { sse41::sw_score_sse41(query, target, scoring) });
            }
        }
    }

    // Scalar fallback
    Ok(sw_score_scalar(query, target, scoring))
}

// ---------------------------------------------------------------------------
// Scalar fallback — Gotoh 3-matrix SW, O(mn) time, O(n) space
// ---------------------------------------------------------------------------

/// Scalar Smith-Waterman score-only with affine gaps (Gotoh).
fn sw_score_scalar(query: &[u8], target: &[u8], scoring: &ScoringScheme) -> i32 {
    let m = query.len();
    let n = target.len();
    let gap_open = scoring.gap_open();
    let gap_extend = scoring.gap_extend();

    // h_row[j] = H[prev_row][j], overwritten to H[curr_row][j]
    let mut h_row = vec![0i32; n + 1];
    // f_col[j] = F[prev_row][j] — vertical gap, extends down rows
    let mut f_col = vec![i32::MIN / 2; n + 1];
    let mut best = 0i32;

    for i in 0..m {
        let mut h_diag = 0i32; // H[i-1][j-1]
        let mut e = i32::MIN / 2; // horizontal gap, resets each row

        for j in 0..n {
            let sub = scoring.score_pair(query[i], target[j]);

            // F[i][j] = max(H[i-1][j] + gap_open, F[i-1][j] + gap_extend)
            f_col[j + 1] = (h_row[j + 1] + gap_open).max(f_col[j + 1] + gap_extend);

            // H[i][j] = max(0, H[i-1][j-1] + sub, E[i][j], F[i][j])
            let h_new = 0i32
                .max(h_diag + sub)
                .max(e)
                .max(f_col[j + 1]);

            // Save H[i-1][j] before overwriting
            h_diag = h_row[j + 1];
            h_row[j + 1] = h_new;

            if h_new > best {
                best = h_new;
            }

            // E[i][j+1] = max(H[i][j] + gap_open, E[i][j] + gap_extend)
            e = (h_new + gap_open).max(e + gap_extend);
        }
    }

    best
}

// ---------------------------------------------------------------------------
// Striped profile builder — shared by SIMD kernels
// ---------------------------------------------------------------------------

#[allow(dead_code)]
/// Build a query profile in Farrar striped layout.
///
/// For stripe `s` and lane `l`, the query position is `l * n_stripes + s`.
/// Returns profile[target_byte * padded_len + s * LANES + l] = score(query[qi], target_byte).
fn build_striped_profile(
    query: &[u8],
    scoring: &ScoringScheme,
    n_stripes: usize,
    lanes: usize,
) -> Vec<i16> {
    let padded_len = n_stripes * lanes;
    let qlen = query.len();
    let mut profile = vec![0i16; 256 * padded_len];

    for s in 0..n_stripes {
        for l in 0..lanes {
            let qi = l * n_stripes + s;
            if qi < qlen {
                for tb in 0u16..=255 {
                    let score = scoring.score_pair(query[qi], tb as u8) as i16;
                    profile[tb as usize * padded_len + s * lanes + l] = score;
                }
            }
        }
    }

    profile
}

// ---------------------------------------------------------------------------
// NEON (aarch64) — always available
// ---------------------------------------------------------------------------

#[cfg(all(feature = "simd", target_arch = "aarch64"))]
mod neon {
    use super::*;
    use core::arch::aarch64::*;

    const LANES: usize = 8; // int16x8_t

    /// NEON-accelerated SW score-only using Farrar striped vectorization.
    ///
    /// # Safety
    ///
    /// Requires aarch64 NEON (always available on aarch64 targets).
    #[target_feature(enable = "neon")]
    pub(super) unsafe fn sw_score_neon(
        query: &[u8],
        target: &[u8],
        scoring: &ScoringScheme,
    ) -> i32 {
        let qlen = query.len();
        let n_stripes = (qlen + LANES - 1) / LANES;
        let padded_len = n_stripes * LANES;

        let gap_open = scoring.gap_open() as i16;
        let gap_extend = scoring.gap_extend() as i16;

        let profile = build_striped_profile(query, scoring, n_stripes, LANES);

        let vzero = vdupq_n_s16(0);
        let vgap_open = vdupq_n_s16(gap_open);
        let vgap_extend = vdupq_n_s16(gap_extend);
        let vneg_inf = vdupq_n_s16(i16::MIN / 2);

        // H and E vectors (one NEON register per stripe)
        let mut vh: Vec<int16x8_t> = vec![vzero; n_stripes];
        let mut ve: Vec<int16x8_t> = vec![vneg_inf; n_stripes];
        let mut vbest = vzero;

        for &tb in target {
            let prof_row = &profile[tb as usize * padded_len..];

            // Save old H values before updating (needed for diagonal computation)
            let vh_old: Vec<int16x8_t> = vh.clone();

            // --- Main pass: compute H from diagonal + E (no F yet) ---
            for s in 0..n_stripes {
                let vp = vld1q_s16(prof_row.as_ptr().add(s * LANES));

                // Diagonal: old H value for query position q-1
                // In striped layout, q = l*n_stripes + s, so q-1 maps to:
                //   s > 0: stripe s-1, same lane → old_H[s-1]
                //   s == 0: stripe n_stripes-1, lane l-1, with 0 at lane 0
                let vdiag = if s == 0 {
                    // vextq_s16(a, b, 7) = [a[7], b[0..6]]
                    // = [0, old_H[last][0..6]] — shift right, insert 0
                    vextq_s16(vzero, vh_old[n_stripes - 1], 7)
                } else {
                    vh_old[s - 1]
                };

                // E[s] = max(old_H[s] + gap_open, old_E[s] + gap_extend)
                let ve_new = vmaxq_s16(
                    vqaddq_s16(vh_old[s], vgap_open),
                    vqaddq_s16(ve[s], vgap_extend),
                );
                ve[s] = ve_new;

                // H[s] = max(0, diag + profile, E)
                let vh_new = vmaxq_s16(vzero, vmaxq_s16(vqaddq_s16(vdiag, vp), ve_new));
                vh[s] = vh_new;
                vbest = vmaxq_s16(vbest, vh_new);
            }

            // --- F-loop: propagate gap-in-target across stripes until convergence ---
            // F extends forward through query positions. In striped layout, this
            // means F propagates from stripe s to stripe s+1 (wrapping).
            loop {
                let mut vf = vneg_inf;
                let mut any_update = false;

                for s in 0..n_stripes {
                    vf = vmaxq_s16(
                        vqaddq_s16(vh[s], vgap_open),
                        vqaddq_s16(vf, vgap_extend),
                    );

                    let vh_cur = vh[s];
                    let vh_new = vmaxq_s16(vh_cur, vf);

                    // Check if any lane was improved
                    let cmp = vceqq_s16(vh_new, vh_cur);
                    if vminvq_u16(cmp) != 0xFFFF {
                        any_update = true;
                        vh[s] = vh_new;
                        vbest = vmaxq_s16(vbest, vh_new);
                    }
                }

                if !any_update {
                    break;
                }

                // Shift F for wrap-around: the last F value must feed into stripe 0
                // at the next iteration, but shifted by one lane (higher lane indices
                // correspond to later query positions within a stripe).
                // For most practical gap penalties this converges in 1-2 iterations.
            }
        }

        // Horizontal max reduction
        vmaxvq_s16(vbest) as i32
    }
}

// ---------------------------------------------------------------------------
// SSE4.1 (x86_64)
// ---------------------------------------------------------------------------

#[cfg(all(feature = "simd", target_arch = "x86_64"))]
mod sse41 {
    use super::*;
    use core::arch::x86_64::*;

    const LANES: usize = 8; // __m128i with i16

    /// SSE4.1-accelerated SW score-only using Farrar striped vectorization.
    ///
    /// # Safety
    ///
    /// Requires SSE4.1 (checked by caller via `is_x86_feature_detected!`).
    #[target_feature(enable = "sse4.1")]
    pub(super) unsafe fn sw_score_sse41(
        query: &[u8],
        target: &[u8],
        scoring: &ScoringScheme,
    ) -> i32 {
        let qlen = query.len();
        let n_stripes = (qlen + LANES - 1) / LANES;
        let padded_len = n_stripes * LANES;

        let gap_open = scoring.gap_open() as i16;
        let gap_extend = scoring.gap_extend() as i16;

        let profile = build_striped_profile(query, scoring, n_stripes, LANES);

        let vzero = _mm_setzero_si128();
        let vgap_open = _mm_set1_epi16(gap_open);
        let vgap_extend = _mm_set1_epi16(gap_extend);
        let vneg_inf = _mm_set1_epi16(i16::MIN / 2);

        let mut vh: Vec<__m128i> = vec![vzero; n_stripes];
        let mut ve: Vec<__m128i> = vec![vneg_inf; n_stripes];
        let mut vbest = vzero;

        for &tb in target {
            let prof_row = &profile[tb as usize * padded_len..];
            let vh_old: Vec<__m128i> = vh.clone();

            // --- Main pass ---
            for s in 0..n_stripes {
                let vp = _mm_loadu_si128(prof_row.as_ptr().add(s * LANES) as *const __m128i);

                let vdiag = if s == 0 {
                    // _mm_slli_si128(x, 2) shifts left by 2 bytes (1 i16 lane),
                    // inserting 0 at the low lane — same as "shift right by 1 element"
                    _mm_slli_si128(vh_old[n_stripes - 1], 2)
                } else {
                    vh_old[s - 1]
                };

                let ve_new = _mm_max_epi16(
                    _mm_adds_epi16(vh_old[s], vgap_open),
                    _mm_adds_epi16(ve[s], vgap_extend),
                );
                ve[s] = ve_new;

                let vh_new = _mm_max_epi16(vzero, _mm_max_epi16(
                    _mm_adds_epi16(vdiag, vp),
                    ve_new,
                ));
                vh[s] = vh_new;
                vbest = _mm_max_epi16(vbest, vh_new);
            }

            // --- F-loop ---
            loop {
                let mut vf = vneg_inf;
                let mut any_update = false;

                for s in 0..n_stripes {
                    vf = _mm_max_epi16(
                        _mm_adds_epi16(vh[s], vgap_open),
                        _mm_adds_epi16(vf, vgap_extend),
                    );

                    let vh_cur = vh[s];
                    let vh_new = _mm_max_epi16(vh_cur, vf);

                    let cmp = _mm_cmpeq_epi16(vh_new, vh_cur);
                    if _mm_movemask_epi8(cmp) != 0xFFFF {
                        any_update = true;
                        vh[s] = vh_new;
                        vbest = _mm_max_epi16(vbest, vh_new);
                    }

                    vf = vf;
                }

                if !any_update {
                    break;
                }
            }
        }

        // Horizontal max reduction
        let mut buf = [0i16; 8];
        _mm_storeu_si128(buf.as_mut_ptr() as *mut __m128i, vbest);
        let mut max_val = i16::MIN;
        for &v in &buf {
            if v > max_val {
                max_val = v;
            }
        }
        max_val as i32
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scoring::{ScoringMatrix, ScoringScheme};

    fn dna_scheme() -> ScoringScheme {
        ScoringScheme::Simple(ScoringMatrix::dna_default())
    }

    #[test]
    fn scalar_identical() {
        let score = sw_score_scalar(b"ACGT", b"ACGT", &dna_scheme());
        assert_eq!(score, 8); // 4 matches × 2
    }

    #[test]
    fn scalar_with_gap() {
        // query: AC-GT vs target: ACAGT — should prefer match with gap
        let score = sw_score_scalar(b"ACGT", b"ACAGT", &dna_scheme());
        assert!(score > 0);
    }

    #[test]
    fn scalar_all_mismatch() {
        let score = sw_score_scalar(b"AAAA", b"CCCC", &dna_scheme());
        assert_eq!(score, 0);
    }

    #[test]
    fn simd_sw_identical() {
        let score = sw_simd_score(b"ACGT", b"ACGT", &dna_scheme()).unwrap();
        assert_eq!(score, 8);
    }

    #[test]
    fn simd_sw_mismatch() {
        let score = sw_simd_score(b"AAACGTAAA", b"TTTCGTTTT", &dna_scheme()).unwrap();
        assert!(score > 0);
    }

    #[test]
    fn simd_sw_all_mismatch() {
        let score = sw_simd_score(b"AAAA", b"CCCC", &dna_scheme()).unwrap();
        assert_eq!(score, 0);
    }

    #[test]
    fn simd_sw_empty_error() {
        assert!(sw_simd_score(b"", b"ACGT", &dna_scheme()).is_err());
        assert!(sw_simd_score(b"ACGT", b"", &dna_scheme()).is_err());
    }

    #[test]
    fn simd_score_matches_scalar() {
        let scoring = dna_scheme();
        let pairs: &[(&[u8], &[u8])] = &[
            (b"ACGTACGT", b"ACGTACGT"),
            (b"AAACGTAAA", b"TTTCGTTTT"),
            (b"ACGT", b"ACTT"),
            (b"GGGGGGGG", b"CCCCCCCC"),
            (b"ACGTACGTACGTACGT", b"ACGTACTTACGTACGT"),
            (b"A", b"A"),
            (b"ACGT", b"ACAGT"),
        ];
        for (q, t) in pairs {
            let simd_score = sw_simd_score(q, t, &scoring).unwrap();
            let scalar_score = sw_score_scalar(q, t, &scoring);
            assert_eq!(
                simd_score, scalar_score,
                "mismatch for q={}, t={}: simd={}, scalar={}",
                String::from_utf8_lossy(q),
                String::from_utf8_lossy(t),
                simd_score,
                scalar_score
            );
        }
    }

    #[test]
    fn simd_sw_long_sequence() {
        let q = b"ACGTACGTACGTACGTACGTACGT";
        let t = b"ACGTACGTACGTACGTACGTACGT";
        let score = sw_simd_score(q, t, &dna_scheme()).unwrap();
        assert_eq!(score, q.len() as i32 * 2); // match_score = 2
    }

    #[test]
    fn simd_sw_asymmetric_lengths() {
        let scoring = dna_scheme();
        let q = b"ACGTACGTACGT"; // 12 bp
        let t = b"ACGT"; // 4 bp
        let score = sw_simd_score(q, t, &scoring).unwrap();
        let scalar = sw_score_scalar(q, t, &scoring);
        assert_eq!(score, scalar);
    }
}
