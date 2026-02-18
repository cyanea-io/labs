//! GPU-accelerated Smith-Waterman alignment for protein sequences.
//!
//! Provides batch protein SW alignment with affine gap penalties and a
//! user-supplied substitution matrix (e.g., BLOSUM62). Supports GPU
//! dispatch (Metal/CUDA) with automatic CPU fallback.

use cyanea_core::{CyaneaError, Result};

/// Result of a single Smith-Waterman alignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SwResult {
    /// Best alignment score.
    pub score: i32,
    /// Query position (0-based) where the best alignment ends.
    pub query_end: usize,
    /// Target position (0-based) where the best alignment ends.
    pub target_end: usize,
}

/// Amino acid to index mapping (24 amino acids + ambiguous + stop).
/// A=0, R=1, N=2, D=3, C=4, Q=5, E=6, G=7, H=8, I=9,
/// L=10, K=11, M=12, F=13, P=14, S=15, T=16, W=17, Y=18, V=19,
/// B=20, Z=21, X=22, *=23
fn aa_to_index(aa: u8) -> usize {
    match aa {
        b'A' | b'a' => 0,
        b'R' | b'r' => 1,
        b'N' | b'n' => 2,
        b'D' | b'd' => 3,
        b'C' | b'c' => 4,
        b'Q' | b'q' => 5,
        b'E' | b'e' => 6,
        b'G' | b'g' => 7,
        b'H' | b'h' => 8,
        b'I' | b'i' => 9,
        b'L' | b'l' => 10,
        b'K' | b'k' => 11,
        b'M' | b'm' => 12,
        b'F' | b'f' => 13,
        b'P' | b'p' => 14,
        b'S' | b's' => 15,
        b'T' | b't' => 16,
        b'W' | b'w' => 17,
        b'Y' | b'y' => 18,
        b'V' | b'v' => 19,
        b'B' | b'b' => 20,
        b'Z' | b'z' => 21,
        b'X' | b'x' => 22,
        b'*' => 23,
        _ => 22, // unknown → X
    }
}

/// Runs batch Smith-Waterman protein alignment.
///
/// Each pair in `pairs` is a `(query, target)` protein sequence. The
/// `substitution_matrix` is a 24×24 flattened row-major matrix (e.g.,
/// BLOSUM62 extended to 24 AA). Gap penalties follow the convention
/// `gap_open + gap_extend` for opening a gap, `gap_extend` for each
/// additional gap position (both should be negative).
///
/// # Errors
///
/// Returns an error if `pairs` is empty or `substitution_matrix` length
/// is not 576 (24×24).
pub fn gpu_smith_waterman_batch(
    pairs: &[(&[u8], &[u8])],
    substitution_matrix: &[i32],
    gap_open: i32,
    gap_extend: i32,
) -> Result<Vec<SwResult>> {
    validate_sw_inputs(pairs, substitution_matrix)?;

    #[cfg(feature = "metal")]
    {
        if let Ok(results) = gpu_sw_metal(pairs, substitution_matrix, gap_open, gap_extend) {
            return Ok(results);
        }
    }

    gpu_sw_cpu(pairs, substitution_matrix, gap_open, gap_extend)
}

/// CPU fallback for batch Smith-Waterman.
pub fn gpu_sw_cpu(
    pairs: &[(&[u8], &[u8])],
    substitution_matrix: &[i32],
    gap_open: i32,
    gap_extend: i32,
) -> Result<Vec<SwResult>> {
    validate_sw_inputs(pairs, substitution_matrix)?;

    let mut results = Vec::with_capacity(pairs.len());
    for &(query, target) in pairs {
        results.push(sw_affine(query, target, substitution_matrix, gap_open, gap_extend));
    }
    Ok(results)
}

/// Single-pair SW with affine gaps.
fn sw_affine(
    query: &[u8],
    target: &[u8],
    submat: &[i32],
    gap_open: i32,
    gap_extend: i32,
) -> SwResult {
    let m = query.len();
    let n = target.len();

    if m == 0 || n == 0 {
        return SwResult {
            score: 0,
            query_end: 0,
            target_end: 0,
        };
    }

    // H[j]: best score ending at column j
    // E[j]: best score ending with gap in query (horizontal)
    // F: best score ending with gap in target (vertical), per-cell
    let mut h_prev = vec![0i32; n + 1];
    let mut e = vec![i32::MIN / 2; n + 1];

    let mut best_score = 0i32;
    let mut best_qi = 0usize;
    let mut best_ti = 0usize;

    for i in 1..=m {
        let qi = aa_to_index(query[i - 1]);
        let mut diag = 0i32; // H[i-1][j-1]
        let mut f = i32::MIN / 2;

        for j in 1..=n {
            let ti = aa_to_index(target[j - 1]);
            let sc = submat[qi * 24 + ti];

            // E: gap in query (from left)
            let e_open = h_prev[j - 1] + gap_open + gap_extend;
            let e_ext = e[j - 1] + gap_extend;
            e[j] = e_open.max(e_ext);

            // F: gap in target (from above)
            let f_open = h_prev[j] + gap_open + gap_extend;
            let f_ext = f + gap_extend;
            f = f_open.max(f_ext);

            let old_h = h_prev[j];
            let mut h = diag + sc;
            h = h.max(e[j]).max(f).max(0);

            diag = old_h;
            h_prev[j] = h;

            if h > best_score {
                best_score = h;
                best_qi = i - 1;
                best_ti = j - 1;
            }
        }
    }

    SwResult {
        score: best_score,
        query_end: best_qi,
        target_end: best_ti,
    }
}

fn validate_sw_inputs(pairs: &[(&[u8], &[u8])], submat: &[i32]) -> Result<()> {
    if pairs.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "smith_waterman: pairs must not be empty".to_string(),
        ));
    }
    if submat.len() != 576 {
        return Err(CyaneaError::InvalidInput(format!(
            "smith_waterman: substitution matrix must have 576 (24×24) elements, got {}",
            submat.len()
        )));
    }
    Ok(())
}

// ── Metal GPU dispatch ────────────────────────────────────────────

#[cfg(feature = "metal")]
fn gpu_sw_metal(
    pairs: &[(&[u8], &[u8])],
    substitution_matrix: &[i32],
    gap_open: i32,
    gap_extend: i32,
) -> Result<Vec<SwResult>> {
    use metal_rs::{CompileOptions, Device, MTLResourceOptions, MTLSize};
    use std::ffi::c_void;

    let device = Device::system_default()
        .ok_or_else(|| CyaneaError::Other("no Metal device available".into()))?;
    let queue = device.new_command_queue();

    let shader_src = include_str!("shaders/sw.metal");
    let opts = CompileOptions::new();
    let lib = device
        .new_library_with_source(shader_src, &opts)
        .map_err(|e| CyaneaError::Other(format!("Metal SW shader compile: {e}")))?;
    let func = lib
        .get_function("smith_waterman_batch", None)
        .map_err(|e| CyaneaError::Other(format!("Metal SW function: {e}")))?;
    let pipeline = device
        .new_compute_pipeline_state_with_function(&func)
        .map_err(|e| CyaneaError::Other(format!("Metal SW pipeline: {e}")))?;

    let n_pairs = pairs.len();

    // Pack sequences into flat arrays with offset/length metadata
    let mut all_queries = Vec::new();
    let mut all_targets = Vec::new();
    let mut query_lengths = Vec::with_capacity(n_pairs);
    let mut target_lengths = Vec::with_capacity(n_pairs);
    let mut query_offsets = Vec::with_capacity(n_pairs);
    let mut target_offsets = Vec::with_capacity(n_pairs);

    for &(q, t) in pairs {
        query_offsets.push(all_queries.len() as u32);
        target_offsets.push(all_targets.len() as u32);
        query_lengths.push(q.len() as u32);
        target_lengths.push(t.len() as u32);
        all_queries.extend_from_slice(q);
        all_targets.extend_from_slice(t);
    }

    // Ensure non-empty buffers
    if all_queries.is_empty() {
        all_queries.push(0);
    }
    if all_targets.is_empty() {
        all_targets.push(0);
    }

    let make_buf = |data: &[u8]| {
        device.new_buffer_with_data(
            data.as_ptr() as *const c_void,
            data.len() as u64,
            MTLResourceOptions::StorageModeShared,
        )
    };
    let make_u32_buf = |data: &[u32]| {
        device.new_buffer_with_data(
            data.as_ptr() as *const c_void,
            (data.len() * 4) as u64,
            MTLResourceOptions::StorageModeShared,
        )
    };
    let make_i32_buf = |data: &[i32]| {
        device.new_buffer_with_data(
            data.as_ptr() as *const c_void,
            (data.len() * 4) as u64,
            MTLResourceOptions::StorageModeShared,
        )
    };
    let make_scalar_i32 = |val: i32| {
        device.new_buffer_with_data(
            &val as *const i32 as *const c_void,
            4,
            MTLResourceOptions::StorageModeShared,
        )
    };

    let q_buf = make_buf(&all_queries);
    let t_buf = make_buf(&all_targets);
    let ql_buf = make_u32_buf(&query_lengths);
    let tl_buf = make_u32_buf(&target_lengths);
    let qo_buf = make_u32_buf(&query_offsets);
    let to_buf = make_u32_buf(&target_offsets);
    let sm_buf = make_i32_buf(substitution_matrix);
    let go_buf = make_scalar_i32(gap_open);
    let ge_buf = make_scalar_i32(gap_extend);
    let scores_buf = device.new_buffer(
        (n_pairs * 4) as u64,
        MTLResourceOptions::StorageModeShared,
    );
    let qe_buf = device.new_buffer(
        (n_pairs * 4) as u64,
        MTLResourceOptions::StorageModeShared,
    );
    let te_buf = device.new_buffer(
        (n_pairs * 4) as u64,
        MTLResourceOptions::StorageModeShared,
    );

    let cmd = queue.new_command_buffer();
    let enc = cmd.new_compute_command_encoder();
    enc.set_compute_pipeline_state(&pipeline);
    enc.set_buffer(0, Some(&q_buf), 0);
    enc.set_buffer(1, Some(&t_buf), 0);
    enc.set_buffer(2, Some(&ql_buf), 0);
    enc.set_buffer(3, Some(&tl_buf), 0);
    enc.set_buffer(4, Some(&qo_buf), 0);
    enc.set_buffer(5, Some(&to_buf), 0);
    enc.set_buffer(6, Some(&sm_buf), 0);
    enc.set_buffer(7, Some(&go_buf), 0);
    enc.set_buffer(8, Some(&ge_buf), 0);
    enc.set_buffer(9, Some(&scores_buf), 0);
    enc.set_buffer(10, Some(&qe_buf), 0);
    enc.set_buffer(11, Some(&te_buf), 0);

    let grid = MTLSize::new(n_pairs as u64, 1, 1);
    let tg = MTLSize::new((n_pairs as u64).min(256).max(1), 1, 1);
    enc.dispatch_threads(grid, tg);
    enc.end_encoding();
    cmd.commit();
    cmd.wait_until_completed();

    // Read results
    let scores_ptr = scores_buf.contents() as *const i32;
    let qe_ptr = qe_buf.contents() as *const u32;
    let te_ptr = te_buf.contents() as *const u32;

    let scores = unsafe { std::slice::from_raw_parts(scores_ptr, n_pairs) };
    let query_ends = unsafe { std::slice::from_raw_parts(qe_ptr, n_pairs) };
    let target_ends = unsafe { std::slice::from_raw_parts(te_ptr, n_pairs) };

    let mut results = Vec::with_capacity(n_pairs);
    for i in 0..n_pairs {
        results.push(SwResult {
            score: scores[i],
            query_end: query_ends[i] as usize,
            target_end: target_ends[i] as usize,
        });
    }
    Ok(results)
}

/// Standard BLOSUM62 matrix extended to 24×24 (includes B, Z, X, *).
/// Row/column order: A R N D C Q E G H I L K M F P S T W Y V B Z X *
#[rustfmt::skip]
pub const BLOSUM62_24: [i32; 576] = [
    // A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sw_identical_proteins() {
        let seq = b"ACDEFGHIKLMNPQRSTVWY";
        let results = gpu_sw_cpu(
            &[(seq.as_slice(), seq.as_slice())],
            &BLOSUM62_24,
            -11,
            -1,
        )
        .unwrap();
        assert_eq!(results.len(), 1);
        // Perfect self-alignment should have positive score
        assert!(results[0].score > 0);
        // Should end at last position
        assert_eq!(results[0].query_end, seq.len() - 1);
        assert_eq!(results[0].target_end, seq.len() - 1);
    }

    #[test]
    fn sw_blosum62_known() {
        // Align "HEAGAWGHEE" vs "PAWHEAE" — well-known example
        let query = b"HEAGAWGHEE";
        let target = b"PAWHEAE";
        let results = gpu_sw_cpu(
            &[(query.as_slice(), target.as_slice())],
            &BLOSUM62_24,
            -11,
            -1,
        )
        .unwrap();
        // Score should be positive (there is a good local alignment)
        assert!(results[0].score > 0);
    }

    #[test]
    fn sw_gap_penalty() {
        // Two identical sequences with a gap inserted
        let query = b"AAAAA";
        let target = b"AAGAAA"; // G inserted
        let results = gpu_sw_cpu(
            &[(query.as_slice(), target.as_slice())],
            &BLOSUM62_24,
            -11,
            -1,
        )
        .unwrap();
        // Score with gap should be less than perfect match
        let perfect = gpu_sw_cpu(
            &[(query.as_slice(), query.as_slice())],
            &BLOSUM62_24,
            -11,
            -1,
        )
        .unwrap();
        assert!(results[0].score < perfect[0].score);
    }

    #[test]
    fn sw_batch_multiple_pairs() {
        let pairs = vec![
            (b"ACDE".as_slice(), b"ACDE".as_slice()),
            (b"WWWW".as_slice(), b"AAAA".as_slice()),
            (b"MFIL".as_slice(), b"MFIL".as_slice()),
        ];
        let results = gpu_sw_cpu(&pairs, &BLOSUM62_24, -11, -1).unwrap();
        assert_eq!(results.len(), 3);
        // Self-alignments should score higher than mismatches
        assert!(results[0].score > results[1].score);
        assert!(results[2].score > results[1].score);
    }

    #[test]
    fn sw_empty_error() {
        assert!(gpu_sw_cpu(&[], &BLOSUM62_24, -11, -1).is_err());
        let seq = b"ACDE";
        assert!(gpu_sw_cpu(&[(seq.as_slice(), seq.as_slice())], &[0; 100], -11, -1).is_err());
    }
}
