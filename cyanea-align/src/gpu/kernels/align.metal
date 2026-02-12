#include <metal_stdlib>
using namespace metal;

// Banded affine-gap alignment kernel for Metal.
// One thread block per sequence pair. Threads cooperate on anti-diagonal wavefront.
//
// Parameters:
//   seqs       - packed sequence data (raw bytes)
//   index      - per-pair (query_offset, query_len, target_offset, target_len) as uint4
//   scores_out - output score per pair (int)
//   end_pos    - output (end_i, end_j) per pair as uint2
//   traceback  - output traceback band (bandwidth cells per query row per pair)
//   match_score, mismatch_score, gap_open, gap_extend - scoring params
//   bandwidth  - band width
//   mode       - 0=Local, 1=Global, 2=SemiGlobal
//   num_pairs  - total pairs

// Traceback encoding:
// Bits [1:0]: H dir (0=diag, 1=up, 2=left, 3=stop)
// Bit [2]: E from open (1) vs extend (0)
// Bit [3]: F from open (1) vs extend (0)

kernel void align_pairs(
    device const uchar* seqs [[buffer(0)]],
    device const uint4* index [[buffer(1)]],
    device int* scores_out [[buffer(2)]],
    device uint2* end_pos [[buffer(3)]],
    device uchar* traceback [[buffer(4)]],
    device const int& match_score [[buffer(5)]],
    device const int& mismatch_score [[buffer(6)]],
    device const int& gap_open [[buffer(7)]],
    device const int& gap_extend [[buffer(8)]],
    device const uint& bandwidth [[buffer(9)]],
    device const uint& mode [[buffer(10)]],
    device const uint& num_pairs [[buffer(11)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= num_pairs) return;

    uint4 idx = index[gid];
    uint q_off = idx.x;
    uint q_len = idx.y;
    uint t_off = idx.z;
    uint t_len = idx.w;

    int bw = int(bandwidth);
    uint tb_offset = gid * (q_len + 1) * bandwidth;

    // DP state: we only keep two rows of the band for H, E, F
    // Since we run one thread per pair, this is a serial DP per pair.
    int best_score = 0;
    uint best_i = 0;
    uint best_j = 0;

    // Allocate band storage in thread-local arrays
    // For bandwidth <= 128, this fits in registers/local memory
    int h_prev[129]; // previous row H values across band
    int h_curr[129]; // current row H values
    int e_band[129]; // E values for current row
    int f_band[129]; // F values for current row

    // Initialize
    for (int b = 0; b < bw + 1; b++) {
        h_prev[b] = 0;
        h_curr[b] = 0;
        e_band[b] = -1000000;
        f_band[b] = -1000000;
    }

    // Global mode: initialize first row/column with gap penalties
    if (mode == 1) {
        for (int b = 0; b < bw / 2 && b <= int(t_len); b++) {
            h_prev[bw / 2 + b] = gap_open + b * gap_extend;
        }
    }

    // Fill DP matrix row by row with banding
    for (uint i = 1; i <= q_len; i++) {
        uchar q_base = seqs[q_off + i - 1];
        int band_center = int(i) * int(t_len) / max(int(q_len), 1);

        // Initialize first element for global
        int h_left = 0;
        if (mode == 1) {
            h_left = gap_open + int(i) * gap_extend;
        }

        for (int b = 0; b < bw; b++) {
            int j = band_center - bw / 2 + b;
            if (j < 0 || j >= int(t_len)) {
                h_curr[b] = 0;
                e_band[b] = -1000000;
                uint tb_idx = tb_offset + i * bandwidth + uint(b);
                if (tb_idx < tb_offset + (q_len + 1) * bandwidth)
                    traceback[tb_idx] = 3; // STOP
                continue;
            }

            uchar t_base = seqs[t_off + uint(j)];
            int sub = (q_base == t_base || (q_base ^ 32) == t_base ||
                       q_base == (t_base ^ 32)) ? match_score : mismatch_score;

            // E: gap in query (insertion)
            int e_open = h_left + gap_open;
            int e_ext = (b > 0) ? e_band[b - 1] + gap_extend : -1000000;
            e_band[b] = max(e_open, e_ext);
            uchar tb_e = (e_open >= e_ext) ? 0x04u : 0x00u;

            // F: gap in target (deletion)
            int f_open = h_prev[b] + gap_open;
            int f_ext = f_band[b] + gap_extend;
            f_band[b] = max(f_open, f_ext);
            uchar tb_f = (f_open >= f_ext) ? 0x08u : 0x00u;

            // H: diagonal
            int h_diag = h_prev[b > 0 ? b - 1 : 0] + sub;
            if (b == 0 && i > 1) {
                // Use previous row's value at adjusted position
                h_diag = h_prev[0] + sub;
            }

            int h_val = max(h_diag, max(e_band[b], f_band[b]));
            uchar tb_h;
            if (h_val == h_diag) tb_h = 0; // DIAG
            else if (h_val == f_band[b]) tb_h = 1; // UP
            else tb_h = 2; // LEFT

            // Local: clamp to 0
            if (mode == 0) {
                if (h_val < 0) {
                    h_val = 0;
                    tb_h = 3; // STOP
                }
            }

            h_curr[b] = h_val;
            h_left = h_val;

            // Store traceback
            uint tb_idx = tb_offset + i * bandwidth + uint(b);
            if (tb_idx < tb_offset + (q_len + 1) * bandwidth)
                traceback[tb_idx] = tb_h | tb_e | tb_f;

            // Track best for local
            if (mode == 0 && h_val > best_score) {
                best_score = h_val;
                best_i = i;
                best_j = uint(j) + 1;
            }
        }

        // Swap rows
        for (int b = 0; b < bw; b++) {
            h_prev[b] = h_curr[b];
        }
    }

    // Determine final score and position based on mode
    if (mode == 1) {
        // Global: score at (m, n)
        int band_center = int(t_len);
        int b = band_center - (int(q_len) * int(t_len) / max(int(q_len), 1)) + bw / 2;
        if (b >= 0 && b < bw) {
            best_score = h_prev[b];
        }
        best_i = q_len;
        best_j = t_len;
    } else if (mode == 2) {
        // Semi-global: best in last row
        best_score = -1000000;
        for (int b = 0; b < bw; b++) {
            if (h_prev[b] > best_score) {
                best_score = h_prev[b];
                best_i = q_len;
                int band_center = int(q_len) * int(t_len) / max(int(q_len), 1);
                best_j = uint(max(0, band_center - bw / 2 + b)) + 1;
            }
        }
    }

    scores_out[gid] = best_score;
    end_pos[gid] = uint2(best_i, best_j);
}
