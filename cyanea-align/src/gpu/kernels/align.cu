// Banded affine-gap alignment kernel for CUDA.
// One thread block per sequence pair (serial within each pair).
// Mirrors the Metal kernel logic for NVIDIA GPUs.

extern "C" {

// Traceback encoding:
// Bits [1:0]: H dir (0=diag, 1=up, 2=left, 3=stop)
// Bit [2]: E from open vs extend
// Bit [3]: F from open vs extend

__global__ void align_pairs(
    const unsigned char* seqs,
    const unsigned int* index, // 4 uints per pair: q_off, q_len, t_off, t_len
    int* scores_out,
    unsigned int* end_pos,     // 2 uints per pair: end_i, end_j
    unsigned char* traceback,
    int match_score,
    int mismatch_score,
    int gap_open,
    int gap_extend,
    unsigned int bandwidth,
    unsigned int mode,         // 0=Local, 1=Global, 2=SemiGlobal
    unsigned int num_pairs
) {
    unsigned int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= num_pairs) return;

    unsigned int q_off = index[gid * 4 + 0];
    unsigned int q_len = index[gid * 4 + 1];
    unsigned int t_off = index[gid * 4 + 2];
    unsigned int t_len = index[gid * 4 + 3];

    int bw = (int)bandwidth;
    unsigned int tb_offset = gid * (q_len + 1) * bandwidth;

    int best_score = 0;
    unsigned int best_i = 0;
    unsigned int best_j = 0;

    // Thread-local band storage
    int h_prev[129];
    int h_curr[129];
    int e_band[129];
    int f_band[129];

    for (int b = 0; b < bw + 1; b++) {
        h_prev[b] = 0;
        h_curr[b] = 0;
        e_band[b] = -1000000;
        f_band[b] = -1000000;
    }

    if (mode == 1) {
        for (int b = 0; b < bw / 2 && b <= (int)t_len; b++) {
            h_prev[bw / 2 + b] = gap_open + b * gap_extend;
        }
    }

    for (unsigned int i = 1; i <= q_len; i++) {
        unsigned char q_base = seqs[q_off + i - 1];
        int band_center = (int)i * (int)t_len / max((int)q_len, 1);
        int h_left = 0;
        if (mode == 1) h_left = gap_open + (int)i * gap_extend;

        for (int b = 0; b < bw; b++) {
            int j = band_center - bw / 2 + b;
            if (j < 0 || j >= (int)t_len) {
                h_curr[b] = 0;
                e_band[b] = -1000000;
                unsigned int tb_idx = tb_offset + i * bandwidth + (unsigned int)b;
                if (tb_idx < tb_offset + (q_len + 1) * bandwidth)
                    traceback[tb_idx] = 3;
                continue;
            }

            unsigned char t_base = seqs[t_off + (unsigned int)j];
            int sub = (q_base == t_base || (q_base ^ 32) == t_base ||
                       q_base == (t_base ^ 32)) ? match_score : mismatch_score;

            int e_open = h_left + gap_open;
            int e_ext = (b > 0) ? e_band[b - 1] + gap_extend : -1000000;
            e_band[b] = (e_open > e_ext) ? e_open : e_ext;

            int f_open = h_prev[b] + gap_open;
            int f_ext = f_band[b] + gap_extend;
            f_band[b] = (f_open > f_ext) ? f_open : f_ext;

            int h_diag = h_prev[b > 0 ? b - 1 : 0] + sub;
            if (b == 0 && i > 1) h_diag = h_prev[0] + sub;

            int h_val = h_diag;
            unsigned char tb_h = 0;
            if (e_band[b] > h_val) { h_val = e_band[b]; tb_h = 2; }
            if (f_band[b] > h_val) { h_val = f_band[b]; tb_h = 1; }

            if (mode == 0 && h_val < 0) { h_val = 0; tb_h = 3; }

            h_curr[b] = h_val;
            h_left = h_val;

            unsigned int tb_idx = tb_offset + i * bandwidth + (unsigned int)b;
            if (tb_idx < tb_offset + (q_len + 1) * bandwidth) {
                unsigned char tb_e = (e_open >= e_ext) ? 0x04 : 0x00;
                unsigned char tb_f = (f_open >= f_ext) ? 0x08 : 0x00;
                traceback[tb_idx] = tb_h | tb_e | tb_f;
            }

            if (mode == 0 && h_val > best_score) {
                best_score = h_val;
                best_i = i;
                best_j = (unsigned int)j + 1;
            }
        }

        for (int b = 0; b < bw; b++) h_prev[b] = h_curr[b];
    }

    if (mode == 1) {
        int band_center = (int)t_len;
        int b = band_center - ((int)q_len * (int)t_len / max((int)q_len, 1)) + bw / 2;
        if (b >= 0 && b < bw) best_score = h_prev[b];
        best_i = q_len;
        best_j = t_len;
    } else if (mode == 2) {
        best_score = -1000000;
        for (int b = 0; b < bw; b++) {
            if (h_prev[b] > best_score) {
                best_score = h_prev[b];
                best_i = q_len;
                int band_center = (int)q_len * (int)t_len / max((int)q_len, 1);
                best_j = (unsigned int)(band_center - bw / 2 + b) + 1;
                if (best_j > t_len) best_j = t_len;
            }
        }
    }

    scores_out[gid] = best_score;
    end_pos[gid * 2 + 0] = best_i;
    end_pos[gid * 2 + 1] = best_j;
}

} // extern "C"
