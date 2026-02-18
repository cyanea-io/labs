#include <metal_stdlib>
using namespace metal;

// GPU Smith-Waterman for protein sequences with affine gap penalties.
// One thread per sequence pair — serial DP within each thread.
// Substitution matrix (24×24) passed as buffer.
// Amino acid index: A=0,R=1,N=2,D=3,C=4,Q=5,E=6,G=7,H=8,I=9,
// L=10,K=11,M=12,F=13,P=14,S=15,T=16,W=17,Y=18,V=19,B=20,Z=21,X=22,*=23

constant int AA_MAP[128] = {
    22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,
    22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,
    22,22,22,22,22,22,22,22,22,22,23,22,22,22,22,22,
    22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,
    22, 0,20, 4, 3, 6,13, 7, 8, 9,22,11,10,12, 2,22,
    14, 5, 1,15,16,22,19,17,22,18,21,22,22,22,22,22,
    22, 0,20, 4, 3, 6,13, 7, 8, 9,22,11,10,12, 2,22,
    14, 5, 1,15,16,22,19,17,22,18,21,22,22,22,22,22
};

kernel void smith_waterman_batch(
    device const uchar* queries [[buffer(0)]],
    device const uchar* targets [[buffer(1)]],
    device const uint* query_lengths [[buffer(2)]],
    device const uint* target_lengths [[buffer(3)]],
    device const uint* query_offsets [[buffer(4)]],
    device const uint* target_offsets [[buffer(5)]],
    device const int* submat [[buffer(6)]],
    device const int& gap_open [[buffer(7)]],
    device const int& gap_extend [[buffer(8)]],
    device int* scores [[buffer(9)]],
    device uint* query_ends [[buffer(10)]],
    device uint* target_ends [[buffer(11)]],
    uint gid [[thread_position_in_grid]]
) {
    uint qlen = query_lengths[gid];
    uint tlen = target_lengths[gid];
    uint qoff = query_offsets[gid];
    uint toff = target_offsets[gid];

    int best = 0;
    uint best_qi = 0, best_ti = 0;

    // Row-based DP: prev_h[j], curr_h[j], e[j] (gap-in-query)
    // f (gap-in-target) computed per-cell
    // Max sequence length per pair: 512
    int prev_h[513];
    int curr_e[513];

    for (uint j = 0; j <= tlen; j++) {
        prev_h[j] = 0;
        curr_e[j] = 0;
    }

    for (uint i = 1; i <= qlen; i++) {
        int qi = AA_MAP[queries[qoff + i - 1] & 0x7F];
        int diag = 0;
        int f = 0;

        for (uint j = 1; j <= tlen; j++) {
            int ti = AA_MAP[targets[toff + j - 1] & 0x7F];
            int sc = submat[qi * 24 + ti];

            // E: gap in query (horizontal extension)
            int eo = prev_h[j] + gap_open + gap_extend;
            int ee = curr_e[j] + gap_extend;
            curr_e[j] = (eo > ee) ? eo : ee;

            // F: gap in target (vertical extension)
            int fo = diag + gap_open + gap_extend;
            int fe = f + gap_extend;
            f = (fo > fe) ? fo : fe;

            // H: best of match, E, F, or 0
            int old_h = prev_h[j];
            int h = old_h + sc; // using prev_h[j] from previous row as diag
            // Wait — diag should be prev_h[j-1] from previous row
            // Let me fix: we track diag = prev_h[j-1] before updating
            h = diag + sc; // diag holds prev_h[j-1]
            if (curr_e[j] > h) h = curr_e[j];
            if (f > h) h = f;
            if (h < 0) h = 0;

            diag = old_h; // save prev_h[j] as next iteration's diag
            prev_h[j] = h;

            if (h > best) {
                best = h;
                best_qi = i - 1;
                best_ti = j - 1;
            }
        }
    }

    scores[gid] = best;
    query_ends[gid] = best_qi;
    target_ends[gid] = best_ti;
}
