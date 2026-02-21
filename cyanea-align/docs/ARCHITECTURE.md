# cyanea-align Architecture

## Module Map

```
cyanea-align
 +-- types.rs              AlignmentMode, CigarOp (full SAM alphabet), AlignmentResult
 +-- cigar.rs              CIGAR parsing, validation, coordinate queries, arithmetic, MD tags
 +-- scoring.rs            ScoringMatrix, SubstitutionMatrix (BLOSUM/PAM), ScoringScheme
 +-- needleman_wunsch.rs   Global alignment (Gotoh 3-matrix affine gaps)
 +-- smith_waterman.rs     Local alignment (Gotoh 3-matrix)
 +-- semi_global.rs        Semi-global alignment (free leading/trailing gaps)
 +-- batch.rs              Batch pairwise alignment dispatcher
 +-- simd.rs               Banded alignment with traceback (all modes)
 +-- simd_sw.rs            [simd] Farrar striped SIMD score-only SW
 +-- msa.rs                Progressive MSA (ClustalW-style guide tree)
 +-- poa.rs                Partial Order Alignment (Lee 2002 DAG)
 +-- minimizers.rs         (w,k)-minimizer extraction with invertible hash
 +-- seed_extend.rs        Seed chaining + banded extension
 +-- wfa.rs                [wfa] Wavefront alignment (O(ns))
 +-- xdrop.rs              X-drop/Z-drop extension alignment
 +-- spliced.rs            Intron-aware DP alignment (GT-AG/GC-AG/AT-AC)
 +-- lcsk.rs               LCSk++ sparse alignment with Fenwick tree
 +-- pair_hmm.rs           Pair HMM: forward + Viterbi (3-state, log-space)
 +-- profile_hmm.rs        Profile HMM: Plan 7, Viterbi/Forward/Backward, E-value
 +-- gpu/
      +-- mod.rs           GPU dispatch, backend selection, CPU fallback
      +-- common.rs        Sequence encoding, traceback reconstruction, pair partitioning
      +-- metal_align.rs   [metal] Metal banded affine-gap aligner
      +-- cuda_align.rs    [cuda] CUDA banded affine-gap aligner
      +-- kernels/
           +-- align.metal  MSL compute shader
           +-- align.cu     CUDA C kernel
```

## Design Decisions

### Private push_cigar Helpers

Each alignment algorithm (NW, SW, semi-global, seed-extend, xdrop, spliced) has its own private `push_cigar` helper function that merges consecutive identical operations during traceback. This avoids a shared mutable abstraction and keeps each module self-contained. The `cigar.rs` module provides the public `collapse_cigar` for post-hoc merging.

### Scoring Matrices

`ScoringScheme` is an enum wrapping either `ScoringMatrix` (simple match/mismatch/gap) or `SubstitutionMatrix` (full NxN lookup). All alignment functions accept `&ScoringScheme` and dispatch internally.

Built-in substitution matrices:
- **BLOSUM**: 30, 45, 62, 80, 90 -- amino acid substitution based on observed blocks
- **PAM**: 40, 120, 200, 250 -- point accepted mutation evolutionary distance

Gap penalties are affine (open + extend) stored in the matrix itself, enabling different penalty schemes per matrix family.

### SIMD: Farrar Striped Layout

`simd_sw.rs` implements the Farrar (2007) striped approach for score-only Smith-Waterman:

- **NEON** (aarch64 / Apple Silicon): 8 x i16 lanes, always available
- **AVX2** (x86_64): 16 x i16 lanes, runtime feature detection
- **SSE4.1** (x86_64): 8 x i16 lanes, fallback from AVX2

The query profile is pre-computed and striped across SIMD lanes. The inner loop processes one target base at a time, updating all query positions in parallel. Score-only mode uses O(n) memory (no traceback matrix).

The banded alignment in `simd.rs` is a separate implementation that does full traceback but restricts computation to a diagonal band.

### GPU Dispatch Architecture

The GPU alignment pipeline:

1. **Pair partitioning**: Pairs where both sequences fit within `max_bandwidth` diagonals go to GPU; others fall back to CPU banded alignment.
2. **Sequence packing**: Query and target sequences are packed into contiguous byte buffers with offsets for efficient GPU memory transfer.
3. **Kernel execution**: One GPU thread per sequence pair. Each thread runs serial banded affine-gap DP (Gotoh 3-matrix) with traceback into shared memory.
4. **Traceback reconstruction**: Host-side code reads the traceback buffer and reconstructs `AlignmentResult` with full CIGAR.

Backend-specific details:
- **Metal**: MSL compute shader compiled at runtime, `StorageModeShared` buffers for CPU/GPU memory sharing (Apple Silicon unified memory).
- **CUDA**: NVRTC-compiled kernel via cudarc driver API. Explicit host-device memory copies.
- **CPU fallback**: Always available. Also used when batch size < `min_pairs_for_gpu` (default 64).

### WFA: Wavefront Algorithm

Implements the Marco-Sola (2021) gap-affine wavefront algorithm:

- Operates in **penalty space**, not score space. The algorithm expands wavefronts for increasing penalty values.
- Three wavefront components: M (match), I (insertion), D (deletion) for affine gaps.
- At each penalty step, wavefronts are extended greedily along matching diagonals (`extend` phase) before expanding to new penalties (`next` phase).
- Time complexity: O(ns) where s is the optimal alignment penalty. For highly similar sequences (low s), this is dramatically faster than O(nm) Needleman-Wunsch.
- The result is converted to standard `AlignmentResult` with CIGAR via traceback through stored wavefronts.

### POA: Lee 2002 DAG Alignment

Partial Order Alignment aligns a new sequence against a directed acyclic graph (DAG) of previously aligned sequences:

1. **Topological sort** the graph nodes.
2. **DP alignment** against the linearized DAG: for each node, the DP row considers all predecessor edges (not just the single diagonal predecessor in linear DP).
3. **Traceback** through the DP matrix to find the optimal path.
4. **Integration**: merge the new sequence into the graph by adding new nodes for insertions and reusing existing nodes for matches.

Consensus is extracted via heaviest-path traversal: at each node, follow the outgoing edge with the highest total weight (sum of sequence memberships).

### Pair HMM: Three-State Log-Space Model

The pair HMM has three states:
- **M** (Match/Mismatch): both sequences advance
- **X** (Insert in seq_a): seq_a advances, seq_b stays
- **Y** (Insert in seq_b): seq_b advances, seq_a stays

All computation is in log-space to avoid floating-point underflow. The `log_sum_exp` helper computes `ln(exp(a) + exp(b))` stably.

### Profile HMM: Plan 7 Architecture

The Plan 7 architecture has 7 transitions per match position:
- M_k -> M_{k+1}, M_k -> I_k, M_k -> D_{k+1}
- I_k -> M_{k+1}, I_k -> I_k
- D_k -> M_{k+1}, D_k -> D_{k+1}

Plus entry/exit probabilities for local alignment mode (uniform Begin->M_k and M_k->End).

Construction from MSA:
1. Classify columns: columns with gap fraction > threshold become insert states.
2. Count emissions and transitions with pseudocounts.
3. Normalize to probabilities and convert to log-space.

E-value calibration: run Viterbi on random sequences, fit a Gumbel distribution to the score distribution via method of moments, then compute E-value as `K * db_size * exp(-lambda * score)`.
