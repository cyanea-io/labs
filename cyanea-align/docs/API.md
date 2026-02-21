# cyanea-align API Reference

Pairwise sequence alignment with affine gap penalties. Supports DNA, RNA, and protein alignment with standard substitution matrices.

## Public API

### Alignment modes and results (`types.rs`)

| Type | Description |
|------|-------------|
| `AlignmentMode` | Enum: `Local`, `Global`, `SemiGlobal` |
| `CigarOp` | Enum: `Match(usize)`, `Mismatch(usize)`, `Insert(usize)`, `Delete(usize)`, `SoftClip(usize)`, `HardClip(usize)`, `RefSkip(usize)`, `Pad(usize)`, `SeqMatch(usize)`, `SeqMismatch(usize)` |
| `AlignmentResult` | Full result: `score`, `query_start/end`, `target_start/end`, `cigar`, `aligned_query/target` |

**AlignmentResult methods:**

| Method | Description |
|--------|-------------|
| `cigar_string() -> String` | Compact CIGAR string (e.g., `3=1X4=`) |
| `identity() -> f64` | Fraction of matching positions |
| `matches() -> usize` | Number of matches |
| `mismatches() -> usize` | Number of mismatches |
| `gaps() -> usize` | Number of gap positions |
| `length() -> usize` | Total alignment length |

### CIGAR string utilities (`cigar.rs`)

Full SAM CIGAR alphabet (9 ops: M, I, D, N, S, H, P, =, X) with parsing, validation, coordinate queries, statistics, arithmetic, and MD tag generation.

**Parsing and formatting:**

| Function | Description |
|----------|-------------|
| `parse_cigar(s) -> Result<Vec<CigarOp>>` | Parse SAM CIGAR string into operations |
| `cigar_string(ops) -> String` | Format operations as compact CIGAR string |
| `validate_cigar(ops) -> Result<()>` | Validate CIGAR operation sequence |

**Coordinate queries:**

| Function | Description |
|----------|-------------|
| `reference_consumed(ops) -> usize` | Number of reference bases consumed (M+D+N+=+X) |
| `query_consumed(ops) -> usize` | Number of query bases consumed (M+I+S+=+X) |
| `alignment_length(ops) -> usize` | Alignment length (M+I+D+N+S+=+X) |
| `query_to_ref(ops, qpos) -> Option<usize>` | Map query position to reference position |
| `ref_to_query(ops, rpos) -> Option<usize>` | Map reference position to query position |

**Statistics:**

| Function | Description |
|----------|-------------|
| `cigar_stats(ops) -> CigarStats` | Match/mismatch/insert/delete/clip counts |
| `identity_from_cigar(ops) -> f64` | Identity fraction from CIGAR operations |

**Arithmetic:**

| Function | Description |
|----------|-------------|
| `merge_cigars(a, b) -> Vec<CigarOp>` | Concatenate two CIGAR strings with run merging |
| `reverse_cigar(ops) -> Vec<CigarOp>` | Reverse operation order |
| `split_cigar(ops, pos) -> (Vec<CigarOp>, Vec<CigarOp>)` | Split at a reference position |
| `collapse_cigar(ops) -> Vec<CigarOp>` | Merge adjacent same-type operations |

**Conversion:**

| Function | Description |
|----------|-------------|
| `alignment_to_cigar(query, target) -> Vec<CigarOp>` | Infer CIGAR from aligned sequences |
| `generate_md_tag(ops, ref_seq) -> String` | Generate SAM MD tag from CIGAR and reference |

### Scoring (`scoring.rs`)

| Type | Description |
|------|-------------|
| `ScoringMatrix` | Simple match/mismatch/gap scores |
| `SubstitutionMatrix` | Full substitution matrices (20x20 for protein) |
| `ScoringScheme` | Enum: `Simple(ScoringMatrix)`, `Substitution(SubstitutionMatrix)` |

**Built-in matrices:**

| Constructor | Description |
|-------------|-------------|
| `ScoringMatrix::dna_default()` | +2/-1/-5/-2 (match/mismatch/gap_open/gap_extend) |
| `SubstitutionMatrix::blosum30()` | BLOSUM30 (gap_open: -14, gap_extend: -4) |
| `SubstitutionMatrix::blosum45()` | BLOSUM45 |
| `SubstitutionMatrix::blosum62()` | BLOSUM62 |
| `SubstitutionMatrix::blosum80()` | BLOSUM80 |
| `SubstitutionMatrix::blosum90()` | BLOSUM90 |
| `SubstitutionMatrix::pam40()` | PAM40 (gap_open: -10, gap_extend: -2) |
| `SubstitutionMatrix::pam120()` | PAM120 (gap_open: -11, gap_extend: -1) |
| `SubstitutionMatrix::pam200()` | PAM200 (gap_open: -11, gap_extend: -1) |
| `SubstitutionMatrix::pam250()` | PAM250 |

### Alignment functions (`needleman_wunsch.rs`, `smith_waterman.rs`, `semi_global.rs`)

| Function | Description |
|----------|-------------|
| `align(query, target, mode, scoring) -> Result<AlignmentResult>` | Dispatcher -- routes to NW, SW, or semi-global |
| `needleman_wunsch(query, target, scoring) -> Result<AlignmentResult>` | Global alignment (Gotoh 3-matrix) |
| `smith_waterman(query, target, scoring) -> Result<AlignmentResult>` | Local alignment (Gotoh 3-matrix) |
| `semi_global(query, target, scoring) -> Result<AlignmentResult>` | Semi-global (free leading/trailing gaps) |

### Batch alignment (`batch.rs`)

| Function | Description |
|----------|-------------|
| `align_batch(pairs, mode, scoring) -> Result<Vec<AlignmentResult>>` | Batch pairwise alignment (all modes) |

### Banded alignment (`simd.rs`)

| Function | Description |
|----------|-------------|
| `banded_nw(query, target, scoring, bandwidth) -> Result<AlignmentResult>` | Banded global alignment |
| `banded_sw(query, target, scoring, bandwidth) -> Result<AlignmentResult>` | Banded local alignment |
| `banded_semi_global(query, target, scoring, bandwidth) -> Result<AlignmentResult>` | Banded semi-global |
| `banded_score_only(query, target, scoring, bandwidth, mode) -> Result<i32>` | Score-only banded |

### SIMD Smith-Waterman (`simd_sw.rs`, `simd` feature)

| Function | Description |
|----------|-------------|
| `sw_simd_score(query, target, scoring) -> Result<i32>` | SIMD-accelerated SW score-only (NEON/AVX2/SSE4.1 + scalar fallback) |

### Multiple sequence alignment (`msa.rs`)

| Type/Function | Description |
|---------------|-------------|
| `MsaResult` | Aligned sequences, column count, conservation score |
| `progressive_msa(sequences, scoring) -> Result<MsaResult>` | ClustalW-style progressive alignment |

### Partial Order Alignment (`poa.rs`)

| Type | Description |
|------|-------------|
| `PoaGraph` | DAG for multiple sequence alignment |
| `PoaScoring` | Scoring parameters (default: match=2, mismatch=-1, gap=-2) |

| Method | Description |
|--------|-------------|
| `PoaGraph::from_sequence(seq) -> Self` | Initialize with first sequence |
| `PoaGraph::add_sequence(seq, scoring) -> Result<()>` | Align and integrate new sequence |
| `PoaGraph::consensus() -> Vec<u8>` | Heaviest-path consensus |
| `PoaGraph::topological_sort() -> Vec<usize>` | Topological ordering |
| `PoaGraph::to_dot() -> String` | DOT format visualization |

### Minimizers (`minimizers.rs`)

| Type | Description |
|------|-------------|
| `Minimizer` | A minimizer hit: `position`, `hash` |

| Function | Description |
|----------|-------------|
| `minimizers(seq, w, k) -> Result<Vec<Minimizer>>` | Extract (w,k)-minimizers from DNA |
| `find_seed_matches(a_mins, b_mins) -> Vec<(usize, usize)>` | Find matching minimizer positions |

### Seed-and-extend (`seed_extend.rs`)

| Type | Description |
|------|-------------|
| `Seed` | Anchor: `query_pos`, `target_pos` |
| `SeedChain` | Chain of colinear seeds with total score |

| Function | Description |
|----------|-------------|
| `chain_seeds(seeds, k, max_gap) -> SeedChain` | DP chaining of colinear seeds |
| `seed_and_extend(query, target, k, w, max_gap, scoring, bandwidth) -> Result<AlignmentResult>` | Full pipeline: minimizers, chain, extend |

### Wavefront alignment (`wfa.rs`, `wfa` feature)

| Function | Description |
|----------|-------------|
| `wfa_align(query, target, scoring) -> Result<AlignmentResult>` | WFA global alignment (O(ns) time) |

### X-drop/Z-drop extension (`xdrop.rs`)

| Type | Description |
|------|-------------|
| `XDropParams` | X-drop configuration: `x_drop` (default 50), `bandwidth` (default 100) |
| `ZDropParams` | Z-drop configuration: `z_drop` (default 200), `bandwidth` (default 500) |
| `ExtensionResult` | Result: `score`, `query_len`, `target_len`, `cigar`, `dropped` |

| Function | Description |
|----------|-------------|
| `x_drop_extend(query, target, scoring, params) -> Result<ExtensionResult>` | X-drop seed extension |
| `z_drop_extend(query, target, scoring, params) -> Result<ExtensionResult>` | Z-drop seed extension |

### Spliced alignment (`spliced.rs`)

| Type | Description |
|------|-------------|
| `SpliceSiteType` | Enum: `GtAg`, `GcAg`, `AtAc`, `NonCanonical` |
| `SplicedAlignParams` | Splice penalties, intron length bounds, max exons |
| `ExonAlignment` | Single exon: query/target coordinates, score, CIGAR |
| `SplicedAlignResult` | Combined alignment + exon list + introns |

| Function | Description |
|----------|-------------|
| `spliced_align(query, target, scoring, params) -> Result<SplicedAlignResult>` | Intron-aware DP alignment |
| `chain_exons(exons, target, params) -> Result<SplicedAlignResult>` | Chain pre-computed exon alignments |

### LCSk++ sparse alignment (`lcsk.rs`)

| Type | Description |
|------|-------------|
| `SparseAlignResult` | Result: `score`, `anchors`, `num_matches` |

| Function | Description |
|----------|-------------|
| `find_kmer_matches(a, b, k) -> Vec<(usize, usize)>` | Find shared k-mer positions |
| `lcsk_plusplus(matches, k) -> (usize, Vec<(usize, usize)>)` | LCSk++ with Fenwick tree |
| `sparse_align(a, b, k) -> Result<SparseAlignResult>` | Full pipeline |

### Pair HMM (`pair_hmm.rs`)

| Type | Description |
|------|-------------|
| `PairHmmParams` | Log-space transition and emission parameters |
| `PairHmmState` | Enum: `Match`, `InsertX`, `InsertY` |
| `PairHmmAlignment` | Result: `score`, `states` |

| Function | Description |
|----------|-------------|
| `pair_hmm_forward(a, b, params) -> Result<f64>` | Forward algorithm (log-likelihood) |
| `pair_hmm_viterbi(a, b, params) -> Result<(f64, PairHmmAlignment)>` | Viterbi decoding |

### Profile HMM (`profile_hmm.rs`)

| Type | Description |
|------|-------------|
| `Alphabet` | Enum: `Dna`, `Rna`, `Protein` |
| `ProfileHmmConfig` | Construction parameters |
| `ProfileHmmState` | Enum: `Match(usize)`, `Insert(usize)`, `Delete(usize)` |
| `ProfileHmmResult` | Result: `score`, `states` |
| `ProfileHmm` | Plan 7 Profile HMM |

| Method | Description |
|--------|-------------|
| `ProfileHmm::from_msa(msa, config) -> Result<Self>` | Build from MSA |
| `ProfileHmm::viterbi(seq) -> Result<ProfileHmmResult>` | Viterbi decoding |
| `ProfileHmm::forward(seq) -> Result<f64>` | Forward algorithm |
| `ProfileHmm::backward(seq) -> Result<f64>` | Backward algorithm |
| `ProfileHmm::calibrate(n_samples, seq_length, seed) -> Result<()>` | Gumbel calibration |
| `ProfileHmm::evalue(score, db_size) -> Result<f64>` | E-value from calibration |

### GPU batch alignment (`gpu/`)

| Type/Function | Description |
|---------------|-------------|
| `GpuBackend` | Enum: `Metal`, `Cuda` |
| `GpuAlignConfig` | Configuration: `max_bandwidth`, `min_pairs_for_gpu` |
| `align_batch_gpu(pairs, mode, scoring) -> Result<Vec<AlignmentResult>>` | Auto-select GPU |
| `align_batch_gpu_with_config(pairs, mode, scoring, config) -> Result<Vec<AlignmentResult>>` | Explicit config |
| `align_batch_on(pairs, mode, scoring, backend, config) -> Result<Vec<AlignmentResult>>` | Specific backend |
| `available_backends() -> Vec<GpuBackend>` | List available backends |

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support |
| `simd` | Yes | NEON/SSE4.1/AVX2 banded and score-only alignment |
| `wasm` | No | WASM target |
| `wfa` | No | Wavefront alignment |
| `serde` | No | Serialization support |
| `cuda` | No | CUDA GPU batch alignment |
| `metal` | No | Metal GPU batch alignment |
| `parallel` | No | Rayon parallelism |

## Dependencies

- `cyanea-core` -- error types
- `metal-rs` (feature = "metal") -- Apple Metal bindings
- `cudarc` (feature = "cuda") -- Safe CUDA driver API

## Tests

321 unit + 8 doc tests across 20+ source files.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 158 | Module declarations, `align()` dispatcher, property tests |
| `types.rs` | 286 | `AlignmentMode`, `CigarOp`, `AlignmentResult` |
| `cigar.rs` | ~600 | CIGAR parsing, validation, arithmetic, MD tags |
| `scoring.rs` | 493 | Scoring matrices (simple + BLOSUM/PAM) |
| `needleman_wunsch.rs` | 270 | Global alignment (Gotoh 3-matrix) |
| `smith_waterman.rs` | 280 | Local alignment (Gotoh 3-matrix) |
| `semi_global.rs` | 313 | Semi-global alignment |
| `batch.rs` | 80 | Batch pairwise alignment |
| `simd.rs` | 443 | Banded alignment (all modes) |
| `simd_sw.rs` | ~400 | SIMD SW score-only (Farrar striped) |
| `msa.rs` | 384 | Progressive multiple sequence alignment |
| `poa.rs` | 766 | Partial Order Alignment (Lee 2002) |
| `minimizers.rs` | ~200 | (w,k)-minimizer extraction |
| `seed_extend.rs` | ~300 | Seed chaining and extension |
| `wfa.rs` | ~500 | Wavefront alignment |
| `xdrop.rs` | ~400 | X-drop/Z-drop extension |
| `spliced.rs` | ~500 | Spliced alignment |
| `lcsk.rs` | 598 | LCSk++ sparse alignment |
| `pair_hmm.rs` | 697 | Pair HMM forward and Viterbi |
| `profile_hmm.rs` | ~900 | Profile HMM (Plan 7) |
| `gpu/mod.rs` | 210 | GPU dispatch, backend selection |
| `gpu/common.rs` | 335 | Sequence encoding, traceback |
| `gpu/metal_align.rs` | 231 | Metal aligner |
| `gpu/cuda_align.rs` | 186 | CUDA aligner |
