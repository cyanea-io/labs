# cyanea-align

Pairwise sequence alignment with affine gap penalties. Supports DNA, RNA, and protein alignment with standard substitution matrices.

## Status: Mostly Complete

Core alignment algorithms (Needleman-Wunsch, Smith-Waterman) are fully implemented with affine gap scoring, BLOSUM/PAM matrices, and batch alignment. SIMD acceleration, multiple sequence alignment, and GPU backends are stubbed for future work.

## Public API

### Alignment modes and results (`types.rs`)

| Type | Description |
|------|-------------|
| `AlignmentMode` | Enum: `Local`, `Global`, `SemiGlobal` |
| `CigarOp` | Enum: `Match`, `Mismatch`, `Insert`, `Delete` |
| `AlignmentResult` | Full result: `score`, `query_start/end`, `target_start/end`, `cigar`, `aligned_query/target` |

**AlignmentResult derived methods:**

| Method | Description |
|--------|-------------|
| `cigar_string() -> String` | Compact CIGAR string (e.g., `3=1X4=`) |
| `identity() -> f64` | Fraction of matching positions |
| `matches() -> usize` | Number of matches |
| `mismatches() -> usize` | Number of mismatches |
| `gaps() -> usize` | Number of gap positions |
| `length() -> usize` | Total alignment length |

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
| `SubstitutionMatrix::blosum62()` | BLOSUM62 |
| `SubstitutionMatrix::blosum45()` | BLOSUM45 |
| `SubstitutionMatrix::blosum80()` | BLOSUM80 |
| `SubstitutionMatrix::pam250()` | PAM250 |

### Alignment functions

| Function | Description |
|----------|-------------|
| `align(query, target, mode, scoring) -> Result<AlignmentResult>` | Dispatcher -- routes to NW or SW based on mode |
| `needleman_wunsch(query, target, scoring) -> Result<AlignmentResult>` | Global/semiglobal alignment |
| `smith_waterman(query, target, scoring) -> Result<AlignmentResult>` | Local alignment |
| `align_batch(pairs, mode, scoring) -> Result<Vec<AlignmentResult>>` | Batch pairwise alignment |

### Planned (stubbed)

| Module | Description |
|--------|-------------|
| `simd` | SIMD-accelerated DP kernels (striped, banded, AVX2/NEON) |
| `msa` | Multiple sequence alignment (progressive, profile HMMs, iterative refinement) |
| `gpu` | GPU-accelerated alignment (CUDA, Metal backends) |

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support |
| `wasm` | No | WASM target |
| `serde` | No | Serialization support |
| `cuda` | No | CUDA GPU backend (future) |
| `metal` | No | Metal GPU backend (future) |

## Dependencies

- `cyanea-core` -- error types

## Tests

42 tests across `scoring.rs`, `needleman_wunsch.rs`, `smith_waterman.rs`, `batch.rs`, and `types.rs`.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 98 | Module declarations, `align()` dispatcher |
| `types.rs` | 286 | `AlignmentMode`, `CigarOp`, `AlignmentResult` |
| `scoring.rs` | 493 | Scoring matrices (simple + BLOSUM/PAM) |
| `needleman_wunsch.rs` | 270 | Global/semiglobal alignment |
| `smith_waterman.rs` | 280 | Local alignment |
| `batch.rs` | 83 | Batch pairwise alignment |
| `simd.rs` | 29 | Stub |
| `msa.rs` | 28 | Stub |
| `gpu.rs` | 45 | Stub |
