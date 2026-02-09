# cyanea-align

Pairwise sequence alignment with affine gap penalties. Supports DNA, RNA, and protein alignment with standard substitution matrices.

## Status: Complete

All three alignment modes (global, local, semi-global) are fully implemented with affine gap scoring (Gotoh 3-matrix), BLOSUM/PAM matrices, batch alignment, and banded variants. GPU backends are stubbed pending hardware SDK availability; true SIMD vectorization deferred pending profile-guided benchmarks (scalar banded kernels serve as baseline).

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
| `align(query, target, mode, scoring) -> Result<AlignmentResult>` | Dispatcher -- routes to NW, SW, or semi-global based on mode |
| `needleman_wunsch(query, target, scoring) -> Result<AlignmentResult>` | Global alignment (Gotoh 3-matrix) |
| `smith_waterman(query, target, scoring) -> Result<AlignmentResult>` | Local alignment (Gotoh 3-matrix) |
| `semi_global(query, target, scoring) -> Result<AlignmentResult>` | Semi-global alignment (free leading/trailing gaps) |
| `align_batch(pairs, mode, scoring) -> Result<Vec<AlignmentResult>>` | Batch pairwise alignment (all modes) |

### Banded alignment (`simd.rs`)

| Function | Description |
|----------|-------------|
| `banded_nw(query, target, scoring, bandwidth) -> Result<AlignmentResult>` | Banded global alignment |
| `banded_sw(query, target, scoring, bandwidth) -> Result<AlignmentResult>` | Banded local alignment |
| `banded_semi_global(query, target, scoring, bandwidth) -> Result<AlignmentResult>` | Banded semi-global alignment |
| `banded_score_only(query, target, scoring, bandwidth, mode) -> Result<i32>` | Score-only banded alignment (all modes) |

### Stubbed (out of scope)

| Module | Description |
|--------|-------------|
| `msa` | Multiple sequence alignment (progressive, profile HMMs, iterative refinement) |
| `gpu` | GPU-accelerated alignment (CUDA, Metal backends -- requires hardware SDK) |

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

75 unit tests + 2 doc tests across all source files.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 97 | Module declarations, `align()` dispatcher |
| `types.rs` | 286 | `AlignmentMode`, `CigarOp`, `AlignmentResult` |
| `scoring.rs` | 493 | Scoring matrices (simple + BLOSUM/PAM) |
| `needleman_wunsch.rs` | 270 | Global alignment (Gotoh 3-matrix) |
| `smith_waterman.rs` | 280 | Local alignment (Gotoh 3-matrix) |
| `semi_global.rs` | 313 | Semi-global alignment (free leading/trailing gaps) |
| `batch.rs` | 80 | Batch pairwise alignment (all modes) |
| `simd.rs` | 443 | Banded alignment (global, local, semi-global) |
| `msa.rs` | 28 | Stub |
| `gpu.rs` | 45 | Stub |
