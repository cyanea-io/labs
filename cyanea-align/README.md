# cyanea-align

> Sequence alignment algorithms for DNA, RNA, and protein sequences.

## What's Inside

- **Global alignment** -- Needleman-Wunsch with Gotoh 3-matrix affine gaps
- **Local alignment** -- Smith-Waterman with affine gap scoring
- **Semi-global alignment** -- free leading/trailing gaps
- **Banded alignment** -- O(n*w) for all three modes
- **Batch alignment** -- align multiple pairs in a single call
- **SIMD** -- NEON (aarch64), SSE4.1/AVX2 (x86_64) acceleration
- **Substitution matrices** -- BLOSUM30/45/62/80, PAM40/120/200/250
- **Multiple sequence alignment** -- ClustalW-style progressive MSA
- **GPU batch alignment** -- Metal and CUDA banded affine-gap DP with automatic CPU fallback
- **Partial Order Alignment** -- Lee 2002 DAG alignment for consensus sequences
- **LCSk++ sparse alignment** -- Fenwick tree-based O(n log n) DP for long sequences
- **Pair HMM** -- three-state (Match/InsertX/InsertY) forward and Viterbi in log-space
- **Profile HMM** -- Plan 7 architecture, Viterbi/Forward/Backward, E-value calibration
- **CIGAR utilities** -- full SAM alphabet (9 ops), parsing, validation, coordinate queries, arithmetic
- **X-drop/Z-drop extension** -- seed extension with early termination
- **Spliced alignment** -- intron-aware alignment for RNA-seq

## Quick Start

```toml
[dependencies]
cyanea-align = "0.1"
```

```rust
use cyanea_align::{smith_waterman, ScoringMatrix};

let scoring = ScoringMatrix::dna_default(); // +2/-1/-5/-2
let result = smith_waterman(b"ACGTACGT", b"ACGTGCGT", &scoring.into()).unwrap();

println!("Score: {}", result.score);
println!("CIGAR: {}", result.cigar_string());
println!("Identity: {:.1}%", result.identity() * 100.0);
```

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support |
| `simd` | Yes | NEON/SSE4.1/AVX2 banded alignment |
| `wasm` | No | WASM target marker |
| `wfa` | No | Wavefront alignment |
| `serde` | No | Serialize/Deserialize derives |
| `cuda` | No | CUDA GPU batch alignment (cudarc + NVRTC) |
| `metal` | No | Metal GPU batch alignment (metal-rs) |
| `parallel` | No | Rayon parallelism |

## Modules

| Module | Description |
|--------|-------------|
| `types` | `AlignmentMode`, `CigarOp`, `AlignmentResult` |
| `scoring` | `ScoringMatrix`, `SubstitutionMatrix` (BLOSUM/PAM) |
| `needleman_wunsch` | Global alignment (Gotoh 3-matrix) |
| `smith_waterman` | Local alignment (Gotoh 3-matrix) |
| `semi_global` | Semi-global alignment |
| `batch` | Batch pairwise alignment |
| `simd` | Banded alignment (all modes) |
| `msa` | Progressive multiple sequence alignment |
| `gpu` | GPU batch alignment dispatch (Metal/CUDA/CPU) |
| `lcsk` | LCSk++ sparse alignment |
| `poa` | Partial Order Alignment (DAG) |
| `pair_hmm` | Pair HMM forward/Viterbi |
| `profile_hmm` | Profile HMM (Plan 7) |
| `cigar` | CIGAR string utilities |

## See Also

- [API Reference (STATUS.md)](docs/STATUS.md)
- [Architecture](../ARCHITECTURE.md)
- [Build Guide](../docs/BUILDING.md)
