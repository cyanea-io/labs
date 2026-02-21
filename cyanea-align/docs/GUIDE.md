# cyanea-align Usage Guide

## Installation

```toml
[dependencies]
cyanea-align = { version = "0.1", features = ["simd"] }
```

Additional features: `wfa` (wavefront alignment), `metal` (Apple GPU), `cuda` (NVIDIA GPU), `parallel` (Rayon).

## Global Alignment (Needleman-Wunsch)

```rust
use cyanea_align::{needleman_wunsch, ScoringMatrix, ScoringScheme};

let scoring = ScoringScheme::Simple(ScoringMatrix::dna_default()); // +2/-1/-5/-2
let result = needleman_wunsch(b"ACGTACGT", b"ACGTTCGT", &scoring).unwrap();

println!("Score: {}", result.score);
println!("CIGAR: {}", result.cigar_string());
println!("Identity: {:.1}%", result.identity() * 100.0);
println!("Aligned query:  {}", String::from_utf8_lossy(&result.aligned_query));
println!("Aligned target: {}", String::from_utf8_lossy(&result.aligned_target));
```

## Local Alignment (Smith-Waterman)

```rust
use cyanea_align::{smith_waterman, ScoringMatrix, ScoringScheme};

let scoring = ScoringScheme::Simple(ScoringMatrix::dna_default());
let result = smith_waterman(b"XXXACGTACGTXXX", b"ACGTACGT", &scoring).unwrap();

println!("Score: {}", result.score);
println!("Query region: {}..{}", result.query_start, result.query_end);
println!("Target region: {}..{}", result.target_start, result.target_end);
```

## Semi-global Alignment

Free leading/trailing gaps -- useful for read mapping where the read may be a subsequence of the reference.

```rust
use cyanea_align::{semi_global, ScoringMatrix, ScoringScheme};

let scoring = ScoringScheme::Simple(ScoringMatrix::dna_default());
let result = semi_global(b"ACGT", b"NNNACGTNNN", &scoring).unwrap();
// Leading/trailing gaps in the query are not penalized
```

## Banded Alignment for Long Sequences

Use banded alignment when sequences are similar and you want O(n*w) instead of O(n*m).

```rust
use cyanea_align::{ScoringMatrix, ScoringScheme};
use cyanea_align::simd::{banded_nw, banded_sw, banded_score_only};
use cyanea_align::AlignmentMode;

let scoring = ScoringScheme::Simple(ScoringMatrix::dna_default());
let bandwidth = 50;

// Banded global alignment
let result = banded_nw(b"ACGT".repeat(100).as_bytes(),
                        b"ACGT".repeat(100).as_bytes(),
                        &scoring, bandwidth).unwrap();

// Score-only (no traceback, faster)
let score = banded_score_only(b"ACGT".repeat(100).as_bytes(),
                               b"ACGT".repeat(100).as_bytes(),
                               &scoring, bandwidth, AlignmentMode::Global).unwrap();
```

## Protein Alignment with Substitution Matrices

```rust
use cyanea_align::{smith_waterman, SubstitutionMatrix, ScoringScheme};

let blosum62 = SubstitutionMatrix::blosum62();
let scoring = ScoringScheme::Substitution(blosum62);

let result = smith_waterman(b"MKWVTFISLLFLFSSAYS", b"MKWVTFISLLFSSAYS", &scoring).unwrap();
println!("Score: {}", result.score);
```

## Multiple Sequence Alignment (MSA)

```rust
use cyanea_align::{ScoringMatrix, ScoringScheme};
use cyanea_align::msa::progressive_msa;

let sequences: Vec<&[u8]> = vec![
    b"ACGTACGT",
    b"ACGTTCGT",
    b"ACGTACG",
    b"ACGTGCGT",
];

let scoring = ScoringScheme::Simple(ScoringMatrix::dna_default());
let msa = progressive_msa(&sequences, &scoring).unwrap();

println!("{} columns, conservation: {:.2}", msa.column_count, msa.conservation);
for seq in &msa.aligned {
    println!("{}", String::from_utf8_lossy(seq));
}
```

## Partial Order Alignment (POA) for Consensus

```rust
use cyanea_align::poa::{PoaGraph, PoaScoring};

let scoring = PoaScoring::default(); // match=2, mismatch=-1, gap=-2

let mut graph = PoaGraph::from_sequence(b"ACGTACGT");
graph.add_sequence(b"ACGTTCGT", &scoring).unwrap();
graph.add_sequence(b"ACGTACGT", &scoring).unwrap();
graph.add_sequence(b"ACGTGCGT", &scoring).unwrap();

let consensus = graph.consensus();
println!("Consensus: {}", String::from_utf8_lossy(&consensus));
println!("Graph: {} nodes, {} sequences", graph.num_nodes(), graph.num_sequences());
```

## Seed-and-Extend with Minimizers

For long sequences where full DP is too expensive.

```rust
use cyanea_align::{ScoringMatrix, ScoringScheme};
use cyanea_align::seed_extend::seed_and_extend;

let scoring = ScoringScheme::Simple(ScoringMatrix::dna_default());
let result = seed_and_extend(
    b"ACGTACGTACGTACGT",  // query
    b"ACGTACGTACGTACGT",  // target
    8,                      // k (k-mer size)
    5,                      // w (minimizer window)
    100,                    // max_gap between seeds
    &scoring,
    50,                     // bandwidth for extension
).unwrap();

println!("Score: {}, Identity: {:.1}%", result.score, result.identity() * 100.0);
```

## CIGAR String Parsing and Manipulation

```rust
use cyanea_align::cigar::{
    parse_cigar, cigar_string, reference_consumed, query_consumed,
    merge_cigars, reverse_cigar, generate_md_tag,
};

// Parse a SAM CIGAR string
let ops = parse_cigar("10M3I4D2S").unwrap();
assert_eq!(cigar_string(&ops), "10M3I4D2S");

// Coordinate queries
assert_eq!(reference_consumed(&ops), 14);  // M + D
assert_eq!(query_consumed(&ops), 15);      // M + I + S

// Arithmetic
let merged = merge_cigars(&parse_cigar("5M").unwrap(), &parse_cigar("3M2I").unwrap());
let reversed = reverse_cigar(&ops);

// MD tag generation
// let md = generate_md_tag(&ops, reference_seq);
```

## Wavefront Alignment (WFA)

O(ns) exact alignment -- dramatically faster than Needleman-Wunsch for similar sequences.

```rust
use cyanea_align::{wfa_align, ScoringMatrix, ScoringScheme};

let scoring = ScoringScheme::Simple(ScoringMatrix::dna_default());
let result = wfa_align(b"ACGTACGT", b"ACGTACGT", &scoring).unwrap();
println!("Score: {}", result.score);
```

## GPU-Accelerated Batch Alignment

```rust
use cyanea_align::{ScoringMatrix, ScoringScheme, AlignmentMode};
use cyanea_align::gpu::{align_batch_gpu, available_backends};

// Check available GPU backends
let backends = available_backends();
println!("Available: {:?}", backends);

let pairs: Vec<(&[u8], &[u8])> = vec![
    (b"ACGTACGT", b"ACGTTCGT"),
    (b"ACGTACGT", b"ACGTGCGT"),
    // ... hundreds more pairs for GPU efficiency
];

let scoring = ScoringScheme::Simple(ScoringMatrix::dna_default());
let results = align_batch_gpu(&pairs, AlignmentMode::Local, &scoring).unwrap();
// Automatically selects Metal/CUDA/CPU based on availability and batch size
```
