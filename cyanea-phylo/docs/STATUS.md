# cyanea-phylo

Phylogenetics: tree data structures, Newick I/O, evolutionary distance models, tree comparison, and tree construction algorithms.

## Status: Mostly Complete

Tree representation, Newick parsing/writing, distance models (p-distance, Jukes-Cantor, Kimura 2-parameter), Robinson-Foulds comparison, and UPGMA/NJ tree construction are implemented. NEXUS format and ancestral reconstruction are stubbed.

## Public API

### Tree types (`tree.rs`)

| Type | Description |
|------|-------------|
| `NodeId` | `usize` identifier |
| `Node` | `id`, `name`, `branch_length`, `children`, `parent` |
| `PhyloTree` | Rooted phylogenetic tree |
| `PreorderIter` | Depth-first pre-order traversal |
| `PostorderIter` | Post-order traversal |

**PhyloTree methods:**

| Method | Description |
|--------|-------------|
| `from_newick(input) -> Result<Self>` | Parse from Newick string |
| `leaf_count() -> usize` | Number of leaves |
| `leaf_names() -> Vec<&str>` | Leaf taxon names |
| Pre-order / post-order iterators | Tree traversal |

### Newick I/O (`newick.rs`)

| Function | Description |
|----------|-------------|
| `parse(input) -> Result<PhyloTree>` | Parse Newick format |
| `write(tree) -> String` | Serialize to Newick format |

### Evolutionary distances (`distance.rs`)

| Type/Function | Description |
|---------------|-------------|
| `DistanceModel` | Enum: `PDistance`, `JukesCantor`, `Kimura2p` |
| `p_distance(a, b) -> Result<f64>` | Raw proportion of differences |
| `jukes_cantor(p) -> Result<f64>` | Jukes-Cantor correction: `-3/4 * ln(1 - 4p/3)` |
| `kimura_2p(transitions, transversions) -> Result<f64>` | Kimura 2-parameter model |
| `sequence_distance_matrix(seqs, model) -> Result<DistanceMatrix>` | Pairwise distance matrix (requires `ml` feature) |

### Tree comparison (`compare.rs`)

| Function | Description |
|----------|-------------|
| `robinson_foulds(t1, t2) -> Result<usize>` | Unweighted Robinson-Foulds distance |
| `robinson_foulds_normalized(t1, t2) -> Result<f64>` | Normalized RF: `RF / (2 * (n - 3))` |
| `branch_score_distance(t1, t2) -> Result<f64>` | Weighted branch score distance |

### Tree construction (`construct.rs`, requires `ml` feature)

| Function | Description |
|----------|-------------|
| `upgma(distances, leaf_names) -> Result<PhyloTree>` | UPGMA clustering |
| `neighbor_joining(distances, leaf_names) -> Result<PhyloTree>` | Neighbor-joining algorithm |

### Planned (stubbed)

| Module | Description |
|--------|-------------|
| `nexus` | NEXUS format parsing and writing |
| `reconstruct` | Ancestral sequence reconstruction (Fitch, Sankoff, ML) |

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support |
| `wasm` | No | WASM target |
| `serde` | No | Serialization support |
| `ml` | No | Enables `cyanea-ml` for distance matrices and tree construction |

## Dependencies

- `cyanea-core` -- error types
- `cyanea-ml` -- distance matrices (optional, `ml` feature)

## Tests

62 tests across 5 source files: tree (15), newick (12), distance (10), compare (10), construct (15).

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 36 | Module declarations, re-exports |
| `tree.rs` | 420 | PhyloTree data structure and traversal |
| `newick.rs` | 285 | Newick parser and writer |
| `distance.rs` | 291 | Evolutionary distance models |
| `compare.rs` | 230 | Robinson-Foulds and branch score comparison |
| `construct.rs` | 425 | UPGMA and neighbor-joining |
| `nexus.rs` | 7 | Stub |
| `reconstruct.rs` | 7 | Stub |
