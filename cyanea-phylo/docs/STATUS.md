# cyanea-phylo

Phylogenetics: tree data structures, Newick I/O, evolutionary distance models, tree comparison, and tree construction algorithms.

## Status: Complete

All phylogenetics functionality is implemented including tree representation, Newick and NEXUS I/O, evolutionary distance models (p-distance, Jukes-Cantor, Kimura 2-parameter), Robinson-Foulds comparison, UPGMA/NJ tree construction, ancestral sequence reconstruction (Fitch and Sankoff parsimony), substitution models (JC69), maximum likelihood tree scoring (Felsenstein pruning), NNI tree search, and bootstrap support estimation.

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

### NEXUS I/O (`nexus.rs`)

| Type/Function | Description |
|---------------|-------------|
| `NexusFile` | `taxa`, `trees: Vec<NamedTree>` |
| `NamedTree` | `name`, `tree: PhyloTree` |
| `parse(input) -> Result<NexusFile>` | Parse NEXUS format |
| `write(taxa, trees) -> String` | Serialize to NEXUS format |

### Ancestral reconstruction (`reconstruct.rs`)

| Type/Function | Description |
|---------------|-------------|
| `AncestralStates` | `states: Vec<u8>`, `n_changes: usize` |
| `CostMatrix` | Cost matrix for Sankoff algorithm |
| `CostMatrix::uniform(n_states) -> Self` | Equal-cost matrix |
| `CostMatrix::from_flat(costs, n_states) -> Result<Self>` | Custom cost matrix |
| `fitch(tree, leaf_states) -> Result<AncestralStates>` | Fitch maximum parsimony |
| `sankoff(tree, leaf_states, cost_matrix) -> Result<AncestralStates>` | Sankoff weighted parsimony |
| `reconstruct_sequences(tree, alignment) -> Result<Vec<AncestralStates>>` | Per-site ancestral reconstruction |

### Substitution models (`models.rs`)

| Type/Function | Description |
|---------------|-------------|
| `NUM_STATES` | Constant: number of nucleotide states (4) |
| `JC69_FREQ` | Constant: JC69 equilibrium frequency (0.25) |
| `nucleotide_index(b) -> Option<usize>` | Map nucleotide byte to index (A=0, C=1, G=2, T/U=3) |
| `jc69_probability(t) -> [[f64; 4]; 4]` | JC69 transition probability matrix for branch length `t` |

### Maximum likelihood (`likelihood.rs`)

| Function | Description |
|----------|-------------|
| `tree_likelihood(tree, sequences, model) -> Result<f64>` | Log-likelihood via Felsenstein's pruning algorithm |
| `nni_search(tree, sequences, model) -> Result<PhyloTree>` | Hill-climbing NNI search for improved tree topology |

### Bootstrap support (`bootstrap.rs`)

| Function | Description |
|----------|-------------|
| `bootstrap_support(sequences, original_tree, tree_builder, n_replicates) -> Result<Vec<f64>>` | Non-parametric bootstrap support values for internal edges |
| `bipartitions(tree) -> Vec<BTreeSet<String>>` | Extract non-trivial bipartitions (splits) from a tree |

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

110 tests across 10 source files.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 44 | Module declarations, re-exports |
| `tree.rs` | 420 | PhyloTree data structure and traversal |
| `newick.rs` | 285 | Newick parser and writer |
| `distance.rs` | 291 | Evolutionary distance models |
| `compare.rs` | 230 | Robinson-Foulds and branch score comparison |
| `construct.rs` | 425 | UPGMA and neighbor-joining |
| `nexus.rs` | 314 | NEXUS format parser and writer |
| `reconstruct.rs` | 454 | Fitch and Sankoff ancestral reconstruction |
| `models.rs` | 186 | JC69 substitution model and nucleotide indexing |
| `likelihood.rs` | 475 | Felsenstein pruning and NNI tree search |
| `bootstrap.rs` | 317 | Bootstrap support and bipartition extraction |
