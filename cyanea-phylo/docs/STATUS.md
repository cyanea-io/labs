# cyanea-phylo

Phylogenetics: tree data structures, Newick I/O, evolutionary distance models, tree comparison, and tree construction algorithms.

## Status: Complete

All phylogenetics functionality is implemented including tree representation, Newick and NEXUS I/O, evolutionary distance models (p-distance, Jukes-Cantor, Kimura 2-parameter), Robinson-Foulds comparison, UPGMA/NJ tree construction, ancestral sequence reconstruction (Fitch/Sankoff parsimony and marginal ML), substitution models (JC69, GTR, HKY85, discrete gamma rate heterogeneity), maximum likelihood tree scoring (Felsenstein pruning with GTR+Gamma support), NNI tree search, bootstrap support estimation, tree rerooting/midpoint rooting, subtree extraction, molecular clock dating, tree drawing coordinates (rectangular/cladogram/radial), and consensus tree construction.

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
| `to_newick() -> String` | Serialize to Newick format |
| `leaf_count() -> usize` | Number of leaves |
| `leaf_names() -> Vec<String>` | Sorted leaf taxon names |
| `subtree_leaf_names(node_id) -> BTreeSet<String>` | Leaf names in subtree |
| `total_branch_length() -> f64` | Sum of all branch lengths |
| `reroot(node_id, position) -> Result<Self>` | Reroot on edge to node |
| `midpoint_root() -> Result<Self>` | Midpoint rerooting |
| `extract_subtree(leaf_names) -> Result<Self>` | Extract subtree with given leaves |
| `subtree_at(node_id) -> Result<Self>` | Extract subtree rooted at node |
| `mrca(a, b) -> Result<NodeId>` | Most recent common ancestor |
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
| `GtrParams` | GTR model: 6 exchangeability rates + 4 base frequencies |
| `GtrParams::new(rates, freqs) -> Result<Self>` | Validated GTR parameter construction |
| `GtrParams::rate_matrix() -> [[f64; 4]; 4]` | Normalized instantaneous rate matrix Q |
| `GtrParams::probability_fn() -> impl Fn(f64) -> [[f64; 4]; 4]` | Eigendecomposed P(t) = exp(Qt) closure |
| `gtr_probability(params, t) -> [[f64; 4]; 4]` | Convenience wrapper for GTR P(t) |
| `hky85_params(kappa, freqs) -> Result<GtrParams>` | HKY85 as GTR special case |
| `GammaRates` | Discrete gamma rate heterogeneity (Yang 1994) |
| `GammaRates::new(alpha, n_categories) -> Result<Self>` | Create gamma rate categories |
| `GammaRates::category_rates() -> Vec<f64>` | Discrete quantile rates averaging to 1.0 |

### Maximum likelihood (`likelihood.rs`)

| Function | Description |
|----------|-------------|
| `tree_likelihood(tree, sequences, model) -> Result<f64>` | Log-likelihood via Felsenstein's pruning algorithm |
| `tree_likelihood_gtr(tree, sequences, model, freqs, gamma) -> Result<f64>` | Generalized log-likelihood with custom model and gamma rates |
| `nni_search(tree, sequences, model) -> Result<PhyloTree>` | Hill-climbing NNI search for improved tree topology |

### Marginal ancestral reconstruction (`marginal.rs`)

| Type/Function | Description |
|---------------|-------------|
| `MarginalPosterior` | `probs: [f64; 4]` posterior over nucleotide states |
| `MarginalReconstruction` | `posteriors: Vec<Vec<MarginalPosterior>>`, `map_states: Vec<Vec<u8>>` |
| `marginal_reconstruct(tree, sequences, model, freqs) -> Result<MarginalReconstruction>` | Two-pass ML marginal ancestral reconstruction |

### Bootstrap support (`bootstrap.rs`)

| Function | Description |
|----------|-------------|
| `bootstrap_support(sequences, original_tree, tree_builder, n_replicates) -> Result<Vec<f64>>` | Non-parametric bootstrap support values for internal edges |
| `bipartitions(tree) -> Vec<BTreeSet<String>>` | Extract non-trivial bipartitions (splits) from a tree |

### Consensus trees (`consensus.rs`)

| Type/Function | Description |
|---------------|-------------|
| `ConsensusType` | Enum: `Strict`, `MajorityRule`, `ExtendedMajorityRule` |
| `SupportedBipartition` | `leaves: BTreeSet<String>`, `support: f64` |
| `consensus_tree(trees, consensus_type) -> Result<PhyloTree>` | Build consensus tree from tree collection |
| `bipartition_frequencies(trees) -> Result<Vec<SupportedBipartition>>` | Bipartition frequencies across trees |

### Molecular clock dating (`dating.rs`)

| Type/Function | Description |
|---------------|-------------|
| `Calibration` | `node_id: NodeId`, `age: f64` |
| `DatingResult` | `node_ages: Vec<f64>`, `rate: f64` |
| `strict_clock(tree, calibration) -> Result<DatingResult>` | Strict molecular clock dating |
| `root_to_tip_regression(tree, tip_dates) -> Result<DatingResult>` | Root-to-tip regression dating |

### Tree drawing (`drawing.rs`)

| Type/Function | Description |
|---------------|-------------|
| `LayoutStyle` | Enum: `Rectangular`, `Cladogram`, `Radial` |
| `NodeCoord` | `node_id`, `x`, `y` |
| `Edge` | `from`, `to` |
| `TreeLayout` | `coords: Vec<NodeCoord>`, `edges: Vec<Edge>` |
| `tree_layout(tree, style, width, height) -> Result<TreeLayout>` | Compute 2D layout coordinates |

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

134 tests across 14 source files.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 50 | Module declarations, re-exports |
| `tree.rs` | 750 | PhyloTree data structure, traversal, rerooting, subtree extraction |
| `newick.rs` | 285 | Newick parser and writer |
| `distance.rs` | 291 | Evolutionary distance models |
| `compare.rs` | 210 | Robinson-Foulds and branch score comparison |
| `construct.rs` | 425 | UPGMA and neighbor-joining |
| `nexus.rs` | 314 | NEXUS format parser and writer |
| `reconstruct.rs` | 454 | Fitch and Sankoff ancestral reconstruction |
| `models.rs` | 630 | JC69, GTR, HKY85 substitution models, discrete gamma rates |
| `likelihood.rs` | 620 | Felsenstein pruning, NNI search, GTR+Gamma likelihood |
| `bootstrap.rs` | 300 | Bootstrap support and bipartition extraction |
| `marginal.rs` | 240 | Marginal ML ancestral reconstruction |
| `dating.rs` | 200 | Molecular clock dating |
| `drawing.rs` | 230 | Tree drawing coordinates (rectangular, cladogram, radial) |
| `consensus.rs` | 380 | Consensus tree construction |
