# API Reference -- cyanea-phylo

Phylogenetics: tree data structures, Newick/NEXUS I/O, evolutionary distance models, tree construction, ancestral reconstruction, substitution models, maximum likelihood inference, Bayesian MCMC, tree search, model selection, protein models, UniFrac, species tree estimation, and simulation.

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

### NEXUS I/O (`nexus.rs`)

| Type/Function | Description |
|---------------|-------------|
| `NexusFile` | `taxa`, `trees: Vec<NamedTree>` |
| `NamedTree` | `name`, `tree: PhyloTree` |
| `parse(input) -> Result<NexusFile>` | Parse NEXUS format |
| `write(taxa, trees) -> String` | Serialize to NEXUS format |

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

### Marginal ancestral reconstruction (`marginal.rs`)

| Type/Function | Description |
|---------------|-------------|
| `MarginalPosterior` | `probs: [f64; 4]` posterior over nucleotide states |
| `MarginalReconstruction` | `posteriors: Vec<Vec<MarginalPosterior>>`, `map_states: Vec<Vec<u8>>` |
| `marginal_reconstruct(tree, sequences, model, freqs) -> Result<MarginalReconstruction>` | Two-pass ML marginal ancestral reconstruction |

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

### Substitution model trait (`subst_model.rs`)

| Type/Function | Description |
|---------------|-------------|
| `trait SubstitutionModel` | Unified interface for N-state models: `n_states()`, `frequencies()`, `rate_matrix()`, `transition_probs(t)`, `n_free_params()` |
| `Jc69Model` | JC69 wrapper implementing `SubstitutionModel` |
| `Hky85Model` | HKY85 wrapper implementing `SubstitutionModel` |
| `GtrModel` | GTR wrapper implementing `SubstitutionModel` |

### Protein substitution models (`protein_models.rs`)

| Type/Function | Description |
|---------------|-------------|
| `AA_STATES` | Constant: number of amino acid states (20) |
| `amino_acid_index(aa) -> Option<usize>` | Map amino acid byte to index (A=0, R=1, ..., V=19) |
| `LgModel` | LG amino acid model (Le & Gascuel 2008) implementing `SubstitutionModel` |
| `WagModel` | WAG amino acid model (Whelan & Goldman 2001) implementing `SubstitutionModel` |
| `JttModel` | JTT amino acid model (Jones, Taylor & Thornton 1992) implementing `SubstitutionModel` |
| `DayhoffModel` | Dayhoff amino acid model (Dayhoff, Schwartz & Orcutt 1978) implementing `SubstitutionModel` |
| `CustomAaModel` | User-provided 20x20 exchangeabilities + equilibrium frequencies |
| `load_aa_model(name) -> Result<Box<dyn SubstitutionModel>>` | Load a named protein model ("LG", "WAG", "JTT", "Dayhoff") |

### Generic likelihood (`generic_likelihood.rs`)

| Function | Description |
|----------|-------------|
| `generic_tree_likelihood(tree, sequences, model, gamma) -> Result<f64>` | Log-likelihood via Felsenstein pruning for any N-state model (DNA or protein) |
| `site_likelihoods(tree, sequences, model, gamma) -> Result<Vec<f64>>` | Per-site log-likelihoods under any substitution model |

### Maximum likelihood (`likelihood.rs`)

| Function | Description |
|----------|-------------|
| `tree_likelihood(tree, sequences, model) -> Result<f64>` | Log-likelihood via Felsenstein's pruning algorithm (4-state) |
| `tree_likelihood_gtr(tree, sequences, model, freqs, gamma) -> Result<f64>` | Generalized log-likelihood with custom model and gamma rates |
| `nni_search(tree, sequences, model) -> Result<PhyloTree>` | Hill-climbing NNI search for improved tree topology |

### Tree search (`tree_search.rs`)

| Type/Function | Description |
|---------------|-------------|
| `AnnealingConfig` | Configuration for simulated annealing: `initial_temp`, `cooling_rate`, `min_temp`, `max_iterations` |
| `spr_move(tree, prune, regraft) -> Result<PhyloTree>` | Subtree pruning and regrafting move |
| `tbr_move(tree, bisect_edge, reconnect_a, reconnect_b) -> Result<PhyloTree>` | Tree bisection and reconnection move |
| `parsimony_ratchet(tree, sequences, model, gamma, iterations) -> Result<PhyloTree>` | Parsimony ratchet with reweighted sites |
| `stochastic_nni(tree, sequences, model, gamma, config) -> Result<PhyloTree>` | Stochastic NNI with simulated annealing acceptance |
| `spr_search(tree, sequences, model, gamma, max_rounds) -> Result<PhyloTree>` | SPR hill-climbing search |
| `nni_search_generic(tree, sequences, model, gamma) -> Result<PhyloTree>` | NNI search using the generic N-state likelihood |

### Model selection (`model_selection.rs`)

| Type/Function | Description |
|---------------|-------------|
| `ModelResult` | `name`, `log_likelihood`, `n_params`, `n_sites`, `aic`, `bic`, `aicc` |
| `ModelSelectionResult` | `models: Vec<ModelResult>`, `best_aic`, `best_bic` |
| `LrtResult` | `stat`, `df`, `p_value`, `reject` |
| `aic(log_likelihood, n_params) -> f64` | Akaike Information Criterion |
| `bic(log_likelihood, n_params, n_sites) -> f64` | Bayesian Information Criterion |
| `aicc(log_likelihood, n_params, n_sites) -> f64` | Corrected AIC for small samples |
| `lrt(ll_null, ll_alt, df) -> LrtResult` | Likelihood ratio test between nested models |
| `model_finder(tree, sequences, models, gamma) -> Result<ModelSelectionResult>` | Compare multiple models and rank by AIC/BIC |

### Bayesian MCMC (`mcmc.rs`)

| Type/Function | Description |
|---------------|-------------|
| `McmcConfig` | `n_generations`, `sample_every`, `burnin`, `proposal_weights` |
| `ProposalWeights` | Relative weights: `nni`, `spr`, `branch_scale`, `branch_slide`, `model_params` |
| `TreePrior` | Enum: `Uniform`, `CoalescentConstant`, `CoalescentExponential`, `BirthDeath` |
| `ClockModel` | Enum: `Strict { rate }`, `UncorrelatedLognormal { mean_rate, stdev }` |
| `McmcSample` | `generation`, `tree`, `log_likelihood`, `log_prior`, `log_posterior` |
| `McmcResult` | `samples: Vec<McmcSample>`, `acceptance_rates` |
| `ConvergenceDiag` | `ess`, `mean_log_likelihood`, `variance_log_likelihood` |
| `PosteriorSummary` | `map_tree`, `mean_branch_lengths`, `node_age_intervals`, `clade_credibilities` |
| `mcmc_sample(tree, sequences, model, prior, clock, config) -> Result<McmcResult>` | Run MCMC sampling for Bayesian inference |
| `convergence_diagnostics(samples) -> ConvergenceDiag` | ESS, mean/variance of log-likelihood |
| `posterior_summary(samples) -> PosteriorSummary` | MAP tree, credible intervals, clade probabilities |

### Simulation (`simulation.rs`)

| Type/Function | Description |
|---------------|-------------|
| `SimulatedAlignment` | `names: Vec<String>`, `sequences: Vec<Vec<u8>>`, `n_substitutions: usize` |
| `simulate_evolution(tree, model, seq_length, seed) -> Result<SimulatedAlignment>` | Sequence evolution along a phylogenetic tree |
| `simulate_coalescent(n_samples, pop_size, seed) -> Result<PhyloTree>` | Coalescent tree simulation (constant population) |
| `simulate_coalescent_growth(n_samples, current_pop, growth_rate, seed) -> Result<PhyloTree>` | Coalescent with exponential population growth |

### Species tree (`species_tree.rs`)

| Type/Function | Description |
|---------------|-------------|
| `ReconciliationResult` | `duplications`, `losses`, `deep_coalescences`, `cost` |
| `ConcordanceFactors` | `gene_cf: Vec<(NodeId, f64)>`, `site_cf: Vec<(NodeId, f64)>` |
| `astral_species_tree(gene_trees) -> Result<PhyloTree>` | ASTRAL-style quartet-based species tree estimation |
| `reconcile(gene_tree, species_tree) -> Result<ReconciliationResult>` | Gene tree / species tree reconciliation |
| `concordance_factors(species_tree, gene_trees, sequences) -> Result<ConcordanceFactors>` | Gene and site concordance factors per branch |

### UniFrac (`unifrac.rs`)

| Type/Function | Description |
|---------------|-------------|
| `UnifracMethod` | Enum: `Unweighted`, `Weighted`, `Generalized(f64)` |
| `UnifracResult` | `distances: Vec<Vec<f64>>`, `sample_names: Vec<String>` |
| `faiths_pd(tree, taxa) -> Result<f64>` | Faith's Phylogenetic Diversity |
| `unweighted_unifrac(tree, sample_a, sample_b) -> Result<f64>` | Unweighted UniFrac distance |
| `weighted_unifrac(tree, sample_a, sample_b) -> Result<f64>` | Weighted UniFrac distance |
| `generalized_unifrac(tree, sample_a, sample_b, alpha) -> Result<f64>` | Generalized UniFrac with parameter alpha |
| `unifrac_matrix(tree, samples, method) -> Result<UnifracResult>` | Pairwise UniFrac distance matrix |

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
| `parallel` | No | Rayon parallelism |

## Dependencies

- `cyanea-core` -- error types
- `cyanea-ml` -- distance matrices (optional, `ml` feature)

## Tests

225 unit tests + 5 doc tests across 24 source files.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 101 | Module declarations, re-exports |
| `tree.rs` | 750 | PhyloTree data structure, traversal, rerooting, subtree extraction |
| `newick.rs` | 285 | Newick parser and writer |
| `nexus.rs` | 314 | NEXUS format parser and writer |
| `distance.rs` | 291 | Evolutionary distance models |
| `compare.rs` | 210 | Robinson-Foulds and branch score comparison |
| `construct.rs` | 425 | UPGMA and neighbor-joining (ml feature) |
| `reconstruct.rs` | 454 | Fitch and Sankoff ancestral reconstruction |
| `marginal.rs` | 240 | Marginal ML ancestral reconstruction |
| `models.rs` | 630 | JC69, GTR, HKY85 substitution models, discrete gamma rates |
| `subst_model.rs` | ~200 | SubstitutionModel trait, Jc69Model, Hky85Model, GtrModel wrappers |
| `protein_models.rs` | ~450 | LG, WAG, JTT, Dayhoff amino acid models |
| `generic_likelihood.rs` | ~200 | N-state Felsenstein pruning (DNA + protein) |
| `likelihood.rs` | 620 | Felsenstein pruning, NNI search, GTR+Gamma likelihood |
| `tree_search.rs` | ~450 | SPR, TBR, parsimony ratchet, stochastic NNI, SPR search |
| `model_selection.rs` | ~200 | AIC, BIC, AICc, LRT, ModelFinder |
| `mcmc.rs` | ~400 | MCMC sampler, Metropolis-Hastings, tree priors, convergence |
| `simulation.rs` | ~400 | Sequence evolution, coalescent simulation |
| `species_tree.rs` | ~300 | ASTRAL species tree, reconciliation, concordance factors |
| `unifrac.rs` | ~350 | UniFrac distances, Faith's PD |
| `bootstrap.rs` | 300 | Bootstrap support and bipartition extraction |
| `consensus.rs` | 380 | Consensus tree construction |
| `dating.rs` | 200 | Molecular clock dating |
| `drawing.rs` | 230 | Tree drawing coordinates (rectangular, cladogram, radial) |
