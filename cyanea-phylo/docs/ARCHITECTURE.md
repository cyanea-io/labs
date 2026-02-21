# Architecture -- cyanea-phylo

Internal design documentation for the phylogenetics crate.

## PhyloTree Data Structure

Arena-allocated rooted tree using `Vec<Node>` with integer indices:

```
PhyloTree {
    nodes: Vec<Node>,     // arena -- all nodes stored in a flat vector
    root: NodeId,         // index of the root node
}

Node {
    id: NodeId,           // index into the arena (== position in Vec)
    parent: Option<NodeId>,
    children: Vec<NodeId>,
    branch_length: Option<f64>,
    name: Option<String>,
}
```

Leaves are nodes with empty `children`. The root has `parent = None`. This arena layout avoids reference-counted pointers and enables O(1) node access by index.

Tree traversal is implemented via iterators:
- `PreorderIter`: uses an explicit stack, visits parent before children
- `PostorderIter`: two-pass (reverse postorder via stack, then reverse), visits children before parent

Rerooting creates a new tree by reversing parent-child relationships along the path from the old root to the new root edge.

## Newick Parser (`newick.rs`)

Recursive descent parser for the Newick format:

- Handles nested parentheses for subtrees
- Parses optional labels after closing parentheses
- Parses optional branch lengths after `:` separator
- Bootstrap values read from internal node labels
- Supports quoted labels and whitespace tolerance
- Round-trip safe: `parse(write(tree))` preserves topology and branch lengths

## NEXUS Parser (`nexus.rs`)

Block-structured format parser:

- Scans for `BEGIN ... END;` blocks (case-insensitive)
- TAXA block: reads `DIMENSIONS NTAX=n` and `TAXLABELS`
- TREES block: reads `TREE name = newick;` entries, delegates to Newick parser
- Ignores unknown blocks (DATA, CHARACTERS, etc.) gracefully

## Distance Matrix

O(n^2) pairwise computation with evolutionary model correction:

- **p-distance**: raw proportion of differing sites (gaps excluded)
- **Jukes-Cantor**: `-3/4 * ln(1 - 4p/3)`, undefined when p >= 0.75
- **Kimura 2-parameter**: separate treatment of transitions and transversions

The `sequence_distance_matrix` function (ml feature) returns a flat `Vec<f64>` of n*(n-1)/2 pairwise distances for use by UPGMA and NJ.

## Felsenstein Pruning (`likelihood.rs`, `generic_likelihood.rs`)

Postorder traversal computing conditional likelihood vectors at each node:

1. Leaf nodes: set partial likelihood to 1.0 for the observed state, 0.0 for others
2. Internal nodes: for each site and each state i, multiply over children: `L_i = product_c(sum_j(P(i->j, t_c) * L_j_c))`
3. Root: sum over states weighted by equilibrium frequencies: `L_site = sum_i(pi_i * L_i_root)`
4. Total log-likelihood: sum of `ln(L_site)` over all sites

The 4-state version (`likelihood.rs`) uses `[[f64; 4]; 4]` arrays for speed. The generic version (`generic_likelihood.rs`) uses `Vec<Vec<f64>>` to support both 4-state (DNA) and 20-state (protein) models.

Discrete gamma rate heterogeneity: evaluates the pruning algorithm once per rate category, averages the site likelihoods, then takes the log.

## NNI / SPR / TBR Tree Rearrangements

**NNI** (Nearest Neighbor Interchange):
- Identifies internal edges with two subtrees on each side
- Swaps one subtree from each side to produce two alternative topologies
- Hill-climbing: accept if log-likelihood improves

**SPR** (Subtree Pruning and Regrafting):
- Detaches a subtree at the prune node, collapses the resulting degree-2 node
- Reattaches the subtree at a new edge (bisecting it with a new internal node)
- Broader search radius than NNI

**TBR** (Tree Bisection and Reconnection):
- Bisects the tree by removing an edge, creating two subtrees
- Reconnects by choosing one edge in each subtree
- Most general rearrangement; subsumes SPR and NNI

`parsimony_ratchet`: iteratively reweights sites, performs NNI, accepts/rejects via parsimony score. Effective at escaping local optima.

`stochastic_nni`: simulated annealing with NNI proposals, temperature-controlled acceptance of worse topologies.

## MCMC Sampler (`mcmc.rs`)

Metropolis-Hastings algorithm with mixed proposal kernel:

1. **Proposal selection**: weighted random choice among NNI, SPR, branch scaling, branch sliding, model parameters
2. **NNI/SPR proposals**: topology changes with Hastings ratio = 1 (symmetric)
3. **Branch scaling**: multiply all branch lengths by exp(U(-delta, delta)), Hastings ratio accounts for the Jacobian
4. **Branch sliding**: perturb a single branch length, symmetric uniform proposal
5. **Acceptance**: `min(1, (L_new * prior_new / L_old * prior_old) * hastings_ratio)`

Tree priors:
- Uniform: flat prior on topologies
- Coalescent (constant): exponential waiting times with rate k*(k-1)/(4*N)
- Coalescent (exponential growth): integrated rate with growth parameter
- Birth-death: speciation/extinction process prior on branch lengths

Clock models:
- Strict: single rate for all branches
- Uncorrelated lognormal: per-branch rates drawn from LogNormal(mean_rate, stdev)

Convergence: ESS estimated via autocorrelation time, mean/variance of log-likelihood trace.

## Coalescent Simulation (`simulation.rs`)

Backward-in-time simulation of gene genealogies:

1. Start with n lineages at time 0
2. Sample exponential waiting time with rate k*(k-1)/(4*N) where k is the current number of lineages
3. Pick two random lineages to coalesce, create a new internal node
4. Set branch lengths from coalescence times
5. Repeat until one lineage remains (the root)

For exponential growth, the waiting time is computed by inverting the cumulative hazard function analytically.

Sequence evolution (`simulate_evolution`) works top-down:
1. Generate root sequence from equilibrium frequencies
2. Pre-order traversal: for each child, compute P(t) from the model, sample child state for each site
3. Collect leaf sequences as the simulated alignment

## UniFrac (`unifrac.rs`)

Phylogenetic beta-diversity metrics computed by tree traversal:

**Unweighted**: postorder propagation of presence/absence flags from leaves to root. Accumulates branch length where exactly one sample is present (unique) vs. either sample (total). Distance = unique/total.

**Weighted**: propagates proportional abundances (normalized to sum=1 per sample). Each branch contributes `bl * |p_A - p_B|` to the numerator, `bl * (p_A + p_B)` to the denominator.

**Generalized**: parameterized by alpha in [0,1]. Alpha=0 gives unweighted-like behavior, alpha=1 gives weighted. Intermediate values control sensitivity to abundant vs. rare taxa.

**Faith's PD**: sum of branch lengths in the minimum spanning subtree connecting observed taxa to the root. Postorder mark-and-sweep of observed taxon ancestors.

## Species Tree (`species_tree.rs`)

**ASTRAL-style estimation**: enumerates quartets across gene trees, counts support for each resolution, builds species tree greedily by majority quartet topology. Practical for up to ~50 taxa (quartet enumeration is O(n^4) but capped at 5000 per gene tree).

**Reconciliation**: maps gene tree nodes onto species tree via LCA mapping. Counts duplications (gene node maps to same species node as a child) and losses (missing lineages in species tree branches).

**Concordance factors**: gene CF is the fraction of gene trees supporting each species tree branch. Site CF is computed from per-site likelihoods under each possible resolution.

## Module Dependencies

```
tree  <--  all modules
newick  <--  nexus, bootstrap, consensus, simulation
models  <--  likelihood, generic_likelihood, mcmc, tree_search, model_selection, simulation
subst_model  <--  protein_models, generic_likelihood, mcmc, tree_search, model_selection, simulation
generic_likelihood  <--  tree_search, model_selection, mcmc
likelihood  <--  tree_search, mcmc
bootstrap  <--  mcmc, species_tree
compare  <--  (standalone)
distance  <--  construct
```

All modules depend on `cyanea-core` for `CyaneaError` and `Result`.
