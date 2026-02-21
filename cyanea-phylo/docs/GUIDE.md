# Usage Guide -- cyanea-phylo

Practical examples for phylogenetic analysis.

## Newick/NEXUS Parsing and Writing

```rust
use cyanea_phylo::{PhyloTree, parse_newick, write_newick};

// Parse a Newick string
let tree = PhyloTree::from_newick("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);").unwrap();
println!("Leaves: {:?}", tree.leaf_names());
println!("Total branch length: {:.2}", tree.total_branch_length());

// Serialize back to Newick
let newick_str = tree.to_newick();
println!("Newick: {}", newick_str);
```

```rust
use cyanea_phylo::nexus;

// Parse NEXUS format
let nexus_text = "#NEXUS\nBEGIN TAXA;\n  DIMENSIONS NTAX=3;\n  TAXLABELS A B C;\nEND;\n\
    BEGIN TREES;\n  TREE t1 = ((A:0.1,B:0.2):0.3,C:0.4);\nEND;";
let nexus_file = nexus::parse(nexus_text).unwrap();
println!("Taxa: {:?}", nexus_file.taxa);
println!("Trees: {}", nexus_file.trees.len());

// Write NEXUS
let output = nexus::write(&nexus_file.taxa, &nexus_file.trees);
```

## Distance Matrix Computation

```rust
use cyanea_phylo::distance::{p_distance, jukes_cantor, kimura_2p};

let seq_a = b"ATCGATCG";
let seq_b = b"ATCAATCG";

let p = p_distance(seq_a, seq_b).unwrap();
println!("p-distance: {:.4}", p);

let jc = jukes_cantor(p).unwrap();
println!("JC distance: {:.4}", jc);

let k2p = kimura_2p(1, 0).unwrap(); // 1 transition, 0 transversions
println!("K2P distance: {:.4}", k2p);
```

## UPGMA and Neighbor-Joining Tree Building

Requires the `ml` feature.

```rust
use cyanea_phylo::construct::{upgma, neighbor_joining};
use cyanea_phylo::distance::{sequence_distance_matrix, DistanceModel};

let sequences: Vec<&[u8]> = vec![b"ATCGATCG", b"ATCAATCG", b"ACCGATCG", b"ATCGATCC"];
let names = vec!["Seq1".into(), "Seq2".into(), "Seq3".into(), "Seq4".into()];

let dm = sequence_distance_matrix(&sequences, DistanceModel::JukesCantor).unwrap();

let upgma_tree = upgma(&dm, &names).unwrap();
println!("UPGMA: {}", upgma_tree.to_newick());

let nj_tree = neighbor_joining(&dm, &names).unwrap();
println!("NJ: {}", nj_tree.to_newick());
```

## Bootstrap Analysis

```rust
use cyanea_phylo::{bootstrap_support, PhyloTree};
use cyanea_phylo::construct::neighbor_joining;
use cyanea_phylo::distance::{sequence_distance_matrix, DistanceModel};

let sequences: Vec<&[u8]> = vec![b"ATCGATCG", b"ATCAATCG", b"ACCGATCG", b"ATCGATCC"];
let tree = PhyloTree::from_newick("((Seq1:0.1,Seq2:0.1):0.05,(Seq3:0.1,Seq4:0.1):0.05);").unwrap();

let builder = |seqs: &[&[u8]]| {
    let dm = sequence_distance_matrix(seqs, DistanceModel::JukesCantor).unwrap();
    let names: Vec<String> = (0..seqs.len()).map(|i| format!("Seq{}", i + 1)).collect();
    neighbor_joining(&dm, &names).unwrap()
};

let support = bootstrap_support(&sequences, &tree, builder, 100).unwrap();
println!("Bootstrap support: {:?}", support);
```

## Maximum Likelihood with NNI Search

```rust
use cyanea_phylo::{PhyloTree, tree_likelihood, nni_search};
use cyanea_phylo::models::jc69_probability;

let tree = PhyloTree::from_newick("((A:0.1,B:0.2):0.3,C:0.4);").unwrap();
let sequences: Vec<&[u8]> = vec![b"ATCGATCG", b"ATCAATCG", b"ACCGATCG"];

let ll = tree_likelihood(&tree, &sequences, jc69_probability).unwrap();
println!("Log-likelihood: {:.4}", ll);

let improved = nni_search(&tree, &sequences, jc69_probability).unwrap();
let ll_improved = tree_likelihood(&improved, &sequences, jc69_probability).unwrap();
println!("Improved log-likelihood: {:.4}", ll_improved);
```

## Substitution Models (JC69, HKY85, GTR+G)

```rust
use cyanea_phylo::models::{GtrParams, GammaRates, hky85_params};
use cyanea_phylo::{tree_likelihood_gtr, PhyloTree};

let tree = PhyloTree::from_newick("((A:0.1,B:0.2):0.3,C:0.4);").unwrap();
let sequences: Vec<&[u8]> = vec![b"ATCGATCG", b"ATCAATCG", b"ACCGATCG"];

// HKY85 with kappa=2.0 and empirical frequencies
let freqs = [0.3, 0.2, 0.2, 0.3];
let hky = hky85_params(2.0, freqs).unwrap();
let prob_fn = hky.probability_fn();

// Add discrete gamma rate heterogeneity (4 categories, alpha=0.5)
let gamma = GammaRates::new(0.5, 4).unwrap();

let ll = tree_likelihood_gtr(&tree, &sequences, &prob_fn, &freqs, Some(&gamma)).unwrap();
println!("HKY85+G4 log-likelihood: {:.4}", ll);
```

## Ancestral Reconstruction

```rust
use cyanea_phylo::{PhyloTree, reconstruct};

let tree = PhyloTree::from_newick("((A:0.1,B:0.2):0.3,C:0.4);").unwrap();
let leaf_states: Vec<u8> = vec![0, 0, 1]; // A=0, B=0, C=1

let result = reconstruct::fitch(&tree, &leaf_states).unwrap();
println!("Ancestral states: {:?}", result.states);
println!("Minimum changes: {}", result.n_changes);
```

For marginal ML reconstruction:

```rust
use cyanea_phylo::{marginal_reconstruct, PhyloTree};
use cyanea_phylo::models::jc69_probability;

let tree = PhyloTree::from_newick("((A:0.1,B:0.2):0.3,C:0.4);").unwrap();
let sequences: Vec<&[u8]> = vec![b"ATCGATCG", b"ATCAATCG", b"ACCGATCG"];

let result = marginal_reconstruct(&tree, &sequences, jc69_probability, &[0.25; 4]).unwrap();
for (node_idx, posteriors) in result.posteriors.iter().enumerate() {
    println!("Node {}: {} sites reconstructed", node_idx, posteriors.len());
}
```

## Tree Comparison (Robinson-Foulds)

```rust
use cyanea_phylo::{PhyloTree, robinson_foulds, robinson_foulds_normalized, branch_score_distance};

let tree1 = PhyloTree::from_newick("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);").unwrap();
let tree2 = PhyloTree::from_newick("((A:0.1,C:0.2):0.3,(B:0.4,D:0.5):0.6);").unwrap();

let rf = robinson_foulds(&tree1, &tree2).unwrap();
let rf_norm = robinson_foulds_normalized(&tree1, &tree2).unwrap();
let bsd = branch_score_distance(&tree1, &tree2).unwrap();

println!("RF distance: {}", rf);
println!("Normalized RF: {:.4}", rf_norm);
println!("Branch score distance: {:.4}", bsd);
```

## Bayesian MCMC

```rust
use cyanea_phylo::{PhyloTree, mcmc_sample, convergence_diagnostics, posterior_summary};
use cyanea_phylo::mcmc::{McmcConfig, TreePrior, ClockModel};
use cyanea_phylo::subst_model::Jc69Model;

let tree = PhyloTree::from_newick("((A:0.1,B:0.2):0.3,C:0.4);").unwrap();
let sequences: Vec<&[u8]> = vec![b"ATCGATCG", b"ATCAATCG", b"ACCGATCG"];
let model = Jc69Model::new();

let config = McmcConfig {
    n_generations: 5000,
    sample_every: 50,
    burnin: 500,
    ..Default::default()
};

let result = mcmc_sample(
    &tree, &sequences, &model,
    &TreePrior::CoalescentConstant { pop_size: 1000.0 },
    &ClockModel::Strict { rate: 0.01 },
    &config,
).unwrap();

let diag = convergence_diagnostics(&result.samples);
println!("ESS: {:?}", diag.ess);

let summary = posterior_summary(&result.samples);
println!("MAP tree: {}", summary.map_tree.to_newick());
```

## UniFrac for Microbiome Analysis

```rust
use std::collections::HashMap;
use cyanea_phylo::{PhyloTree, faiths_pd, unweighted_unifrac, weighted_unifrac, unifrac_matrix};
use cyanea_phylo::unifrac::UnifracMethod;

let tree = PhyloTree::from_newick("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);").unwrap();

// Faith's Phylogenetic Diversity
let taxa: std::collections::HashSet<String> = ["A", "B", "C"].iter().map(|s| s.to_string()).collect();
let pd = faiths_pd(&tree, &taxa).unwrap();
println!("Faith's PD: {:.4}", pd);

// Pairwise UniFrac
let sample_a: HashMap<String, f64> = [("A", 10.0), ("B", 5.0)].iter().map(|&(k, v)| (k.into(), v)).collect();
let sample_b: HashMap<String, f64> = [("C", 8.0), ("D", 3.0)].iter().map(|&(k, v)| (k.into(), v)).collect();

let uw = unweighted_unifrac(&tree, &sample_a, &sample_b).unwrap();
let ww = weighted_unifrac(&tree, &sample_a, &sample_b).unwrap();
println!("Unweighted UniFrac: {:.4}", uw);
println!("Weighted UniFrac: {:.4}", ww);
```

## Species Tree Estimation

```rust
use cyanea_phylo::{PhyloTree, astral_species_tree, reconcile, concordance_factors};

let gene_tree1 = PhyloTree::from_newick("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);").unwrap();
let gene_tree2 = PhyloTree::from_newick("((A:0.1,C:0.2):0.3,(B:0.4,D:0.5):0.6);").unwrap();
let gene_trees = vec![gene_tree1.clone(), gene_tree2.clone()];

// ASTRAL-style species tree
let species_tree = astral_species_tree(&gene_trees).unwrap();
println!("Species tree: {}", species_tree.to_newick());

// Gene-species reconciliation
let recon = reconcile(&gene_tree1, &species_tree).unwrap();
println!("Duplications: {}, Losses: {}", recon.duplications.len(), recon.losses);
```

## Simulation

```rust
use cyanea_phylo::{PhyloTree, simulate_evolution, simulate_coalescent};
use cyanea_phylo::subst_model::Jc69Model;

// Simulate a coalescent tree
let tree = simulate_coalescent(20, 10000.0, 42).unwrap();
println!("Coalescent tree with {} leaves", tree.leaf_count());

// Simulate sequence evolution along the tree
let model = Jc69Model::new();
let alignment = simulate_evolution(&tree, &model, 500, 42).unwrap();
println!("Simulated {} sequences of length {}",
    alignment.names.len(), alignment.sequences[0].len());
println!("Total substitutions: {}", alignment.n_substitutions);
```

## Model Selection

```rust
use cyanea_phylo::{PhyloTree, model_finder, aic, bic, lrt};
use cyanea_phylo::subst_model::{Jc69Model, Hky85Model, SubstitutionModel};
use cyanea_phylo::models::GammaRates;

let tree = PhyloTree::from_newick("((A:0.1,B:0.2):0.3,C:0.4);").unwrap();
let sequences: Vec<&[u8]> = vec![b"ATCGATCG", b"ATCAATCG", b"ACCGATCG"];

let jc = Jc69Model::new();
let hky = Hky85Model::new(2.0, [0.3, 0.2, 0.2, 0.3]).unwrap();

let candidates: Vec<(&str, &dyn SubstitutionModel)> = vec![
    ("JC69", &jc),
    ("HKY85", &hky),
];

let result = model_finder(&tree, &sequences, &candidates, None).unwrap();
println!("Best model (AIC): {}", result.best_aic);
println!("Best model (BIC): {}", result.best_bic);
```
