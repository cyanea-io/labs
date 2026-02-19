//! Phylogenetics and trees for the Cyanea bioinformatics ecosystem.
//!
//! Provides phylogenetic tree data structures, Newick format I/O,
//! evolutionary distance models, distance-based tree construction,
//! tree comparison metrics, substitution models, maximum likelihood
//! inference, and bootstrap support.
//!
//! # Quick start
//!
//! ```
//! use cyanea_phylo::PhyloTree;
//!
//! let tree = PhyloTree::from_newick("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);").unwrap();
//! assert_eq!(tree.leaf_count(), 4);
//! assert_eq!(tree.leaf_names(), vec!["A", "B", "C", "D"]);
//! ```

pub mod bootstrap;
pub mod compare;
pub mod consensus;
pub mod dating;
pub mod distance;
pub mod drawing;
pub mod generic_likelihood;
pub mod likelihood;
pub mod marginal;
pub mod mcmc;
pub mod model_selection;
pub mod models;
pub mod newick;
pub mod nexus;
pub mod protein_models;
pub mod reconstruct;
pub mod species_tree;
pub mod subst_model;
pub mod tree;
pub mod tree_search;
pub mod unifrac;

#[cfg(feature = "ml")]
pub mod construct;

pub use bootstrap::{bipartitions, bootstrap_support};
pub use compare::{branch_score_distance, robinson_foulds, robinson_foulds_normalized};
pub use consensus::{bipartition_frequencies, consensus_tree, ConsensusType, SupportedBipartition};
pub use dating::{root_to_tip_regression, strict_clock, Calibration, DatingResult};
pub use distance::{jukes_cantor, kimura_2p, p_distance};
pub use drawing::{tree_layout, Edge, LayoutStyle, NodeCoord, TreeLayout};
pub use likelihood::{nni_search, tree_likelihood, tree_likelihood_gtr};
pub use marginal::{marginal_reconstruct, MarginalPosterior, MarginalReconstruction};
pub use models::{gtr_probability, hky85_params, jc69_probability, nucleotide_index, GammaRates, GtrParams};
pub use newick::{parse as parse_newick, write as write_newick};
pub use tree::{Node, NodeId, PhyloTree};
pub use unifrac::{
    faiths_pd, generalized_unifrac, unweighted_unifrac, weighted_unifrac, unifrac_matrix,
    UnifracMethod, UnifracResult,
};

// Re-export substitution model trait and wrappers
pub use subst_model::{SubstitutionModel, Jc69Model, Hky85Model, GtrModel};

// Re-export protein models
pub use protein_models::{
    amino_acid_index, load_aa_model, CustomAaModel, DayhoffModel, JttModel, LgModel, WagModel,
    AA_STATES,
};

// Re-export generic likelihood
pub use generic_likelihood::{generic_tree_likelihood, site_likelihoods};

// Re-export tree search
pub use tree_search::{
    nni_search_generic, parsimony_ratchet, spr_move, spr_search, stochastic_nni, tbr_move,
    AnnealingConfig,
};

// Re-export model selection
pub use model_selection::{aic, aicc, bic, lrt, model_finder, LrtResult, ModelResult, ModelSelectionResult};

// Re-export MCMC
pub use mcmc::{
    convergence_diagnostics, mcmc_sample, posterior_summary, ClockModel, ConvergenceDiag,
    McmcConfig, McmcResult, McmcSample, PosteriorSummary, ProposalWeights, TreePrior,
};

// Re-export species tree
pub use species_tree::{
    astral_species_tree, concordance_factors, reconcile, ConcordanceFactors,
    ReconciliationResult,
};

#[cfg(feature = "ml")]
pub use construct::{neighbor_joining, upgma};

#[cfg(feature = "ml")]
pub use distance::{sequence_distance_matrix, DistanceModel};
