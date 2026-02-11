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
pub mod distance;
pub mod likelihood;
pub mod models;
pub mod newick;
pub mod nexus;
pub mod reconstruct;
pub mod tree;

#[cfg(feature = "ml")]
pub mod construct;

pub use bootstrap::{bipartitions, bootstrap_support};
pub use compare::{branch_score_distance, robinson_foulds, robinson_foulds_normalized};
pub use distance::{jukes_cantor, kimura_2p, p_distance};
pub use likelihood::{nni_search, tree_likelihood};
pub use models::{jc69_probability, nucleotide_index};
pub use newick::{parse as parse_newick, write as write_newick};
pub use tree::{Node, NodeId, PhyloTree};

#[cfg(feature = "ml")]
pub use construct::{neighbor_joining, upgma};

#[cfg(feature = "ml")]
pub use distance::{sequence_distance_matrix, DistanceModel};
