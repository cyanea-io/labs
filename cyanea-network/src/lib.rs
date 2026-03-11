//! Network and pathway biology for the Cyanea bioinformatics ecosystem.
//!
//! Provides graph data structures, network topology metrics, community detection,
//! PPI network analysis, gene regulatory network inference, and pathway analysis.
//!
//! # Modules
//!
//! - [`graph`] — Core graph (directed/undirected, weighted, attributed)
//! - [`topology`] — Centrality, clustering, shortest paths, components
//! - [`community`] — Louvain modularity optimization, label propagation
//! - [`ppi`] — Protein-protein interaction networks, propagation, hub/bottleneck
//! - [`grn`] — Gene regulatory network inference (correlation, MI, CLR)
//! - [`pathway`] — Pathway representation, topology scoring, crosstalk, GMT I/O
//! - [`formats`] — GraphML, SIF, GEXF parsing/writing
//!
//! # Example
//!
//! ```
//! use cyanea_network::{Graph, graph::GraphType};
//!
//! let mut g = Graph::new(GraphType::Undirected);
//! g.add_node("TP53", "TP53").unwrap();
//! g.add_node("MDM2", "MDM2").unwrap();
//! g.add_node("BAX", "BAX").unwrap();
//! g.add_edge("TP53", "MDM2", 0.95).unwrap();
//! g.add_edge("TP53", "BAX", 0.85).unwrap();
//!
//! assert_eq!(g.node_count(), 3);
//! assert_eq!(g.edge_count(), 2);
//! assert_eq!(g.degree("TP53").unwrap(), 2);
//! ```

pub mod community;
pub mod formats;
pub mod graph;
pub mod grn;
pub mod pathway;
pub mod ppi;
pub mod topology;

// Re-export core graph types
pub use graph::{Edge, Graph, GraphType, Node};

// Re-export topology metrics
pub use topology::{
    average_clustering, betweenness_centrality, bfs_distances, closeness_centrality,
    clustering_coefficient, connected_components, degree_centrality, diameter, pagerank,
    shortest_path, CentralityResult,
};

// Re-export community detection
pub use community::{label_propagation, louvain, modularity, CommunityResult};

// Re-export PPI analysis
pub use ppi::{
    build_ppi_network, hub_bottleneck_analysis, neighborhood_similarity, network_propagation,
    predict_interactions, HubBottleneckResult, PpiEvidence,
};

// Re-export GRN inference
pub use grn::{clr, infer_grn, CorrelationMethod, ExpressionMatrix, GrnResult};

// Re-export pathway analysis
pub use pathway::{
    parse_gmt, pathway_crosstalk, score_pathways_by_topology, write_gmt, CrosstalkResult,
    Pathway, PathwayDatabase, PathwayScore, PathwaySource,
};

// Re-export format I/O
pub use formats::{parse_graphml, parse_sif, write_gexf, write_graphml, write_sif};
