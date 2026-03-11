# cyanea-network Usage Guide

## Overview

`cyanea-network` provides graph data structures and algorithms for biological network analysis: building and querying graphs, computing topology metrics, detecting communities, analyzing protein-protein interaction networks, inferring gene regulatory networks, and scoring pathways.

## 1. Building Graphs

### Create a simple undirected graph

```rust
use cyanea_network::{Graph, graph::GraphType};

let mut g = Graph::new(GraphType::Undirected);
g.add_node("TP53", "TP53").unwrap();
g.add_node("MDM2", "MDM2").unwrap();
g.add_node("BAX", "BAX").unwrap();
g.add_edge("TP53", "MDM2", 0.95).unwrap();
g.add_edge("TP53", "BAX", 0.85).unwrap();

assert_eq!(g.node_count(), 3);
assert_eq!(g.edge_count(), 2);
assert_eq!(g.degree("TP53").unwrap(), 2);
```

### Nodes and edges with attributes

```rust
use cyanea_network::graph::{Node, Edge, GraphType};
use cyanea_network::Graph;

let mut g = Graph::new(GraphType::Directed);

let node = Node::new("TP53", "TP53")
    .with_attr("type", "tumor_suppressor")
    .with_attr("organism", "human");
g.add_node_struct(node).unwrap();
g.add_node("MDM2", "MDM2").unwrap();

let edge = Edge::weighted("TP53", "MDM2", 0.95)
    .with_type("activation")
    .with_attr("source", "STRING");
g.add_edge_struct(edge).unwrap();
```

### Subgraphs and properties

```rust
use cyanea_network::{Graph, graph::GraphType};

let mut g = Graph::new(GraphType::Undirected);
for id in &["A", "B", "C", "D"] {
    g.add_node(id, id).unwrap();
}
g.add_edge("A", "B", 1.0).unwrap();
g.add_edge("B", "C", 1.0).unwrap();
g.add_edge("C", "D", 1.0).unwrap();
g.add_edge("A", "D", 1.0).unwrap();

println!("Density: {:.3}", g.density());

let sub = g.subgraph(&["A", "B", "C"]).unwrap();
assert_eq!(sub.node_count(), 3);
assert_eq!(sub.edge_count(), 2); // A-B and B-C

let (ids, matrix) = g.adjacency_matrix();
println!("Adjacency matrix: {:?}", matrix);
```

## 2. Topology Analysis

### Centrality measures

```rust
use cyanea_network::{
    Graph, graph::GraphType,
    degree_centrality, closeness_centrality, betweenness_centrality, pagerank,
};

let mut g = Graph::new(GraphType::Undirected);
for id in &["A", "B", "C", "D"] {
    g.add_node(id, id).unwrap();
}
g.add_edge("A", "B", 1.0).unwrap();
g.add_edge("B", "C", 1.0).unwrap();
g.add_edge("C", "D", 1.0).unwrap();

// Degree centrality
let dc = degree_centrality(&g).unwrap();
let top = dc.top_k(2);
println!("Top 2 by degree: {:?}", top);

// Closeness centrality
let cc = closeness_centrality(&g).unwrap();
println!("Closeness of B: {:.3}", cc.scores["B"]);

// Betweenness centrality
let bc = betweenness_centrality(&g).unwrap();
println!("Betweenness of B: {:.3}", bc.scores["B"]);

// PageRank
let pr = pagerank(&g, 0.85, 100, 1e-6).unwrap();
for (node, score) in pr.top_k(3) {
    println!("{}: {:.4}", node, score);
}
```

### Clustering, paths, and components

```rust
use cyanea_network::{
    Graph, graph::GraphType,
    clustering_coefficient, average_clustering,
    shortest_path, bfs_distances, diameter, connected_components,
};

let mut g = Graph::new(GraphType::Undirected);
for id in &["A", "B", "C"] {
    g.add_node(id, id).unwrap();
}
g.add_edge("A", "B", 1.0).unwrap();
g.add_edge("B", "C", 1.0).unwrap();
g.add_edge("A", "C", 1.0).unwrap();

// Clustering coefficient (triangle = 1.0)
let cc = clustering_coefficient(&g, "A").unwrap();
println!("Clustering(A): {:.2}", cc);

let avg = average_clustering(&g).unwrap();
println!("Average clustering: {:.2}", avg);

// Shortest paths
let path = shortest_path(&g, "A", "C").unwrap().unwrap();
println!("Path A->C: {:?}", path);

let dists = bfs_distances(&g, "A");
println!("Distances from A: {:?}", dists);

// Diameter
let d = diameter(&g).unwrap();
println!("Diameter: {:?}", d);

// Connected components
let comps = connected_components(&g);
println!("Components: {}", comps.len());
```

## 3. Community Detection

### Louvain modularity optimization

```rust
use cyanea_network::{Graph, graph::GraphType, louvain};

// Two cliques connected by a bridge
let mut g = Graph::new(GraphType::Undirected);
for id in &["A", "B", "C", "D", "E", "F"] {
    g.add_node(id, id).unwrap();
}
// Clique 1
g.add_edge("A", "B", 1.0).unwrap();
g.add_edge("B", "C", 1.0).unwrap();
g.add_edge("A", "C", 1.0).unwrap();
// Clique 2
g.add_edge("D", "E", 1.0).unwrap();
g.add_edge("E", "F", 1.0).unwrap();
g.add_edge("D", "F", 1.0).unwrap();
// Bridge
g.add_edge("C", "D", 1.0).unwrap();

let result = louvain(&g, 1.0).unwrap();
println!("Communities: {}, Modularity: {:.3}",
    result.num_communities, result.modularity);

for (community, size) in result.community_sizes() {
    let members = result.members(community);
    println!("  Community {}: {} members — {:?}", community, size, members);
}
```

### Label propagation

```rust
use cyanea_network::{Graph, graph::GraphType, label_propagation, modularity};

let mut g = Graph::new(GraphType::Undirected);
// ... build graph ...
# for id in &["A", "B", "C"] { g.add_node(id, id).unwrap(); }
# g.add_edge("A", "B", 1.0).unwrap();
# g.add_edge("B", "C", 1.0).unwrap();

let result = label_propagation(&g, 100).unwrap();
println!("Communities: {}, Modularity: {:.3}",
    result.num_communities, result.modularity);
```

## 4. PPI Network Analysis

### Build a PPI network from evidence

```rust
use cyanea_network::ppi::{PpiEvidence, build_ppi_network};

let mut ev1 = PpiEvidence::new("TP53", "MDM2");
ev1.experimental = 0.95;
ev1.database = 0.9;

let mut ev2 = PpiEvidence::new("TP53", "BAX");
ev2.experimental = 0.8;

let mut ev3 = PpiEvidence::new("LOW1", "LOW2");
ev3.textmining = 0.1; // Low confidence — will be filtered

let graph = build_ppi_network(&[ev1, ev2, ev3], 0.5).unwrap();
println!("Nodes: {}, Edges: {}", graph.node_count(), graph.edge_count());
// LOW1-LOW2 excluded (combined score 0.1 < 0.5)
```

### Network propagation (guilt-by-association)

```rust
use cyanea_network::{Graph, graph::GraphType};
use cyanea_network::ppi::network_propagation;
use std::collections::HashMap;

let mut g = Graph::new(GraphType::Undirected);
for id in &["BRCA1", "TP53", "MDM2", "BAX", "BCL2"] {
    g.add_node(id, id).unwrap();
}
g.add_edge("BRCA1", "TP53", 0.9).unwrap();
g.add_edge("TP53", "MDM2", 0.95).unwrap();
g.add_edge("TP53", "BAX", 0.85).unwrap();
g.add_edge("BAX", "BCL2", 0.9).unwrap();

// Seed BRCA1, propagate through network
let seeds: HashMap<String, f64> =
    vec![("BRCA1".into(), 1.0)].into_iter().collect();
let scores = network_propagation(&g, &seeds, 0.3, 100, 1e-8).unwrap();

let mut ranked: Vec<_> = scores.iter().collect();
ranked.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap());
for (gene, score) in &ranked {
    println!("{}: {:.4}", gene, score);
}
```

### Hub and bottleneck analysis

```rust
use cyanea_network::{Graph, graph::GraphType};
use cyanea_network::ppi::hub_bottleneck_analysis;

let mut g = Graph::new(GraphType::Undirected);
// ... build PPI network ...
# g.add_node("center", "center").unwrap();
# for i in 0..5 { let id = format!("n{}", i); g.add_node(&id, &id).unwrap(); g.add_edge("center", &id, 1.0).unwrap(); }

let result = hub_bottleneck_analysis(&g, 0.8, 0.8).unwrap();
println!("Hubs: {:?}", result.hubs);
println!("Bottlenecks: {:?}", result.bottlenecks);
println!("Hub-bottlenecks: {:?}", result.hub_bottlenecks);
```

### Link prediction

```rust
use cyanea_network::{Graph, graph::GraphType};
use cyanea_network::ppi::{neighborhood_similarity, predict_interactions};

let mut g = Graph::new(GraphType::Undirected);
for id in &["A", "B", "C", "D"] {
    g.add_node(id, id).unwrap();
}
g.add_edge("A", "C", 1.0).unwrap();
g.add_edge("A", "D", 1.0).unwrap();
g.add_edge("B", "C", 1.0).unwrap();
g.add_edge("B", "D", 1.0).unwrap();

// Jaccard similarity between neighborhoods
let sim = neighborhood_similarity(&g, "A", "B").unwrap();
println!("Neighborhood similarity A-B: {:.3}", sim);

// Top predicted interactions
let preds = predict_interactions(&g, 5).unwrap();
for (a, b, score) in &preds {
    println!("Predicted: {} - {} (Jaccard={:.3})", a, b, score);
}
```

## 5. GRN Inference

### Correlation-based inference

```rust
use cyanea_network::grn::{ExpressionMatrix, infer_grn, CorrelationMethod};

let expr = ExpressionMatrix::new(
    vec!["TF1".into(), "TF2".into(), "G1".into(), "G2".into(), "G3".into()],
    (0..10).map(|i| format!("S{}", i)).collect(),
    vec![
        vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0], // TF1
        vec![10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0], // TF2
        vec![1.1, 2.2, 2.8, 4.1, 5.2, 5.9, 7.1, 8.0, 8.8, 10.2], // G1
        vec![0.9, 1.8, 3.2, 3.9, 4.8, 6.1, 7.2, 7.9, 9.1, 9.8],  // G2
        vec![9.8, 9.1, 7.9, 7.2, 6.1, 4.8, 3.9, 3.2, 1.8, 0.9],  // G3
    ],
).unwrap();

// TF indices: TF1=0, TF2=1
let result = infer_grn(&expr, &[0, 1], CorrelationMethod::Pearson, 0.9).unwrap();

println!("Method: {}", result.method);
println!("Edges: {}", result.network.edge_count());
for ((tf, target), score) in &result.scores {
    println!("  {} -> {}: {:.3}", tf, target, score);
}
```

### CLR (Context Likelihood of Relatedness)

```rust
use cyanea_network::grn::{ExpressionMatrix, clr};

// ... build expression matrix as above ...
# let expr = cyanea_network::grn::ExpressionMatrix::new(
#     vec!["TF1".into(), "TF2".into(), "G1".into(), "G2".into(), "G3".into()],
#     (0..10).map(|i| format!("S{}", i)).collect(),
#     vec![
#         vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
#         vec![10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0],
#         vec![1.1, 2.2, 2.8, 4.1, 5.2, 5.9, 7.1, 8.0, 8.8, 10.2],
#         vec![0.9, 1.8, 3.2, 3.9, 4.8, 6.1, 7.2, 7.9, 9.1, 9.8],
#         vec![9.8, 9.1, 7.9, 7.2, 6.1, 4.8, 3.9, 3.2, 1.8, 0.9],
#     ],
# ).unwrap();

let result = clr(&expr, &[0, 1], 1.0).unwrap();
println!("CLR edges: {}", result.network.edge_count());
for ((tf, target), score) in &result.scores {
    println!("  {} -> {}: CLR={:.3}", tf, target, score);
}
```

## 6. Pathway Analysis

### Load pathways and score by topology

```rust
use cyanea_network::{Graph, graph::GraphType};
use cyanea_network::pathway::{
    Pathway, PathwayDatabase, PathwaySource,
    score_pathways_by_topology, pathway_crosstalk,
};

// Build a PPI network
let mut g = Graph::new(GraphType::Undirected);
for id in &["CDK1", "CDK2", "CCNB1", "TP53", "BAX", "BCL2", "CASP3"] {
    g.add_node(id, id).unwrap();
}
g.add_edge("CDK1", "CDK2", 0.9).unwrap();
g.add_edge("CDK1", "CCNB1", 0.95).unwrap();
g.add_edge("CDK2", "TP53", 0.8).unwrap();
g.add_edge("TP53", "BAX", 0.85).unwrap();
g.add_edge("BAX", "BCL2", 0.9).unwrap();
g.add_edge("BAX", "CASP3", 0.8).unwrap();

// Build pathway database
let mut db = PathwayDatabase::new();
db.add(Pathway::new("PW001", "Cell Cycle", PathwaySource::KEGG,
    vec!["CDK1".into(), "CDK2".into(), "CCNB1".into(), "TP53".into()]));
db.add(Pathway::new("PW002", "Apoptosis", PathwaySource::KEGG,
    vec!["TP53".into(), "BAX".into(), "BCL2".into(), "CASP3".into()]));

// Topology scoring
let scores = score_pathways_by_topology(&g, &db).unwrap();
for s in &scores {
    println!("{} ({}): topology={:.3}, density={:.3}, centrality={:.3}",
        s.pathway_name, s.pathway_id,
        s.topology_score, s.internal_density, s.avg_centrality);
}

// Pathway crosstalk
let crosstalk = pathway_crosstalk(&g, &db, 1).unwrap();
for ct in &crosstalk {
    println!("{} <-> {}: shared={}, Jaccard={:.3}, cross-edges={}",
        ct.pathway_a, ct.pathway_b,
        ct.shared_members, ct.jaccard, ct.cross_edges);
}
```

### GMT format I/O

```rust
use cyanea_network::pathway::{parse_gmt, write_gmt};

let gmt_text = "PW001\tCell Cycle\tCDK1\tCDK2\tCCNB1\n\
                PW002\tApoptosis\tTP53\tBAX\tBCL2\n";

let db = parse_gmt(gmt_text).unwrap();
println!("Loaded {} pathways", db.len());

let pw = &db.pathways["PW001"];
println!("{}: {} members", pw.name, pw.size());

// Write back to GMT
let output = write_gmt(&db);
println!("{}", output);
```

## 7. Format I/O

### SIF (Cytoscape)

```rust
use cyanea_network::formats::{parse_sif, write_sif};

let sif = "TP53\tpp\tMDM2\nTP53\tpp\tBAX\nMDM2\tpp\tMDMX\n";
let g = parse_sif(sif).unwrap();
println!("Nodes: {}, Edges: {}", g.node_count(), g.edge_count());

let output = write_sif(&g);
println!("{}", output);
```

### GraphML

```rust
use cyanea_network::formats::{parse_graphml, write_graphml};

let graphml = r#"<?xml version="1.0" encoding="UTF-8"?>
<graphml>
  <graph edgedefault="undirected">
    <node id="TP53"/>
    <node id="MDM2"/>
    <edge source="TP53" target="MDM2"/>
  </graph>
</graphml>"#;

let g = parse_graphml(graphml).unwrap();
println!("Nodes: {}, Type: {:?}", g.node_count(), g.graph_type);

// Write to GraphML
let xml = write_graphml(&g);
```

### GEXF (Gephi)

```rust
use cyanea_network::{Graph, graph::GraphType};
use cyanea_network::formats::write_gexf;

let mut g = Graph::new(GraphType::Directed);
g.add_node("TF1", "TF1").unwrap();
g.add_node("G1", "Gene 1").unwrap();
g.add_edge("TF1", "G1", 0.9).unwrap();

let gexf = write_gexf(&g);
println!("{}", gexf);
```

## 8. Complete Pipeline

```rust
use cyanea_network::{
    Graph, graph::GraphType, louvain, betweenness_centrality,
};
use cyanea_network::ppi::{PpiEvidence, build_ppi_network, hub_bottleneck_analysis};
use cyanea_network::pathway::{
    PathwayDatabase, Pathway, PathwaySource,
    parse_gmt, score_pathways_by_topology,
};
use cyanea_network::formats::write_graphml;

// 1. Build PPI network from evidence
let evidence = vec![
    { let mut e = PpiEvidence::new("TP53", "MDM2"); e.experimental = 0.95; e },
    { let mut e = PpiEvidence::new("TP53", "BAX"); e.experimental = 0.8; e.database = 0.7; e },
    { let mut e = PpiEvidence::new("BAX", "BCL2"); e.experimental = 0.9; e },
    { let mut e = PpiEvidence::new("MDM2", "MDMX"); e.experimental = 0.85; e },
];
let graph = build_ppi_network(&evidence, 0.5).unwrap();

// 2. Identify hubs and bottlenecks
let hb = hub_bottleneck_analysis(&graph, 0.75, 0.75).unwrap();

// 3. Detect communities
let communities = louvain(&graph, 1.0).unwrap();

// 4. Score pathways
let mut db = PathwayDatabase::new();
db.add(Pathway::new("PW001", "Apoptosis", PathwaySource::KEGG,
    vec!["TP53".into(), "BAX".into(), "BCL2".into()]));
let pathway_scores = score_pathways_by_topology(&graph, &db).unwrap();

// 5. Export
let graphml = write_graphml(&graph);
```
