# cyanea-network API Reference

## Types

### `Graph`
Adjacency-list graph for biological networks. Supports directed and undirected graphs with weighted edges and attributes.

| Field | Type | Description |
|-------|------|-------------|
| `graph_type` | `GraphType` | Directed or undirected |
| `nodes` | `HashMap<String, Node>` | Node ID to Node |
| `edges` | `Vec<Edge>` | All edges (canonical storage) |
| `attributes` | `HashMap<String, String>` | Graph-level attributes |

**Methods:**

| Method | Signature | Description |
|--------|-----------|-------------|
| `new` | `(GraphType) -> Self` | Create empty graph |
| `node_count` | `() -> usize` | Number of nodes |
| `edge_count` | `() -> usize` | Number of edges (undirected counted once) |
| `add_node` | `(id, label) -> Result<()>` | Add node by ID and label |
| `add_node_struct` | `(Node) -> Result<()>` | Add node from struct |
| `add_edge` | `(source, target, weight) -> Result<()>` | Add weighted edge |
| `add_edge_struct` | `(Edge) -> Result<()>` | Add edge from struct |
| `neighbors` | `(node_id) -> Result<Vec<(&str, f64)>>` | Neighbors as (id, weight) pairs |
| `degree` | `(node_id) -> Result<usize>` | Number of edges |
| `in_degree` | `(node_id) -> Result<usize>` | In-degree (directed) or degree (undirected) |
| `out_degree` | `(node_id) -> Result<usize>` | Out-degree (directed) or degree (undirected) |
| `has_edge` | `(source, target) -> bool` | Check if edge exists |
| `edge_weight` | `(source, target) -> Option<f64>` | Get edge weight |
| `node_ids` | `() -> Vec<&str>` | All node IDs |
| `get_node` | `(id) -> Option<&Node>` | Get node by ID |
| `remove_node` | `(id) -> Result<()>` | Remove node and its edges |
| `subgraph` | `(node_ids) -> Result<Graph>` | Extract induced subgraph |
| `density` | `() -> f64` | Graph density |
| `adjacency_matrix` | `() -> (Vec<String>, Vec<Vec<f64>>)` | Sorted-ID adjacency matrix |

### `Node`
| Field | Type | Description |
|-------|------|-------------|
| `id` | `String` | Unique identifier |
| `label` | `String` | Human-readable label (e.g., gene symbol) |
| `attributes` | `HashMap<String, String>` | Key-value attributes |

**Methods:**

| Method | Signature | Description |
|--------|-----------|-------------|
| `new` | `(id, label) -> Self` | Create node |
| `with_attr` | `(self, key, value) -> Self` | Builder: add attribute |

### `Edge`
| Field | Type | Description |
|-------|------|-------------|
| `source` | `String` | Source node ID |
| `target` | `String` | Target node ID |
| `weight` | `f64` | Edge weight (default 1.0) |
| `edge_type` | `Option<String>` | Interaction type (e.g., "activation", "binding") |
| `attributes` | `HashMap<String, String>` | Key-value attributes |

**Methods:**

| Method | Signature | Description |
|--------|-----------|-------------|
| `new` | `(source, target) -> Self` | Create unweighted edge |
| `weighted` | `(source, target, weight) -> Self` | Create weighted edge |
| `with_type` | `(self, edge_type) -> Self` | Builder: set interaction type |
| `with_attr` | `(self, key, value) -> Self` | Builder: add attribute |

### `GraphType`
`Directed`, `Undirected`

### `CentralityResult`
| Field | Type | Description |
|-------|------|-------------|
| `scores` | `HashMap<String, f64>` | Node ID to centrality score |

**Methods:**

| Method | Signature | Description |
|--------|-----------|-------------|
| `top_k` | `(k) -> Vec<(&str, f64)>` | Top-k nodes by score |

### `CommunityResult`
| Field | Type | Description |
|-------|------|-------------|
| `assignments` | `HashMap<String, usize>` | Node ID to community label |
| `modularity` | `f64` | Modularity score (Q) |
| `num_communities` | `usize` | Number of communities |

**Methods:**

| Method | Signature | Description |
|--------|-----------|-------------|
| `members` | `(community) -> Vec<&str>` | Nodes in a community |
| `community_sizes` | `() -> Vec<(usize, usize)>` | (community_id, size) sorted by size desc |

### `PpiEvidence`
STRING-style evidence channels for PPI scoring.

| Field | Type | Description |
|-------|------|-------------|
| `protein_a` | `String` | Protein A identifier |
| `protein_b` | `String` | Protein B identifier |
| `experimental` | `f64` | Experimental evidence (0.0-1.0) |
| `database` | `f64` | Database/curated evidence |
| `textmining` | `f64` | Text-mining evidence |
| `coexpression` | `f64` | Co-expression evidence |
| `neighborhood` | `f64` | Genomic context evidence |
| `cooccurrence` | `f64` | Phylogenetic profiles evidence |

**Methods:**

| Method | Signature | Description |
|--------|-----------|-------------|
| `new` | `(protein_a, protein_b) -> Self` | Create with all scores at 0.0 |
| `combined_score` | `() -> f64` | STRING noisy-OR combined score |

### `HubBottleneckResult`
| Field | Type | Description |
|-------|------|-------------|
| `degrees` | `HashMap<String, usize>` | Node degrees |
| `betweenness` | `HashMap<String, f64>` | Betweenness centrality scores |
| `hubs` | `Vec<String>` | High-degree nodes |
| `bottlenecks` | `Vec<String>` | High-betweenness nodes |
| `hub_bottlenecks` | `Vec<String>` | Nodes that are both hubs and bottlenecks |

### `ExpressionMatrix`
Gene-by-sample expression matrix for GRN inference.

| Field | Type | Description |
|-------|------|-------------|
| `genes` | `Vec<String>` | Gene identifiers (row labels) |
| `samples` | `Vec<String>` | Sample identifiers (column labels) |
| `data` | `Vec<Vec<f64>>` | Expression values (genes x samples) |

**Methods:**

| Method | Signature | Description |
|--------|-----------|-------------|
| `new` | `(genes, samples, data) -> Result<Self>` | Create with validation |
| `num_genes` | `() -> usize` | Number of genes |
| `num_samples` | `() -> usize` | Number of samples |
| `gene_expression` | `(gene_idx) -> &[f64]` | Expression vector for a gene |

### `CorrelationMethod`
`Pearson`, `Spearman`, `MutualInformation`

### `GrnResult`
| Field | Type | Description |
|-------|------|-------------|
| `network` | `Graph` | Inferred regulatory network (directed) |
| `scores` | `HashMap<(String, String), f64>` | Edge scores: (TF, target) to score |
| `method` | `String` | Method used |

### `Pathway`
| Field | Type | Description |
|-------|------|-------------|
| `id` | `String` | Pathway ID (e.g., "hsa04110", "R-HSA-69278") |
| `name` | `String` | Pathway name |
| `source` | `PathwaySource` | Source database |
| `members` | `Vec<String>` | Gene/protein members |
| `category` | `Option<String>` | Top-level classification |
| `organism` | `Option<String>` | Organism |

**Methods:**

| Method | Signature | Description |
|--------|-----------|-------------|
| `new` | `(id, name, source, members) -> Self` | Create pathway |
| `size` | `() -> usize` | Number of members |
| `contains` | `(gene) -> bool` | Check membership |

### `PathwaySource`
`KEGG`, `Reactome`, `WikiPathways`, `BioCyc`, `Custom`

### `PathwayDatabase`
| Field | Type | Description |
|-------|------|-------------|
| `pathways` | `HashMap<String, Pathway>` | Pathways indexed by ID |

**Methods:**

| Method | Signature | Description |
|--------|-----------|-------------|
| `new` | `() -> Self` | Create empty database |
| `add` | `(Pathway)` | Add a pathway |
| `pathways_for_gene` | `(gene) -> Vec<&Pathway>` | Find pathways containing a gene |
| `len` | `() -> usize` | Number of pathways |
| `is_empty` | `() -> bool` | Check if empty |

### `PathwayScore`
| Field | Type | Description |
|-------|------|-------------|
| `pathway_id` | `String` | Pathway ID |
| `pathway_name` | `String` | Pathway name |
| `members_in_network` | `usize` | Members present in the network |
| `total_members` | `usize` | Total pathway members |
| `avg_centrality` | `f64` | Average degree centrality of members |
| `internal_density` | `f64` | Fraction of internal edges |
| `topology_score` | `f64` | Combined score (0.5 * density + 0.5 * centrality) |

### `CrosstalkResult`
| Field | Type | Description |
|-------|------|-------------|
| `pathway_a` | `String` | Pathway A ID |
| `pathway_b` | `String` | Pathway B ID |
| `shared_members` | `usize` | Shared gene count |
| `jaccard` | `f64` | Jaccard similarity of member sets |
| `cross_edges` | `usize` | Cross-edges in the network |

## Functions

### Topology Metrics

| Function | Signature | Description |
|----------|-----------|-------------|
| `degree_centrality` | `(graph) -> Result<CentralityResult>` | C_d(v) = deg(v) / (n-1) |
| `closeness_centrality` | `(graph) -> Result<CentralityResult>` | Wasserman-Faust normalized closeness |
| `betweenness_centrality` | `(graph) -> Result<CentralityResult>` | Brandes' algorithm, normalized |
| `pagerank` | `(graph, damping, max_iter, tolerance) -> Result<CentralityResult>` | PageRank via power iteration |
| `clustering_coefficient` | `(graph, node_id) -> Result<f64>` | Local clustering coefficient |
| `average_clustering` | `(graph) -> Result<f64>` | Mean clustering over all nodes |
| `bfs_distances` | `(graph, source) -> HashMap<String, usize>` | BFS shortest-path distances |
| `shortest_path` | `(graph, source, target) -> Result<Option<Vec<String>>>` | BFS shortest path |
| `diameter` | `(graph) -> Result<Option<usize>>` | Longest shortest path (None if disconnected) |
| `connected_components` | `(graph) -> Vec<HashSet<String>>` | Connected components |

### Community Detection

| Function | Signature | Description |
|----------|-----------|-------------|
| `modularity` | `(graph, assignments) -> f64` | Newman-Girvan modularity Q |
| `louvain` | `(graph, resolution) -> Result<CommunityResult>` | Louvain modularity optimization |
| `label_propagation` | `(graph, max_iter) -> Result<CommunityResult>` | Label propagation |

### PPI Network Analysis

| Function | Signature | Description |
|----------|-----------|-------------|
| `build_ppi_network` | `(evidence, threshold) -> Result<Graph>` | Build PPI graph from evidence with confidence filter |
| `network_propagation` | `(graph, seeds, restart_prob, max_iter, tolerance) -> Result<HashMap<String, f64>>` | Random walk with restart |
| `hub_bottleneck_analysis` | `(graph, degree_pct, betweenness_pct) -> Result<HubBottleneckResult>` | Identify hubs and bottlenecks |
| `neighborhood_similarity` | `(graph, node_a, node_b) -> Result<f64>` | Jaccard similarity of neighborhoods |
| `predict_interactions` | `(graph, top_k) -> Result<Vec<(String, String, f64)>>` | Link prediction by neighborhood similarity |

### GRN Inference

| Function | Signature | Description |
|----------|-----------|-------------|
| `infer_grn` | `(expr, tf_genes, method, threshold) -> Result<GrnResult>` | Pairwise correlation/MI-based GRN |
| `clr` | `(expr, tf_genes, threshold) -> Result<GrnResult>` | Context Likelihood of Relatedness |

### Pathway Analysis

| Function | Signature | Description |
|----------|-----------|-------------|
| `score_pathways_by_topology` | `(network, db) -> Result<Vec<PathwayScore>>` | Score pathways by network topology |
| `pathway_crosstalk` | `(network, db, min_overlap) -> Result<Vec<CrosstalkResult>>` | Pairwise pathway overlap and connectivity |
| `parse_gmt` | `(content) -> Result<PathwayDatabase>` | Parse GMT format |
| `write_gmt` | `(db) -> String` | Write GMT format |

### Format I/O

| Function | Signature | Description |
|----------|-----------|-------------|
| `parse_sif` | `(content) -> Result<Graph>` | Parse SIF (Cytoscape) |
| `write_sif` | `(graph) -> String` | Write SIF |
| `parse_graphml` | `(content) -> Result<Graph>` | Parse simplified GraphML |
| `write_graphml` | `(graph) -> String` | Write GraphML |
| `write_gexf` | `(graph) -> String` | Write GEXF (Gephi) |
