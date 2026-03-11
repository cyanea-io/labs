# cyanea-network Architecture

## Module Map

```
cyanea-network
 +-- graph.rs       Graph, Node, Edge, GraphType — adjacency-list core
 +-- topology.rs    Centrality, clustering, shortest paths, components
 +-- community.rs   Louvain modularity, label propagation
 +-- ppi.rs         PPI networks, STRING scoring, propagation, hub/bottleneck
 +-- grn.rs         GRN inference: correlation, MI, CLR
 +-- pathway.rs     Pathway representation, topology scoring, crosstalk, GMT I/O
 +-- formats.rs     GraphML, SIF, GEXF parsing and writing
```

## Design Decisions

### Adjacency-List Graph

The `Graph` struct uses a `HashMap<String, Vec<(String, usize)>>` adjacency list where each entry maps a node ID to a list of (neighbor_id, edge_index) pairs. Edges are stored canonically in a `Vec<Edge>`. For undirected graphs, each edge is entered into the adjacency list in both directions but stored only once in the edge vector, so `edge_count()` returns the true edge count.

String-keyed nodes are used instead of integer indices for direct compatibility with gene symbols, protein IDs, and pathway identifiers. The `adjacency_matrix()` method sorts node IDs lexicographically for deterministic matrix output.

Node removal rebuilds the entire adjacency list from the edge vector because edge indices shift after filtering. This is O(E) but keeps the data structure consistent without tombstones.

### Brandes Betweenness Centrality

Betweenness centrality uses Brandes' algorithm (2001), which computes centrality for all nodes in O(VE) time for unweighted graphs. For each source node s:

1. BFS from s, recording shortest-path counts (sigma) and predecessors
2. Stack-based back-propagation accumulates dependency values (delta)
3. Delta values are added to the global centrality scores

Normalization divides by (n-1)(n-2) for directed graphs and (n-1)(n-2)/2 for undirected, following the standard convention.

### Closeness Centrality

Uses Wasserman-Faust normalization for disconnected graphs: `C_c(v) = (reachable / sum_dist) * (reachable / (n-1))`. This avoids the undefined-closeness problem for unreachable nodes by scaling the raw closeness by the fraction of the graph that is reachable.

### PageRank

Iterative power method with configurable damping factor, max iterations, and convergence tolerance. Dangling nodes (zero out-degree) distribute their rank mass uniformly to all nodes. Convergence is checked as the maximum absolute change across all node scores per iteration.

### Louvain Community Detection

Implements the Louvain algorithm for modularity optimization:

1. **Initialization**: Each node starts in its own community.
2. **Local moves**: For each node, compute the modularity gain of moving it to each neighboring community. The gain formula is: `delta_Q = w_ic - resolution * k_i * sigma_tot_c / 2m`, where `w_ic` is the edge weight to community c, `k_i` is the node's weighted degree, and `sigma_tot_c` is the community's total weighted degree.
3. **Convergence**: Repeat passes until no node moves (max 20 passes).
4. **Renumbering**: Community labels are made contiguous (0, 1, 2, ...).

The resolution parameter (default 1.0) controls granularity: values above 1.0 produce smaller communities, below 1.0 produces larger ones. The current implementation performs single-level optimization without the aggregation phase of multi-level Louvain.

### Newman-Girvan Modularity

Q = (1/2m) * sum_ij [A_ij - k_i*k_j/(2m)] * delta(c_i, c_j)

Computed via the full adjacency matrix. Used both as a standalone metric and internally by Louvain and label propagation to score their final partitions.

### Label Propagation

Each node adopts the label with the highest total edge weight among its neighbors. Weighted voting ensures that strong connections have more influence. The algorithm is non-deterministic due to iteration order but fast (typically converges in under 10 iterations). Max iteration cap prevents infinite loops in oscillating cases.

### STRING Noisy-OR Scoring

PPI evidence channels are combined using the noisy-OR model from STRING:

`combined = 1 - product(1 - s_i)` for all non-zero channels.

This treats each evidence channel as an independent source of support. Only channels with positive scores contribute. Scores are clamped to [0, 1] before combination.

### Network Propagation (Random Walk with Restart)

Implements RWR for guilt-by-association analysis:

1. Seed scores are normalized to sum to 1
2. Row-normalize the adjacency matrix (transition matrix)
3. At each step: `score_new = (1-r) * W^T * score + r * seed`, where r is the restart probability
4. Iterate until max absolute change < tolerance

Higher restart probability (0.3-0.5) keeps scores closer to seed nodes; lower values allow broader propagation. The column-vector multiply `W^T * score` propagates scores from neighbors.

### Hub-Bottleneck Analysis

Classifies nodes using percentile thresholds on two metrics:

- **Hubs**: Nodes with degree above the specified degree percentile
- **Bottlenecks**: Nodes with betweenness centrality above the specified betweenness percentile
- **Hub-bottlenecks**: Intersection of the two sets

Betweenness is computed via the Brandes algorithm. These classifications are biologically meaningful: hub-bottleneck proteins tend to be essential genes and drug targets.

### GRN Inference

Two inference approaches, both producing directed graphs from TF to target:

1. **Pairwise correlation/MI**: Computes Pearson, Spearman (rank-based Pearson), or mutual information between each TF and all genes. Edges above the threshold are included with |correlation| or MI as weight.

2. **CLR (Context Likelihood of Relatedness)**: Enhances MI-based inference by converting raw MI values to z-scores relative to the background distribution for each gene. The CLR score is `sqrt(z_tf^2 + z_target^2)`, which highlights interactions that are strong relative to both the TF's and the target's typical MI values. This reduces false positives from highly variable genes.

### Mutual Information Estimation

Uses histogram-based estimation with `sqrt(n)` bins. For each pair of variables, a 2D histogram of joint counts plus 1D marginals is computed, then MI = sum(p_xy * ln(p_xy / (p_x * p_y))). Result is clamped to non-negative. This estimator is simple and fast but may underestimate MI for small sample sizes.

### Pathway Topology Scoring

Each pathway is scored by how coherently its members form a module in the network:

- **Internal density**: Fraction of possible edges between pathway members that actually exist
- **Average centrality**: Mean degree centrality of pathway members, normalized by (n-1)
- **Combined score**: 0.5 * density + 0.5 * centrality

Pathways with fewer than 2 members in the network receive a score of 0. Results are sorted by topology score descending.

### Pathway Crosstalk

For each pair of pathways, computes:

- **Shared members**: Intersection of member sets
- **Jaccard similarity**: |intersection| / |union|
- **Cross-edges**: Network edges from pathway A members to pathway B members (excluding shared members)

Only pairs with at least `min_overlap` shared members are returned.

### GMT Format

The Gene Matrix Transposed format is tab-delimited with one pathway per line: `id\tname\tgene1\tgene2\t...`. Parsing skips lines with fewer than 3 fields. All pathways loaded from GMT get `PathwaySource::Custom`. Writing sorts pathways by ID for deterministic output.

### Format I/O

Three network exchange formats are supported:

- **SIF**: Cytoscape Simple Interaction Format. Line-based: `source\ttype\ttarget1\ttarget2`. Supports comments (#), isolated nodes, multi-target lines. Parsed as undirected.
- **GraphML**: XML-based. Supports `<key>` attribute definitions, `<data>` elements for node/edge attributes, and `edgedefault` for graph direction. Parsing is line-based (not a full XML parser) for simplicity.
- **GEXF**: Gephi exchange format. Write-only. Produces valid GEXF 1.3 XML with nodes, edges, and weights.

## Dependencies

```
cyanea-core (errors, traits)
```

No external dependencies beyond `cyanea-core` and `thiserror`. All graph algorithms, scoring functions, statistical methods, and format parsers are implemented from first principles.

## Testing Strategy

- **Unit tests**: 87 tests inline in each module covering all public functions
- **Doc tests**: 2 doc tests on `Graph` demonstrating basic usage
- **Total**: 89 tests across 7 source files
- Centrality values validated against known graph topologies (triangle, line, star)
- Louvain validated on two-clique benchmark (detects correct partition)
- STRING combined score validated with manual noisy-OR calculations
- SIF and GraphML roundtrip parsing validated
- GMT roundtrip validated
- CLR validated against raw MI (produces non-empty results on correlated data)
