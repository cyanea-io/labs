# cyanea-phylo

> Phylogenetics: tree data structures, I/O, distance models, tree construction, ancestral reconstruction, and maximum likelihood inference.

## What's Inside

- **Tree representation** -- rooted `PhyloTree` with pre-order/post-order traversal, rerooting, subtree extraction, MRCA queries
- **Newick I/O** -- parse and serialize Newick format
- **NEXUS I/O** -- parse and serialize NEXUS format with taxa blocks
- **Evolutionary distances** -- p-distance, Jukes-Cantor, Kimura 2-parameter
- **Tree construction** -- UPGMA and Neighbor-Joining (requires `ml` feature)
- **Tree comparison** -- Robinson-Foulds (unweighted and normalized), branch score distance
- **Ancestral reconstruction** -- Fitch (maximum parsimony), Sankoff (weighted parsimony), marginal ML
- **Substitution models** -- JC69, HKY85, GTR with eigendecomposed P(t), discrete gamma rate heterogeneity
- **Maximum likelihood** -- Felsenstein pruning, GTR+Gamma, NNI tree search
- **Bootstrap support** -- non-parametric bootstrap with bipartition matching
- **Consensus trees** -- strict, majority-rule, extended majority-rule
- **Molecular dating** -- strict clock, root-to-tip regression
- **Tree drawing** -- 2D layout coordinates (rectangular, cladogram, radial)

## Quick Start

```toml
[dependencies]
cyanea-phylo = { version = "0.1", features = ["ml"] }
```

```rust
use cyanea_phylo::{PhyloTree, compare::robinson_foulds};

let tree1 = PhyloTree::from_newick("((A:0.1,B:0.2):0.3,C:0.4);").unwrap();
let tree2 = PhyloTree::from_newick("((A:0.1,C:0.2):0.3,B:0.4);").unwrap();

println!("Leaves: {:?}", tree1.leaf_names());
println!("RF distance: {}", robinson_foulds(&tree1, &tree2).unwrap());
```

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support |
| `wasm` | No | WASM target marker |
| `serde` | No | Serialize/Deserialize derives |
| `ml` | No | cyanea-ml for distance matrices and tree construction |
| `parallel` | No | Rayon parallelism |

## Modules

| Module | Description |
|--------|-------------|
| `tree` | `PhyloTree`, `Node`, traversal iterators |
| `newick` | Newick parser and writer |
| `nexus` | NEXUS parser and writer |
| `distance` | p-distance, Jukes-Cantor, Kimura 2-parameter |
| `compare` | Robinson-Foulds, branch score distance |
| `construct` | UPGMA, Neighbor-Joining (feature-gated) |
| `reconstruct` | Fitch, Sankoff ancestral reconstruction |
| `marginal` | Marginal ML ancestral reconstruction |
| `models` | JC69, HKY85, GTR, discrete gamma rates |
| `likelihood` | Felsenstein pruning, NNI search |
| `bootstrap` | Bootstrap support, bipartition extraction |
| `consensus` | Strict/majority-rule/extended consensus trees |
| `dating` | Strict clock, root-to-tip regression |
| `drawing` | Tree layout coordinates (rectangular/cladogram/radial) |

## See Also

- [API Reference (STATUS.md)](docs/STATUS.md)
- [Architecture](../ARCHITECTURE.md)
- [Build Guide](../docs/BUILDING.md)
