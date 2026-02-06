//! Core phylogenetic tree data structures.
//!
//! Uses arena-style storage: nodes live in a flat `Vec<Node>` and are
//! referenced by `NodeId` (a `usize` index).

use cyanea_core::{CyaneaError, Result, Summarizable};

/// Index into the tree's node arena.
pub type NodeId = usize;

/// A single node in a phylogenetic tree.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Node {
    /// Index of this node in the arena.
    pub id: NodeId,
    /// Parent node (None for root).
    pub parent: Option<NodeId>,
    /// Child nodes.
    pub children: Vec<NodeId>,
    /// Branch length from this node to its parent.
    pub branch_length: Option<f64>,
    /// Taxon or clade label.
    pub name: Option<String>,
}

impl Node {
    /// True if this node has no children.
    pub fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }

    /// True if this node has no parent.
    pub fn is_root(&self) -> bool {
        self.parent.is_none()
    }

    /// True if this node has both parent and children.
    pub fn is_internal(&self) -> bool {
        !self.is_leaf() && !self.is_root()
    }
}

/// A rooted phylogenetic tree stored as an arena of nodes.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PhyloTree {
    nodes: Vec<Node>,
    root: NodeId,
}

impl PhyloTree {
    /// Create a new tree with a single unnamed root node.
    pub fn new() -> Self {
        let root = Node {
            id: 0,
            parent: None,
            children: Vec::new(),
            branch_length: None,
            name: None,
        };
        Self {
            nodes: vec![root],
            root: 0,
        }
    }

    /// Create a tree from pre-built nodes and a root index.
    ///
    /// This is used by the Newick parser and tree construction algorithms.
    pub fn from_nodes(nodes: Vec<Node>, root: NodeId) -> Result<Self> {
        if nodes.is_empty() {
            return Err(CyaneaError::InvalidInput("empty node list".into()));
        }
        if root >= nodes.len() {
            return Err(CyaneaError::InvalidInput(format!(
                "root index {} out of range ({})",
                root,
                nodes.len()
            )));
        }
        Ok(Self { nodes, root })
    }

    /// Add a child to `parent` and return its `NodeId`.
    pub fn add_child(
        &mut self,
        parent: NodeId,
        name: Option<String>,
        branch_length: Option<f64>,
    ) -> Result<NodeId> {
        if parent >= self.nodes.len() {
            return Err(CyaneaError::InvalidInput(format!(
                "parent index {} out of range ({})",
                parent,
                self.nodes.len()
            )));
        }
        let id = self.nodes.len();
        self.nodes.push(Node {
            id,
            parent: Some(parent),
            children: Vec::new(),
            branch_length,
            name,
        });
        self.nodes[parent].children.push(id);
        Ok(id)
    }

    /// Access a node by id.
    pub fn get_node(&self, id: NodeId) -> Option<&Node> {
        self.nodes.get(id)
    }

    /// Mutable access to a node by id.
    pub fn get_node_mut(&mut self, id: NodeId) -> Option<&mut Node> {
        self.nodes.get_mut(id)
    }

    /// The root node id.
    pub fn root(&self) -> NodeId {
        self.root
    }

    /// Total number of nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Number of leaf nodes.
    pub fn leaf_count(&self) -> usize {
        self.nodes.iter().filter(|n| n.is_leaf()).count()
    }

    /// All leaf node ids.
    pub fn leaves(&self) -> Vec<NodeId> {
        self.nodes
            .iter()
            .filter(|n| n.is_leaf())
            .map(|n| n.id)
            .collect()
    }

    /// All internal (non-leaf, non-root) node ids.
    pub fn internal_nodes(&self) -> Vec<NodeId> {
        self.nodes
            .iter()
            .filter(|n| n.is_internal())
            .map(|n| n.id)
            .collect()
    }

    /// Maximum depth from root to any leaf.
    pub fn depth(&self) -> usize {
        fn dfs(tree: &PhyloTree, node: NodeId) -> usize {
            let children = &tree.nodes[node].children;
            if children.is_empty() {
                return 0;
            }
            children.iter().map(|&c| 1 + dfs(tree, c)).max().unwrap_or(0)
        }
        dfs(self, self.root)
    }

    /// Pre-order (parent before children) traversal yielding node ids.
    pub fn iter_preorder(&self) -> PreorderIter<'_> {
        PreorderIter {
            tree: self,
            stack: vec![self.root],
        }
    }

    /// Post-order (children before parent) traversal yielding node ids.
    pub fn iter_postorder(&self) -> PostorderIter {
        // Build postorder sequence by reversing a modified preorder
        // (visit right children first, then reverse the whole thing).
        let mut result = Vec::new();
        let mut stack = vec![self.root];
        while let Some(id) = stack.pop() {
            result.push(id);
            for &child in &self.nodes[id].children {
                stack.push(child);
            }
        }
        result.reverse();
        PostorderIter {
            sequence: result,
            pos: 0,
        }
    }

    /// Most recent common ancestor of two nodes.
    pub fn mrca(&self, a: NodeId, b: NodeId) -> Result<NodeId> {
        if a >= self.nodes.len() || b >= self.nodes.len() {
            return Err(CyaneaError::InvalidInput("node id out of range".into()));
        }
        // Collect ancestors of a (including a itself)
        let mut ancestors_a = Vec::new();
        let mut cur = a;
        loop {
            ancestors_a.push(cur);
            match self.nodes[cur].parent {
                Some(p) => cur = p,
                None => break,
            }
        }
        // Walk up from b and find first match
        cur = b;
        loop {
            if ancestors_a.contains(&cur) {
                return Ok(cur);
            }
            match self.nodes[cur].parent {
                Some(p) => cur = p,
                None => break,
            }
        }
        // Should always find root as common ancestor
        Ok(self.root)
    }

    /// Sorted list of leaf names (leaves without names are excluded).
    pub fn leaf_names(&self) -> Vec<String> {
        let mut names: Vec<String> = self
            .nodes
            .iter()
            .filter(|n| n.is_leaf())
            .filter_map(|n| n.name.clone())
            .collect();
        names.sort();
        names
    }

    /// Parse a Newick format string into a tree.
    pub fn from_newick(input: &str) -> Result<Self> {
        crate::newick::parse(input)
    }

    /// Serialize the tree to a Newick format string.
    pub fn to_newick(&self) -> String {
        crate::newick::write(self)
    }
}

impl Default for PhyloTree {
    fn default() -> Self {
        Self::new()
    }
}

impl Summarizable for PhyloTree {
    fn summary(&self) -> String {
        let leaves = self.leaf_count();
        let internal = self.node_count() - leaves;
        format!(
            "PhyloTree: {} nodes ({} leaves, {} internal)",
            self.node_count(),
            leaves,
            internal
        )
    }
}

/// Pre-order iterator over node ids.
pub struct PreorderIter<'a> {
    tree: &'a PhyloTree,
    stack: Vec<NodeId>,
}

impl<'a> Iterator for PreorderIter<'a> {
    type Item = NodeId;

    fn next(&mut self) -> Option<Self::Item> {
        let id = self.stack.pop()?;
        // Push children in reverse order so leftmost is visited first.
        for &child in self.tree.nodes[id].children.iter().rev() {
            self.stack.push(child);
        }
        Some(id)
    }
}

/// Post-order iterator over node ids.
pub struct PostorderIter {
    sequence: Vec<NodeId>,
    pos: usize,
}

impl Iterator for PostorderIter {
    type Item = NodeId;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos < self.sequence.len() {
            let id = self.sequence[self.pos];
            self.pos += 1;
            Some(id)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_tree() -> PhyloTree {
        // ((A:0.1,B:0.2)AB:0.3,(C:0.4,D:0.5)CD:0.6)root;
        let mut tree = PhyloTree::new();
        tree.get_node_mut(0).unwrap().name = Some("root".into());
        let ab = tree.add_child(0, Some("AB".into()), Some(0.3)).unwrap();
        let cd = tree.add_child(0, Some("CD".into()), Some(0.6)).unwrap();
        tree.add_child(ab, Some("A".into()), Some(0.1)).unwrap();
        tree.add_child(ab, Some("B".into()), Some(0.2)).unwrap();
        tree.add_child(cd, Some("C".into()), Some(0.4)).unwrap();
        tree.add_child(cd, Some("D".into()), Some(0.5)).unwrap();
        tree
    }

    #[test]
    fn new_tree_has_single_root() {
        let tree = PhyloTree::new();
        assert_eq!(tree.node_count(), 1);
        assert_eq!(tree.leaf_count(), 1);
        assert!(tree.get_node(0).unwrap().is_root());
    }

    #[test]
    fn add_child_works() {
        let mut tree = PhyloTree::new();
        let c1 = tree.add_child(0, Some("A".into()), Some(1.0)).unwrap();
        assert_eq!(c1, 1);
        assert_eq!(tree.node_count(), 2);
        assert_eq!(tree.get_node(c1).unwrap().parent, Some(0));
        assert_eq!(tree.get_node(0).unwrap().children, vec![1]);
    }

    #[test]
    fn add_child_invalid_parent() {
        let mut tree = PhyloTree::new();
        assert!(tree.add_child(99, None, None).is_err());
    }

    #[test]
    fn leaf_and_internal_counts() {
        let tree = sample_tree();
        assert_eq!(tree.node_count(), 7);
        assert_eq!(tree.leaf_count(), 4);
        assert_eq!(tree.leaves().len(), 4);
        assert_eq!(tree.internal_nodes().len(), 2); // AB and CD (root is not "internal")
    }

    #[test]
    fn depth_of_sample_tree() {
        let tree = sample_tree();
        assert_eq!(tree.depth(), 2);
    }

    #[test]
    fn depth_of_single_node() {
        let tree = PhyloTree::new();
        assert_eq!(tree.depth(), 0);
    }

    #[test]
    fn preorder_traversal() {
        let tree = sample_tree();
        let order: Vec<NodeId> = tree.iter_preorder().collect();
        // root(0), AB(1), A(3), B(4), CD(2), C(5), D(6)
        assert_eq!(order, vec![0, 1, 3, 4, 2, 5, 6]);
    }

    #[test]
    fn postorder_traversal() {
        let tree = sample_tree();
        let order: Vec<NodeId> = tree.iter_postorder().collect();
        // A(3), B(4), AB(1), C(5), D(6), CD(2), root(0)
        assert_eq!(order, vec![3, 4, 1, 5, 6, 2, 0]);
    }

    #[test]
    fn mrca_siblings() {
        let tree = sample_tree();
        // A(3) and B(4) are siblings under AB(1)
        assert_eq!(tree.mrca(3, 4).unwrap(), 1);
    }

    #[test]
    fn mrca_cousins() {
        let tree = sample_tree();
        // A(3) and C(5) share root(0) as MRCA
        assert_eq!(tree.mrca(3, 5).unwrap(), 0);
    }

    #[test]
    fn mrca_parent_child() {
        let tree = sample_tree();
        // AB(1) and A(3): MRCA is AB(1)
        assert_eq!(tree.mrca(1, 3).unwrap(), 1);
    }

    #[test]
    fn mrca_same_node() {
        let tree = sample_tree();
        assert_eq!(tree.mrca(3, 3).unwrap(), 3);
    }

    #[test]
    fn leaf_names_sorted() {
        let tree = sample_tree();
        assert_eq!(tree.leaf_names(), vec!["A", "B", "C", "D"]);
    }

    #[test]
    fn summary_format() {
        let tree = sample_tree();
        assert_eq!(tree.summary(), "PhyloTree: 7 nodes (4 leaves, 3 internal)");
    }
}
