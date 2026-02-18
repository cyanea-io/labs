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

    /// Read-only access to the full node arena.
    pub fn nodes(&self) -> &[Node] {
        &self.nodes
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

    /// Collect all leaf names in the subtree rooted at `node_id`.
    pub fn subtree_leaf_names(&self, node_id: NodeId) -> std::collections::BTreeSet<String> {
        let mut leaves = std::collections::BTreeSet::new();
        let mut stack = vec![node_id];
        while let Some(id) = stack.pop() {
            let node = &self.nodes[id];
            if node.is_leaf() {
                if let Some(ref name) = node.name {
                    leaves.insert(name.clone());
                }
            } else {
                for &child in &node.children {
                    stack.push(child);
                }
            }
        }
        leaves
    }

    /// Sum of all branch lengths in the tree.
    pub fn total_branch_length(&self) -> f64 {
        self.nodes
            .iter()
            .filter_map(|n| n.branch_length)
            .sum()
    }

    /// Reroot the tree on the edge leading to `node_id`.
    ///
    /// `position` is the fraction (0.0–1.0) along the branch from `node_id`'s
    /// parent toward `node_id` where the new root is placed. Defaults to 0.5
    /// (midpoint of the edge).
    pub fn reroot(&self, node_id: NodeId, position: Option<f64>) -> Result<Self> {
        if node_id >= self.nodes.len() {
            return Err(CyaneaError::InvalidInput("node id out of range".into()));
        }
        if node_id == self.root {
            return Ok(self.clone());
        }
        let parent_id = self.nodes[node_id]
            .parent
            .ok_or_else(|| CyaneaError::InvalidInput("cannot reroot at root".into()))?;

        let frac = position.unwrap_or(0.5).clamp(0.0, 1.0);
        let edge_len = self.nodes[node_id].branch_length.unwrap_or(0.0);
        let len_toward_node = edge_len * frac;
        let len_toward_parent = edge_len - len_toward_node;

        // Collect path from node_id's parent up to root (we will reverse edges along this path).
        let mut path = Vec::new();
        let mut cur = parent_id;
        loop {
            path.push(cur);
            match self.nodes[cur].parent {
                Some(p) => cur = p,
                None => break,
            }
        }

        // Build new tree: clone all nodes first.
        let mut nodes: Vec<Node> = self.nodes.iter().cloned().collect();

        // Create new root node.
        let new_root_id = nodes.len();
        nodes.push(Node {
            id: new_root_id,
            parent: None,
            children: vec![node_id, parent_id],
            branch_length: None,
            name: None,
        });

        // node_id becomes child of new root.
        nodes[node_id].parent = Some(new_root_id);
        nodes[node_id].branch_length = Some(len_toward_node);

        // Reverse edges along the path from parent_id to old root.
        // parent_id becomes the other child of new root.
        nodes[parent_id].children.retain(|&c| c != node_id);
        nodes[parent_id].parent = Some(new_root_id);
        // The branch from parent_id to new_root gets the remaining length.
        let old_parent_bl = nodes[parent_id].branch_length;
        nodes[parent_id].branch_length = Some(len_toward_parent);

        // Reverse remaining edges along the path (parent_id → grandparent → ... → old root).
        let mut prev_id = parent_id;
        let mut prev_old_bl = old_parent_bl;
        for &cur_id in &path[1..] {
            // cur_id was the parent of prev_id. Now prev_id becomes parent of cur_id.
            nodes[cur_id].children.retain(|&c| c != prev_id);
            nodes[cur_id].parent = Some(prev_id);
            nodes[prev_id].children.push(cur_id);
            // Branch length: cur_id takes the old branch length that prev_id had to cur_id.
            let old_bl = nodes[cur_id].branch_length;
            nodes[cur_id].branch_length = prev_old_bl;
            prev_old_bl = old_bl;
            prev_id = cur_id;
        }

        // The old root: if it now has degree 1, collapse it.
        let old_root = *path.last().unwrap();
        if nodes[old_root].children.len() == 1 {
            let only_child = nodes[old_root].children[0];
            let grandparent = nodes[old_root].parent;
            // Merge branch lengths.
            let merged_bl = nodes[old_root].branch_length.unwrap_or(0.0)
                + nodes[only_child].branch_length.unwrap_or(0.0);
            nodes[only_child].parent = grandparent;
            nodes[only_child].branch_length = Some(merged_bl);
            if let Some(gp) = grandparent {
                for c in &mut nodes[gp].children {
                    if *c == old_root {
                        *c = only_child;
                    }
                }
            }
            // Mark old root as removed (orphan with no children).
            nodes[old_root].parent = None;
            nodes[old_root].children.clear();
        }

        // Remap to remove orphan nodes and fix IDs.
        Self::compact_nodes(nodes, new_root_id)
    }

    /// Reroot at the midpoint of the longest root-to-root path.
    pub fn midpoint_root(&self) -> Result<Self> {
        let leaves = self.leaves();
        if leaves.len() < 2 {
            return Ok(self.clone());
        }

        // Find the two most distant leaves via two-pass BFS.
        let dists_from_first = self.distances_from(leaves[0]);
        let farthest1 = *leaves
            .iter()
            .max_by(|&&a, &&b| {
                dists_from_first[a]
                    .partial_cmp(&dists_from_first[b])
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap();
        let dists_from_farthest = self.distances_from(farthest1);
        let farthest2 = *leaves
            .iter()
            .max_by(|&&a, &&b| {
                dists_from_farthest[a]
                    .partial_cmp(&dists_from_farthest[b])
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap();
        let leaf_a = farthest1;
        let leaf_b = farthest2;
        let max_dist = dists_from_farthest[farthest2];

        if max_dist == 0.0 {
            return Ok(self.clone());
        }

        // Find the path from leaf_a to leaf_b through MRCA.
        let lca = self.mrca(leaf_a, leaf_b)?;

        // Walk from leaf_a toward LCA, accumulating distance, find where midpoint falls.
        let half = max_dist / 2.0;
        let target = self.find_midpoint_on_path(leaf_a, leaf_b, lca, half);

        self.reroot(target.0, Some(target.1))
    }

    /// Compute distances from a source node to all other nodes.
    fn distances_from(&self, source: NodeId) -> Vec<f64> {
        let n = self.nodes.len();
        let mut dist = vec![f64::MAX; n];
        dist[source] = 0.0;

        // Build adjacency: for each node, neighbors are parent + children.
        let mut visited = vec![false; n];
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(source);
        visited[source] = true;

        while let Some(u) = queue.pop_front() {
            let node = &self.nodes[u];
            // Children
            for &c in &node.children {
                if !visited[c] {
                    visited[c] = true;
                    let bl = self.nodes[c].branch_length.unwrap_or(0.0);
                    dist[c] = dist[u] + bl;
                    queue.push_back(c);
                }
            }
            // Parent
            if let Some(p) = node.parent {
                if !visited[p] {
                    visited[p] = true;
                    let bl = node.branch_length.unwrap_or(0.0);
                    dist[p] = dist[u] + bl;
                    queue.push_back(p);
                }
            }
        }
        dist
    }

    /// Find the edge containing the midpoint on the path between two leaves.
    /// Returns (node_id, fraction) for use with reroot.
    fn find_midpoint_on_path(
        &self,
        leaf_a: NodeId,
        leaf_b: NodeId,
        lca: NodeId,
        target_dist: f64,
    ) -> (NodeId, f64) {
        // Build path from leaf_a to lca.
        let mut path_a = Vec::new();
        let mut cur = leaf_a;
        while cur != lca {
            path_a.push(cur);
            cur = self.nodes[cur].parent.unwrap();
        }
        path_a.push(lca);

        // Build path from leaf_b to lca.
        let mut path_b = Vec::new();
        cur = leaf_b;
        while cur != lca {
            path_b.push(cur);
            cur = self.nodes[cur].parent.unwrap();
        }
        // Full path: leaf_a → ... → lca → ... → leaf_b
        let mut full_path = path_a;
        path_b.pop(); // Remove lca (already in path_a)
        path_b.reverse();
        full_path.extend(path_b);

        // Walk along path accumulating distance.
        let mut acc = 0.0;
        for i in 0..full_path.len() - 1 {
            let node = full_path[i];
            let next = full_path[i + 1];
            // Determine the edge length between node and next.
            let edge_len = if self.nodes[next].parent == Some(node) {
                self.nodes[next].branch_length.unwrap_or(0.0)
            } else {
                self.nodes[node].branch_length.unwrap_or(0.0)
            };
            if acc + edge_len >= target_dist {
                // Midpoint is on this edge.
                let remainder = target_dist - acc;
                // The edge is child→parent. We want to reroot on the child's edge.
                if self.nodes[next].parent == Some(node) {
                    // next is child of node, edge is next's branch
                    let frac = remainder / edge_len.max(1e-15);
                    return (next, frac);
                } else {
                    // node is child of next, edge is node's branch
                    let frac = 1.0 - remainder / edge_len.max(1e-15);
                    return (node, frac);
                }
            }
            acc += edge_len;
        }

        // Fallback: reroot at last node.
        (*full_path.last().unwrap(), 0.5)
    }

    /// Extract a subtree containing only the specified leaves.
    ///
    /// Collapses internal nodes with a single child, summing branch lengths.
    pub fn extract_subtree(&self, leaf_names: &[&str]) -> Result<Self> {
        let target: std::collections::BTreeSet<String> =
            leaf_names.iter().map(|s| s.to_string()).collect();

        // Verify all requested leaves exist.
        let all_leaves = self.leaf_names();
        for name in &target {
            if !all_leaves.contains(name) {
                return Err(CyaneaError::InvalidInput(format!(
                    "leaf '{}' not found in tree",
                    name
                )));
            }
        }
        if target.len() < 2 {
            return Err(CyaneaError::InvalidInput(
                "need at least 2 leaves for subtree extraction".into(),
            ));
        }

        // Mark nodes to keep: target leaves + all their ancestors.
        let mut keep = vec![false; self.nodes.len()];
        for id in 0..self.nodes.len() {
            let node = &self.nodes[id];
            if node.is_leaf() {
                if let Some(ref name) = node.name {
                    if target.contains(name) {
                        // Mark this leaf and all ancestors.
                        let mut cur = id;
                        loop {
                            keep[cur] = true;
                            match self.nodes[cur].parent {
                                Some(p) => cur = p,
                                None => break,
                            }
                        }
                    }
                }
            }
        }

        // Build new tree from kept nodes, collapsing degree-2 internals.
        self.build_filtered_tree(&keep)
    }

    /// Extract the subtree rooted at `node_id`.
    pub fn subtree_at(&self, node_id: NodeId) -> Result<Self> {
        if node_id >= self.nodes.len() {
            return Err(CyaneaError::InvalidInput("node id out of range".into()));
        }

        // DFS clone from node_id.
        let mut new_nodes = Vec::new();
        let mut old_to_new = std::collections::HashMap::new();

        let mut stack = vec![(node_id, None::<NodeId>)];
        while let Some((old_id, new_parent)) = stack.pop() {
            let old_node = &self.nodes[old_id];
            let new_id = new_nodes.len();
            old_to_new.insert(old_id, new_id);

            new_nodes.push(Node {
                id: new_id,
                parent: new_parent,
                children: Vec::new(),
                branch_length: if old_id == node_id {
                    None
                } else {
                    old_node.branch_length
                },
                name: old_node.name.clone(),
            });

            if let Some(np) = new_parent {
                new_nodes[np].children.push(new_id);
            }

            // Push children in reverse for consistent ordering.
            for &child in old_node.children.iter().rev() {
                stack.push((child, Some(new_id)));
            }
        }

        PhyloTree::from_nodes(new_nodes, 0)
    }

    /// Build a filtered tree from nodes marked to keep.
    fn build_filtered_tree(&self, keep: &[bool]) -> Result<Self> {
        // First, build the filtered topology preserving only kept nodes.
        // Then collapse degree-2 internal nodes.
        let mut new_nodes = Vec::new();
        let mut old_to_new = std::collections::HashMap::new();

        // Process in preorder to ensure parents are created before children.
        for old_id in self.iter_preorder() {
            if !keep[old_id] {
                continue;
            }
            let old_node = &self.nodes[old_id];
            let new_id = new_nodes.len();
            old_to_new.insert(old_id, new_id);

            // Find the kept parent.
            let new_parent = if old_id == self.root {
                None
            } else {
                let mut cur = old_node.parent;
                // Walk up to find the nearest kept ancestor.
                while let Some(p) = cur {
                    if keep[p] {
                        break;
                    }
                    cur = self.nodes[p].parent;
                }
                cur.and_then(|p| old_to_new.get(&p).copied())
            };

            // Compute accumulated branch length to kept parent.
            let branch_length = if old_id == self.root {
                old_node.branch_length
            } else {
                let mut bl = old_node.branch_length.unwrap_or(0.0);
                let mut cur = old_node.parent;
                while let Some(p) = cur {
                    if keep[p] {
                        break;
                    }
                    bl += self.nodes[p].branch_length.unwrap_or(0.0);
                    cur = self.nodes[p].parent;
                }
                Some(bl)
            };

            new_nodes.push(Node {
                id: new_id,
                parent: new_parent,
                children: Vec::new(),
                branch_length,
                name: old_node.name.clone(),
            });

            if let Some(np) = new_parent {
                new_nodes[np].children.push(new_id);
            }
        }

        // Collapse degree-2 internal nodes (non-root with exactly 1 child).
        loop {
            let mut collapsed = false;
            for i in 0..new_nodes.len() {
                if new_nodes[i].children.len() == 1
                    && new_nodes[i].parent.is_some()
                    && new_nodes[i].name.is_none()
                {
                    let child = new_nodes[i].children[0];
                    let parent = new_nodes[i].parent.unwrap();
                    // Merge branch lengths.
                    let merged = new_nodes[i].branch_length.unwrap_or(0.0)
                        + new_nodes[child].branch_length.unwrap_or(0.0);
                    new_nodes[child].branch_length = Some(merged);
                    new_nodes[child].parent = Some(parent);
                    // Replace i with child in parent's children.
                    for c in &mut new_nodes[parent].children {
                        if *c == i {
                            *c = child;
                        }
                    }
                    // Mark i as removed.
                    new_nodes[i].parent = None;
                    new_nodes[i].children.clear();
                    collapsed = true;
                    break;
                }
            }
            if !collapsed {
                break;
            }
        }

        // Also collapse root if it has degree 1.
        let mut root_id = 0;
        while new_nodes[root_id].children.len() == 1 && new_nodes[root_id].name.is_none() {
            let child = new_nodes[root_id].children[0];
            new_nodes[root_id].parent = None;
            new_nodes[root_id].children.clear();
            new_nodes[child].parent = None;
            root_id = child;
        }

        // Compact: remap IDs.
        Self::compact_nodes(new_nodes, root_id)
    }

    /// Remove orphan nodes and remap IDs contiguously.
    fn compact_nodes(nodes: Vec<Node>, root: NodeId) -> Result<Self> {
        // Identify live nodes (reachable from root).
        let mut live = vec![false; nodes.len()];
        let mut stack = vec![root];
        while let Some(id) = stack.pop() {
            if live[id] {
                continue;
            }
            live[id] = true;
            for &c in &nodes[id].children {
                stack.push(c);
            }
        }

        // Build mapping from old to new IDs.
        let mut old_to_new = vec![0usize; nodes.len()];
        let mut new_id = 0;
        for (old_id, &is_live) in live.iter().enumerate() {
            if is_live {
                old_to_new[old_id] = new_id;
                new_id += 1;
            }
        }

        let mut new_nodes = Vec::with_capacity(new_id);
        for (old_id, node) in nodes.iter().enumerate() {
            if !live[old_id] {
                continue;
            }
            new_nodes.push(Node {
                id: old_to_new[old_id],
                parent: node.parent.map(|p| old_to_new[p]),
                children: node.children.iter().map(|&c| old_to_new[c]).collect(),
                branch_length: node.branch_length,
                name: node.name.clone(),
            });
        }

        // Fix root parent.
        let new_root = old_to_new[root];
        new_nodes[new_root].parent = None;

        PhyloTree::from_nodes(new_nodes, new_root)
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

    #[test]
    fn subtree_leaf_names_works() {
        let tree = sample_tree();
        let names = tree.subtree_leaf_names(1); // AB clade
        assert_eq!(
            names,
            ["A", "B"].iter().map(|s| s.to_string()).collect::<std::collections::BTreeSet<_>>()
        );
    }

    #[test]
    fn total_branch_length_works() {
        let tree = sample_tree();
        // 0.1 + 0.2 + 0.3 + 0.4 + 0.5 + 0.6 = 2.1
        assert!((tree.total_branch_length() - 2.1).abs() < 1e-12);
    }

    #[test]
    fn reroot_preserves_leaf_set() {
        let tree = PhyloTree::from_newick("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);").unwrap();
        let rerooted = tree.reroot(3, Some(0.5)).unwrap(); // Reroot at leaf A's edge
        let mut orig = tree.leaf_names();
        let mut reroot = rerooted.leaf_names();
        orig.sort();
        reroot.sort();
        assert_eq!(orig, reroot);
    }

    #[test]
    fn reroot_preserves_total_branch_length() {
        let tree = PhyloTree::from_newick("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);").unwrap();
        let orig_bl = tree.total_branch_length();
        let rerooted = tree.reroot(3, Some(0.5)).unwrap();
        let new_bl = rerooted.total_branch_length();
        assert!(
            (orig_bl - new_bl).abs() < 1e-10,
            "branch length changed: {} -> {}",
            orig_bl,
            new_bl
        );
    }

    #[test]
    fn midpoint_root_balanced() {
        // Balanced tree: midpoint should keep it balanced.
        let tree =
            PhyloTree::from_newick("((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0);").unwrap();
        let rooted = tree.midpoint_root().unwrap();
        assert_eq!(rooted.leaf_count(), 4);
        let orig_bl = tree.total_branch_length();
        let new_bl = rooted.total_branch_length();
        assert!(
            (orig_bl - new_bl).abs() < 1e-10,
            "branch length changed: {} -> {}",
            orig_bl,
            new_bl
        );
    }

    #[test]
    fn extract_subtree_correct_leaves() {
        let tree =
            PhyloTree::from_newick("(((A:0.1,B:0.2):0.3,C:0.4):0.5,(D:0.6,E:0.7):0.8);")
                .unwrap();
        let sub = tree.extract_subtree(&["A", "B", "C"]).unwrap();
        let mut names = sub.leaf_names();
        names.sort();
        assert_eq!(names, vec!["A", "B", "C"]);
    }

    #[test]
    fn extract_subtree_preserves_topology() {
        let tree =
            PhyloTree::from_newick("(((A:0.1,B:0.2):0.3,C:0.4):0.5,(D:0.6,E:0.7):0.8);")
                .unwrap();
        let sub = tree.extract_subtree(&["D", "E"]).unwrap();
        assert_eq!(sub.leaf_count(), 2);
        let mut names = sub.leaf_names();
        names.sort();
        assert_eq!(names, vec!["D", "E"]);
    }

    #[test]
    fn subtree_at_extracts_clade() {
        let tree = sample_tree();
        // Node 2 is CD clade
        let sub = tree.subtree_at(2).unwrap();
        assert_eq!(sub.leaf_count(), 2);
        let mut names = sub.leaf_names();
        names.sort();
        assert_eq!(names, vec!["C", "D"]);
        // Root of subtree should not have a branch length
        assert!(sub.get_node(sub.root()).unwrap().branch_length.is_none());
    }
}
