//! Gene tree / species tree methods: ASTRAL, reconciliation, concordance factors.
//!
//! Provides species tree estimation from gene trees, gene-tree/species-tree
//! reconciliation, and concordance factor computation.

use std::collections::{BTreeSet, HashMap, HashSet};

use crate::bootstrap::bipartitions;
use crate::tree::{NodeId, PhyloTree};
use cyanea_core::{CyaneaError, Result};

/// Result of gene tree / species tree reconciliation.
#[derive(Debug, Clone)]
pub struct ReconciliationResult {
    /// Node IDs in the gene tree identified as duplications.
    pub duplications: Vec<NodeId>,
    /// Number of inferred gene losses.
    pub losses: usize,
    /// Number of deep coalescences.
    pub deep_coalescences: usize,
    /// Total reconciliation cost (duplications + losses).
    pub cost: f64,
}

/// Gene and site concordance factors per branch.
#[derive(Debug, Clone)]
pub struct ConcordanceFactors {
    /// Per-branch gene concordance: (species_tree_node_id, fraction).
    pub gene_cf: Vec<(NodeId, f64)>,
    /// Per-branch site concordance: (species_tree_node_id, fraction).
    pub site_cf: Vec<(NodeId, f64)>,
}

/// Estimate a species tree from gene trees using quartet-based ASTRAL approach.
///
/// For each possible internal branch resolution, counts the number of
/// quartets supporting it across all gene trees. Builds the species tree
/// greedily by choosing the resolution with highest quartet support.
///
/// Practical for up to ~50 taxa.
pub fn astral_species_tree(gene_trees: &[PhyloTree]) -> Result<PhyloTree> {
    if gene_trees.is_empty() {
        return Err(CyaneaError::InvalidInput(
            "no gene trees provided".into(),
        ));
    }

    // Collect all taxa across gene trees.
    let mut all_taxa: BTreeSet<String> = BTreeSet::new();
    for gt in gene_trees {
        for name in gt.leaf_names() {
            all_taxa.insert(name);
        }
    }
    let taxa: Vec<String> = all_taxa.into_iter().collect();
    let n = taxa.len();

    if n < 4 {
        // With < 4 taxa, just use the first gene tree topology.
        return Ok(gene_trees[0].clone());
    }

    // Count quartet topologies across all gene trees.
    // For taxa (a,b,c,d), three possible unrooted topologies:
    //   ab|cd, ac|bd, ad|bc
    // Count which topology each gene tree supports.
    let mut quartet_counts: HashMap<(usize, usize, usize, usize), [usize; 3]> = HashMap::new();

    for gt in gene_trees {
        let gt_leaves: HashSet<String> = gt.leaf_names().into_iter().collect();

        // For efficiency, limit to a subset of quartets for larger datasets.
        let max_quartets = 5000;
        let mut count = 0;

        for i in 0..n {
            if count >= max_quartets {
                break;
            }
            for j in (i + 1)..n {
                for k in (j + 1)..n {
                    for l in (k + 1)..n {
                        if count >= max_quartets {
                            break;
                        }

                        // Check all 4 taxa exist in this gene tree.
                        if !gt_leaves.contains(&taxa[i])
                            || !gt_leaves.contains(&taxa[j])
                            || !gt_leaves.contains(&taxa[k])
                            || !gt_leaves.contains(&taxa[l])
                        {
                            continue;
                        }

                        // Determine which quartet topology this gene tree supports.
                        let topo = quartet_topology(gt, &taxa[i], &taxa[j], &taxa[k], &taxa[l]);
                        let entry = quartet_counts.entry((i, j, k, l)).or_insert([0; 3]);
                        entry[topo] += 1;
                        count += 1;
                    }
                }
            }
        }
    }

    // Build species tree greedily using bipartition scoring.
    // Score each potential bipartition by how many quartets support it.
    let mut bipartition_scores: HashMap<BTreeSet<String>, f64> = HashMap::new();

    for gt in gene_trees {
        let bps = bipartitions(gt);
        for bp in bps {
            // Only include bipartitions with taxa we know about.
            let filtered: BTreeSet<String> = bp
                .into_iter()
                .filter(|t| taxa.contains(t))
                .collect();
            if filtered.len() >= 2 && filtered.len() <= n - 2 {
                *bipartition_scores.entry(filtered).or_insert(0.0) += 1.0;
            }
        }
    }

    // Sort bipartitions by score (greedy).
    let mut scored: Vec<(BTreeSet<String>, f64)> = bipartition_scores.into_iter().collect();
    scored.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    // Build tree by adding compatible bipartitions.
    let mut accepted: Vec<BTreeSet<String>> = Vec::new();
    for (bp, _score) in &scored {
        // Check compatibility: two bipartitions are compatible if one is a subset
        // of the other, or they are disjoint, or their complements are disjoint.
        let compatible = accepted.iter().all(|existing| {
            is_compatible(bp, existing, &taxa)
        });

        if compatible {
            accepted.push(bp.clone());
        }

        // Max n-3 internal bipartitions for n taxa.
        if accepted.len() >= n - 3 {
            break;
        }
    }

    // Build tree from accepted bipartitions.
    build_tree_from_bipartitions(&taxa, &accepted)
}

/// Reconcile a gene tree with a species tree.
///
/// Uses LCA mapping: each gene tree node is mapped to the most recent
/// common ancestor of its leaves in the species tree. Gene tree nodes
/// mapped above their children's mappings indicate duplications.
pub fn reconcile(
    gene_tree: &PhyloTree,
    species_tree: &PhyloTree,
    gene_to_species: &HashMap<String, String>,
) -> Result<ReconciliationResult> {
    let gt_n = gene_tree.node_count();
    let _st_leaves: HashSet<String> = species_tree.leaf_names().into_iter().collect();

    // Map gene tree leaves to species tree leaves.
    let mut gt_to_st: Vec<Option<NodeId>> = vec![None; gt_n];
    for &leaf_id in &gene_tree.leaves() {
        if let Some(name) = gene_tree.get_node(leaf_id).and_then(|n| n.name.as_ref()) {
            let species_name = gene_to_species
                .get(name)
                .cloned()
                .unwrap_or_else(|| name.clone());

            // Find corresponding leaf in species tree.
            for &st_leaf in &species_tree.leaves() {
                if let Some(st_name) = species_tree.get_node(st_leaf).and_then(|n| n.name.as_ref())
                {
                    if *st_name == species_name {
                        gt_to_st[leaf_id] = Some(st_leaf);
                        break;
                    }
                }
            }
        }
    }

    // LCA mapping: post-order, map internal nodes to LCA of children's mappings.
    for id in gene_tree.iter_postorder() {
        let node = gene_tree.get_node(id).unwrap();
        if node.is_leaf() {
            continue;
        }

        let child_mappings: Vec<NodeId> = node
            .children
            .iter()
            .filter_map(|&c| gt_to_st[c])
            .collect();

        if child_mappings.is_empty() {
            continue;
        }

        // Compute LCA of all child mappings.
        let mut lca = child_mappings[0];
        for &m in &child_mappings[1..] {
            lca = species_tree.mrca(lca, m)?;
        }
        gt_to_st[id] = Some(lca);
    }

    // Identify duplications and count losses.
    let mut duplications = Vec::new();
    let mut losses = 0;
    let mut deep_coalescences = 0;

    for id in gene_tree.iter_postorder() {
        let node = gene_tree.get_node(id).unwrap();
        if node.is_leaf() || gt_to_st[id].is_none() {
            continue;
        }

        let my_mapping = gt_to_st[id].unwrap();

        // Check if any child maps to the same species tree node (= duplication).
        let is_dup = node.children.iter().any(|&c| gt_to_st[c] == Some(my_mapping));
        if is_dup {
            duplications.push(id);
        }

        // Count losses: number of species tree edges between parent mapping
        // and child mapping.
        for &child_id in &node.children {
            if let Some(child_mapping) = gt_to_st[child_id] {
                let depth_diff = count_edges_between(species_tree, my_mapping, child_mapping);
                if depth_diff > 1 {
                    losses += depth_diff - 1;
                    deep_coalescences += depth_diff - 1;
                }
            }
        }
    }

    let cost = duplications.len() as f64 + losses as f64;

    Ok(ReconciliationResult {
        duplications,
        losses,
        deep_coalescences,
        cost,
    })
}

/// Compute gene and site concordance factors.
///
/// Gene CF: for each branch in the species tree, fraction of gene trees
/// containing the same bipartition.
///
/// Site CF: for each branch, fraction of sites supporting that topology
/// (requires sequences).
pub fn concordance_factors(
    species_tree: &PhyloTree,
    gene_trees: &[PhyloTree],
    sequences: Option<&[&[u8]]>,
) -> Result<ConcordanceFactors> {
    let n_gene_trees = gene_trees.len();

    // Gene concordance factor.
    let mut gene_cf = Vec::new();

    // For each internal node in species tree, compute gCF.
    for node_id in species_tree.iter_preorder() {
        let node = species_tree.get_node(node_id).unwrap();
        if node.is_leaf() || node.is_root() {
            continue;
        }

        let sp_bp = species_tree.subtree_leaf_names(node_id);
        if sp_bp.len() <= 1 {
            continue;
        }

        let mut concordant = 0;
        for gt in gene_trees {
            let gt_bps = bipartitions(gt);
            if gt_bps.contains(&sp_bp) {
                concordant += 1;
            }
        }

        let gcf = if n_gene_trees > 0 {
            concordant as f64 / n_gene_trees as f64
        } else {
            0.0
        };
        gene_cf.push((node_id, gcf));
    }

    // Site concordance factor (if sequences provided).
    let mut site_cf = Vec::new();
    if let Some(seqs) = sequences {
        if !seqs.is_empty() {
            let seq_len = seqs[0].len();
            let sp_leaves = species_tree.leaf_names();

            for node_id in species_tree.iter_preorder() {
                let node = species_tree.get_node(node_id).unwrap();
                if node.is_leaf() || node.is_root() {
                    continue;
                }

                let sp_bp = species_tree.subtree_leaf_names(node_id);
                if sp_bp.len() <= 1 || sp_bp.len() >= sp_leaves.len() - 1 {
                    continue;
                }

                // For each site, check if it supports this bipartition.
                // A site supports the split if taxa in the bipartition share
                // the same state and differ from taxa outside.
                let mut supporting = 0;
                let mut decisive = 0;

                for site in 0..seq_len {
                    // Get states for taxa in and out of the bipartition.
                    let mut in_states = HashSet::new();
                    let mut out_states = HashSet::new();

                    for (i, name) in sp_leaves.iter().enumerate() {
                        if i < seqs.len() {
                            let state = seqs[i][site];
                            if sp_bp.contains(name) {
                                in_states.insert(state);
                            } else {
                                out_states.insert(state);
                            }
                        }
                    }

                    // Decisive: both partitions have informative states.
                    if !in_states.is_empty() && !out_states.is_empty() {
                        decisive += 1;
                        // Supporting: in-states and out-states are disjoint.
                        if in_states.is_disjoint(&out_states) {
                            supporting += 1;
                        }
                    }
                }

                let scf = if decisive > 0 {
                    supporting as f64 / decisive as f64
                } else {
                    0.0
                };
                site_cf.push((node_id, scf));
            }
        }
    }

    Ok(ConcordanceFactors { gene_cf, site_cf })
}

// ── Helpers ────────────────────────────────────────────────────────────────

/// Determine quartet topology for taxa (a,b,c,d) in a tree.
/// Returns 0 for ab|cd, 1 for ac|bd, 2 for ad|bc.
fn quartet_topology(
    tree: &PhyloTree,
    a: &str,
    b: &str,
    c: &str,
    d: &str,
) -> usize {
    // Find leaf nodes for each taxon.
    let find_leaf = |name: &str| -> Option<NodeId> {
        for &id in &tree.leaves() {
            if let Some(n) = tree.get_node(id).and_then(|n| n.name.as_ref()) {
                if n == name {
                    return Some(id);
                }
            }
        }
        None
    };

    let (la, lb, lc, ld) = match (find_leaf(a), find_leaf(b), find_leaf(c), find_leaf(d)) {
        (Some(a), Some(b), Some(c), Some(d)) => (a, b, c, d),
        _ => return 0, // Default if taxa not found
    };

    // Compare MRCA depths. The topology is determined by which pair
    // has the most recent common ancestor.
    let mrca_ab = tree.mrca(la, lb).unwrap_or(tree.root());
    let mrca_cd = tree.mrca(lc, ld).unwrap_or(tree.root());
    let mrca_ac = tree.mrca(la, lc).unwrap_or(tree.root());
    let mrca_bd = tree.mrca(lb, ld).unwrap_or(tree.root());
    let mrca_ad = tree.mrca(la, ld).unwrap_or(tree.root());
    let mrca_bc = tree.mrca(lb, lc).unwrap_or(tree.root());

    // The topology where the two MRCAs are lowest (most specific) wins.
    let depth_ab_cd = node_depth(tree, mrca_ab) + node_depth(tree, mrca_cd);
    let depth_ac_bd = node_depth(tree, mrca_ac) + node_depth(tree, mrca_bd);
    let depth_ad_bc = node_depth(tree, mrca_ad) + node_depth(tree, mrca_bc);

    // Higher depth sum = more resolved = this is the topology.
    if depth_ab_cd >= depth_ac_bd && depth_ab_cd >= depth_ad_bc {
        0 // ab|cd
    } else if depth_ac_bd >= depth_ad_bc {
        1 // ac|bd
    } else {
        2 // ad|bc
    }
}

/// Compute depth (distance from root) of a node.
fn node_depth(tree: &PhyloTree, node_id: NodeId) -> usize {
    let mut depth = 0;
    let mut cur = node_id;
    while let Some(parent) = tree.get_node(cur).and_then(|n| n.parent) {
        depth += 1;
        cur = parent;
    }
    depth
}

/// Count edges between two nodes in a tree.
fn count_edges_between(tree: &PhyloTree, a: NodeId, b: NodeId) -> usize {
    // Walk from a to root, collecting ancestors.
    let mut ancestors_a = Vec::new();
    let mut cur = a;
    loop {
        ancestors_a.push(cur);
        match tree.get_node(cur).and_then(|n| n.parent) {
            Some(p) => cur = p,
            None => break,
        }
    }

    // Walk from b upward until we find a common ancestor.
    cur = b;
    let mut depth_b = 0;
    loop {
        if let Some(pos) = ancestors_a.iter().position(|&x| x == cur) {
            return pos + depth_b;
        }
        match tree.get_node(cur).and_then(|n| n.parent) {
            Some(p) => {
                cur = p;
                depth_b += 1;
            }
            None => break,
        }
    }
    0
}

/// Check if two bipartitions are compatible.
fn is_compatible(bp1: &BTreeSet<String>, bp2: &BTreeSet<String>, all_taxa: &[String]) -> bool {
    let comp1: BTreeSet<String> = all_taxa
        .iter()
        .filter(|t| !bp1.contains(*t))
        .cloned()
        .collect();
    let comp2: BTreeSet<String> = all_taxa
        .iter()
        .filter(|t| !bp2.contains(*t))
        .cloned()
        .collect();

    // Compatible if one of the four intersections is empty.
    bp1.intersection(bp2).next().is_none()
        || bp1.intersection(&comp2).next().is_none()
        || comp1.intersection(bp2).next().is_none()
        || comp1.intersection(&comp2).next().is_none()
}

/// Build a tree from a set of compatible bipartitions.
fn build_tree_from_bipartitions(
    taxa: &[String],
    bipartitions: &[BTreeSet<String>],
) -> Result<PhyloTree> {
    // Start with a star tree, then refine by adding bipartitions.
    let mut tree = PhyloTree::new();
    for taxon in taxa {
        tree.add_child(0, Some(taxon.clone()), Some(0.1))?;
    }

    // For each bipartition, create an internal node grouping those taxa.
    for bp in bipartitions {
        if bp.len() < 2 || bp.len() >= taxa.len() - 1 {
            continue;
        }

        // Find leaves in the bipartition that are currently children of root.
        let root = tree.root();
        let root_children: Vec<NodeId> = tree.get_node(root).unwrap().children.clone();

        let bp_children: Vec<NodeId> = root_children
            .iter()
            .filter(|&&c| {
                let names = tree.subtree_leaf_names(c);
                names.is_subset(bp) && !names.is_empty()
            })
            .copied()
            .collect();

        if bp_children.len() < 2 {
            continue;
        }

        // Create new internal node.
        let new_internal = tree.add_child(root, None, Some(0.1))?;

        // Move children to new internal.
        for &child_id in &bp_children {
            // Update parent pointer.
            if let Some(node) = tree.get_node_mut(child_id) {
                node.parent = Some(new_internal);
            }
            // Add to new internal's children.
            if let Some(internal) = tree.get_node_mut(new_internal) {
                internal.children.push(child_id);
            }
        }

        // Remove moved children from root.
        if let Some(root_node) = tree.get_node_mut(root) {
            root_node.children.retain(|c| !bp_children.contains(c));
        }
    }

    Ok(tree)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn astral_identical_gene_trees() {
        let t1 = PhyloTree::from_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();
        let t2 = PhyloTree::from_newick("((A:0.2,B:0.2):0.2,(C:0.2,D:0.2):0.2);").unwrap();
        let t3 = PhyloTree::from_newick("((A:0.3,B:0.3):0.3,(C:0.3,D:0.3):0.3);").unwrap();

        let result = astral_species_tree(&[t1, t2, t3]).unwrap();
        assert_eq!(result.leaf_count(), 4);
        let mut names = result.leaf_names();
        names.sort();
        assert_eq!(names, vec!["A", "B", "C", "D"]);
    }

    #[test]
    fn reconciliation_matching_trees() {
        let gene = PhyloTree::from_newick("((a:0.1,b:0.1):0.1,(c:0.1,d:0.1):0.1);").unwrap();
        let species =
            PhyloTree::from_newick("((a:0.1,b:0.1):0.1,(c:0.1,d:0.1):0.1);").unwrap();
        let mapping = HashMap::new(); // Identity mapping (names match)

        let result = reconcile(&gene, &species, &mapping).unwrap();
        assert_eq!(result.duplications.len(), 0, "no duplications expected");
        assert_eq!(result.losses, 0, "no losses expected");
    }

    #[test]
    fn reconciliation_detects_duplication() {
        // Gene tree has (a1,a2) which both map to species 'a'.
        let gene = PhyloTree::from_newick("((a1:0.1,a2:0.1):0.1,b:0.1);").unwrap();
        let species = PhyloTree::from_newick("(a:0.1,b:0.1);").unwrap();

        let mut mapping = HashMap::new();
        mapping.insert("a1".to_string(), "a".to_string());
        mapping.insert("a2".to_string(), "a".to_string());
        mapping.insert("b".to_string(), "b".to_string());

        let result = reconcile(&gene, &species, &mapping).unwrap();
        assert!(
            result.duplications.len() > 0,
            "should detect duplication"
        );
    }

    #[test]
    fn gcf_all_agree() {
        let species =
            PhyloTree::from_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();
        let gt1 = PhyloTree::from_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();
        let gt2 = PhyloTree::from_newick("((A:0.2,B:0.2):0.2,(C:0.2,D:0.2):0.2);").unwrap();

        let cf = concordance_factors(&species, &[gt1, gt2], None).unwrap();
        for &(_, gcf) in &cf.gene_cf {
            assert!(
                (gcf - 1.0).abs() < 1e-10,
                "gCF should be 1.0 when all agree, got {}",
                gcf
            );
        }
    }

    #[test]
    fn gcf_with_discordance() {
        let species =
            PhyloTree::from_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();
        let gt_agree =
            PhyloTree::from_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();
        let gt_discord =
            PhyloTree::from_newick("((A:0.1,C:0.1):0.1,(B:0.1,D:0.1):0.1);").unwrap();

        let cf = concordance_factors(&species, &[gt_agree, gt_discord], None).unwrap();
        // At least one branch should have gCF < 1.0
        let has_discordant = cf.gene_cf.iter().any(|&(_, gcf)| gcf < 1.0);
        assert!(has_discordant, "should have discordant gCF values");
    }

    #[test]
    fn scf_computed() {
        let species =
            PhyloTree::from_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();
        let seqs: Vec<Vec<u8>> = vec![
            b"AAACCCAAA".to_vec(),
            b"AAACCCAAA".to_vec(),
            b"CCCAAACCC".to_vec(),
            b"CCCAAACCC".to_vec(),
        ];
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();

        let cf = concordance_factors(&species, &[], Some(&refs)).unwrap();
        assert!(!cf.site_cf.is_empty(), "should have site CF values");
        for &(_, scf) in &cf.site_cf {
            assert!(
                (0.0..=1.0).contains(&scf),
                "sCF should be in [0,1], got {}",
                scf
            );
        }
    }

    #[test]
    fn empty_gene_tree_list() {
        assert!(astral_species_tree(&[]).is_err());
    }

    #[test]
    fn astral_five_taxa() {
        let gt1 = PhyloTree::from_newick(
            "(((A:0.1,B:0.1):0.1,C:0.1):0.1,(D:0.1,E:0.1):0.1);",
        )
        .unwrap();
        let gt2 = PhyloTree::from_newick(
            "(((A:0.1,B:0.1):0.1,C:0.1):0.1,(D:0.1,E:0.1):0.1);",
        )
        .unwrap();

        let result = astral_species_tree(&[gt1, gt2]).unwrap();
        assert_eq!(result.leaf_count(), 5);
    }
}
