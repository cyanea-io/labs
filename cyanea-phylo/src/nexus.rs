//! NEXUS format parser and writer.
//!
//! Supports the core NEXUS blocks used in phylogenetics:
//!
//! - **TAXA block** — taxon names and count
//! - **TREES block** — phylogenetic trees in Newick format, with optional
//!   translate tables that map numeric labels to taxon names
//!
//! # Example
//!
//! ```
//! use cyanea_phylo::nexus;
//!
//! let input = r#"#NEXUS
//! BEGIN TAXA;
//!   DIMENSIONS NTAX=3;
//!   TAXLABELS A B C;
//! END;
//! BEGIN TREES;
//!   TREE tree1 = ((A,B),C);
//! END;
//! "#;
//!
//! let nexus = nexus::parse(input).unwrap();
//! assert_eq!(nexus.taxa, vec!["A", "B", "C"]);
//! assert_eq!(nexus.trees.len(), 1);
//! ```

use cyanea_core::{CyaneaError, Result};

use crate::tree::PhyloTree;

/// A parsed NEXUS file.
#[derive(Debug, Clone)]
pub struct NexusFile {
    /// Taxon names from the TAXA block.
    pub taxa: Vec<String>,
    /// Named trees from the TREES block.
    pub trees: Vec<NamedTree>,
}

/// A single named tree from a TREES block.
#[derive(Debug, Clone)]
pub struct NamedTree {
    /// Tree name (e.g. `"tree1"`).
    pub name: String,
    /// The phylogenetic tree.
    pub tree: PhyloTree,
}

/// Parse a NEXUS format string.
///
/// # Errors
///
/// Returns an error if the input is not valid NEXUS or contains unparseable trees.
pub fn parse(input: &str) -> Result<NexusFile> {
    let trimmed = input.trim();
    if !trimmed.starts_with("#NEXUS") && !trimmed.starts_with("#nexus") {
        return Err(CyaneaError::Parse("not a NEXUS file (missing #NEXUS header)".into()));
    }

    let mut taxa = Vec::new();
    let mut trees = Vec::new();
    let mut translate_table: Vec<(String, String)> = Vec::new();

    let blocks = extract_blocks(trimmed);

    for (block_type, block_content) in &blocks {
        match block_type.to_uppercase().as_str() {
            "TAXA" => {
                taxa = parse_taxa_block(block_content);
            }
            "TREES" => {
                let (trans, named_trees) = parse_trees_block(block_content)?;
                translate_table = trans;
                trees = named_trees;
            }
            _ => {
                // Ignore unknown blocks
            }
        }
    }

    // Apply translate table to tree leaf names
    if !translate_table.is_empty() {
        for nt in &mut trees {
            apply_translate(&mut nt.tree, &translate_table);
        }
    }

    Ok(NexusFile { taxa, trees })
}

/// Write a NEXUS file from taxa and trees.
pub fn write(taxa: &[String], trees: &[(&str, &PhyloTree)]) -> String {
    let mut buf = String::from("#NEXUS\n\n");

    if !taxa.is_empty() {
        buf.push_str("BEGIN TAXA;\n");
        buf.push_str(&format!("  DIMENSIONS NTAX={};\n", taxa.len()));
        buf.push_str("  TAXLABELS");
        for t in taxa {
            buf.push(' ');
            buf.push_str(t);
        }
        buf.push_str(";\n");
        buf.push_str("END;\n\n");
    }

    if !trees.is_empty() {
        buf.push_str("BEGIN TREES;\n");
        for (name, tree) in trees {
            let newick = crate::newick::write(tree);
            buf.push_str(&format!("  TREE {} = {};\n", name, newick));
        }
        buf.push_str("END;\n");
    }

    buf
}

fn extract_blocks(input: &str) -> Vec<(String, String)> {
    let mut blocks = Vec::new();
    let upper = input.to_uppercase();
    let mut pos = 0;

    while let Some(begin_pos) = upper[pos..].find("BEGIN ") {
        let abs_begin = pos + begin_pos + 6;
        if let Some(semi) = upper[abs_begin..].find(';') {
            let block_type = input[abs_begin..abs_begin + semi].trim().to_string();
            let content_start = abs_begin + semi + 1;

            if let Some(end_pos) = upper[content_start..].find("END;") {
                let content = input[content_start..content_start + end_pos].to_string();
                blocks.push((block_type, content));
                pos = content_start + end_pos + 4;
            } else {
                break;
            }
        } else {
            break;
        }
    }

    blocks
}

fn parse_taxa_block(content: &str) -> Vec<String> {
    let upper = content.to_uppercase();
    let mut taxa = Vec::new();

    if let Some(tl_pos) = upper.find("TAXLABELS") {
        let after = &content[tl_pos + 9..];
        if let Some(semi_pos) = after.find(';') {
            let labels_str = &after[..semi_pos];
            for token in labels_str.split_whitespace() {
                let t = token.trim_matches(|c: char| c == '\'' || c == '"');
                if !t.is_empty() {
                    taxa.push(t.to_string());
                }
            }
        }
    }

    taxa
}

fn parse_trees_block(content: &str) -> Result<(Vec<(String, String)>, Vec<NamedTree>)> {
    let upper = content.to_uppercase();
    let mut translate_table = Vec::new();
    let mut trees = Vec::new();

    // Parse TRANSLATE if present
    if let Some(trans_pos) = upper.find("TRANSLATE") {
        let after = &content[trans_pos + 9..];
        if let Some(semi_pos) = after.find(';') {
            let trans_str = &after[..semi_pos];
            for entry in trans_str.split(',') {
                let parts: Vec<&str> = entry.split_whitespace().collect();
                if parts.len() >= 2 {
                    let num = parts[0].trim().to_string();
                    let name = parts[1].trim().trim_matches(|c: char| c == '\'' || c == '"').to_string();
                    translate_table.push((num, name));
                }
            }
        }
    }

    // Parse TREE statements
    let mut search_pos = 0;
    while let Some(tree_pos) = upper[search_pos..].find("TREE ") {
        // Skip TRANSLATE keyword match
        let abs_pos = search_pos + tree_pos;
        if abs_pos > 0 && upper[..abs_pos].ends_with("TRANS") {
            search_pos = abs_pos + 5;
            continue;
        }

        let after = &content[abs_pos + 5..];
        if let Some(eq_pos) = after.find('=') {
            let name = after[..eq_pos].trim().to_string();
            // Skip optional [&R] or [&U] annotations
            let newick_part = &after[eq_pos + 1..];
            if let Some(semi_pos) = newick_part.find(';') {
                let mut newick = newick_part[..semi_pos].trim().to_string();
                // Strip [&...] annotations
                while let Some(bracket_start) = newick.find('[') {
                    if let Some(bracket_end) = newick[bracket_start..].find(']') {
                        newick = format!(
                            "{}{}",
                            &newick[..bracket_start],
                            &newick[bracket_start + bracket_end + 1..]
                        );
                    } else {
                        break;
                    }
                }
                // Add semicolon for the Newick parser
                newick.push(';');
                let tree = crate::newick::parse(&newick)?;
                trees.push(NamedTree {
                    name: name.trim_start_matches('*').trim().to_string(),
                    tree,
                });
            }
        }
        search_pos = abs_pos + 5;
    }

    Ok((translate_table, trees))
}

fn apply_translate(tree: &mut PhyloTree, table: &[(String, String)]) {
    for id in 0..tree.node_count() {
        if let Some(node) = tree.get_node_mut(id) {
            if let Some(ref name) = node.name {
                for (num, taxon) in table {
                    if name == num {
                        node.name = Some(taxon.clone());
                        break;
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_basic_nexus() {
        let input = "#NEXUS\nBEGIN TAXA;\n  DIMENSIONS NTAX=3;\n  TAXLABELS A B C;\nEND;\nBEGIN TREES;\n  TREE t1 = ((A,B),C);\nEND;\n";
        let nexus = parse(input).unwrap();
        assert_eq!(nexus.taxa, vec!["A", "B", "C"]);
        assert_eq!(nexus.trees.len(), 1);
        assert_eq!(nexus.trees[0].name, "t1");
        assert_eq!(nexus.trees[0].tree.leaf_count(), 3);
    }

    #[test]
    fn parse_missing_header() {
        assert!(parse("BEGIN TAXA; END;").is_err());
    }

    #[test]
    fn parse_with_branch_lengths() {
        let input = "#NEXUS\nBEGIN TREES;\n  TREE t1 = ((A:0.1,B:0.2):0.3,C:0.4);\nEND;\n";
        let nexus = parse(input).unwrap();
        assert_eq!(nexus.trees[0].tree.leaf_count(), 3);
    }

    #[test]
    fn parse_with_translate() {
        let input = "#NEXUS\nBEGIN TREES;\n  TRANSLATE\n    1 Alpha,\n    2 Beta,\n    3 Gamma;\n  TREE t1 = ((1,2),3);\nEND;\n";
        let nexus = parse(input).unwrap();
        let names = nexus.trees[0].tree.leaf_names();
        assert!(names.contains(&"Alpha".to_string()));
        assert!(names.contains(&"Beta".to_string()));
        assert!(names.contains(&"Gamma".to_string()));
    }

    #[test]
    fn parse_multiple_trees() {
        let input = "#NEXUS\nBEGIN TREES;\n  TREE t1 = (A,B);\n  TREE t2 = (C,D);\nEND;\n";
        let nexus = parse(input).unwrap();
        assert_eq!(nexus.trees.len(), 2);
    }

    #[test]
    fn parse_taxa_only() {
        let input = "#NEXUS\nBEGIN TAXA;\n  DIMENSIONS NTAX=2;\n  TAXLABELS X Y;\nEND;\n";
        let nexus = parse(input).unwrap();
        assert_eq!(nexus.taxa, vec!["X", "Y"]);
        assert!(nexus.trees.is_empty());
    }

    #[test]
    fn write_roundtrip() {
        let tree = PhyloTree::from_newick("((A,B),C);").unwrap();
        let taxa = vec!["A".into(), "B".into(), "C".into()];
        let output = write(&taxa, &[("t1", &tree)]);
        assert!(output.contains("#NEXUS"));
        assert!(output.contains("NTAX=3"));
        assert!(output.contains("TREE t1"));
    }

    #[test]
    fn write_empty() {
        let output = write(&[], &[]);
        assert!(output.starts_with("#NEXUS"));
    }
}
