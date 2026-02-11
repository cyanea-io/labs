//! Newick format parser and writer.
//!
//! Supports the standard Newick grammar:
//! ```text
//! tree     = subtree ';'
//! subtree  = '(' children ')' label | label
//! children = subtree (',' subtree)*
//! label    = name? (':' length)?
//! ```

use crate::tree::{Node, NodeId, PhyloTree};
use cyanea_core::{CyaneaError, Result};

/// Parse a Newick format string into a `PhyloTree`.
pub fn parse(input: &str) -> Result<PhyloTree> {
    let bytes = input.as_bytes();
    let mut parser = Parser::new(bytes);
    let (nodes, root) = parser.parse_tree()?;
    PhyloTree::from_nodes(nodes, root)
}

/// Serialize a `PhyloTree` to a Newick format string.
pub fn write(tree: &PhyloTree) -> String {
    let mut buf = String::new();
    write_subtree(tree, tree.root(), &mut buf);
    buf.push(';');
    buf
}

fn write_subtree(tree: &PhyloTree, id: NodeId, buf: &mut String) {
    let node = tree.get_node(id).unwrap();
    if !node.children.is_empty() {
        buf.push('(');
        for (i, &child) in node.children.iter().enumerate() {
            if i > 0 {
                buf.push(',');
            }
            write_subtree(tree, child, buf);
        }
        buf.push(')');
    }
    if let Some(ref name) = node.name {
        buf.push_str(name);
    }
    if let Some(len) = node.branch_length {
        buf.push(':');
        // Use enough precision but strip trailing zeros
        let s = format!("{:.10}", len);
        let s = s.trim_end_matches('0');
        let s = s.trim_end_matches('.');
        buf.push_str(s);
    }
}

struct Parser<'a> {
    input: &'a [u8],
    pos: usize,
    nodes: Vec<Node>,
}

impl<'a> Parser<'a> {
    fn new(input: &'a [u8]) -> Self {
        Self {
            input,
            pos: 0,
            nodes: Vec::new(),
        }
    }

    fn parse_tree(&mut self) -> Result<(Vec<Node>, NodeId)> {
        self.skip_whitespace();
        let root = self.parse_subtree(None)?;
        self.skip_whitespace();
        if self.pos >= self.input.len() || self.input[self.pos] != b';' {
            return Err(CyaneaError::Parse("expected ';' at end of Newick string".into()));
        }
        self.pos += 1;
        Ok((std::mem::take(&mut self.nodes), root))
    }

    fn parse_subtree(&mut self, parent: Option<NodeId>) -> Result<NodeId> {
        self.skip_whitespace();
        let id = self.alloc_node(parent);

        if self.peek() == Some(b'(') {
            self.pos += 1; // consume '('
            // Parse children
            let first_child = self.parse_subtree(Some(id))?;
            self.nodes[id].children.push(first_child);

            loop {
                self.skip_whitespace();
                if self.peek() == Some(b',') {
                    self.pos += 1;
                    let child = self.parse_subtree(Some(id))?;
                    self.nodes[id].children.push(child);
                } else {
                    break;
                }
            }
            self.skip_whitespace();
            if self.peek() != Some(b')') {
                return Err(CyaneaError::Parse("expected ')' in Newick string".into()));
            }
            self.pos += 1; // consume ')'
        }

        // Parse optional label (name:length)
        self.parse_label(id)?;
        Ok(id)
    }

    fn parse_label(&mut self, id: NodeId) -> Result<()> {
        self.skip_whitespace();
        // Parse name (everything until ':', ',', ')', ';', '(' or whitespace)
        let name = self.parse_name();
        if !name.is_empty() {
            self.nodes[id].name = Some(name);
        }
        // Parse optional branch length
        self.skip_whitespace();
        if self.peek() == Some(b':') {
            self.pos += 1;
            self.skip_whitespace();
            let len_str = self.parse_float_str();
            if len_str.is_empty() {
                return Err(CyaneaError::Parse("expected number after ':'".into()));
            }
            let len: f64 = len_str.parse().map_err(|_| {
                CyaneaError::Parse(format!("invalid branch length: '{}'", len_str))
            })?;
            self.nodes[id].branch_length = Some(len);
        }
        Ok(())
    }

    fn parse_name(&mut self) -> String {
        let start = self.pos;
        while self.pos < self.input.len() {
            match self.input[self.pos] {
                b':' | b',' | b')' | b'(' | b';' => break,
                b' ' | b'\t' | b'\n' | b'\r' => break,
                _ => self.pos += 1,
            }
        }
        String::from_utf8_lossy(&self.input[start..self.pos]).into_owned()
    }

    fn parse_float_str(&mut self) -> String {
        let start = self.pos;
        while self.pos < self.input.len() {
            match self.input[self.pos] {
                b'0'..=b'9' | b'.' | b'-' | b'+' | b'e' | b'E' => self.pos += 1,
                _ => break,
            }
        }
        String::from_utf8_lossy(&self.input[start..self.pos]).into_owned()
    }

    fn alloc_node(&mut self, parent: Option<NodeId>) -> NodeId {
        let id = self.nodes.len();
        self.nodes.push(Node {
            id,
            parent,
            children: Vec::new(),
            branch_length: None,
            name: None,
        });
        id
    }

    fn peek(&self) -> Option<u8> {
        self.input.get(self.pos).copied()
    }

    fn skip_whitespace(&mut self) {
        while self.pos < self.input.len() {
            match self.input[self.pos] {
                b' ' | b'\t' | b'\n' | b'\r' => self.pos += 1,
                _ => break,
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_simple_pair() {
        let tree = parse("(A,B);").unwrap();
        assert_eq!(tree.node_count(), 3);
        assert_eq!(tree.leaf_count(), 2);
        assert_eq!(tree.leaf_names(), vec!["A", "B"]);
    }

    #[test]
    fn parse_with_branch_lengths() {
        let tree = parse("(A:0.1,B:0.2):0.0;").unwrap();
        assert_eq!(tree.node_count(), 3);
        let root = tree.get_node(tree.root()).unwrap();
        assert_eq!(root.branch_length, Some(0.0));
        // Find leaf A
        let leaves = tree.leaves();
        let a = leaves
            .iter()
            .find(|&&id| tree.get_node(id).unwrap().name.as_deref() == Some("A"))
            .unwrap();
        assert_eq!(tree.get_node(*a).unwrap().branch_length, Some(0.1));
    }

    #[test]
    fn parse_nested() {
        let tree = parse("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);").unwrap();
        assert_eq!(tree.node_count(), 7);
        assert_eq!(tree.leaf_count(), 4);
        assert_eq!(tree.leaf_names(), vec!["A", "B", "C", "D"]);
    }

    #[test]
    fn parse_internal_names() {
        let tree = parse("((A,B)AB,(C,D)CD)root;").unwrap();
        let root = tree.get_node(tree.root()).unwrap();
        assert_eq!(root.name.as_deref(), Some("root"));
    }

    #[test]
    fn parse_single_leaf() {
        let tree = parse("A:1.5;").unwrap();
        assert_eq!(tree.node_count(), 1);
        assert_eq!(tree.leaf_count(), 1);
        let root = tree.get_node(tree.root()).unwrap();
        assert_eq!(root.name.as_deref(), Some("A"));
        assert_eq!(root.branch_length, Some(1.5));
    }

    #[test]
    fn parse_whitespace() {
        let tree = parse("  ( A : 0.1 , B : 0.2 ) ; ").unwrap();
        assert_eq!(tree.node_count(), 3);
        assert_eq!(tree.leaf_count(), 2);
    }

    #[test]
    fn parse_error_unbalanced_parens() {
        assert!(parse("((A,B);").is_err());
    }

    #[test]
    fn parse_error_missing_semicolon() {
        assert!(parse("(A,B)").is_err());
    }

    #[test]
    fn parse_error_bad_float() {
        assert!(parse("(A:abc,B);").is_err());
    }

    #[test]
    fn write_simple() {
        let tree = parse("(A,B);").unwrap();
        let s = write(&tree);
        assert_eq!(s, "(A,B);");
    }

    #[test]
    fn roundtrip_with_lengths() {
        let input = "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);";
        let tree = parse(input).unwrap();
        let output = write(&tree);
        // Re-parse to verify structural equivalence
        let tree2 = parse(&output).unwrap();
        assert_eq!(tree.node_count(), tree2.node_count());
        assert_eq!(tree.leaf_names(), tree2.leaf_names());
    }

    #[test]
    fn roundtrip_no_lengths() {
        let input = "((A,B),(C,D));";
        let tree = parse(input).unwrap();
        let output = write(&tree);
        assert_eq!(output, input);
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use proptest::prelude::*;

    /// Strategy for generating leaf names (simple alphanumeric, no special chars)
    fn leaf_name() -> impl Strategy<Value = String> {
        "[A-Za-z][A-Za-z0-9]{0,5}"
    }

    /// Strategy for a simple Newick tree with 2-6 leaves (no branch lengths)
    fn simple_newick() -> impl Strategy<Value = String> {
        proptest::collection::vec(leaf_name(), 2..=6)
            .prop_map(|leaves| {
                // Build a simple caterpillar tree
                if leaves.len() == 2 {
                    return format!("({},{});", leaves[0], leaves[1]);
                }
                let mut s = format!("({},{}", leaves[0], leaves[1]);
                for leaf in &leaves[2..] {
                    s = format!("({},{})", s, leaf);
                }
                s.push(';');
                s
            })
    }

    proptest! {
        #[test]
        fn newick_roundtrip_preserves_leaf_names(newick in simple_newick()) {
            if let Ok(tree) = parse(&newick) {
                let output = write(&tree);
                let tree2 = parse(&output).unwrap();
                let mut names1 = tree.leaf_names();
                let mut names2 = tree2.leaf_names();
                names1.sort();
                names2.sort();
                prop_assert_eq!(names1, names2);
            }
        }

        #[test]
        fn parse_newick_does_not_panic(s in "\\PC{0,100}") {
            let _ = parse(&s);
        }

        #[test]
        fn node_count_ge_leaf_count(newick in simple_newick()) {
            if let Ok(tree) = parse(&newick) {
                prop_assert!(tree.node_count() >= tree.leaf_count(),
                    "node_count={} < leaf_count={}", tree.node_count(), tree.leaf_count());
            }
        }
    }
}
