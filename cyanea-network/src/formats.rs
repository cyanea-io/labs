//! Network file format parsing and writing.
//!
//! Supports GraphML, SIF (Simple Interaction Format), and GEXF for
//! exchanging biological network data.

use crate::graph::{Edge, Graph, GraphType, Node};
use cyanea_core::{CyaneaError, Result};
use std::collections::HashMap;

// ──────────────────────────── SIF ────────────────────────────

/// Parse SIF (Simple Interaction Format).
///
/// Each line: `nodeA interaction_type nodeB [nodeC ...]`
///
/// SIF is commonly used by Cytoscape.
pub fn parse_sif(content: &str) -> Result<Graph> {
    let mut graph = Graph::new(GraphType::Undirected);
    let mut seen_nodes = std::collections::HashSet::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 3 {
            // Single node line
            if parts.len() == 1 && !seen_nodes.contains(parts[0]) {
                seen_nodes.insert(parts[0].to_string());
                graph.add_node(parts[0], parts[0])?;
            }
            continue;
        }

        let source = parts[0];
        let interaction = parts[1];

        if !seen_nodes.contains(source) {
            seen_nodes.insert(source.to_string());
            graph.add_node(source, source)?;
        }

        for &target in &parts[2..] {
            if !seen_nodes.contains(target) {
                seen_nodes.insert(target.to_string());
                graph.add_node(target, target)?;
            }

            let edge = Edge::new(source, target).with_type(interaction);
            graph.add_edge_struct(edge)?;
        }
    }

    Ok(graph)
}

/// Write a graph in SIF format.
pub fn write_sif(graph: &Graph) -> String {
    let mut lines = Vec::new();

    // Nodes with edges
    let mut nodes_with_edges = std::collections::HashSet::new();
    for edge in &graph.edges {
        let edge_type = edge.edge_type.as_deref().unwrap_or("interacts_with");
        lines.push(format!("{}\t{}\t{}", edge.source, edge_type, edge.target));
        nodes_with_edges.insert(edge.source.clone());
        nodes_with_edges.insert(edge.target.clone());
    }

    // Isolated nodes
    for id in graph.node_ids() {
        if !nodes_with_edges.contains(id) {
            lines.push(id.to_string());
        }
    }

    lines.join("\n")
}

// ──────────────────────────── GraphML ────────────────────────────

/// Parse a simplified GraphML format.
///
/// Supports `<node>`, `<edge>`, `<data>` elements for node/edge attributes,
/// and graph-level `edgedefault` for directed/undirected.
pub fn parse_graphml(content: &str) -> Result<Graph> {
    // Determine graph type
    let graph_type = if content.contains("edgedefault=\"directed\"") {
        GraphType::Directed
    } else {
        GraphType::Undirected
    };

    let mut graph = Graph::new(graph_type);
    let mut key_names: HashMap<String, String> = HashMap::new();

    // Parse key definitions
    for line in content.lines() {
        let line = line.trim();
        if line.starts_with("<key") {
            if let (Some(id), Some(name)) = (extract_attr(line, "id"), extract_attr(line, "attr.name")) {
                key_names.insert(id, name);
            }
        }
    }

    // Parse nodes
    let mut in_node = false;
    let mut current_node_id = String::new();
    let mut current_attrs: HashMap<String, String> = HashMap::new();
    let mut node_label = String::new();

    for line in content.lines() {
        let line = line.trim();

        if line.starts_with("<node") && !line.starts_with("</node") {
            in_node = true;
            current_node_id = extract_attr(line, "id").unwrap_or_default();
            node_label = current_node_id.clone();
            current_attrs.clear();

            if line.contains("/>") {
                // Self-closing node
                let mut node = Node::new(&current_node_id, &node_label);
                node.attributes = current_attrs.clone();
                graph.add_node_struct(node)?;
                in_node = false;
            }
        } else if line.starts_with("</node") {
            let mut node = Node::new(&current_node_id, &node_label);
            node.attributes = current_attrs.clone();
            graph.add_node_struct(node)?;
            in_node = false;
        } else if in_node && line.starts_with("<data") {
            if let Some(key) = extract_attr(line, "key") {
                let value = extract_text_content(line);
                let attr_name = key_names.get(&key).unwrap_or(&key).clone();
                if attr_name == "label" {
                    node_label = value.clone();
                }
                current_attrs.insert(attr_name, value);
            }
        }
    }

    // Parse edges
    for line in content.lines() {
        let line = line.trim();
        if line.starts_with("<edge") {
            let source = extract_attr(line, "source").unwrap_or_default();
            let target = extract_attr(line, "target").unwrap_or_default();

            if !source.is_empty() && !target.is_empty() {
                let weight_str = extract_attr(line, "weight");
                let weight = weight_str
                    .and_then(|s| s.parse::<f64>().ok())
                    .unwrap_or(1.0);

                graph.add_edge(&source, &target, weight)?;
            }
        }
    }

    Ok(graph)
}

/// Write a graph in GraphML format.
pub fn write_graphml(graph: &Graph) -> String {
    let mut lines = Vec::new();

    lines.push(r#"<?xml version="1.0" encoding="UTF-8"?>"#.to_string());
    lines.push(r#"<graphml xmlns="http://graphml.graphstruct.org/graphml">"#.to_string());

    // Key definitions
    lines.push(r#"  <key id="label" for="node" attr.name="label" attr.type="string"/>"#.to_string());
    lines.push(r#"  <key id="weight" for="edge" attr.name="weight" attr.type="double"/>"#.to_string());

    // Collect all node attribute keys
    let mut node_attr_keys = std::collections::HashSet::new();
    for node in graph.nodes.values() {
        for key in node.attributes.keys() {
            if key != "label" {
                node_attr_keys.insert(key.clone());
            }
        }
    }
    for key in &node_attr_keys {
        lines.push(format!(
            r#"  <key id="{}" for="node" attr.name="{}" attr.type="string"/>"#,
            key, key
        ));
    }

    let edge_default = match graph.graph_type {
        GraphType::Directed => "directed",
        GraphType::Undirected => "undirected",
    };
    lines.push(format!(r#"  <graph edgedefault="{}">"#, edge_default));

    // Nodes
    let mut node_ids: Vec<&String> = graph.nodes.keys().collect();
    node_ids.sort();
    for id in node_ids {
        let node = &graph.nodes[id];
        lines.push(format!(r#"    <node id="{}">"#, node.id));
        lines.push(format!(
            r#"      <data key="label">{}</data>"#,
            node.label
        ));
        for (key, value) in &node.attributes {
            lines.push(format!(r#"      <data key="{}">{}</data>"#, key, value));
        }
        lines.push("    </node>".to_string());
    }

    // Edges
    for (i, edge) in graph.edges.iter().enumerate() {
        lines.push(format!(
            r#"    <edge id="e{}" source="{}" target="{}">"#,
            i, edge.source, edge.target
        ));
        lines.push(format!(
            r#"      <data key="weight">{}</data>"#,
            edge.weight
        ));
        lines.push("    </edge>".to_string());
    }

    lines.push("  </graph>".to_string());
    lines.push("</graphml>".to_string());

    lines.join("\n")
}

// ──────────────────────────── GEXF ────────────────────────────

/// Write a graph in GEXF format (Gephi exchange format).
pub fn write_gexf(graph: &Graph) -> String {
    let mut lines = Vec::new();

    lines.push(r#"<?xml version="1.0" encoding="UTF-8"?>"#.to_string());
    lines.push(r#"<gexf xmlns="http://gexf.net/1.3" version="1.3">"#.to_string());

    let edge_type = match graph.graph_type {
        GraphType::Directed => "directed",
        GraphType::Undirected => "undirected",
    };
    lines.push(format!(r#"  <graph defaultedgetype="{}">"#, edge_type));

    // Nodes
    lines.push("    <nodes>".to_string());
    let mut node_ids: Vec<&String> = graph.nodes.keys().collect();
    node_ids.sort();
    for id in node_ids {
        let node = &graph.nodes[id];
        lines.push(format!(
            r#"      <node id="{}" label="{}" />"#,
            node.id, node.label
        ));
    }
    lines.push("    </nodes>".to_string());

    // Edges
    lines.push("    <edges>".to_string());
    for (i, edge) in graph.edges.iter().enumerate() {
        lines.push(format!(
            r#"      <edge id="{}" source="{}" target="{}" weight="{}" />"#,
            i, edge.source, edge.target, edge.weight
        ));
    }
    lines.push("    </edges>".to_string());

    lines.push("  </graph>".to_string());
    lines.push("</gexf>".to_string());

    lines.join("\n")
}

// ──────────────────────────── Helpers ────────────────────────────

/// Extract an XML attribute value.
fn extract_attr(line: &str, attr: &str) -> Option<String> {
    let pattern = format!("{}=\"", attr);
    if let Some(start) = line.find(&pattern) {
        let value_start = start + pattern.len();
        if let Some(end) = line[value_start..].find('"') {
            return Some(line[value_start..value_start + end].to_string());
        }
    }
    None
}

/// Extract text content between XML tags.
fn extract_text_content(line: &str) -> String {
    if let Some(start) = line.find('>') {
        if let Some(end) = line[start + 1..].find('<') {
            return line[start + 1..start + 1 + end].to_string();
        }
    }
    String::new()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_sif() {
        let sif = "TP53\tpp\tMDM2\nTP53\tpp\tBAX\nMDM2\tpp\tMDMX\n";
        let g = parse_sif(sif).unwrap();
        assert_eq!(g.node_count(), 4);
        assert_eq!(g.edge_count(), 3);
        assert!(g.has_edge("TP53", "MDM2"));
    }

    #[test]
    fn test_parse_sif_multiple_targets() {
        let sif = "TP53\tpp\tMDM2\tBAX\tBCL2\n";
        let g = parse_sif(sif).unwrap();
        assert_eq!(g.node_count(), 4);
        assert_eq!(g.edge_count(), 3);
    }

    #[test]
    fn test_parse_sif_isolated() {
        let sif = "ORPHAN\nTP53\tpp\tMDM2\n";
        let g = parse_sif(sif).unwrap();
        assert_eq!(g.node_count(), 3);
        assert!(g.get_node("ORPHAN").is_some());
    }

    #[test]
    fn test_parse_sif_comments() {
        let sif = "# comment\nA\tpp\tB\n";
        let g = parse_sif(sif).unwrap();
        assert_eq!(g.node_count(), 2);
    }

    #[test]
    fn test_write_sif() {
        let mut g = Graph::new(GraphType::Undirected);
        g.add_node("A", "A").unwrap();
        g.add_node("B", "B").unwrap();
        let edge = Edge::new("A", "B").with_type("pp");
        g.add_edge_struct(edge).unwrap();
        let sif = write_sif(&g);
        assert!(sif.contains("A\tpp\tB"));
    }

    #[test]
    fn test_sif_roundtrip() {
        let sif = "A\tpp\tB\nB\tpp\tC\n";
        let g = parse_sif(sif).unwrap();
        let output = write_sif(&g);
        let g2 = parse_sif(&output).unwrap();
        assert_eq!(g.node_count(), g2.node_count());
        assert_eq!(g.edge_count(), g2.edge_count());
    }

    #[test]
    fn test_parse_graphml() {
        let graphml = r#"<?xml version="1.0" encoding="UTF-8"?>
<graphml>
  <graph edgedefault="undirected">
    <node id="A"/>
    <node id="B"/>
    <node id="C"/>
    <edge source="A" target="B"/>
    <edge source="B" target="C"/>
  </graph>
</graphml>"#;
        let g = parse_graphml(graphml).unwrap();
        assert_eq!(g.node_count(), 3);
        assert_eq!(g.edge_count(), 2);
        assert_eq!(g.graph_type, GraphType::Undirected);
    }

    #[test]
    fn test_parse_graphml_directed() {
        let graphml = r#"<graphml>
  <graph edgedefault="directed">
    <node id="TF1"/>
    <node id="G1"/>
    <edge source="TF1" target="G1"/>
  </graph>
</graphml>"#;
        let g = parse_graphml(graphml).unwrap();
        assert_eq!(g.graph_type, GraphType::Directed);
        assert!(g.has_edge("TF1", "G1"));
        assert!(!g.has_edge("G1", "TF1"));
    }

    #[test]
    fn test_write_graphml() {
        let mut g = Graph::new(GraphType::Undirected);
        g.add_node("A", "Gene A").unwrap();
        g.add_node("B", "Gene B").unwrap();
        g.add_edge("A", "B", 0.95).unwrap();
        let xml = write_graphml(&g);
        assert!(xml.contains("edgedefault=\"undirected\""));
        assert!(xml.contains("node id=\"A\""));
        assert!(xml.contains("source=\"A\""));
        assert!(xml.contains("0.95"));
    }

    #[test]
    fn test_graphml_roundtrip() {
        let mut g = Graph::new(GraphType::Undirected);
        g.add_node("X", "X").unwrap();
        g.add_node("Y", "Y").unwrap();
        g.add_edge("X", "Y", 1.0).unwrap();
        let xml = write_graphml(&g);
        let g2 = parse_graphml(&xml).unwrap();
        assert_eq!(g.node_count(), g2.node_count());
        assert_eq!(g.edge_count(), g2.edge_count());
    }

    #[test]
    fn test_write_gexf() {
        let mut g = Graph::new(GraphType::Directed);
        g.add_node("TF1", "TF1").unwrap();
        g.add_node("G1", "G1").unwrap();
        g.add_edge("TF1", "G1", 0.9).unwrap();
        let gexf = write_gexf(&g);
        assert!(gexf.contains("defaultedgetype=\"directed\""));
        assert!(gexf.contains("node id=\"G1\""));
        assert!(gexf.contains("source=\"TF1\""));
    }

    #[test]
    fn test_extract_attr() {
        assert_eq!(
            extract_attr(r#"<node id="test" label="foo">"#, "id"),
            Some("test".into())
        );
        assert_eq!(extract_attr("<node>", "id"), None);
    }
}
