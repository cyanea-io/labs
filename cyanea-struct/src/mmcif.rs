//! mmCIF/PDBx format parser.
//!
//! Parses the mmCIF text format (key-value pairs and `loop_` tables) to extract
//! macromolecular structure data. Only `_atom_site.` loop records are currently
//! used to build a [`Structure`].
//!
//! # Example
//!
//! ```
//! use cyanea_struct::mmcif::parse_mmcif;
//!
//! let mmcif_text = "\
//! data_1TST
//! #
//! loop_
//! _atom_site.group_PDB
//! _atom_site.id
//! _atom_site.type_symbol
//! _atom_site.label_atom_id
//! _atom_site.label_comp_id
//! _atom_site.label_asym_id
//! _atom_site.label_seq_id
//! _atom_site.Cartn_x
//! _atom_site.Cartn_y
//! _atom_site.Cartn_z
//! _atom_site.occupancy
//! _atom_site.B_iso_or_equiv
//! ATOM 1 N N ALA A 1 1.000 2.000 3.000 1.00 10.00
//! ATOM 2 C CA ALA A 1 2.000 2.000 3.000 1.00 12.00
//! #
//! ";
//!
//! let structure = parse_mmcif(mmcif_text).unwrap();
//! assert_eq!(structure.id, "1TST");
//! assert_eq!(structure.atom_count(), 2);
//! ```

use cyanea_core::{CyaneaError, Result};

use crate::types::{Atom, Chain, Point3D, Residue, Structure};

use alloc::collections::BTreeMap;
use alloc::string::{String, ToString};
use alloc::vec::Vec;

/// Parse an mmCIF-format string into a [`Structure`].
///
/// Extracts atom coordinates from the `_atom_site.` loop table and groups them
/// into chains and residues. The structure ID is taken from the `data_` block
/// header.
///
/// # Errors
///
/// Returns an error if no `_atom_site.` loop is found, if required fields are
/// missing, or if coordinate values cannot be parsed.
pub fn parse_mmcif(input: &str) -> Result<Structure> {
    let lines: Vec<&str> = input.lines().collect();

    // Extract structure ID from data_ block
    let id = extract_data_id(&lines);

    // Find and parse _atom_site. loop
    let atom_records = find_and_parse_atom_site_loop(&lines)?;

    if atom_records.is_empty() {
        return Err(CyaneaError::Parse(
            "no _atom_site. records found in mmCIF".into(),
        ));
    }

    // Build Structure from atom records
    build_structure_from_records(&id, &atom_records)
}

/// Parse a `loop_` construct into a list of records (field-name to value maps).
///
/// The input `lines` should start right after the `loop_` keyword line, with
/// field names (`_category.field`) followed by data rows. Parsing stops at
/// a line starting with `#`, `loop_`, `data_`, or at a new category prefix.
pub fn parse_mmcif_loop(lines: &[&str]) -> Vec<BTreeMap<String, String>> {
    let mut field_names: Vec<String> = Vec::new();
    let mut records: Vec<BTreeMap<String, String>> = Vec::new();
    let mut reading_fields = true;

    for line in lines {
        let trimmed = line.trim();

        // Skip empty lines and comments
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with('#') {
            if !field_names.is_empty() && !reading_fields {
                break; // End of data block
            }
            continue;
        }

        // Stop at new structural keywords
        if trimmed.starts_with("loop_") || trimmed.starts_with("data_") {
            break;
        }

        if reading_fields {
            if trimmed.starts_with('_') {
                field_names.push(trimmed.to_string());
                continue;
            } else {
                // Transition to data rows
                reading_fields = false;
            }
        }

        // Parse data row
        if !reading_fields {
            let tokens = tokenize_mmcif_line(trimmed);
            if tokens.len() < field_names.len() {
                // Possibly a continuation or malformed; skip
                continue;
            }
            let mut record = BTreeMap::new();
            for (i, name) in field_names.iter().enumerate() {
                record.insert(name.clone(), tokens[i].clone());
            }
            records.push(record);
        }
    }

    records
}

/// Extract the structure ID from a `data_XXXX` line.
fn extract_data_id(lines: &[&str]) -> String {
    for line in lines {
        let trimmed = line.trim();
        if trimmed.starts_with("data_") {
            return trimmed[5..].to_string();
        }
    }
    String::from("UNKN")
}

/// Locate the `_atom_site.` loop and parse it into records.
fn find_and_parse_atom_site_loop(
    lines: &[&str],
) -> Result<Vec<BTreeMap<String, String>>> {
    let mut i = 0;
    while i < lines.len() {
        let trimmed = lines[i].trim();

        // Look for loop_ followed by _atom_site. fields
        if trimmed == "loop_" {
            // Check if next non-empty/non-comment line starts with _atom_site.
            let mut j = i + 1;
            while j < lines.len() {
                let next = lines[j].trim();
                if next.is_empty() || next.starts_with('#') {
                    j += 1;
                    continue;
                }
                if next.starts_with("_atom_site.") {
                    // Found it — parse from j onward
                    return Ok(parse_mmcif_loop(&lines[j..]));
                }
                break;
            }
        }

        i += 1;
    }

    Err(CyaneaError::Parse(
        "no _atom_site. loop found in mmCIF data".into(),
    ))
}

/// Tokenize a single mmCIF data line, respecting single- and double-quoted strings.
fn tokenize_mmcif_line(line: &str) -> Vec<String> {
    let mut tokens = Vec::new();
    let chars: Vec<char> = line.chars().collect();
    let len = chars.len();
    let mut i = 0;

    while i < len {
        // Skip whitespace
        if chars[i].is_whitespace() {
            i += 1;
            continue;
        }

        // Quoted string (single or double)
        if chars[i] == '\'' || chars[i] == '"' {
            let quote = chars[i];
            i += 1;
            let start = i;
            while i < len && chars[i] != quote {
                i += 1;
            }
            let token: String = chars[start..i].iter().collect();
            tokens.push(token);
            if i < len {
                i += 1; // skip closing quote
            }
            continue;
        }

        // Unquoted token
        let start = i;
        while i < len && !chars[i].is_whitespace() {
            i += 1;
        }
        let token: String = chars[start..i].iter().collect();
        tokens.push(token);
    }

    tokens
}

/// Build a Structure from parsed _atom_site. records.
fn build_structure_from_records(
    id: &str,
    records: &[BTreeMap<String, String>],
) -> Result<Structure> {
    // Group atoms by chain, then by residue
    // We use a Vec of (chain_id, Vec<(res_key, Vec<Atom>)>) to maintain order
    let mut chain_order: Vec<char> = Vec::new();
    let mut chain_residues: BTreeMap<char, Vec<(i32, String, Vec<Atom>)>> = BTreeMap::new();

    // Track current residue per chain to group atoms
    let mut current_res: BTreeMap<char, (i32, String)> = BTreeMap::new();

    for (idx, record) in records.iter().enumerate() {
        // Skip HETATM if needed, but we'll include them like the PDB parser does
        let group = record
            .get("_atom_site.group_PDB")
            .map(|s| s.as_str())
            .unwrap_or("ATOM");
        let is_hetatm = group == "HETATM";

        let serial = get_field_u32(record, "_atom_site.id", idx)?;
        let atom_name = get_field_str(record, "_atom_site.label_atom_id", idx)?;
        let comp_id = get_field_str(record, "_atom_site.label_comp_id", idx)?;
        let asym_id = get_field_str(record, "_atom_site.label_asym_id", idx)?;
        let seq_id = get_field_i32(record, "_atom_site.label_seq_id", idx)?;

        let x = get_field_f64(record, "_atom_site.Cartn_x", idx)?;
        let y = get_field_f64(record, "_atom_site.Cartn_y", idx)?;
        let z = get_field_f64(record, "_atom_site.Cartn_z", idx)?;

        let occupancy = get_field_f64_opt(record, "_atom_site.occupancy").unwrap_or(1.0);
        let b_factor = get_field_f64_opt(record, "_atom_site.B_iso_or_equiv").unwrap_or(0.0);

        let element = record
            .get("_atom_site.type_symbol")
            .filter(|s| s.as_str() != "." && s.as_str() != "?")
            .map(|s| s.clone());

        // Use first character of asym_id as chain ID
        let chain_id = asym_id.chars().next().unwrap_or('A');

        let atom = Atom {
            serial,
            name: atom_name.clone(),
            alt_loc: record
                .get("_atom_site.label_alt_id")
                .and_then(|s| {
                    let c = s.chars().next()?;
                    if c == '.' || c == '?' {
                        None
                    } else {
                        Some(c)
                    }
                }),
            coords: Point3D::new(x, y, z),
            occupancy,
            temp_factor: b_factor,
            element,
            charge: None,
            is_hetatm,
        };

        // Check if we need a new chain or a new residue
        if !chain_residues.contains_key(&chain_id) {
            chain_order.push(chain_id);
            chain_residues.insert(chain_id, Vec::new());
        }

        let residues = chain_residues.get_mut(&chain_id).unwrap();

        let need_new_residue = match current_res.get(&chain_id) {
            Some((cur_seq, cur_name)) => *cur_seq != seq_id || *cur_name != comp_id,
            None => true,
        };

        if need_new_residue {
            residues.push((seq_id, comp_id.clone(), Vec::new()));
            current_res.insert(chain_id, (seq_id, comp_id));
        }

        // Add atom to the last residue of this chain
        if let Some(last) = residues.last_mut() {
            last.2.push(atom);
        }
    }

    // Build chains in order
    let mut chains = Vec::new();
    for chain_id in &chain_order {
        if let Some(residue_data) = chain_residues.remove(chain_id) {
            let residues: Vec<Residue> = residue_data
                .into_iter()
                .map(|(seq_num, name, atoms)| Residue {
                    name,
                    seq_num,
                    i_code: None,
                    atoms,
                })
                .collect();
            if !residues.is_empty() {
                chains.push(Chain::new(*chain_id, residues));
            }
        }
    }

    Ok(Structure {
        id: id.to_string(),
        chains,
    })
}

// ---- Field extraction helpers ----

fn get_field_str(
    record: &BTreeMap<String, String>,
    field: &str,
    row: usize,
) -> Result<String> {
    record
        .get(field)
        .filter(|s| s.as_str() != "." && s.as_str() != "?")
        .map(|s| s.clone())
        .ok_or_else(|| {
            CyaneaError::Parse(alloc::format!(
                "missing field {} in _atom_site row {}",
                field,
                row
            ))
        })
}

fn get_field_u32(
    record: &BTreeMap<String, String>,
    field: &str,
    row: usize,
) -> Result<u32> {
    let s = get_field_str(record, field, row)?;
    s.parse::<u32>().map_err(|e| {
        CyaneaError::Parse(alloc::format!(
            "bad {} value '{}' in row {}: {}",
            field,
            s,
            row,
            e
        ))
    })
}

fn get_field_i32(
    record: &BTreeMap<String, String>,
    field: &str,
    row: usize,
) -> Result<i32> {
    let s = get_field_str(record, field, row)?;
    s.parse::<i32>().map_err(|e| {
        CyaneaError::Parse(alloc::format!(
            "bad {} value '{}' in row {}: {}",
            field,
            s,
            row,
            e
        ))
    })
}

fn get_field_f64(
    record: &BTreeMap<String, String>,
    field: &str,
    row: usize,
) -> Result<f64> {
    let s = get_field_str(record, field, row)?;
    s.parse::<f64>().map_err(|e| {
        CyaneaError::Parse(alloc::format!(
            "bad {} value '{}' in row {}: {}",
            field,
            s,
            row,
            e
        ))
    })
}

fn get_field_f64_opt(
    record: &BTreeMap<String, String>,
    field: &str,
) -> Option<f64> {
    record
        .get(field)
        .filter(|s| s.as_str() != "." && s.as_str() != "?")
        .and_then(|s| s.parse::<f64>().ok())
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    fn minimal_mmcif() -> &'static str {
        "\
data_1TST
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 N N ALA A 1 1.000 2.000 3.000 1.00 10.00
ATOM 2 C CA ALA A 1 2.000 2.000 3.000 1.00 12.00
ATOM 3 C C ALA A 1 3.000 2.000 3.000 1.00 11.00
ATOM 4 O O ALA A 1 3.000 3.000 3.000 1.00 13.00
ATOM 5 N N GLY A 2 4.000 2.000 3.000 1.00 9.00
ATOM 6 C CA GLY A 2 5.000 2.000 3.000 1.00 8.00
ATOM 7 C C GLY A 2 6.000 2.000 3.000 1.00 10.00
ATOM 8 O O GLY A 2 6.000 3.000 3.000 1.00 11.00
#
"
    }

    #[test]
    fn parse_minimal_mmcif() {
        let s = parse_mmcif(minimal_mmcif()).unwrap();
        assert_eq!(s.id, "1TST");
        assert_eq!(s.chain_count(), 1);
        assert_eq!(s.residue_count(), 2);
        assert_eq!(s.atom_count(), 8);
    }

    #[test]
    fn parse_coordinates() {
        let s = parse_mmcif(minimal_mmcif()).unwrap();
        let chain = s.get_chain('A').unwrap();
        let first_atom = &chain.residues[0].atoms[0];
        assert!((first_atom.coords.x - 1.0).abs() < 1e-10);
        assert!((first_atom.coords.y - 2.0).abs() < 1e-10);
        assert!((first_atom.coords.z - 3.0).abs() < 1e-10);
    }

    #[test]
    fn parse_b_factors() {
        let s = parse_mmcif(minimal_mmcif()).unwrap();
        let chain = s.get_chain('A').unwrap();
        let first_atom = &chain.residues[0].atoms[0];
        assert!((first_atom.temp_factor - 10.0).abs() < 1e-10);
        let ca_atom = &chain.residues[0].atoms[1];
        assert!((ca_atom.temp_factor - 12.0).abs() < 1e-10);
    }

    #[test]
    fn parse_residue_info() {
        let s = parse_mmcif(minimal_mmcif()).unwrap();
        let chain = s.get_chain('A').unwrap();
        assert_eq!(chain.residues[0].name, "ALA");
        assert_eq!(chain.residues[0].seq_num, 1);
        assert_eq!(chain.residues[1].name, "GLY");
        assert_eq!(chain.residues[1].seq_num, 2);
    }

    #[test]
    fn parse_atom_names() {
        let s = parse_mmcif(minimal_mmcif()).unwrap();
        let chain = s.get_chain('A').unwrap();
        let names: Vec<&str> = chain.residues[0]
            .atoms
            .iter()
            .map(|a| a.name.as_str())
            .collect();
        assert_eq!(names, vec!["N", "CA", "C", "O"]);
    }

    #[test]
    fn parse_multi_chain() {
        let input = "\
data_2CHN
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C CA ALA A 1 1.000 2.000 3.000 1.00 10.00
ATOM 2 C CA GLY A 2 4.000 5.000 6.000 1.00 11.00
ATOM 3 C CA VAL B 1 7.000 8.000 9.000 1.00 12.00
#
";
        let s = parse_mmcif(input).unwrap();
        assert_eq!(s.id, "2CHN");
        assert_eq!(s.chain_count(), 2);
        assert_eq!(s.get_chain('A').unwrap().residue_count(), 2);
        assert_eq!(s.get_chain('B').unwrap().residue_count(), 1);
    }

    #[test]
    fn parse_hetatm_records() {
        let input = "\
data_HET
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C CA ALA A 1 1.000 2.000 3.000 1.00 10.00
HETATM 2 O O HOH A 2 4.000 5.000 6.000 1.00 15.00
#
";
        let s = parse_mmcif(input).unwrap();
        let chain = s.get_chain('A').unwrap();
        assert!(!chain.residues[0].atoms[0].is_hetatm);
        assert!(chain.residues[1].atoms[0].is_hetatm);
    }

    #[test]
    fn missing_data_block() {
        let input = "\
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 N N ALA A 1 1.000 2.000 3.000 1.00 10.00
#
";
        let s = parse_mmcif(input).unwrap();
        assert_eq!(s.id, "UNKN");
        assert_eq!(s.atom_count(), 1);
    }

    #[test]
    fn no_atom_site_loop_error() {
        let input = "\
data_EMPTY
#
loop_
_entity.id
_entity.type
1 polymer
#
";
        assert!(parse_mmcif(input).is_err());
    }

    #[test]
    fn quoted_values() {
        let tokens = tokenize_mmcif_line("ATOM 1 C 'C A' ALA A 1 1.0 2.0 3.0 1.00 10.00");
        assert_eq!(tokens[3], "C A");
    }

    #[test]
    fn parse_loop_helper() {
        let lines = vec![
            "_test.field_a",
            "_test.field_b",
            "_test.field_c",
            "val1 val2 val3",
            "val4 val5 val6",
            "#",
        ];
        let records = parse_mmcif_loop(&lines);
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].get("_test.field_a").unwrap(), "val1");
        assert_eq!(records[0].get("_test.field_b").unwrap(), "val2");
        assert_eq!(records[1].get("_test.field_c").unwrap(), "val6");
    }

    #[test]
    fn element_symbol_extracted() {
        let s = parse_mmcif(minimal_mmcif()).unwrap();
        let chain = s.get_chain('A').unwrap();
        assert_eq!(
            chain.residues[0].atoms[0].element.as_deref(),
            Some("N")
        );
        assert_eq!(
            chain.residues[0].atoms[1].element.as_deref(),
            Some("C")
        );
    }
}
