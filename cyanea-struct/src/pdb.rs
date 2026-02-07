//! PDB format parser.
//!
//! Parses ATOM, HETATM, TER, HEADER, and MODEL/ENDMDL records from PDB-format
//! text. Only the first MODEL is returned for multi-model (NMR) files.

use cyanea_core::{CyaneaError, Result};

use crate::types::{Atom, Chain, Point3D, Residue, Structure};

use alloc::string::{String, ToString};
use alloc::vec::Vec;

/// Parse a PDB-format string into a [`Structure`].
///
/// # Errors
///
/// Returns an error if no ATOM records are found or if an ATOM record is
/// malformed (wrong column widths, unparseable coordinates).
pub fn parse_pdb(input: &str) -> Result<Structure> {
    let mut id = String::from("UNKN");
    let mut current_chain_id: Option<char> = None;
    let mut current_residue_key: Option<(i32, Option<char>, String)> = None;
    let mut chains: Vec<Chain> = Vec::new();
    let mut residues: Vec<Residue> = Vec::new();
    let mut atoms: Vec<Atom> = Vec::new();
    let mut atom_count = 0u32;
    let mut in_first_model = true;

    for line in input.lines() {
        if line.starts_with("ENDMDL") {
            break; // only first model
        }
        if line.starts_with("MODEL") {
            if !in_first_model {
                break;
            }
            in_first_model = false;
            continue;
        }

        if line.starts_with("HEADER") && line.len() >= 66 {
            let pdb_id = line[62..66].trim();
            if !pdb_id.is_empty() {
                id = pdb_id.into();
            }
        }

        let is_atom = line.starts_with("ATOM  ");
        let is_hetatm = line.starts_with("HETATM");

        if is_atom || is_hetatm {
            let atom = parse_atom_record(line, is_hetatm)?;
            atom_count += 1;

            let chain_id = parse_chain_id(line);
            let seq_num = parse_residue_seq(line)?;
            let i_code = parse_insertion_code(line);
            let res_name = parse_residue_name(line);

            let res_key = (seq_num, i_code, res_name.clone());

            // Chain changed?
            if current_chain_id != Some(chain_id) {
                // Flush current residue
                if !atoms.is_empty() {
                    if let Some((sn, ic, rn)) = current_residue_key.take() {
                        residues.push(Residue {
                            name: rn,
                            seq_num: sn,
                            i_code: ic,
                            atoms: core::mem::take(&mut atoms),
                        });
                    }
                }
                // Flush current chain
                if let Some(cid) = current_chain_id {
                    if !residues.is_empty() {
                        chains.push(Chain::new(cid, core::mem::take(&mut residues)));
                    }
                }
                current_chain_id = Some(chain_id);
                current_residue_key = Some(res_key);
            } else if current_residue_key.as_ref() != Some(&res_key) {
                // Residue changed within same chain
                if !atoms.is_empty() {
                    if let Some((sn, ic, rn)) = current_residue_key.take() {
                        residues.push(Residue {
                            name: rn,
                            seq_num: sn,
                            i_code: ic,
                            atoms: core::mem::take(&mut atoms),
                        });
                    }
                }
                current_residue_key = Some(res_key);
            }

            atoms.push(atom);
        }

        if line.starts_with("TER") {
            // Flush residue + chain
            if !atoms.is_empty() {
                if let Some((sn, ic, rn)) = current_residue_key.take() {
                    residues.push(Residue {
                        name: rn,
                        seq_num: sn,
                        i_code: ic,
                        atoms: core::mem::take(&mut atoms),
                    });
                }
            }
            if let Some(cid) = current_chain_id.take() {
                if !residues.is_empty() {
                    chains.push(Chain::new(cid, core::mem::take(&mut residues)));
                }
            }
        }
    }

    // Flush remaining
    if !atoms.is_empty() {
        if let Some((sn, ic, rn)) = current_residue_key.take() {
            residues.push(Residue {
                name: rn,
                seq_num: sn,
                i_code: ic,
                atoms,
            });
        }
    }
    if let Some(cid) = current_chain_id {
        if !residues.is_empty() {
            chains.push(Chain::new(cid, residues));
        }
    }

    if atom_count == 0 {
        return Err(CyaneaError::Parse("no ATOM records found".into()));
    }

    Ok(Structure { id, chains })
}

/// Parse a PDB file from disk.
#[cfg(feature = "std")]
pub fn parse_pdb_file(path: impl AsRef<::std::path::Path>) -> Result<Structure> {
    let contents = ::std::fs::read_to_string(path)?;
    parse_pdb(&contents)
}

fn parse_atom_record(line: &str, is_hetatm: bool) -> Result<Atom> {
    // PDB format is fixed-width columns. We need at least 54 chars for coords.
    if line.len() < 54 {
        return Err(CyaneaError::Parse(alloc::format!(
            "ATOM record too short ({} chars): {}",
            line.len(),
            line
        )));
    }

    let serial = safe_slice(line, 6, 11)
        .trim()
        .parse::<u32>()
        .map_err(|e| CyaneaError::Parse(alloc::format!("bad atom serial: {}", e)))?;

    let name = safe_slice(line, 12, 16).to_string();

    let alt_loc = {
        let c = safe_slice(line, 16, 17).chars().next().unwrap_or(' ');
        if c == ' ' {
            None
        } else {
            Some(c)
        }
    };

    let x = safe_slice(line, 30, 38)
        .trim()
        .parse::<f64>()
        .map_err(|e| CyaneaError::Parse(alloc::format!("bad x coordinate: {}", e)))?;
    let y = safe_slice(line, 38, 46)
        .trim()
        .parse::<f64>()
        .map_err(|e| CyaneaError::Parse(alloc::format!("bad y coordinate: {}", e)))?;
    let z = safe_slice(line, 46, 54)
        .trim()
        .parse::<f64>()
        .map_err(|e| CyaneaError::Parse(alloc::format!("bad z coordinate: {}", e)))?;

    let occupancy = if line.len() >= 60 {
        safe_slice(line, 54, 60).trim().parse::<f64>().unwrap_or(1.0)
    } else {
        1.0
    };

    let temp_factor = if line.len() >= 66 {
        safe_slice(line, 60, 66).trim().parse::<f64>().unwrap_or(0.0)
    } else {
        0.0
    };

    let element = if line.len() >= 78 {
        let e = safe_slice(line, 76, 78).trim();
        if e.is_empty() {
            None
        } else {
            Some(e.to_string())
        }
    } else {
        None
    };

    let charge = if line.len() >= 80 {
        let ch = safe_slice(line, 78, 80).trim();
        parse_pdb_charge(ch)
    } else {
        None
    };

    Ok(Atom {
        serial,
        name,
        alt_loc,
        coords: Point3D::new(x, y, z),
        occupancy,
        temp_factor,
        element,
        charge,
        is_hetatm,
    })
}

fn parse_chain_id(line: &str) -> char {
    safe_slice(line, 21, 22).chars().next().unwrap_or(' ')
}

fn parse_residue_seq(line: &str) -> Result<i32> {
    safe_slice(line, 22, 26)
        .trim()
        .parse::<i32>()
        .map_err(|e| CyaneaError::Parse(alloc::format!("bad residue seq number: {}", e)))
}

fn parse_insertion_code(line: &str) -> Option<char> {
    let c = safe_slice(line, 26, 27).chars().next().unwrap_or(' ');
    if c == ' ' {
        None
    } else {
        Some(c)
    }
}

fn parse_residue_name(line: &str) -> String {
    safe_slice(line, 17, 20).trim().to_string()
}

fn parse_pdb_charge(s: &str) -> Option<i8> {
    if s.is_empty() {
        return None;
    }
    // PDB charge format: "2+" or "1-" or "+2" etc.
    let s = s.trim();
    if s.is_empty() {
        return None;
    }
    // Try "2+" / "2-" format
    if s.len() == 2 {
        let chars: Vec<char> = s.chars().collect();
        if chars[1] == '+' || chars[1] == '-' {
            if let Some(d) = chars[0].to_digit(10) {
                let sign: i8 = if chars[1] == '+' { 1 } else { -1 };
                return Some(sign * d as i8);
            }
        }
        if chars[0] == '+' || chars[0] == '-' {
            if let Some(d) = chars[1].to_digit(10) {
                let sign: i8 = if chars[0] == '+' { 1 } else { -1 };
                return Some(sign * d as i8);
            }
        }
    }
    None
}

/// Safe substring that handles short lines gracefully.
fn safe_slice(s: &str, start: usize, end: usize) -> &str {
    let bytes = s.as_bytes();
    let len = bytes.len();
    if start >= len {
        return "";
    }
    let actual_end = end.min(len);
    // Safety: PDB files are ASCII, so byte boundaries = char boundaries
    &s[start..actual_end]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn minimal_pdb() -> &'static str {
        "\
HEADER                                                        1CRN\n\
ATOM      1  N   THR A   1       2.464   9.901  13.546  1.00 10.00           N\n\
ATOM      2  CA  THR A   1       2.135  10.226  12.120  1.00 10.00           C\n\
ATOM      3  C   THR A   1       3.427  10.018  11.354  1.00 10.00           C\n\
ATOM      4  O   THR A   1       3.426  10.335  10.184  1.00 10.00           O\n\
ATOM      5  N   ILE A   2       4.462   9.470  11.952  1.00 10.00           N\n\
ATOM      6  CA  ILE A   2       5.735   9.197  11.275  1.00 10.00           C\n\
TER       7      ILE A   2\n\
END\n"
    }

    #[test]
    fn parse_single_chain() {
        let s = parse_pdb(minimal_pdb()).unwrap();
        assert_eq!(s.id, "1CRN");
        assert_eq!(s.chain_count(), 1);
        assert_eq!(s.get_chain('A').unwrap().residue_count(), 2);
        assert_eq!(s.atom_count(), 6);
    }

    #[test]
    fn parse_multi_chain() {
        let input = "\
ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C\n\
TER\n\
ATOM      2  CA  ALA B   1       4.000   5.000   6.000  1.00  0.00           C\n\
TER\n\
END\n";
        let s = parse_pdb(input).unwrap();
        assert_eq!(s.chain_count(), 2);
        assert_eq!(s.get_chain('A').unwrap().atom_count(), 1);
        assert_eq!(s.get_chain('B').unwrap().atom_count(), 1);
    }

    #[test]
    fn parse_hetatm() {
        let input = "\
ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C\n\
HETATM    2  O   HOH A   2       4.000   5.000   6.000  1.00  0.00           O\n\
END\n";
        let s = parse_pdb(input).unwrap();
        let chain = s.get_chain('A').unwrap();
        assert_eq!(chain.residue_count(), 2);
        let water = &chain.residues[1];
        assert_eq!(water.name, "HOH");
        assert!(water.atoms[0].is_hetatm);
    }

    #[test]
    fn parse_insertion_codes() {
        // Use properly column-formatted PDB lines.
        // Col  1- 6: record type
        // Col  7-11: serial
        // Col 13-16: atom name
        // Col 18-20: res name
        // Col 22:    chain ID
        // Col 23-26: seq num (right-justified)
        // Col 27:    iCode
        // Col 31-38: x  (8.3f)
        // Col 39-46: y  (8.3f)
        // Col 47-54: z  (8.3f)
        //       1         2         3         4         5         6         7
        // 01234567890123456789012345678901234567890123456789012345678901234567890123456789
        let input = "\
ATOM      1  CA  ALA A  10       1.000   2.000   3.000  1.00  0.00           C\n\
ATOM      2  CA  ALA A  10A      4.000   5.000   6.000  1.00  0.00           C\n\
END\n";
        let s = parse_pdb(input).unwrap();
        let chain = s.get_chain('A').unwrap();
        assert_eq!(chain.residue_count(), 2);
        assert_eq!(chain.residues[0].i_code, None);
        assert_eq!(chain.residues[1].i_code, Some('A'));
    }

    #[test]
    fn parse_malformed_record() {
        let input = "ATOM   BAD\n";
        assert!(parse_pdb(input).is_err());
    }
}
