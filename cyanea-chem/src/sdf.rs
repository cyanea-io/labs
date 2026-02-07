//! SDF and MOL V2000 file format parser.

use cyanea_core::{CyaneaError, Result};

use crate::element::element_by_symbol;
use crate::molecule::{Bond, BondOrder, MolAtom, Molecule};

/// Parse a MOL V2000 block into a `Molecule`.
pub fn parse_mol_v2000(input: &str) -> Result<Molecule> {
    let lines: Vec<&str> = input.lines().collect();

    if lines.len() < 4 {
        return Err(CyaneaError::Parse("MOL block too short".into()));
    }

    // Header: line 0 = molecule name, 1 = program/timestamp, 2 = comment
    let name = lines[0].trim().to_string();

    // Counts line (line 3): aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
    let counts_line = lines[3];
    if counts_line.len() < 6 {
        return Err(CyaneaError::Parse("counts line too short".into()));
    }

    let num_atoms: usize = counts_line[0..3]
        .trim()
        .parse()
        .map_err(|_| CyaneaError::Parse("invalid atom count".into()))?;
    let num_bonds: usize = counts_line[3..6]
        .trim()
        .parse()
        .map_err(|_| CyaneaError::Parse("invalid bond count".into()))?;

    let atom_start = 4;
    let bond_start = atom_start + num_atoms;

    if lines.len() < bond_start + num_bonds {
        return Err(CyaneaError::Parse("MOL block truncated".into()));
    }

    // Parse atom block
    let mut atoms = Vec::with_capacity(num_atoms);
    for i in 0..num_atoms {
        let line = lines[atom_start + i];
        let atom = parse_atom_line(line)?;
        atoms.push(atom);
    }

    // Parse bond block
    let mut bonds = Vec::with_capacity(num_bonds);
    for i in 0..num_bonds {
        let line = lines[bond_start + i];
        let bond = parse_bond_line(line)?;
        bonds.push(bond);
    }

    // Parse properties block for charges (M  CHG)
    for i in (bond_start + num_bonds)..lines.len() {
        let line = lines[i];
        if line.starts_with("M  END") {
            break;
        }
        if line.starts_with("M  CHG") {
            parse_charge_line(line, &mut atoms)?;
        }
    }

    Ok(Molecule::new(name, atoms, bonds))
}

/// Parse a multi-molecule SDF string, returning results for each molecule.
pub fn parse_sdf(input: &str) -> Vec<Result<Molecule>> {
    let mut results = Vec::new();

    // Split on $$$$ delimiter
    for block in input.split("$$$$") {
        let trimmed = block.trim();
        if trimmed.is_empty() {
            continue;
        }
        results.push(parse_mol_v2000(trimmed));
    }

    results
}

/// Parse an SDF file from disk.
#[cfg(feature = "std")]
pub fn parse_sdf_file(path: impl AsRef<std::path::Path>) -> Result<Vec<Molecule>> {
    let content = std::fs::read_to_string(path.as_ref()).map_err(|e| {
        CyaneaError::Io(std::io::Error::new(
            e.kind(),
            format!("{}: {}", path.as_ref().display(), e),
        ))
    })?;
    parse_sdf(&content)
        .into_iter()
        .collect::<Result<Vec<_>>>()
}

fn parse_atom_line(line: &str) -> Result<MolAtom> {
    // V2000 atom line: xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
    // Columns: x(0..10) y(10..20) z(20..30) _ symbol(31..34) ...
    if line.len() < 34 {
        return Err(CyaneaError::Parse(format!(
            "atom line too short: '{line}'"
        )));
    }

    let symbol = line[31..34].trim();
    let elem = element_by_symbol(symbol).ok_or_else(|| {
        CyaneaError::Parse(format!("unknown element '{symbol}' in MOL atom block"))
    })?;

    // Charge from atom line (column 36-38), old-style: 0=none, 1=+3, 2=+2, 3=+1, 4=doublet, 5=-1, 6=-2, 7=-3
    let charge = if line.len() >= 39 {
        match line[36..39].trim().parse::<u8>() {
            Ok(0) | Err(_) => 0i8,
            Ok(1) => 3,
            Ok(2) => 2,
            Ok(3) => 1,
            Ok(5) => -1,
            Ok(6) => -2,
            Ok(7) => -3,
            _ => 0,
        }
    } else {
        0
    };

    Ok(MolAtom {
        atomic_number: elem.atomic_number,
        formal_charge: charge,
        isotope: None,
        is_aromatic: false,
        implicit_hydrogens: 0, // V2000 doesn't store implicit H — computed from valence if needed
    })
}

fn parse_bond_line(line: &str) -> Result<Bond> {
    // V2000 bond line: 111222tttsssxxxrrrccc
    // Columns: atom1(0..3) atom2(3..6) type(6..9) ...
    if line.len() < 9 {
        return Err(CyaneaError::Parse(format!(
            "bond line too short: '{line}'"
        )));
    }

    let a1: usize = line[0..3]
        .trim()
        .parse()
        .map_err(|_| CyaneaError::Parse("invalid bond atom1".into()))?;
    let a2: usize = line[3..6]
        .trim()
        .parse()
        .map_err(|_| CyaneaError::Parse("invalid bond atom2".into()))?;
    let bond_type: u8 = line[6..9]
        .trim()
        .parse()
        .map_err(|_| CyaneaError::Parse("invalid bond type".into()))?;

    let order = match bond_type {
        1 => BondOrder::Single,
        2 => BondOrder::Double,
        3 => BondOrder::Triple,
        4 => BondOrder::Aromatic,
        _ => BondOrder::Single,
    };

    // V2000 atom indices are 1-based
    Ok(Bond {
        atom1: a1 - 1,
        atom2: a2 - 1,
        order,
        is_aromatic: bond_type == 4,
    })
}

fn parse_charge_line(line: &str, atoms: &mut [MolAtom]) -> Result<()> {
    // M  CHG  n  aaa vvv  aaa vvv ...
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 4 {
        return Ok(());
    }
    let count: usize = parts[2]
        .parse()
        .map_err(|_| CyaneaError::Parse("invalid charge count".into()))?;

    for i in 0..count {
        let idx_pos = 3 + i * 2;
        let val_pos = 4 + i * 2;
        if val_pos >= parts.len() {
            break;
        }
        let atom_idx: usize = parts[idx_pos]
            .parse::<usize>()
            .map_err(|_| CyaneaError::Parse("invalid charge atom index".into()))?
            - 1;
        let charge: i8 = parts[val_pos]
            .parse()
            .map_err(|_| CyaneaError::Parse("invalid charge value".into()))?;
        if atom_idx < atoms.len() {
            atoms[atom_idx].formal_charge = charge;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use cyanea_core::Annotated;

    fn minimal_mol() -> &'static str {
        "\
Methane
     RDKit          3D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END"
    }

    #[test]
    fn parse_minimal_mol() {
        let mol = parse_mol_v2000(minimal_mol()).unwrap();
        assert_eq!(mol.name(), "Methane");
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(mol.bond_count(), 0);
        assert_eq!(mol.atoms[0].atomic_number, 6);
    }

    #[test]
    fn parse_mol_with_charges() {
        let mol_str = "\
Charged
     RDKit          3D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  CHG  2   1   1   2  -1
M  END";
        let mol = parse_mol_v2000(mol_str).unwrap();
        assert_eq!(mol.atoms[0].formal_charge, 1);
        assert_eq!(mol.atoms[1].formal_charge, -1);
    }

    #[test]
    fn parse_multi_molecule_sdf() {
        let sdf = format!("{}\n$$$$\n{}\n$$$$\n", minimal_mol(), minimal_mol());
        let results = parse_sdf(&sdf);
        assert_eq!(results.len(), 2);
        assert!(results[0].is_ok());
        assert!(results[1].is_ok());
    }

    #[test]
    fn malformed_mol_error() {
        assert!(parse_mol_v2000("too\nshort").is_err());
        assert!(parse_mol_v2000("name\nprog\ncomment\nabc  0").is_err());
    }
}
