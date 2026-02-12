//! SDF and MOL V2000/V3000 file format parser.

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

/// Parse a MOL V3000 (Enhanced) block into a `Molecule`.
///
/// The V3000 format uses `M  V30` line prefixes and free-format fields
/// instead of fixed-width columns. Atom and bond indices are 1-based.
pub fn parse_mol_v3000(input: &str) -> Result<Molecule> {
    let lines: Vec<&str> = input.lines().collect();

    if lines.len() < 4 {
        return Err(CyaneaError::Parse("V3000 MOL block too short".into()));
    }

    let name = lines[0].trim().to_string();

    // Find COUNTS line
    let counts_line = lines
        .iter()
        .find(|l| l.contains("M  V30 COUNTS"))
        .ok_or_else(|| CyaneaError::Parse("V3000: missing COUNTS line".into()))?;
    let counts_parts: Vec<&str> = counts_line.split_whitespace().collect();
    // Expected: M V30 COUNTS natoms nbonds ...
    if counts_parts.len() < 5 {
        return Err(CyaneaError::Parse("V3000: COUNTS line too short".into()));
    }
    let num_atoms: usize = counts_parts[3]
        .parse()
        .map_err(|_| CyaneaError::Parse("V3000: invalid atom count".into()))?;
    let num_bonds: usize = counts_parts[4]
        .parse()
        .map_err(|_| CyaneaError::Parse("V3000: invalid bond count".into()))?;

    // Find ATOM block
    let atom_begin = lines
        .iter()
        .position(|l| l.contains("M  V30 BEGIN ATOM"))
        .ok_or_else(|| CyaneaError::Parse("V3000: missing BEGIN ATOM".into()))?;
    let atom_end = lines
        .iter()
        .position(|l| l.contains("M  V30 END ATOM"))
        .ok_or_else(|| CyaneaError::Parse("V3000: missing END ATOM".into()))?;

    let mut atoms = Vec::with_capacity(num_atoms);
    for line in &lines[atom_begin + 1..atom_end] {
        let atom = parse_v3000_atom_line(line)?;
        atoms.push(atom);
    }

    if atoms.len() != num_atoms {
        return Err(CyaneaError::Parse(format!(
            "V3000: expected {num_atoms} atoms, found {}",
            atoms.len()
        )));
    }

    // Find BOND block (optional — molecules can have zero bonds)
    let mut bonds = Vec::with_capacity(num_bonds);
    if let Some(bond_begin) = lines.iter().position(|l| l.contains("M  V30 BEGIN BOND")) {
        let bond_end = lines
            .iter()
            .position(|l| l.contains("M  V30 END BOND"))
            .ok_or_else(|| CyaneaError::Parse("V3000: missing END BOND".into()))?;

        for line in &lines[bond_begin + 1..bond_end] {
            let bond = parse_v3000_bond_line(line)?;
            bonds.push(bond);
        }

        if bonds.len() != num_bonds {
            return Err(CyaneaError::Parse(format!(
                "V3000: expected {num_bonds} bonds, found {}",
                bonds.len()
            )));
        }
    } else if num_bonds != 0 {
        return Err(CyaneaError::Parse(
            "V3000: COUNTS declares bonds but no BOND block found".into(),
        ));
    }

    Ok(Molecule::new(name, atoms, bonds))
}

/// Parse a V3000 atom line: `M  V30 idx symbol x y z aamap [CHG=val] [RAD=val] ...`
fn parse_v3000_atom_line(line: &str) -> Result<MolAtom> {
    let content = line
        .trim()
        .strip_prefix("M  V30 ")
        .or_else(|| line.trim().strip_prefix("M  V30\t"))
        .ok_or_else(|| CyaneaError::Parse(format!("V3000: atom line missing prefix: '{line}'")))?;

    let parts: Vec<&str> = content.split_whitespace().collect();
    // idx(0) symbol(1) x(2) y(3) z(4) aamap(5) [keyword=value ...]
    if parts.len() < 6 {
        return Err(CyaneaError::Parse(format!(
            "V3000: atom line too short: '{line}'"
        )));
    }

    let symbol = parts[1];
    let elem = element_by_symbol(symbol).ok_or_else(|| {
        CyaneaError::Parse(format!("V3000: unknown element '{symbol}'"))
    })?;

    // Parse optional keyword=value pairs
    let mut charge: i8 = 0;
    for part in &parts[6..] {
        if let Some(val) = part.strip_prefix("CHG=") {
            charge = val
                .parse()
                .map_err(|_| CyaneaError::Parse(format!("V3000: invalid CHG value '{val}'")))?;
        }
    }

    Ok(MolAtom {
        atomic_number: elem.atomic_number,
        formal_charge: charge,
        isotope: None,
        is_aromatic: false,
        implicit_hydrogens: 0,
        chirality: crate::molecule::Chirality::None,
    })
}

/// Parse a V3000 bond line: `M  V30 idx type atom1 atom2 [CFG=val] ...`
fn parse_v3000_bond_line(line: &str) -> Result<Bond> {
    let content = line
        .trim()
        .strip_prefix("M  V30 ")
        .or_else(|| line.trim().strip_prefix("M  V30\t"))
        .ok_or_else(|| CyaneaError::Parse(format!("V3000: bond line missing prefix: '{line}'")))?;

    let parts: Vec<&str> = content.split_whitespace().collect();
    // idx(0) type(1) atom1(2) atom2(3) [keyword=value ...]
    if parts.len() < 4 {
        return Err(CyaneaError::Parse(format!(
            "V3000: bond line too short: '{line}'"
        )));
    }

    let bond_type: u8 = parts[1]
        .parse()
        .map_err(|_| CyaneaError::Parse("V3000: invalid bond type".into()))?;
    let a1: usize = parts[2]
        .parse()
        .map_err(|_| CyaneaError::Parse("V3000: invalid bond atom1".into()))?;
    let a2: usize = parts[3]
        .parse()
        .map_err(|_| CyaneaError::Parse("V3000: invalid bond atom2".into()))?;

    if a1 == 0 || a2 == 0 {
        return Err(CyaneaError::Parse(
            "V3000: bond atom indices must be >= 1".into(),
        ));
    }

    let order = match bond_type {
        1 => BondOrder::Single,
        2 => BondOrder::Double,
        3 => BondOrder::Triple,
        4 => BondOrder::Aromatic,
        _ => BondOrder::Single,
    };

    // V3000 atom indices are 1-based
    Ok(Bond {
        atom1: a1 - 1,
        atom2: a2 - 1,
        order,
        is_aromatic: bond_type == 4,
        stereo: crate::molecule::BondStereo::None,
    })
}

/// Detect whether a MOL block uses V3000 format.
fn is_v3000(block: &str) -> bool {
    block.lines().take(5).any(|line| line.contains("V3000"))
        || block.contains("M  V30 BEGIN CTAB")
}

/// Parse a multi-molecule SDF string, returning results for each molecule.
///
/// Auto-detects V2000 vs V3000 format for each molecule block.
pub fn parse_sdf(input: &str) -> Vec<Result<Molecule>> {
    let mut results = Vec::new();

    // Split on $$$$ delimiter
    for block in input.split("$$$$") {
        let trimmed = block.trim();
        if trimmed.is_empty() {
            continue;
        }
        if is_v3000(trimmed) {
            results.push(parse_mol_v3000(trimmed));
        } else {
            results.push(parse_mol_v2000(trimmed));
        }
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
        chirality: crate::molecule::Chirality::None,
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
        stereo: crate::molecule::BondStereo::None,
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

    // ---- V3000 tests ----

    fn minimal_v3000_methane() -> &'static str {
        "\
Methane
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0  0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 END CTAB
M  END"
    }

    #[test]
    fn parse_v3000_minimal_methane() {
        let mol = parse_mol_v3000(minimal_v3000_methane()).unwrap();
        assert_eq!(mol.name(), "Methane");
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(mol.bond_count(), 0);
        assert_eq!(mol.atoms[0].atomic_number, 6);
        assert_eq!(mol.atoms[0].formal_charge, 0);
    }

    #[test]
    fn parse_v3000_with_charges() {
        let mol_str = "\
ChargedMol
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0  0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N 0.0000 0.0000 0.0000 0 CHG=1
M  V30 2 O 1.0000 0.0000 0.0000 0 CHG=-1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 END CTAB
M  END";
        let mol = parse_mol_v3000(mol_str).unwrap();
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.atoms[0].atomic_number, 7); // N
        assert_eq!(mol.atoms[0].formal_charge, 1);
        assert_eq!(mol.atoms[1].atomic_number, 8); // O
        assert_eq!(mol.atoms[1].formal_charge, -1);
    }

    #[test]
    fn parse_v3000_with_bonds() {
        let mol_str = "\
Water
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0  0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 O 0.0000 0.0000 0.0000 0
M  V30 2 H 0.7572 0.5866 0.0000 0
M  V30 3 H -0.7572 0.5866 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 END BOND
M  V30 END CTAB
M  END";
        let mol = parse_mol_v3000(mol_str).unwrap();
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 2);
        // bond 1: O-H (atoms 0,1)
        assert_eq!(mol.bonds[0].atom1, 0);
        assert_eq!(mol.bonds[0].atom2, 1);
        assert_eq!(mol.bonds[0].order, BondOrder::Single);
        // bond 2: O-H (atoms 0,2)
        assert_eq!(mol.bonds[1].atom1, 0);
        assert_eq!(mol.bonds[1].atom2, 2);
    }

    #[test]
    fn parse_sdf_autodetect_v3000() {
        let sdf = format!("{}\n$$$$\n", minimal_v3000_methane());
        let results = parse_sdf(&sdf);
        assert_eq!(results.len(), 1);
        let mol = results[0].as_ref().unwrap();
        assert_eq!(mol.name(), "Methane");
        assert_eq!(mol.atom_count(), 1);
    }

    #[test]
    fn malformed_v3000_error() {
        // Missing BEGIN ATOM
        let bad = "\
Bad
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0  0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 END CTAB
M  END";
        assert!(parse_mol_v3000(bad).is_err());

        // Too short
        assert!(parse_mol_v3000("too\nshort").is_err());

        // Atom count mismatch
        let mismatch = "\
Mismatch
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0  0 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 END CTAB
M  END";
        assert!(parse_mol_v3000(mismatch).is_err());
    }

    #[test]
    fn parse_mixed_sdf_v2000_and_v3000() {
        let v2000_block = minimal_mol();
        let v3000_block = minimal_v3000_methane();
        let sdf = format!("{v2000_block}\n$$$$\n{v3000_block}\n$$$$\n");
        let results = parse_sdf(&sdf);
        assert_eq!(results.len(), 2);
        let mol1 = results[0].as_ref().unwrap();
        let mol2 = results[1].as_ref().unwrap();
        assert_eq!(mol1.name(), "Methane");
        assert_eq!(mol2.name(), "Methane");
        assert_eq!(mol1.atom_count(), 1);
        assert_eq!(mol2.atom_count(), 1);
    }
}
