//! SMILES string parser.

use std::collections::BTreeMap;

use cyanea_core::{CyaneaError, Result};

use crate::element::element_by_symbol;
use crate::molecule::{Bond, BondOrder, MolAtom, Molecule};

/// Parse a SMILES string into a `Molecule`.
pub fn parse_smiles(smiles: &str) -> Result<Molecule> {
    parse_smiles_named(smiles, "")
}

/// Parse a SMILES string into a `Molecule` with a given name.
pub fn parse_smiles_named(smiles: &str, name: &str) -> Result<Molecule> {
    let mut parser = SmilesParser::new(smiles);
    parser.parse()?;
    parser.resolve_ring_closures()?;
    parser.compute_implicit_hydrogens();
    Ok(Molecule::new(name.to_string(), parser.atoms, parser.bonds))
}

struct SmilesParser<'a> {
    input: &'a [u8],
    pos: usize,
    atoms: Vec<MolAtom>,
    bonds: Vec<Bond>,
    /// ring_closures[digit] = (atom_idx, Option<BondOrder>)
    ring_closures: BTreeMap<u16, (usize, Option<BondOrder>)>,
    /// Stack of atom indices for branch handling
    stack: Vec<usize>,
    /// Index of the previous atom (for bonding)
    prev_atom: Option<usize>,
    /// Pending bond order for the next bond
    pending_bond: Option<BondOrder>,
}

impl<'a> SmilesParser<'a> {
    fn new(input: &'a str) -> Self {
        SmilesParser {
            input: input.as_bytes(),
            pos: 0,
            atoms: Vec::new(),
            bonds: Vec::new(),
            ring_closures: BTreeMap::new(),
            stack: Vec::new(),
            prev_atom: None,
            pending_bond: None,
        }
    }

    fn peek(&self) -> Option<u8> {
        self.input.get(self.pos).copied()
    }

    fn advance(&mut self) -> Option<u8> {
        let ch = self.input.get(self.pos).copied();
        if ch.is_some() {
            self.pos += 1;
        }
        ch
    }

    fn parse(&mut self) -> Result<()> {
        while self.pos < self.input.len() {
            match self.peek() {
                Some(b'(') => {
                    self.advance();
                    if let Some(prev) = self.prev_atom {
                        self.stack.push(prev);
                    }
                }
                Some(b')') => {
                    self.advance();
                    self.prev_atom = self.stack.pop();
                    self.pending_bond = None;
                }
                Some(b'-') => {
                    self.advance();
                    self.pending_bond = Some(BondOrder::Single);
                }
                Some(b'=') => {
                    self.advance();
                    self.pending_bond = Some(BondOrder::Double);
                }
                Some(b'#') => {
                    self.advance();
                    self.pending_bond = Some(BondOrder::Triple);
                }
                Some(b':') => {
                    self.advance();
                    self.pending_bond = Some(BondOrder::Aromatic);
                }
                Some(b'/') | Some(b'\\') => {
                    // Stereo bond markers — consume and ignore
                    self.advance();
                }
                Some(b'%') => {
                    self.advance();
                    let ring_num = self.parse_two_digit_ring()?;
                    self.handle_ring_closure(ring_num)?;
                }
                Some(b'[') => {
                    self.parse_bracket_atom()?;
                }
                Some(ch) if ch.is_ascii_digit() => {
                    self.advance();
                    let ring_num = (ch - b'0') as u16;
                    self.handle_ring_closure(ring_num)?;
                }
                Some(ch) if is_organic_atom_start(ch) => {
                    self.parse_organic_atom()?;
                }
                Some(b'.') => {
                    // Disconnected fragments
                    self.advance();
                    self.prev_atom = None;
                    self.pending_bond = None;
                }
                Some(ch) => {
                    return Err(CyaneaError::Parse(format!(
                        "unexpected character '{}' at position {}",
                        ch as char, self.pos
                    )));
                }
                None => break,
            }
        }
        Ok(())
    }

    fn parse_organic_atom(&mut self) -> Result<()> {
        let ch = self.advance().unwrap();
        let is_aromatic = ch.is_ascii_lowercase();
        let upper = ch.to_ascii_uppercase();

        let symbol = match upper {
            b'B' => {
                if !is_aromatic && self.peek() == Some(b'r') {
                    self.advance();
                    "Br"
                } else {
                    "B"
                }
            }
            b'C' => {
                if !is_aromatic && self.peek() == Some(b'l') {
                    self.advance();
                    "Cl"
                } else {
                    "C"
                }
            }
            b'N' => "N",
            b'O' => "O",
            b'P' => "P",
            b'S' => {
                if !is_aromatic && self.peek() == Some(b'i') {
                    self.advance();
                    "Si"
                } else if !is_aromatic && self.peek() == Some(b'e') {
                    self.advance();
                    "Se"
                } else {
                    "S"
                }
            }
            b'F' => "F",
            b'I' => "I",
            _ => {
                return Err(CyaneaError::Parse(format!(
                    "unknown organic atom '{}'",
                    upper as char
                )));
            }
        };

        let elem = element_by_symbol(symbol).ok_or_else(|| {
            CyaneaError::Parse(format!("unknown element '{symbol}'"))
        })?;

        let atom = MolAtom {
            atomic_number: elem.atomic_number,
            formal_charge: 0,
            isotope: None,
            is_aromatic,
            implicit_hydrogens: 0, // computed later
        };

        let atom_idx = self.atoms.len();
        self.atoms.push(atom);
        self.add_bond_to_prev(atom_idx)?;
        self.prev_atom = Some(atom_idx);
        Ok(())
    }

    fn parse_bracket_atom(&mut self) -> Result<()> {
        self.advance(); // consume '['

        // Optional isotope
        let isotope = self.parse_optional_number();

        // Atom symbol
        let ch = self.advance().ok_or_else(|| {
            CyaneaError::Parse("unexpected end of SMILES in bracket atom".into())
        })?;

        let is_aromatic = ch.is_ascii_lowercase();
        let upper = ch.to_ascii_uppercase();

        // Try two-letter symbol first
        let symbol = if let Some(next) = self.peek() {
            if next.is_ascii_lowercase() && next != b'@' {
                let two_letter = format!("{}{}", upper as char, next as char);
                if element_by_symbol(&two_letter).is_some() {
                    self.advance();
                    two_letter
                } else {
                    String::from(upper as char)
                }
            } else {
                String::from(upper as char)
            }
        } else {
            String::from(upper as char)
        };

        let elem = element_by_symbol(&symbol).ok_or_else(|| {
            CyaneaError::Parse(format!("unknown element '{symbol}'"))
        })?;

        // Skip stereochemistry markers
        while self.peek() == Some(b'@') {
            self.advance();
        }

        // Optional hydrogen count
        let mut explicit_h = 0u8;
        if self.peek() == Some(b'H') {
            self.advance();
            if let Some(d) = self.peek() {
                if d.is_ascii_digit() {
                    self.advance();
                    explicit_h = d - b'0';
                } else {
                    explicit_h = 1;
                }
            } else {
                explicit_h = 1;
            }
        }

        // Optional charge
        let mut charge: i8 = 0;
        match self.peek() {
            Some(b'+') => {
                self.advance();
                charge = if let Some(d) = self.peek() {
                    if d.is_ascii_digit() {
                        self.advance();
                        (d - b'0') as i8
                    } else if self.peek() == Some(b'+') {
                        // Count consecutive +
                        let mut c = 1i8;
                        while self.peek() == Some(b'+') {
                            self.advance();
                            c += 1;
                        }
                        c
                    } else {
                        1
                    }
                } else {
                    1
                };
            }
            Some(b'-') => {
                self.advance();
                charge = if let Some(d) = self.peek() {
                    if d.is_ascii_digit() {
                        self.advance();
                        -((d - b'0') as i8)
                    } else if self.peek() == Some(b'-') {
                        let mut c = -1i8;
                        while self.peek() == Some(b'-') {
                            self.advance();
                            c -= 1;
                        }
                        c
                    } else {
                        -1
                    }
                } else {
                    -1
                };
            }
            _ => {}
        }

        // Closing bracket
        if self.advance() != Some(b']') {
            return Err(CyaneaError::Parse("expected ']' in bracket atom".into()));
        }

        let atom = MolAtom {
            atomic_number: elem.atomic_number,
            formal_charge: charge,
            isotope: isotope.map(|n| n as u16),
            is_aromatic,
            implicit_hydrogens: explicit_h, // bracket atoms specify H explicitly
        };

        let atom_idx = self.atoms.len();
        self.atoms.push(atom);
        self.add_bond_to_prev(atom_idx)?;
        self.prev_atom = Some(atom_idx);
        Ok(())
    }

    fn parse_optional_number(&mut self) -> Option<u32> {
        let mut n: u32 = 0;
        let mut found = false;
        while let Some(ch) = self.peek() {
            if ch.is_ascii_digit() {
                self.advance();
                n = n * 10 + (ch - b'0') as u32;
                found = true;
            } else {
                break;
            }
        }
        if found { Some(n) } else { None }
    }

    fn parse_two_digit_ring(&mut self) -> Result<u16> {
        let d1 = self.advance().ok_or_else(|| {
            CyaneaError::Parse("expected digit after '%'".into())
        })?;
        let d2 = self.advance().ok_or_else(|| {
            CyaneaError::Parse("expected second digit after '%'".into())
        })?;
        if !d1.is_ascii_digit() || !d2.is_ascii_digit() {
            return Err(CyaneaError::Parse("invalid ring closure number after '%'".into()));
        }
        Ok((d1 - b'0') as u16 * 10 + (d2 - b'0') as u16)
    }

    fn handle_ring_closure(&mut self, ring_num: u16) -> Result<()> {
        let current = self.prev_atom.ok_or_else(|| {
            CyaneaError::Parse("ring closure without preceding atom".into())
        })?;

        if let Some((open_atom, open_bond)) = self.ring_closures.remove(&ring_num) {
            // Close the ring
            let order = self.pending_bond.or(open_bond).unwrap_or(BondOrder::Single);
            let is_aromatic = self.atoms[open_atom].is_aromatic && self.atoms[current].is_aromatic;
            let order = if is_aromatic && order == BondOrder::Single {
                BondOrder::Aromatic
            } else {
                order
            };
            self.bonds.push(Bond {
                atom1: open_atom,
                atom2: current,
                order,
                is_aromatic,
            });
            self.pending_bond = None;
        } else {
            // Open a ring closure
            self.ring_closures.insert(ring_num, (current, self.pending_bond.take()));
        }
        Ok(())
    }

    fn add_bond_to_prev(&mut self, atom_idx: usize) -> Result<()> {
        if let Some(prev) = self.prev_atom {
            let both_aromatic = self.atoms[prev].is_aromatic && self.atoms[atom_idx].is_aromatic;
            let order = self.pending_bond.take().unwrap_or_else(|| {
                if both_aromatic {
                    BondOrder::Aromatic
                } else {
                    BondOrder::Single
                }
            });
            let is_aromatic = both_aromatic && order == BondOrder::Aromatic;
            self.bonds.push(Bond {
                atom1: prev,
                atom2: atom_idx,
                order,
                is_aromatic,
            });
        }
        self.pending_bond = None;
        Ok(())
    }

    fn resolve_ring_closures(&self) -> Result<()> {
        if !self.ring_closures.is_empty() {
            let open: Vec<_> = self.ring_closures.keys().collect();
            return Err(CyaneaError::Parse(format!(
                "unmatched ring closure(s): {:?}",
                open
            )));
        }
        if !self.stack.is_empty() {
            return Err(CyaneaError::Parse(format!(
                "{} unmatched '(' in SMILES",
                self.stack.len()
            )));
        }
        Ok(())
    }

    fn compute_implicit_hydrogens(&mut self) {
        // For organic subset atoms only (bracket atoms already have explicit H)
        for i in 0..self.atoms.len() {
            let atom = &self.atoms[i];
            // Bracket atoms have H already set — identified by having isotope or charge
            if atom.isotope.is_some() || atom.formal_charge != 0 {
                continue;
            }

            let degree = self.bond_degree(i);
            let is_aromatic = atom.is_aromatic;
            let target = self.target_valence(atom.atomic_number);

            if let Some(target) = target {
                // For aromatic atoms: 1 electron goes to the pi system,
                // remaining (target - 1) are available for sigma bonds + implicit H.
                // For non-aromatic atoms: all target valence electrons are for bonds + H.
                let available = if is_aromatic {
                    target.saturating_sub(1)
                } else {
                    target
                };
                // Implicit H = available sigma capacity - actual sigma bonds (degree)
                // For non-aromatic atoms with double/triple bonds, use the full bond-order sum
                let used = if is_aromatic {
                    degree // aromatic bonds count as 1 sigma each
                } else {
                    self.bond_order_sum(i)
                };
                if available > used {
                    self.atoms[i].implicit_hydrogens = (available - used) as u8;
                }
            }
        }
    }

    /// Number of explicit bonds (graph degree) for an atom.
    fn bond_degree(&self, atom_idx: usize) -> usize {
        self.bonds
            .iter()
            .filter(|b| b.atom1 == atom_idx || b.atom2 == atom_idx)
            .count()
    }

    /// Sum of bond orders (integer-rounded) for non-aromatic valence calculation.
    fn bond_order_sum(&self, atom_idx: usize) -> usize {
        let mut v = 0.0f64;
        for bond in &self.bonds {
            if bond.atom1 == atom_idx || bond.atom2 == atom_idx {
                v += bond.order.as_f64();
            }
        }
        v.round() as usize
    }

    fn target_valence(&self, atomic_number: u8) -> Option<usize> {
        match atomic_number {
            5 => Some(3),    // B
            6 => Some(4),    // C
            7 => Some(3),    // N
            8 => Some(2),    // O
            15 => Some(3),   // P
            16 => Some(2),   // S
            9 => Some(1),    // F
            17 => Some(1),   // Cl
            35 => Some(1),   // Br
            53 => Some(1),   // I
            _ => None,
        }
    }
}

fn is_organic_atom_start(ch: u8) -> bool {
    matches!(
        ch,
        b'B' | b'C' | b'N' | b'O' | b'P' | b'S' | b'F' | b'I'
            | b'b' | b'c' | b'n' | b'o' | b'p' | b's'
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_methane() {
        let mol = parse_smiles("C").unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(mol.bond_count(), 0);
        assert_eq!(mol.atoms[0].atomic_number, 6);
        assert_eq!(mol.atoms[0].implicit_hydrogens, 4);
    }

    #[test]
    fn parse_ethanol() {
        let mol = parse_smiles("CCO").unwrap();
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 2);
        // C has 3H, C has 2H, O has 1H
        assert_eq!(mol.atoms[0].implicit_hydrogens, 3);
        assert_eq!(mol.atoms[1].implicit_hydrogens, 2);
        assert_eq!(mol.atoms[2].implicit_hydrogens, 1);
    }

    #[test]
    fn parse_benzene() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 6); // 5 chain + 1 ring closure
        for atom in &mol.atoms {
            assert!(atom.is_aromatic);
            assert_eq!(atom.implicit_hydrogens, 1);
        }
    }

    #[test]
    fn parse_branching() {
        // Isobutane: CC(C)C
        let mol = parse_smiles("CC(C)C").unwrap();
        assert_eq!(mol.atom_count(), 4);
        assert_eq!(mol.bond_count(), 3);
        assert_eq!(mol.degree(1), 3); // Central carbon
    }

    #[test]
    fn parse_double_bond() {
        // Ethene: C=C
        let mol = parse_smiles("C=C").unwrap();
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 1);
        assert_eq!(mol.bonds[0].order, BondOrder::Double);
        assert_eq!(mol.atoms[0].implicit_hydrogens, 2);
        assert_eq!(mol.atoms[1].implicit_hydrogens, 2);
    }

    #[test]
    fn parse_bracket_atom() {
        // Ammonium: [NH4+]
        let mol = parse_smiles("[NH4+]").unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(mol.atoms[0].atomic_number, 7);
        assert_eq!(mol.atoms[0].formal_charge, 1);
        assert_eq!(mol.atoms[0].implicit_hydrogens, 4);
    }

    #[test]
    fn parse_two_digit_ring_closure() {
        // Ring closure with %nn notation
        let mol = parse_smiles("C%10CCCCCCCCC%10").unwrap();
        assert_eq!(mol.atom_count(), 10);
        // 9 chain bonds + 1 ring closure
        assert_eq!(mol.bond_count(), 10);
    }

    #[test]
    fn invalid_smiles_error() {
        assert!(parse_smiles("C(").is_err());
        assert!(parse_smiles("C1CC").is_err()); // unmatched ring closure
        assert!(parse_smiles("[").is_err());
    }
}

#[cfg(test)]
mod proptests {
    use super::*;
    use crate::properties::molecular_formula;
    use proptest::prelude::*;

    /// Strategy for valid simple SMILES: chains of organic subset atoms
    fn simple_smiles() -> impl Strategy<Value = String> {
        let atoms = prop_oneof![
            Just("C"),
            Just("N"),
            Just("O"),
            Just("S"),
            Just("c"),
            Just("n"),
            Just("o"),
        ];
        proptest::collection::vec(atoms, 1..=20)
            .prop_map(|parts| parts.join(""))
    }

    proptest! {
        #[test]
        fn parse_smiles_does_not_panic(s in "\\PC{0,100}") {
            let _ = parse_smiles(&s);
        }

        #[test]
        fn formula_is_deterministic(smi in simple_smiles()) {
            if let Ok(mol) = parse_smiles(&smi) {
                let f1 = molecular_formula(&mol);
                let f2 = molecular_formula(&mol);
                prop_assert_eq!(f1, f2);
            }
        }

        #[test]
        fn atom_count_positive_on_success(smi in simple_smiles()) {
            if let Ok(mol) = parse_smiles(&smi) {
                prop_assert!(mol.atom_count() > 0);
            }
        }
    }
}
