//! SMARTS pattern parsing and substructure matching.
//!
//! SMARTS extends SMILES with atom/bond query primitives and logical operators
//! for flexible substructure searching.

use cyanea_core::{CyaneaError, Result};

use crate::molecule::{BondOrder, Molecule};
use crate::ring;
use crate::substructure::SubstructureMatch;

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/// A primitive atom query in SMARTS.
#[derive(Debug, Clone, PartialEq)]
pub enum AtomPrimitive {
    /// Match by atomic number (`#6` = carbon).
    AtomicNum(u8),
    /// Match aromatic atoms (`a`).
    Aromatic,
    /// Match aliphatic atoms (`A`).
    Aliphatic,
    /// Match by graph degree (`D2`).
    Degree(u8),
    /// Match by explicit hydrogen count (`H1`).
    HCount(u8),
    /// Match by total hydrogen count including implicit (`h1`).
    TotalHCount(u8),
    /// Match by formal charge (`+1`, `-2`).
    Charge(i8),
    /// Match atoms in any ring (`R`).
    RingMember,
    /// Match atoms in a ring of specific size (`r5`).
    RingSize(u8),
    /// Match by connectivity / total degree including implicit H (`X3`).
    Connectivity(u8),
    /// Match by total valence (`v3`).
    Valence(u8),
    /// Wildcard — matches any atom (`*`).
    Wildcard,
}

/// A logical atom expression in SMARTS.
#[derive(Debug, Clone, PartialEq)]
pub enum AtomExpr {
    /// A single primitive test.
    Prim(AtomPrimitive),
    /// High-precedence AND (juxtaposition or `&`).
    And(Vec<AtomExpr>),
    /// OR (`,`).
    Or(Vec<AtomExpr>),
    /// Negation (`!`).
    Not(Box<AtomExpr>),
    /// Recursive SMARTS (`$(...)`).
    Recursive(SmartsPattern),
}

/// A bond expression in SMARTS.
#[derive(Debug, Clone, PartialEq)]
pub enum BondExpr {
    Single,
    Double,
    Triple,
    Aromatic,
    Ring,
    Any,
    Not(Box<BondExpr>),
    And(Vec<BondExpr>),
    Or(Vec<BondExpr>),
}

/// A single atom in a SMARTS pattern.
#[derive(Debug, Clone, PartialEq)]
pub struct SmartsAtom {
    pub expr: AtomExpr,
}

/// A bond between two atoms in a SMARTS pattern.
#[derive(Debug, Clone, PartialEq)]
pub struct SmartsBond {
    pub atom1: usize,
    pub atom2: usize,
    pub expr: BondExpr,
}

/// A parsed SMARTS pattern.
#[derive(Debug, Clone, PartialEq)]
pub struct SmartsPattern {
    pub atoms: Vec<SmartsAtom>,
    pub bonds: Vec<SmartsBond>,
    adjacency: Vec<Vec<(usize, usize)>>,
}

impl SmartsPattern {
    fn new(atoms: Vec<SmartsAtom>, bonds: Vec<SmartsBond>) -> Self {
        let mut adjacency = vec![Vec::new(); atoms.len()];
        for (bi, bond) in bonds.iter().enumerate() {
            adjacency[bond.atom1].push((bond.atom2, bi));
            adjacency[bond.atom2].push((bond.atom1, bi));
        }
        SmartsPattern { atoms, bonds, adjacency }
    }
}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

struct SmartsParser<'a> {
    input: &'a [u8],
    pos: usize,
    atoms: Vec<SmartsAtom>,
    bonds: Vec<SmartsBond>,
    stack: Vec<usize>,
    prev_atom: Option<usize>,
    pending_bond: Option<BondExpr>,
    ring_closures: std::collections::BTreeMap<u16, (usize, Option<BondExpr>)>,
}

impl<'a> SmartsParser<'a> {
    fn new(input: &'a str) -> Self {
        SmartsParser {
            input: input.as_bytes(),
            pos: 0,
            atoms: Vec::new(),
            bonds: Vec::new(),
            stack: Vec::new(),
            prev_atom: None,
            pending_bond: None,
            ring_closures: std::collections::BTreeMap::new(),
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
                    self.pending_bond = Some(BondExpr::Single);
                }
                Some(b'=') => {
                    self.advance();
                    self.pending_bond = Some(BondExpr::Double);
                }
                Some(b'#') if !self.is_atom_num_context() => {
                    self.advance();
                    self.pending_bond = Some(BondExpr::Triple);
                }
                Some(b':') => {
                    self.advance();
                    self.pending_bond = Some(BondExpr::Aromatic);
                }
                Some(b'~') => {
                    self.advance();
                    self.pending_bond = Some(BondExpr::Any);
                }
                Some(b'@') => {
                    self.advance();
                    self.pending_bond = Some(BondExpr::Ring);
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
                Some(b'*') => {
                    self.advance();
                    let atom = SmartsAtom { expr: AtomExpr::Prim(AtomPrimitive::Wildcard) };
                    let idx = self.atoms.len();
                    self.atoms.push(atom);
                    self.add_bond_to_prev(idx);
                    self.prev_atom = Some(idx);
                }
                Some(b'.') => {
                    self.advance();
                    self.prev_atom = None;
                    self.pending_bond = None;
                }
                Some(ch) if is_organic_smarts(ch) => {
                    self.parse_organic_atom()?;
                }
                Some(ch) => {
                    return Err(CyaneaError::Parse(format!(
                        "unexpected character '{}' at position {} in SMARTS",
                        ch as char, self.pos
                    )));
                }
                None => break,
            }
        }
        self.resolve_ring_closures()?;
        Ok(())
    }

    /// Check if `#` is in an atom-number context (inside brackets after `[`).
    fn is_atom_num_context(&self) -> bool {
        // `#` as triple bond only outside brackets — but we handle bracket atoms
        // separately, so this is always a bond context.
        false
    }

    fn parse_organic_atom(&mut self) -> Result<()> {
        let ch = self.advance().unwrap();
        let is_aromatic = ch.is_ascii_lowercase();
        let upper = ch.to_ascii_uppercase();

        // Two-letter organic subset
        let atomic_num = match upper {
            b'B' => {
                if !is_aromatic && self.peek() == Some(b'r') {
                    self.advance();
                    35 // Br
                } else {
                    5 // B
                }
            }
            b'C' => {
                if !is_aromatic && self.peek() == Some(b'l') {
                    self.advance();
                    17 // Cl
                } else {
                    6 // C
                }
            }
            b'N' => 7,
            b'O' => 8,
            b'P' => 15,
            b'S' => {
                if !is_aromatic && self.peek() == Some(b'i') {
                    self.advance();
                    14 // Si
                } else if !is_aromatic && self.peek() == Some(b'e') {
                    self.advance();
                    34 // Se
                } else {
                    16 // S
                }
            }
            b'F' => 9,
            b'I' => 53,
            _ => {
                return Err(CyaneaError::Parse(format!(
                    "unknown organic atom '{}' in SMARTS",
                    upper as char
                )));
            }
        };

        let expr = if is_aromatic {
            AtomExpr::And(vec![
                AtomExpr::Prim(AtomPrimitive::AtomicNum(atomic_num)),
                AtomExpr::Prim(AtomPrimitive::Aromatic),
            ])
        } else {
            AtomExpr::Prim(AtomPrimitive::AtomicNum(atomic_num))
        };

        let atom = SmartsAtom { expr };
        let idx = self.atoms.len();
        self.atoms.push(atom);
        self.add_bond_to_prev(idx);
        self.prev_atom = Some(idx);
        Ok(())
    }

    fn parse_bracket_atom(&mut self) -> Result<()> {
        self.advance(); // consume '['

        // Check for recursive SMARTS: $( ... )
        if self.peek() == Some(b'$') {
            self.advance(); // consume '$'
            if self.advance() != Some(b'(') {
                return Err(CyaneaError::Parse("expected '(' after '$' in SMARTS".into()));
            }
            // Find matching ')'
            let start = self.pos;
            let mut depth = 1i32;
            while self.pos < self.input.len() && depth > 0 {
                match self.input[self.pos] {
                    b'(' => depth += 1,
                    b')' => depth -= 1,
                    _ => {}
                }
                if depth > 0 {
                    self.pos += 1;
                }
            }
            if depth != 0 {
                return Err(CyaneaError::Parse("unmatched '(' in recursive SMARTS".into()));
            }
            let inner = std::str::from_utf8(&self.input[start..self.pos])
                .map_err(|_| CyaneaError::Parse("invalid UTF-8 in recursive SMARTS".into()))?;
            self.pos += 1; // consume ')'

            let sub_pattern = parse_smarts(inner)?;

            if self.advance() != Some(b']') {
                return Err(CyaneaError::Parse("expected ']' after recursive SMARTS".into()));
            }

            let atom = SmartsAtom { expr: AtomExpr::Recursive(sub_pattern) };
            let idx = self.atoms.len();
            self.atoms.push(atom);
            self.add_bond_to_prev(idx);
            self.prev_atom = Some(idx);
            return Ok(());
        }

        // Parse atom expression with logical operators
        let expr = self.parse_atom_or()?;

        if self.advance() != Some(b']') {
            return Err(CyaneaError::Parse("expected ']' in SMARTS bracket atom".into()));
        }

        let atom = SmartsAtom { expr };
        let idx = self.atoms.len();
        self.atoms.push(atom);
        self.add_bond_to_prev(idx);
        self.prev_atom = Some(idx);
        Ok(())
    }

    // Atom expression parsing with operator precedence:
    //   or_expr  = and_expr (',' and_expr)*
    //   and_expr = not_expr (('&' | implicit_and) not_expr)*
    //   not_expr = '!'? primitive

    fn parse_atom_or(&mut self) -> Result<AtomExpr> {
        let mut terms = vec![self.parse_atom_and()?];
        while self.peek() == Some(b',') {
            self.advance();
            terms.push(self.parse_atom_and()?);
        }
        if terms.len() == 1 {
            Ok(terms.remove(0))
        } else {
            Ok(AtomExpr::Or(terms))
        }
    }

    fn parse_atom_and(&mut self) -> Result<AtomExpr> {
        let mut terms = vec![self.parse_atom_not()?];
        loop {
            match self.peek() {
                Some(b'&') => {
                    self.advance();
                    terms.push(self.parse_atom_not()?);
                }
                // Implicit AND: next char starts a new primitive (not ], not ,)
                Some(ch) if ch != b']' && ch != b',' && is_atom_prim_start(ch) => {
                    terms.push(self.parse_atom_not()?);
                }
                _ => break,
            }
        }
        if terms.len() == 1 {
            Ok(terms.remove(0))
        } else {
            Ok(AtomExpr::And(terms))
        }
    }

    fn parse_atom_not(&mut self) -> Result<AtomExpr> {
        if self.peek() == Some(b'!') {
            self.advance();
            let inner = self.parse_atom_primitive()?;
            Ok(AtomExpr::Not(Box::new(inner)))
        } else {
            self.parse_atom_primitive()
        }
    }

    fn parse_atom_primitive(&mut self) -> Result<AtomExpr> {
        match self.peek() {
            Some(b'#') => {
                self.advance();
                let n = self.parse_number()? as u8;
                Ok(AtomExpr::Prim(AtomPrimitive::AtomicNum(n)))
            }
            Some(b'a') => {
                self.advance();
                Ok(AtomExpr::Prim(AtomPrimitive::Aromatic))
            }
            Some(b'A') => {
                self.advance();
                // Check if followed by lowercase letter -> element symbol like Al, As, Ag
                if let Some(next) = self.peek() {
                    if next == b'l' || next == b's' || next == b'g' || next == b'r' || next == b'u' {
                        // This is an element symbol, not the aliphatic primitive
                        // But in SMARTS context, plain 'A' means aliphatic
                        // Element symbols inside brackets need to be distinguished
                        // We treat standalone 'A' as aliphatic
                    }
                }
                Ok(AtomExpr::Prim(AtomPrimitive::Aliphatic))
            }
            Some(b'D') => {
                self.advance();
                let n = self.parse_optional_digit().unwrap_or(1);
                Ok(AtomExpr::Prim(AtomPrimitive::Degree(n)))
            }
            Some(b'H') => {
                self.advance();
                let n = self.parse_optional_digit().unwrap_or(1);
                Ok(AtomExpr::Prim(AtomPrimitive::HCount(n)))
            }
            Some(b'h') => {
                self.advance();
                let n = self.parse_optional_digit().unwrap_or(1);
                Ok(AtomExpr::Prim(AtomPrimitive::TotalHCount(n)))
            }
            Some(b'R') => {
                self.advance();
                Ok(AtomExpr::Prim(AtomPrimitive::RingMember))
            }
            Some(b'r') => {
                self.advance();
                if let Some(n) = self.parse_optional_digit() {
                    Ok(AtomExpr::Prim(AtomPrimitive::RingSize(n)))
                } else {
                    Ok(AtomExpr::Prim(AtomPrimitive::RingMember))
                }
            }
            Some(b'X') => {
                self.advance();
                let n = self.parse_optional_digit().unwrap_or(1);
                Ok(AtomExpr::Prim(AtomPrimitive::Connectivity(n)))
            }
            Some(b'v') => {
                self.advance();
                let n = self.parse_optional_digit().unwrap_or(1);
                Ok(AtomExpr::Prim(AtomPrimitive::Valence(n)))
            }
            Some(b'*') => {
                self.advance();
                Ok(AtomExpr::Prim(AtomPrimitive::Wildcard))
            }
            Some(b'+') => {
                self.advance();
                let n = self.parse_optional_digit().unwrap_or(1);
                Ok(AtomExpr::Prim(AtomPrimitive::Charge(n as i8)))
            }
            Some(b'-') => {
                self.advance();
                let n = self.parse_optional_digit().unwrap_or(1);
                Ok(AtomExpr::Prim(AtomPrimitive::Charge(-(n as i8))))
            }
            // Organic element symbols inside brackets
            Some(ch) if ch.is_ascii_uppercase() || ch.is_ascii_lowercase() => {
                let atomic_num = self.parse_element_symbol()?;
                Ok(AtomExpr::Prim(AtomPrimitive::AtomicNum(atomic_num)))
            }
            Some(ch) => Err(CyaneaError::Parse(format!(
                "unexpected '{}' in SMARTS atom expression at position {}",
                ch as char, self.pos
            ))),
            None => Err(CyaneaError::Parse("unexpected end of SMARTS".into())),
        }
    }

    fn parse_element_symbol(&mut self) -> Result<u8> {
        let ch = self.advance().unwrap();
        let is_aromatic = ch.is_ascii_lowercase();
        let upper = ch.to_ascii_uppercase();

        // Try two-letter elements
        let (symbol, atomic_num) = if let Some(next) = self.peek() {
            if next.is_ascii_lowercase() && next != b',' && next != b'&' && next != b']' && next != b'!' {
                let two = format!("{}{}", upper as char, next as char);
                if let Some(num) = symbol_to_atomic_num(&two) {
                    self.advance();
                    (two, num)
                } else {
                    let one = format!("{}", upper as char);
                    let num = symbol_to_atomic_num(&one).ok_or_else(|| {
                        CyaneaError::Parse(format!("unknown element '{}' in SMARTS", upper as char))
                    })?;
                    (one, num)
                }
            } else {
                let one = format!("{}", upper as char);
                let num = symbol_to_atomic_num(&one).ok_or_else(|| {
                    CyaneaError::Parse(format!("unknown element '{}' in SMARTS", upper as char))
                })?;
                (one, num)
            }
        } else {
            let one = format!("{}", upper as char);
            let num = symbol_to_atomic_num(&one).ok_or_else(|| {
                CyaneaError::Parse(format!("unknown element '{}' in SMARTS", upper as char))
            })?;
            (one, num)
        };

        let _ = symbol;

        if is_aromatic {
            // Aromatic element in SMARTS atom context — still returns atomic number.
            // The aromaticity check is handled separately if needed.
        }

        Ok(atomic_num)
    }

    fn parse_number(&mut self) -> Result<u32> {
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
        if !found {
            return Err(CyaneaError::Parse(format!(
                "expected number at position {} in SMARTS",
                self.pos
            )));
        }
        Ok(n)
    }

    fn parse_optional_digit(&mut self) -> Option<u8> {
        if let Some(ch) = self.peek() {
            if ch.is_ascii_digit() {
                self.advance();
                return Some(ch - b'0');
            }
        }
        None
    }

    fn parse_two_digit_ring(&mut self) -> Result<u16> {
        let d1 = self.advance().ok_or_else(|| {
            CyaneaError::Parse("expected digit after '%' in SMARTS".into())
        })?;
        let d2 = self.advance().ok_or_else(|| {
            CyaneaError::Parse("expected second digit after '%' in SMARTS".into())
        })?;
        if !d1.is_ascii_digit() || !d2.is_ascii_digit() {
            return Err(CyaneaError::Parse("invalid ring closure after '%' in SMARTS".into()));
        }
        Ok((d1 - b'0') as u16 * 10 + (d2 - b'0') as u16)
    }

    fn handle_ring_closure(&mut self, ring_num: u16) -> Result<()> {
        let current = self.prev_atom.ok_or_else(|| {
            CyaneaError::Parse("ring closure without preceding atom in SMARTS".into())
        })?;

        if let Some((open_atom, open_bond)) = self.ring_closures.remove(&ring_num) {
            let expr = self.pending_bond.take().or(open_bond).unwrap_or(BondExpr::Any);
            self.bonds.push(SmartsBond { atom1: open_atom, atom2: current, expr });
        } else {
            self.ring_closures.insert(ring_num, (current, self.pending_bond.take()));
        }
        Ok(())
    }

    fn add_bond_to_prev(&mut self, atom_idx: usize) {
        if let Some(prev) = self.prev_atom {
            let expr = self.pending_bond.take().unwrap_or(BondExpr::Any);
            self.bonds.push(SmartsBond { atom1: prev, atom2: atom_idx, expr });
        }
        self.pending_bond = None;
    }

    fn resolve_ring_closures(&self) -> Result<()> {
        if !self.ring_closures.is_empty() {
            let open: Vec<_> = self.ring_closures.keys().collect();
            return Err(CyaneaError::Parse(format!(
                "unmatched ring closure(s) in SMARTS: {:?}",
                open
            )));
        }
        Ok(())
    }
}

fn is_organic_smarts(ch: u8) -> bool {
    matches!(
        ch,
        b'B' | b'C' | b'N' | b'O' | b'P' | b'S' | b'F' | b'I'
            | b'b' | b'c' | b'n' | b'o' | b'p' | b's'
    )
}

fn is_atom_prim_start(ch: u8) -> bool {
    matches!(
        ch,
        b'#' | b'a' | b'A' | b'D' | b'H' | b'h' | b'R' | b'r'
            | b'X' | b'v' | b'*' | b'+' | b'-' | b'!'
    ) || ch.is_ascii_alphabetic()
}

fn symbol_to_atomic_num(symbol: &str) -> Option<u8> {
    match symbol {
        "H" => Some(1), "He" => Some(2), "Li" => Some(3), "Be" => Some(4),
        "B" => Some(5), "C" => Some(6), "N" => Some(7), "O" => Some(8),
        "F" => Some(9), "Ne" => Some(10), "Na" => Some(11), "Mg" => Some(12),
        "Al" => Some(13), "Si" => Some(14), "P" => Some(15), "S" => Some(16),
        "Cl" => Some(17), "Ar" => Some(18), "K" => Some(19), "Ca" => Some(20),
        "Sc" => Some(21), "Ti" => Some(22), "V" => Some(23), "Cr" => Some(24),
        "Mn" => Some(25), "Fe" => Some(26), "Co" => Some(27), "Ni" => Some(28),
        "Cu" => Some(29), "Zn" => Some(30), "Ga" => Some(31), "Ge" => Some(32),
        "As" => Some(33), "Se" => Some(34), "Br" => Some(35), "Kr" => Some(36),
        "Rb" => Some(37), "Sr" => Some(38), "I" => Some(53),
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Parse a SMARTS string into a `SmartsPattern`.
pub fn parse_smarts(smarts: &str) -> Result<SmartsPattern> {
    if smarts.is_empty() {
        return Err(CyaneaError::Parse("empty SMARTS string".into()));
    }
    let mut parser = SmartsParser::new(smarts);
    parser.parse()?;
    Ok(SmartsPattern::new(parser.atoms, parser.bonds))
}

/// Check if any substructure of `target` matches the SMARTS `pattern`.
pub fn smarts_match(target: &Molecule, pattern: &SmartsPattern) -> bool {
    let rings = ring::find_sssr(target);
    let ring_membership = build_ring_membership(target, &rings);
    let ring_sizes = build_ring_sizes(target, &rings);
    let mut state = SmartsVf2State::new(target, pattern, &rings, &ring_membership, &ring_sizes);
    state.search(true);
    !state.matches.is_empty()
}

/// Find all substructure matches of SMARTS `pattern` in `target`.
pub fn smarts_find_all(target: &Molecule, pattern: &SmartsPattern) -> Vec<SubstructureMatch> {
    let rings = ring::find_sssr(target);
    let ring_membership = build_ring_membership(target, &rings);
    let ring_sizes = build_ring_sizes(target, &rings);
    let mut state = SmartsVf2State::new(target, pattern, &rings, &ring_membership, &ring_sizes);
    state.search(false);
    state.matches
}

// ---------------------------------------------------------------------------
// Ring helpers
// ---------------------------------------------------------------------------

fn build_ring_membership(mol: &Molecule, rings: &[Vec<usize>]) -> Vec<bool> {
    let mut member = vec![false; mol.atom_count()];
    for ring in rings {
        for &idx in ring {
            member[idx] = true;
        }
    }
    member
}

fn build_ring_sizes(mol: &Molecule, rings: &[Vec<usize>]) -> Vec<Vec<u8>> {
    let mut sizes = vec![Vec::new(); mol.atom_count()];
    for ring in rings {
        let sz = ring.len() as u8;
        for &idx in ring {
            if !sizes[idx].contains(&sz) {
                sizes[idx].push(sz);
            }
        }
    }
    sizes
}

fn build_ring_bond_set(mol: &Molecule, rings: &[Vec<usize>]) -> std::collections::HashSet<usize> {
    let mut set = std::collections::HashSet::new();
    for ring in rings {
        for i in 0..ring.len() {
            let a1 = ring[i];
            let a2 = ring[(i + 1) % ring.len()];
            for &(neighbor, bond_idx) in &mol.adjacency[a1] {
                if neighbor == a2 {
                    set.insert(bond_idx);
                }
            }
        }
    }
    set
}

// ---------------------------------------------------------------------------
// Atom/bond expression evaluation
// ---------------------------------------------------------------------------

fn eval_atom_expr(
    expr: &AtomExpr,
    mol: &Molecule,
    atom_idx: usize,
    ring_membership: &[bool],
    ring_sizes: &[Vec<u8>],
    rings: &[Vec<usize>],
) -> bool {
    match expr {
        AtomExpr::Prim(prim) => eval_atom_prim(prim, mol, atom_idx, ring_membership, ring_sizes),
        AtomExpr::And(terms) => terms.iter().all(|t| eval_atom_expr(t, mol, atom_idx, ring_membership, ring_sizes, rings)),
        AtomExpr::Or(terms) => terms.iter().any(|t| eval_atom_expr(t, mol, atom_idx, ring_membership, ring_sizes, rings)),
        AtomExpr::Not(inner) => !eval_atom_expr(inner, mol, atom_idx, ring_membership, ring_sizes, rings),
        AtomExpr::Recursive(sub_pattern) => {
            // Check if any match of sub_pattern includes this atom as its first atom
            let sub_matches = smarts_find_all(mol, sub_pattern);
            sub_matches.iter().any(|m| {
                m.atom_mapping.first().map(|&(_, t)| t) == Some(atom_idx)
            })
        }
    }
}

fn eval_atom_prim(
    prim: &AtomPrimitive,
    mol: &Molecule,
    atom_idx: usize,
    ring_membership: &[bool],
    ring_sizes: &[Vec<u8>],
) -> bool {
    let atom = &mol.atoms[atom_idx];
    match prim {
        AtomPrimitive::AtomicNum(n) => atom.atomic_number == *n,
        AtomPrimitive::Aromatic => atom.is_aromatic,
        AtomPrimitive::Aliphatic => !atom.is_aromatic,
        AtomPrimitive::Degree(d) => mol.degree(atom_idx) == *d as usize,
        AtomPrimitive::HCount(h) => atom.implicit_hydrogens == *h,
        AtomPrimitive::TotalHCount(h) => {
            let explicit_h = mol.adjacency[atom_idx]
                .iter()
                .filter(|&&(n, _)| mol.atoms[n].atomic_number == 1)
                .count() as u8;
            atom.implicit_hydrogens + explicit_h == *h
        }
        AtomPrimitive::Charge(c) => atom.formal_charge == *c,
        AtomPrimitive::RingMember => ring_membership[atom_idx],
        AtomPrimitive::RingSize(s) => ring_sizes[atom_idx].contains(s),
        AtomPrimitive::Connectivity(x) => {
            (mol.degree(atom_idx) + atom.implicit_hydrogens as usize) == *x as usize
        }
        AtomPrimitive::Valence(v) => {
            let bond_order_sum: f64 = mol.adjacency[atom_idx]
                .iter()
                .map(|&(_, bi)| mol.bonds[bi].order.as_f64())
                .sum();
            (bond_order_sum.round() as u8 + atom.implicit_hydrogens) == *v
        }
        AtomPrimitive::Wildcard => true,
    }
}

fn eval_bond_expr(
    expr: &BondExpr,
    mol: &Molecule,
    bond_idx: usize,
    ring_bond_set: &std::collections::HashSet<usize>,
) -> bool {
    let bond = &mol.bonds[bond_idx];
    match expr {
        BondExpr::Single => bond.order == BondOrder::Single && !bond.is_aromatic,
        BondExpr::Double => bond.order == BondOrder::Double,
        BondExpr::Triple => bond.order == BondOrder::Triple,
        BondExpr::Aromatic => bond.is_aromatic || bond.order == BondOrder::Aromatic,
        BondExpr::Ring => ring_bond_set.contains(&bond_idx),
        BondExpr::Any => true,
        BondExpr::Not(inner) => !eval_bond_expr(inner, mol, bond_idx, ring_bond_set),
        BondExpr::And(terms) => terms.iter().all(|t| eval_bond_expr(t, mol, bond_idx, ring_bond_set)),
        BondExpr::Or(terms) => terms.iter().any(|t| eval_bond_expr(t, mol, bond_idx, ring_bond_set)),
    }
}

// ---------------------------------------------------------------------------
// VF2 subgraph isomorphism for SMARTS
// ---------------------------------------------------------------------------

struct SmartsVf2State<'a> {
    target: &'a Molecule,
    pattern: &'a SmartsPattern,
    rings: &'a [Vec<usize>],
    ring_membership: &'a [bool],
    ring_sizes: &'a [Vec<u8>],
    ring_bond_set: std::collections::HashSet<usize>,
    core_target: Vec<Option<usize>>,
    core_pattern: Vec<Option<usize>>,
    matches: Vec<SubstructureMatch>,
}

impl<'a> SmartsVf2State<'a> {
    fn new(
        target: &'a Molecule,
        pattern: &'a SmartsPattern,
        rings: &'a [Vec<usize>],
        ring_membership: &'a [bool],
        ring_sizes: &'a [Vec<u8>],
    ) -> Self {
        let ring_bond_set = build_ring_bond_set(target, rings);
        SmartsVf2State {
            target,
            pattern,
            rings,
            ring_membership,
            ring_sizes,
            ring_bond_set,
            core_target: vec![None; target.atom_count()],
            core_pattern: vec![None; pattern.atoms.len()],
            matches: Vec::new(),
        }
    }

    fn search(&mut self, early_exit: bool) {
        if self.pattern.atoms.is_empty() {
            return;
        }
        if self.pattern.atoms.len() > self.target.atom_count() {
            return;
        }
        self.match_recursive(0, early_exit);
    }

    fn match_recursive(&mut self, depth: usize, early_exit: bool) {
        if early_exit && !self.matches.is_empty() {
            return;
        }

        if depth == self.pattern.atoms.len() {
            let mapping: Vec<(usize, usize)> = self
                .core_pattern
                .iter()
                .enumerate()
                .map(|(p, t)| (p, t.unwrap()))
                .collect();
            self.matches.push(SubstructureMatch { atom_mapping: mapping });
            return;
        }

        let pattern_atom = depth;
        let candidates = self.find_candidates(pattern_atom);

        for target_atom in candidates {
            if self.core_target[target_atom].is_some() {
                continue;
            }

            if self.is_feasible(pattern_atom, target_atom) {
                self.core_pattern[pattern_atom] = Some(target_atom);
                self.core_target[target_atom] = Some(pattern_atom);

                self.match_recursive(depth + 1, early_exit);

                self.core_pattern[pattern_atom] = None;
                self.core_target[target_atom] = None;

                if early_exit && !self.matches.is_empty() {
                    return;
                }
            }
        }
    }

    fn find_candidates(&self, pattern_atom: usize) -> Vec<usize> {
        let mut candidates: Option<Vec<usize>> = None;

        for &(p_neighbor, _) in &self.pattern.adjacency[pattern_atom] {
            if let Some(t_mapped) = self.core_pattern[p_neighbor] {
                let t_neighbors: Vec<usize> = self
                    .target
                    .adjacency[t_mapped]
                    .iter()
                    .map(|&(n, _)| n)
                    .filter(|&n| self.core_target[n].is_none())
                    .collect();

                candidates = Some(match candidates {
                    None => t_neighbors,
                    Some(existing) => existing
                        .into_iter()
                        .filter(|n| t_neighbors.contains(n))
                        .collect(),
                });
            }
        }

        match candidates {
            Some(c) => c,
            None => (0..self.target.atom_count())
                .filter(|&i| self.core_target[i].is_none())
                .collect(),
        }
    }

    fn is_feasible(&self, pattern_atom: usize, target_atom: usize) -> bool {
        // Check atom expression
        if !eval_atom_expr(
            &self.pattern.atoms[pattern_atom].expr,
            self.target,
            target_atom,
            self.ring_membership,
            self.ring_sizes,
            self.rings,
        ) {
            return false;
        }

        // Check bonds to already-mapped neighbors
        for &(p_neighbor, p_bond_idx) in &self.pattern.adjacency[pattern_atom] {
            if let Some(t_mapped) = self.core_pattern[p_neighbor] {
                // Find the bond in target between target_atom and t_mapped
                let t_bond_idx = self.target.adjacency[target_atom]
                    .iter()
                    .find(|&&(n, _)| n == t_mapped)
                    .map(|&(_, bi)| bi);

                match t_bond_idx {
                    None => return false,
                    Some(tbi) => {
                        if !eval_bond_expr(
                            &self.pattern.bonds[p_bond_idx].expr,
                            self.target,
                            tbi,
                            &self.ring_bond_set,
                        ) {
                            return false;
                        }
                    }
                }
            }
        }

        true
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::smiles::parse_smiles;

    #[test]
    fn atomic_num_matches_carbon() {
        let pattern = parse_smarts("[#6]").unwrap();
        let mol = parse_smiles("C").unwrap();
        assert!(smarts_match(&mol, &pattern));
    }

    #[test]
    fn atomic_num_excludes_nitrogen() {
        let pattern = parse_smarts("[#6]").unwrap();
        let mol = parse_smiles("N").unwrap();
        assert!(!smarts_match(&mol, &pattern));
    }

    #[test]
    fn not_carbon() {
        let pattern = parse_smarts("[!#6]").unwrap();
        let mol_n = parse_smiles("N").unwrap();
        assert!(smarts_match(&mol_n, &pattern));
        let mol_c = parse_smiles("C").unwrap();
        // C has only one atom (#6), so [!#6] should not match
        // But C also has implicit H — implicit H are not graph atoms though
        assert!(!smarts_match(&mol_c, &pattern));
    }

    #[test]
    fn aromatic_matches() {
        let pattern = parse_smarts("[a]").unwrap();
        let benzene = parse_smiles("c1ccccc1").unwrap();
        assert!(smarts_match(&benzene, &pattern));
        let ethane = parse_smiles("CC").unwrap();
        assert!(!smarts_match(&ethane, &pattern));
    }

    #[test]
    fn degree_matches() {
        let pattern = parse_smarts("[D3]").unwrap();
        // Isobutane center carbon has degree 3
        let isobutane = parse_smiles("CC(C)C").unwrap();
        assert!(smarts_match(&isobutane, &pattern));
    }

    #[test]
    fn ring_member_matches() {
        let pattern = parse_smarts("[R]").unwrap();
        let benzene = parse_smiles("c1ccccc1").unwrap();
        assert!(smarts_match(&benzene, &pattern));
        let ethane = parse_smiles("CC").unwrap();
        assert!(!smarts_match(&ethane, &pattern));
    }

    #[test]
    fn ring_size_matches() {
        let pattern = parse_smarts("[r6]").unwrap();
        let benzene = parse_smiles("c1ccccc1").unwrap();
        assert!(smarts_match(&benzene, &pattern));
        let cyclopentane = parse_smiles("C1CCCC1").unwrap();
        assert!(!smarts_match(&cyclopentane, &pattern));
    }

    #[test]
    fn recursive_smarts_hydroxyl() {
        let pattern = parse_smarts("[$([OH])]").unwrap();
        let phenol = parse_smiles("Oc1ccccc1").unwrap();
        assert!(smarts_match(&phenol, &pattern));
    }

    #[test]
    fn or_combination() {
        // Match nitrogen or oxygen
        let pattern = parse_smarts("[#7,#8]").unwrap();
        let mol_n = parse_smiles("N").unwrap();
        assert!(smarts_match(&mol_n, &pattern));
        let mol_o = parse_smiles("O").unwrap();
        assert!(smarts_match(&mol_o, &pattern));
        let mol_c = parse_smiles("C").unwrap();
        assert!(!smarts_match(&mol_c, &pattern));
    }

    #[test]
    fn benzene_pattern_matches() {
        let pattern = parse_smarts("c1ccccc1").unwrap();
        let benzene = parse_smiles("c1ccccc1").unwrap();
        assert!(smarts_match(&benzene, &pattern));
    }

    #[test]
    fn wildcard_matches_any() {
        let pattern = parse_smarts("*").unwrap();
        let mol = parse_smiles("C").unwrap();
        assert!(smarts_match(&mol, &pattern));
        let mol2 = parse_smiles("N").unwrap();
        assert!(smarts_match(&mol2, &pattern));
    }

    #[test]
    fn invalid_smarts_error() {
        assert!(parse_smarts("").is_err());
        assert!(parse_smarts("[").is_err());
    }

    #[test]
    fn find_all_returns_mappings() {
        let pattern = parse_smarts("[#6]").unwrap();
        let mol = parse_smiles("CCO").unwrap();
        let matches = smarts_find_all(&mol, &pattern);
        assert_eq!(matches.len(), 2); // Two carbons
    }

    #[test]
    fn chain_pattern() {
        let pattern = parse_smarts("[#6][#8]").unwrap();
        let mol = parse_smiles("CCO").unwrap();
        assert!(smarts_match(&mol, &pattern));
        let mol2 = parse_smiles("CCC").unwrap();
        assert!(!smarts_match(&mol2, &pattern));
    }
}
