//! Secondary structure assignment using DSSP-like algorithms.
//!
//! Provides both a simplified distance-based DSSP (the original
//! [`assign_secondary_structure`]) and a full hydrogen-bond energy-based DSSP
//! ([`dssp`]) that classifies residues into the standard 8-state scheme.

use cyanea_core::{CyaneaError, Result, Summarizable};

use crate::geometry::{angle_points, dihedral_points};
use crate::types::{Chain, Point3D, Residue};

#[allow(unused_imports)]
use alloc::format;
use alloc::string::String;
use alloc::vec;
use alloc::vec::Vec;

/// Secondary structure classification for a single residue (simplified 4-state).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum SecondaryStructure {
    Helix,
    Sheet,
    Turn,
    Coil,
}

impl SecondaryStructure {
    /// Single-character DSSP code.
    pub fn code(&self) -> char {
        match self {
            SecondaryStructure::Helix => 'H',
            SecondaryStructure::Sheet => 'E',
            SecondaryStructure::Turn => 'T',
            SecondaryStructure::Coil => 'C',
        }
    }
}

/// Full 8-state DSSP classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum DsspState {
    /// Alpha-helix (i -> i+4 H-bond pattern, >= 4 consecutive).
    H,
    /// 3_10-helix (i -> i+3 H-bond pattern).
    G,
    /// Pi-helix (i -> i+5 H-bond pattern).
    I,
    /// Extended strand in beta-sheet.
    E,
    /// Isolated beta-bridge residue.
    B,
    /// Hydrogen-bonded turn.
    T,
    /// Bend (CA angle > 70 degrees).
    S,
    /// Coil / loop (none of the above).
    C,
}

impl DsspState {
    /// Single-character DSSP code.
    pub fn code(&self) -> char {
        match self {
            DsspState::H => 'H',
            DsspState::G => 'G',
            DsspState::I => 'I',
            DsspState::E => 'E',
            DsspState::B => 'B',
            DsspState::T => 'T',
            DsspState::S => 'S',
            DsspState::C => 'C',
        }
    }

    /// Convert to the simplified 4-state classification.
    pub fn to_simplified(&self) -> SecondaryStructure {
        match self {
            DsspState::H | DsspState::G | DsspState::I => SecondaryStructure::Helix,
            DsspState::E | DsspState::B => SecondaryStructure::Sheet,
            DsspState::T | DsspState::S => SecondaryStructure::Turn,
            DsspState::C => SecondaryStructure::Coil,
        }
    }
}

/// Full DSSP assignment result for a chain.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DsspAssignment {
    /// Chain identifier.
    pub chain_id: char,
    /// One 8-state assignment per residue.
    pub states: Vec<DsspState>,
}

impl DsspAssignment {
    /// Convert to a string of single-character DSSP codes.
    pub fn to_string_codes(&self) -> String {
        self.states.iter().map(|s| s.code()).collect()
    }

    /// Convert to the simplified 4-state assignment.
    pub fn to_simplified(&self) -> SecondaryStructureAssignment {
        SecondaryStructureAssignment {
            chain_id: self.chain_id,
            assignments: self.states.iter().map(|s| s.to_simplified()).collect(),
        }
    }

    /// Count of each DSSP state.
    pub fn counts(&self) -> DsspCounts {
        let mut counts = DsspCounts::default();
        for s in &self.states {
            match s {
                DsspState::H => counts.h += 1,
                DsspState::G => counts.g += 1,
                DsspState::I => counts.i += 1,
                DsspState::E => counts.e += 1,
                DsspState::B => counts.b += 1,
                DsspState::T => counts.t += 1,
                DsspState::S => counts.s += 1,
                DsspState::C => counts.c += 1,
            }
        }
        counts
    }
}

impl Summarizable for DsspAssignment {
    fn summary(&self) -> String {
        let c = self.counts();
        format!(
            "Chain {} DSSP: {} residue(s) — H:{} G:{} I:{} E:{} B:{} T:{} S:{} C:{}",
            self.chain_id,
            self.states.len(),
            c.h, c.g, c.i, c.e, c.b, c.t, c.s, c.c,
        )
    }
}

/// Counts of each DSSP state.
#[derive(Debug, Clone, Default)]
pub struct DsspCounts {
    pub h: usize,
    pub g: usize,
    pub i: usize,
    pub e: usize,
    pub b: usize,
    pub t: usize,
    pub s: usize,
    pub c: usize,
}

/// Secondary structure assignment for a chain.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SecondaryStructureAssignment {
    /// Chain identifier.
    pub chain_id: char,
    /// One assignment per residue.
    pub assignments: Vec<SecondaryStructure>,
}

impl SecondaryStructureAssignment {
    /// Count of each secondary structure type.
    pub fn counts(&self) -> (usize, usize, usize, usize) {
        let mut h = 0;
        let mut e = 0;
        let mut t = 0;
        let mut c = 0;
        for ss in &self.assignments {
            match ss {
                SecondaryStructure::Helix => h += 1,
                SecondaryStructure::Sheet => e += 1,
                SecondaryStructure::Turn => t += 1,
                SecondaryStructure::Coil => c += 1,
            }
        }
        (h, e, t, c)
    }

    /// Fraction of residues in helix.
    pub fn helix_fraction(&self) -> f64 {
        if self.assignments.is_empty() {
            return 0.0;
        }
        let (h, _, _, _) = self.counts();
        h as f64 / self.assignments.len() as f64
    }

    /// Fraction of residues in sheet.
    pub fn sheet_fraction(&self) -> f64 {
        if self.assignments.is_empty() {
            return 0.0;
        }
        let (_, e, _, _) = self.counts();
        e as f64 / self.assignments.len() as f64
    }
}

impl Summarizable for SecondaryStructureAssignment {
    fn summary(&self) -> String {
        let (h, e, t, c) = self.counts();
        format!(
            "Chain {} SS: {} residue(s) — H:{} E:{} T:{} C:{}",
            self.chain_id,
            self.assignments.len(),
            h,
            e,
            t,
            c,
        )
    }
}

/// Assign secondary structure to each residue in a chain.
///
/// Uses a simplified DSSP-like algorithm:
/// - **Helix**: backbone O(i) to N(i+4) distance < 3.5 Å for at least 3 consecutive residues
/// - **Turn**: backbone O(i) to N(i+3) distance < 3.5 Å
/// - **Sheet**: detected via antiparallel beta-bridge patterns (O(i)→N(j) and O(j)→N(i) both < 3.5 Å with |i-j| > 4)
/// - **Coil**: everything else
pub fn assign_secondary_structure(chain: &Chain) -> Result<SecondaryStructureAssignment> {
    let n = chain.residues.len();
    if n == 0 {
        return Err(CyaneaError::InvalidInput(
            "cannot assign SS to empty chain".into(),
        ));
    }

    let mut assignments = vec![SecondaryStructure::Coil; n];

    // Precompute backbone atom positions (N, O) for each residue
    let backbone: Vec<Option<(Point3D, Point3D)>> = chain
        .residues
        .iter()
        .map(|r| {
            let n_atom = r.get_atom("N")?;
            let o_atom = r.get_atom("O")?;
            Some((n_atom.coords, o_atom.coords))
        })
        .collect();

    let hbond_cutoff = 3.5;

    // Check i→i+4 H-bonds (helix pattern)
    let mut helix_hbond = vec![false; n];
    for i in 0..n.saturating_sub(4) {
        if let (Some((_, o_i)), Some((n_i4, _))) = (&backbone[i], &backbone[i + 4]) {
            if o_i.distance_to(n_i4) < hbond_cutoff {
                helix_hbond[i] = true;
            }
        }
    }

    // Helix: 3+ consecutive i→i+4 H-bonds
    for i in 0..n.saturating_sub(6) {
        if helix_hbond[i] && helix_hbond[i + 1] && helix_hbond[i + 2] {
            // Mark residues i through i+6 (the full helix span) as Helix
            for j in i..=(i + 5).min(n - 1) {
                assignments[j] = SecondaryStructure::Helix;
            }
        }
    }

    // Check antiparallel beta-bridges (sheet pattern)
    for i in 0..n {
        for j in (i + 5)..n {
            if let (Some((n_i, o_i)), Some((n_j, o_j))) = (&backbone[i], &backbone[j]) {
                let oi_nj = o_i.distance_to(n_j);
                let oj_ni = o_j.distance_to(n_i);
                if oi_nj < hbond_cutoff && oj_ni < hbond_cutoff {
                    if assignments[i] == SecondaryStructure::Coil {
                        assignments[i] = SecondaryStructure::Sheet;
                    }
                    if assignments[j] == SecondaryStructure::Coil {
                        assignments[j] = SecondaryStructure::Sheet;
                    }
                }
            }
        }
    }

    // Check i→i+3 (turn pattern)
    for i in 0..n.saturating_sub(3) {
        if let (Some((_, o_i)), Some((n_i3, _))) = (&backbone[i], &backbone[i + 3]) {
            if o_i.distance_to(n_i3) < hbond_cutoff {
                // Only assign turn if not already helix/sheet
                for j in i..=(i + 3).min(n - 1) {
                    if assignments[j] == SecondaryStructure::Coil {
                        assignments[j] = SecondaryStructure::Turn;
                    }
                }
            }
        }
    }

    Ok(SecondaryStructureAssignment {
        chain_id: chain.id,
        assignments,
    })
}

/// Full DSSP secondary structure assignment using hydrogen-bond energies.
///
/// Implements the Kabsch & Sander (1983) algorithm:
///
/// 1. Identify backbone atoms (N, CA, C, O) for each residue
/// 2. Calculate NH->CO hydrogen bond energies using the electrostatic model:
///    `E = 0.084 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN) * 332` kcal/mol
///    A bond exists if E < -0.5 kcal/mol
/// 3. Classify into 8 states based on H-bond patterns:
///    - H: alpha-helix (i->i+4 pattern, >= 4 consecutive)
///    - G: 3_10-helix (i->i+3 pattern)
///    - I: pi-helix (i->i+5 pattern)
///    - E: extended strand (parallel or antiparallel H-bonds to distant residues)
///    - B: isolated beta-bridge
///    - T: hydrogen-bonded turn
///    - S: bend (CA angle > 70 degrees)
///    - C: coil (everything else)
///
/// # Errors
///
/// Returns an error if the chain is empty.
pub fn dssp(chain: &Chain) -> Result<DsspAssignment> {
    let n = chain.residues.len();
    if n == 0 {
        return Err(CyaneaError::InvalidInput(
            "cannot assign DSSP to empty chain".into(),
        ));
    }

    // Step 1: Extract backbone atom positions
    let backbone = extract_backbone(chain);

    // Step 2: Compute hydrogen bond energy matrix
    // hbond_energy[i][j] stores the energy of the NH(i) -> CO(j) hydrogen bond
    let hbond = compute_hbond_matrix(&backbone, n);

    let mut states = vec![DsspState::C; n];

    // Step 3a: Identify helix patterns (must be done before bridges)
    // Check n-turns of size 3, 4, 5
    let turn3 = find_turns(&hbond, n, 3);
    let turn4 = find_turns(&hbond, n, 4);
    let turn5 = find_turns(&hbond, n, 5);

    // Alpha-helix (H): 4+ consecutive i->i+4 H-bonds
    assign_helix(&turn4, n, 4, DsspState::H, &mut states);
    // 3_10-helix (G): 3+ consecutive i->i+3 H-bonds
    assign_helix(&turn3, n, 3, DsspState::G, &mut states);
    // Pi-helix (I): 5+ consecutive i->i+5 H-bonds
    assign_helix(&turn5, n, 5, DsspState::I, &mut states);

    // Step 3b: Identify beta-bridges and strands
    let bridges = find_bridges(&hbond, n);
    assign_bridges_and_strands(&bridges, n, &mut states);

    // Step 3c: Assign turns (only to residues not already assigned helix/strand)
    assign_turns(&turn3, &turn4, &turn5, n, &mut states);

    // Step 3d: Assign bends
    assign_bends(chain, &mut states);

    Ok(DsspAssignment {
        chain_id: chain.id,
        states,
    })
}

/// Backbone atom coordinates for a single residue.
#[allow(dead_code)]
struct BackboneAtoms {
    n: Option<Point3D>,
    ca: Option<Point3D>,
    c: Option<Point3D>,
    o: Option<Point3D>,
}

/// Extract backbone atom positions for all residues.
fn extract_backbone(chain: &Chain) -> Vec<BackboneAtoms> {
    chain
        .residues
        .iter()
        .map(|r| BackboneAtoms {
            n: r.get_atom("N").map(|a| a.coords),
            ca: r.get_atom("CA").map(|a| a.coords),
            c: r.get_atom("C").map(|a| a.coords),
            o: r.get_atom("O").map(|a| a.coords),
        })
        .collect()
}

/// Compute the NH(i) -> CO(j) hydrogen bond energy in kcal/mol.
///
/// Uses the Kabsch-Sander electrostatic model:
/// E = 0.084 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN) * 332
///
/// The H atom position is estimated as being along the N-C_prev bond direction,
/// 1.0 A from N.
fn hbond_energy(donor: &BackboneAtoms, acceptor: &BackboneAtoms, prev_c: Option<&Point3D>) -> f64 {
    // Need: donor N and H (estimated), acceptor C and O
    let n_pos = match donor.n {
        Some(p) => p,
        None => return 0.0,
    };
    let o_pos = match acceptor.o {
        Some(p) => p,
        None => return 0.0,
    };
    let c_pos = match acceptor.c {
        Some(p) => p,
        None => return 0.0,
    };

    // Estimate H position: H is on the N, opposite to the C(i-1)->N direction
    // H = N + (N - C_prev).normalize() * 1.0
    let h_pos = match prev_c {
        Some(c_prev) => {
            let nc_dir = n_pos.sub(c_prev).normalize();
            n_pos.add(&nc_dir.scale(1.0))
        }
        None => return 0.0, // Can't estimate H for first residue
    };

    let r_on = n_pos.distance_to(&o_pos);
    let r_ch = h_pos.distance_to(&c_pos);
    let r_oh = h_pos.distance_to(&o_pos);
    let r_cn = n_pos.distance_to(&c_pos);

    // Avoid division by zero
    if r_on < 0.5 || r_ch < 0.5 || r_oh < 0.5 || r_cn < 0.5 {
        return 0.0;
    }

    0.084 * (1.0 / r_on + 1.0 / r_ch - 1.0 / r_oh - 1.0 / r_cn) * 332.0
}

/// Compute the hydrogen bond energy matrix.
/// Returns a flat Vec where hbond[i * n + j] = energy of NH(i) -> CO(j).
/// Only stores energies below the threshold.
fn compute_hbond_matrix(backbone: &[BackboneAtoms], n: usize) -> Vec<f64> {
    let mut hbond = vec![0.0_f64; n * n];
    let hbond_threshold = -0.5;

    for i in 1..n {
        // Donor is residue i (NH), need C of residue i-1 for H estimation
        let prev_c = backbone[i - 1].c.as_ref();

        for j in 0..n {
            if i == j {
                continue;
            }
            // Skip if donor and acceptor are too close in sequence for meaningful bond
            // (i, i+1) and (i, i-1) are peptide bonds, not H-bonds
            if (i as isize - j as isize).unsigned_abs() < 2 {
                continue;
            }

            let energy = hbond_energy(&backbone[i], &backbone[j], prev_c);
            if energy < hbond_threshold {
                hbond[i * n + j] = energy;
            }
        }
    }

    hbond
}

/// Find n-turns: residue i has a turn of size `turn_size` if NH(i+turn_size) -> CO(i) H-bond exists.
fn find_turns(hbond: &[f64], n: usize, turn_size: usize) -> Vec<bool> {
    let mut turns = vec![false; n];
    for i in 0..n.saturating_sub(turn_size) {
        let donor = i + turn_size;
        if donor < n && hbond[donor * n + i] < -0.5 {
            turns[i] = true;
        }
    }
    turns
}

/// Assign helix states based on consecutive turns.
fn assign_helix(
    turns: &[bool],
    n: usize,
    turn_size: usize,
    state: DsspState,
    states: &mut [DsspState],
) {
    // A helix requires at least 2 consecutive turns of the same type
    // (meaning the pattern i->i+k and (i+1)->(i+1+k) both exist)
    // The helix then spans from i+1 to i+turn_size for each such pair.
    let min_consecutive = 2;
    let mut consecutive = 0;

    for i in 0..n.saturating_sub(turn_size) {
        if turns[i] {
            consecutive += 1;
            if consecutive >= min_consecutive {
                // Mark the helix residues for this span
                // The helix starts at i - (consecutive - min_consecutive) + 1
                // and extends to i + turn_size
                let start = (i + 2).saturating_sub(consecutive);
                let end = (i + turn_size).min(n - 1);
                for j in start..=end {
                    // Only assign if current state is Coil or same type
                    // H takes precedence over G, G over I
                    if states[j] == DsspState::C
                        || (state == DsspState::H && (states[j] == DsspState::G || states[j] == DsspState::I))
                        || (state == DsspState::G && states[j] == DsspState::I)
                    {
                        states[j] = state;
                    }
                }
            }
        } else {
            consecutive = 0;
        }
    }
}

/// A beta-bridge between two residues.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
struct BetaBridge {
    i: usize,
    j: usize,
    is_parallel: bool,
}

/// Find all beta-bridges in the chain.
fn find_bridges(hbond: &[f64], n: usize) -> Vec<BetaBridge> {
    let mut bridges = Vec::new();
    let threshold = -0.5;

    for i in 1..n.saturating_sub(1) {
        for j in (i + 3)..n.saturating_sub(1) {
            // Antiparallel bridge: NH(i)->CO(j) and NH(j)->CO(i)
            let anti = hbond[i * n + j] < threshold && hbond[j * n + i] < threshold;
            // Alternative antiparallel: NH(i)->CO(j-1) and NH(j+1)->CO(i)
            // or NH(i+1)->CO(j) and NH(j)->CO(i-1)
            let anti_alt1 = j > 0
                && i + 1 < n
                && hbond[i * n + (j - 1)] < threshold
                && (j + 1) < n
                && hbond[(j + 1) * n + i] < threshold;
            let anti_alt2 = i > 0
                && (i + 1) < n
                && hbond[(i + 1) * n + j] < threshold
                && hbond[j * n + (i - 1)] < threshold;

            if anti || anti_alt1 || anti_alt2 {
                bridges.push(BetaBridge {
                    i,
                    j,
                    is_parallel: false,
                });
                continue;
            }

            // Parallel bridge: NH(i)->CO(j-1) and NH(j)->CO(i)  (but different pattern)
            // Pattern 1: NH(i-1)->CO(j) and NH(j)->CO(i)
            let para1 = i > 0
                && hbond[(i - 1) * n + j] < threshold
                && hbond[j * n + i] < threshold;
            // Pattern 2: NH(i)->CO(j) and NH(j+1)->CO(i)
            let para2 = j + 1 < n
                && hbond[i * n + j] < threshold
                && hbond[(j + 1) * n + i] < threshold;

            if para1 || para2 {
                bridges.push(BetaBridge {
                    i,
                    j,
                    is_parallel: true,
                });
            }
        }
    }

    bridges
}

/// Assign E (extended strand) and B (isolated bridge) states.
fn assign_bridges_and_strands(bridges: &[BetaBridge], _n: usize, states: &mut [DsspState]) {
    // First pass: mark all bridge residues
    for bridge in bridges {
        // Only override coil
        if states[bridge.i] == DsspState::C {
            states[bridge.i] = DsspState::B;
        }
        if states[bridge.j] == DsspState::C {
            states[bridge.j] = DsspState::B;
        }
    }

    // Second pass: if two adjacent residues both have bridges, they form a strand (E)
    // A strand = consecutive bridge residues
    // Check for consecutive bridge residues
    let bridge_residues: Vec<bool> = states.iter().map(|s| *s == DsspState::B).collect();
    let n = bridge_residues.len();

    for i in 0..n {
        if !bridge_residues[i] {
            continue;
        }
        // Check if this residue has an adjacent bridge residue
        let has_neighbor = (i > 0 && bridge_residues[i - 1])
            || (i + 1 < n && bridge_residues[i + 1]);
        if has_neighbor {
            states[i] = DsspState::E;
        }
    }

    // Upgrade remaining B neighbors of E to E
    let strand_residues: Vec<bool> = states.iter().map(|s| *s == DsspState::E).collect();
    for i in 0..n {
        if states[i] == DsspState::B {
            let has_e_neighbor = (i > 0 && strand_residues[i - 1])
                || (i + 1 < n && strand_residues[i + 1]);
            if has_e_neighbor {
                states[i] = DsspState::E;
            }
        }
    }
}

/// Assign T (turn) to residues involved in H-bonded turns but not helices/strands.
fn assign_turns(
    turn3: &[bool],
    turn4: &[bool],
    turn5: &[bool],
    n: usize,
    states: &mut [DsspState],
) {
    for i in 0..n {
        // Mark residues i+1 to i+turn_size-1 as Turn if there's a turn at i
        for (turns, turn_size) in [(turn3, 3usize), (turn4, 4usize), (turn5, 5usize)] {
            if i < turns.len() && turns[i] {
                let start = i + 1;
                let end = (i + turn_size).min(n);
                for j in start..end {
                    if states[j] == DsspState::C {
                        states[j] = DsspState::T;
                    }
                }
            }
        }
    }
}

/// Assign S (bend) where the angle between direction vectors
/// CA(i-2)->CA(i) and CA(i)->CA(i+2) exceeds 70 degrees.
///
/// `angle_points(p1, p2, p3)` returns the vertex angle at p2, which is
/// `180 - angle_between_direction_vectors`. So a bend exists when the
/// vertex angle < 110 degrees.
fn assign_bends(chain: &Chain, states: &mut [DsspState]) {
    let n = chain.residues.len();
    if n < 5 {
        return;
    }

    for i in 2..n.saturating_sub(2) {
        if states[i] != DsspState::C {
            continue;
        }

        let ca_prev2 = chain.residues[i - 2].get_atom("CA");
        let ca_curr = chain.residues[i].get_atom("CA");
        let ca_next2 = chain.residues[i + 2].get_atom("CA");

        if let (Some(p1), Some(p2), Some(p3)) = (ca_prev2, ca_curr, ca_next2) {
            // Vertex angle at CA(i) = 180 - angle_between_direction_vectors
            // DSSP bend: direction angle > 70° means vertex angle < 110°
            let vertex_angle = angle_points(&p1.coords, &p2.coords, &p3.coords);
            if vertex_angle < 110.0 {
                states[i] = DsspState::S;
            }
        }
    }
}

/// Compute backbone dihedral angles (phi, psi) for a residue given its neighbors.
///
/// - `prev`: previous residue (needed for phi: C(i-1)–N(i)–CA(i)–C(i))
/// - `curr`: current residue
/// - `next`: next residue (needed for psi: N(i)–CA(i)–C(i)–N(i+1))
///
/// Returns `Some((phi, psi))` in degrees if all required atoms are present.
pub fn backbone_dihedrals(
    prev: Option<&Residue>,
    curr: &Residue,
    next: Option<&Residue>,
) -> Option<(f64, f64)> {
    let n_i = curr.get_atom("N")?.coords;
    let ca_i = curr.get_atom("CA")?.coords;
    let c_i = curr.get_atom("C")?.coords;

    let phi = if let Some(prev_res) = prev {
        let c_prev = prev_res.get_atom("C")?.coords;
        dihedral_points(&c_prev, &n_i, &ca_i, &c_i)
    } else {
        return None;
    };

    let psi = if let Some(next_res) = next {
        let n_next = next_res.get_atom("N")?.coords;
        dihedral_points(&n_i, &ca_i, &c_i, &n_next)
    } else {
        return None;
    };

    Some((phi, psi))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{Atom, Chain, Residue};

    fn make_atom(name: &str, x: f64, y: f64, z: f64) -> Atom {
        Atom {
            serial: 1,
            name: name.into(),
            alt_loc: None,
            coords: Point3D::new(x, y, z),
            occupancy: 1.0,
            temp_factor: 0.0,
            element: None,
            charge: None,
            is_hetatm: false,
        }
    }

    fn make_helix_chain() -> Chain {
        // Build a synthetic alpha-helix with proper H-bond geometry.
        // Alpha helix: 3.6 residues/turn, rise 1.5Å/residue, radius 2.3Å.
        // Key: O(i) must be < 3.5Å from N(i+4).
        //
        // Strategy: first place all N atoms on the helix, then place O atoms
        // offset toward N(i+4) so the H-bond distance criterion is met.
        let n_residues = 12;
        let rise = 1.5_f64; // Å per residue
        let r = 2.3_f64; // helix radius
        let angle_per_res = 100.0_f64.to_radians(); // ~100° per residue

        // Compute N positions first
        let n_positions: Vec<Point3D> = (0..n_residues)
            .map(|i| {
                let angle = i as f64 * angle_per_res;
                Point3D::new(r * angle.cos(), r * angle.sin(), i as f64 * rise)
            })
            .collect();

        let mut residues = Vec::new();
        for i in 0..n_residues {
            let n_pos = n_positions[i];
            let angle = i as f64 * angle_per_res;
            let z = i as f64 * rise;

            let ca_pos = Point3D::new(
                r * (angle + 0.3).cos(),
                r * (angle + 0.3).sin(),
                z + 0.4,
            );
            let c_pos = Point3D::new(
                r * (angle + 0.6).cos(),
                r * (angle + 0.6).sin(),
                z + 0.8,
            );

            // Place O near N(i+4) so that dist(O_i, N_{i+4}) < 3.5Å.
            // We place O 2.5Å from C along the C→N(i+4) direction, giving
            // an O→N(i+4) distance of ~2.7Å (well within the 3.5Å cutoff).
            let o_pos = if i + 4 < n_residues {
                let target = n_positions[i + 4];
                let dir = target.sub(&c_pos);
                let d = dir.norm();
                if d > 1e-10 {
                    c_pos.add(&dir.scale(2.5 / d))
                } else {
                    Point3D::new(c_pos.x, c_pos.y + 1.24, c_pos.z)
                }
            } else {
                Point3D::new(c_pos.x, c_pos.y + 1.24, c_pos.z)
            };

            residues.push(Residue {
                name: "ALA".into(),
                seq_num: i as i32 + 1,
                i_code: None,
                atoms: vec![
                    make_atom("N", n_pos.x, n_pos.y, n_pos.z),
                    make_atom("CA", ca_pos.x, ca_pos.y, ca_pos.z),
                    make_atom("C", c_pos.x, c_pos.y, c_pos.z),
                    make_atom("O", o_pos.x, o_pos.y, o_pos.z),
                ],
            });
        }
        Chain::new('A', residues)
    }

    #[test]
    fn short_chain_all_coil() {
        let chain = Chain::new(
            'A',
            vec![
                Residue {
                    name: "ALA".into(),
                    seq_num: 1,
                    i_code: None,
                    atoms: vec![
                        make_atom("N", 0.0, 0.0, 0.0),
                        make_atom("CA", 1.5, 0.0, 0.0),
                        make_atom("C", 3.0, 0.0, 0.0),
                        make_atom("O", 3.5, 1.0, 0.0),
                    ],
                },
                Residue {
                    name: "GLY".into(),
                    seq_num: 2,
                    i_code: None,
                    atoms: vec![
                        make_atom("N", 4.0, 0.0, 0.0),
                        make_atom("CA", 5.5, 0.0, 0.0),
                        make_atom("C", 7.0, 0.0, 0.0),
                        make_atom("O", 7.5, 1.0, 0.0),
                    ],
                },
            ],
        );
        let ss = assign_secondary_structure(&chain).unwrap();
        assert_eq!(ss.assignments.len(), 2);
        // Too short for helix/sheet, should be coil or turn
        for a in &ss.assignments {
            assert!(
                *a == SecondaryStructure::Coil || *a == SecondaryStructure::Turn,
                "short chain should be coil/turn"
            );
        }
    }

    #[test]
    fn helix_detection() {
        let chain = make_helix_chain();
        let ss = assign_secondary_structure(&chain).unwrap();
        // At least some residues should be helix
        let (h, _, _, _) = ss.counts();
        assert!(
            h > 0,
            "helix chain should have helical residues, got: {:?}",
            ss.assignments
        );
    }

    #[test]
    fn sheet_detection() {
        // Build two strands running antiparallel with close O-N contacts.
        let mut residues = Vec::new();
        // Strand 1: residues 1-5, running along +x
        for i in 0..5 {
            let x = i as f64 * 3.5;
            residues.push(Residue {
                name: "ALA".into(),
                seq_num: i + 1,
                i_code: None,
                atoms: vec![
                    make_atom("N", x, 0.0, 0.0),
                    make_atom("CA", x + 1.0, 0.0, 0.0),
                    make_atom("C", x + 2.0, 0.0, 0.0),
                    make_atom("O", x + 2.0, 1.0, 0.0),
                ],
            });
        }
        // Spacer residues 6-10 (far away)
        for i in 5..10 {
            residues.push(Residue {
                name: "ALA".into(),
                seq_num: i + 1,
                i_code: None,
                atoms: vec![
                    make_atom("N", 50.0, 50.0, i as f64 * 10.0),
                    make_atom("CA", 51.0, 50.0, i as f64 * 10.0),
                    make_atom("C", 52.0, 50.0, i as f64 * 10.0),
                    make_atom("O", 52.0, 51.0, i as f64 * 10.0),
                ],
            });
        }
        // Strand 2: residues 11-15, antiparallel (running along -x)
        // O of strand1[i] close to N of strand2[4-i] and vice versa
        for i in 0..5 {
            let x = (4 - i) as f64 * 3.5;
            residues.push(Residue {
                name: "ALA".into(),
                seq_num: (i + 11) as i32,
                i_code: None,
                atoms: vec![
                    make_atom("N", x + 2.0, 2.5, 0.0), // Close to O of strand 1
                    make_atom("CA", x + 1.0, 2.5, 0.0),
                    make_atom("C", x, 2.5, 0.0),
                    make_atom("O", x, 1.5, 0.0), // Close to N of strand 1
                ],
            });
        }

        let chain = Chain::new('A', residues);
        let ss = assign_secondary_structure(&chain).unwrap();
        let (_, e, _, _) = ss.counts();
        assert!(e > 0, "should detect sheet residues, got: {:?}", ss.assignments);
    }

    #[test]
    fn backbone_dihedral_computation() {
        let prev = Residue {
            name: "ALA".into(),
            seq_num: 1,
            i_code: None,
            atoms: vec![
                make_atom("N", -1.0, 0.0, 0.0),
                make_atom("CA", 0.0, 0.0, 0.0),
                make_atom("C", 1.0, 0.0, 0.0),
                make_atom("O", 1.0, 1.0, 0.0),
            ],
        };
        let curr = Residue {
            name: "GLY".into(),
            seq_num: 2,
            i_code: None,
            atoms: vec![
                make_atom("N", 2.0, 0.0, 0.0),
                make_atom("CA", 3.0, 0.0, 0.0),
                make_atom("C", 4.0, 0.0, 0.0),
                make_atom("O", 4.0, 1.0, 0.0),
            ],
        };
        let next = Residue {
            name: "ALA".into(),
            seq_num: 3,
            i_code: None,
            atoms: vec![
                make_atom("N", 5.0, 0.0, 0.0),
                make_atom("CA", 6.0, 0.0, 0.0),
                make_atom("C", 7.0, 0.0, 0.0),
                make_atom("O", 7.0, 1.0, 0.0),
            ],
        };
        let result = backbone_dihedrals(Some(&prev), &curr, Some(&next));
        assert!(result.is_some(), "should compute phi/psi");
    }

    // ---- Full DSSP tests ----

    #[test]
    fn dssp_empty_chain_error() {
        let chain = Chain::new('A', vec![]);
        assert!(dssp(&chain).is_err());
    }

    #[test]
    fn dssp_short_chain_coil() {
        let chain = Chain::new(
            'A',
            vec![
                Residue {
                    name: "ALA".into(),
                    seq_num: 1,
                    i_code: None,
                    atoms: vec![
                        make_atom("N", 0.0, 0.0, 0.0),
                        make_atom("CA", 1.5, 0.0, 0.0),
                        make_atom("C", 3.0, 0.0, 0.0),
                        make_atom("O", 3.5, 1.0, 0.0),
                    ],
                },
                Residue {
                    name: "GLY".into(),
                    seq_num: 2,
                    i_code: None,
                    atoms: vec![
                        make_atom("N", 4.0, 0.0, 0.0),
                        make_atom("CA", 5.5, 0.0, 0.0),
                        make_atom("C", 7.0, 0.0, 0.0),
                        make_atom("O", 7.5, 1.0, 0.0),
                    ],
                },
            ],
        );
        let result = dssp(&chain).unwrap();
        assert_eq!(result.states.len(), 2);
        // Short chain without H-bonds should not be helix or strand
        for s in &result.states {
            assert!(
                *s != DsspState::H && *s != DsspState::E,
                "short chain should not be helix/strand"
            );
        }
    }

    #[test]
    fn dssp_helix_chain() {
        let chain = make_helix_chain();
        let result = dssp(&chain).unwrap();
        assert_eq!(result.states.len(), 12);
        // The DSSP should detect some helical or turn structure
        // (energy-based DSSP on synthetic helix may classify as H, G, or T)
        let has_structure = result.states.iter().any(|s| *s != DsspState::C);
        assert!(
            has_structure,
            "helix chain should have some structure assigned, got: {}",
            result.to_string_codes()
        );
    }

    #[test]
    fn dssp_state_codes() {
        assert_eq!(DsspState::H.code(), 'H');
        assert_eq!(DsspState::G.code(), 'G');
        assert_eq!(DsspState::I.code(), 'I');
        assert_eq!(DsspState::E.code(), 'E');
        assert_eq!(DsspState::B.code(), 'B');
        assert_eq!(DsspState::T.code(), 'T');
        assert_eq!(DsspState::S.code(), 'S');
        assert_eq!(DsspState::C.code(), 'C');
    }

    #[test]
    fn dssp_to_simplified_mapping() {
        assert_eq!(DsspState::H.to_simplified(), SecondaryStructure::Helix);
        assert_eq!(DsspState::G.to_simplified(), SecondaryStructure::Helix);
        assert_eq!(DsspState::I.to_simplified(), SecondaryStructure::Helix);
        assert_eq!(DsspState::E.to_simplified(), SecondaryStructure::Sheet);
        assert_eq!(DsspState::B.to_simplified(), SecondaryStructure::Sheet);
        assert_eq!(DsspState::T.to_simplified(), SecondaryStructure::Turn);
        assert_eq!(DsspState::S.to_simplified(), SecondaryStructure::Turn);
        assert_eq!(DsspState::C.to_simplified(), SecondaryStructure::Coil);
    }

    #[test]
    fn dssp_assignment_to_simplified() {
        let chain = make_helix_chain();
        let dssp_result = dssp(&chain).unwrap();
        let simplified = dssp_result.to_simplified();
        assert_eq!(simplified.chain_id, 'A');
        assert_eq!(simplified.assignments.len(), 12);
    }

    #[test]
    fn dssp_string_codes() {
        let assignment = DsspAssignment {
            chain_id: 'A',
            states: vec![DsspState::H, DsspState::H, DsspState::C, DsspState::E],
        };
        assert_eq!(assignment.to_string_codes(), "HHCE");
    }

    #[test]
    fn dssp_counts() {
        let assignment = DsspAssignment {
            chain_id: 'A',
            states: vec![
                DsspState::H,
                DsspState::H,
                DsspState::G,
                DsspState::E,
                DsspState::B,
                DsspState::T,
                DsspState::S,
                DsspState::C,
            ],
        };
        let counts = assignment.counts();
        assert_eq!(counts.h, 2);
        assert_eq!(counts.g, 1);
        assert_eq!(counts.i, 0);
        assert_eq!(counts.e, 1);
        assert_eq!(counts.b, 1);
        assert_eq!(counts.t, 1);
        assert_eq!(counts.s, 1);
        assert_eq!(counts.c, 1);
    }

    #[test]
    fn dssp_bend_detection() {
        // Build a chain with a sharp bend at residue 3 (CA angle > 70 degrees)
        // Place CAs: 0=(0,0,0), 1=(3.8,0,0), 2=(7.6,0,0), 3=(7.6,3.8,0), 4=(7.6,7.6,0)
        // Angle at CA(2): CA(0)->CA(2)->CA(4)
        // Vector CA(0)->CA(2) = (7.6,0,0), Vector CA(4)->CA(2) = (0,-7.6,0) => 90 degrees > 70
        let mut residues = Vec::new();
        let positions = vec![
            (0.0, 0.0, 0.0),
            (3.8, 0.0, 0.0),
            (7.6, 0.0, 0.0),
            (7.6, 3.8, 0.0),
            (7.6, 7.6, 0.0),
        ];
        for (i, (x, y, z)) in positions.iter().enumerate() {
            residues.push(Residue {
                name: "ALA".into(),
                seq_num: i as i32 + 1,
                i_code: None,
                atoms: vec![
                    make_atom("N", x - 0.5, *y, *z),
                    make_atom("CA", *x, *y, *z),
                    make_atom("C", x + 0.5, *y, *z),
                    make_atom("O", x + 0.5, y + 1.0, *z),
                ],
            });
        }
        let chain = Chain::new('A', residues);
        let result = dssp(&chain).unwrap();
        // Residue at index 2 should be detected as a bend (S)
        assert_eq!(
            result.states[2],
            DsspState::S,
            "residue at bend should be S, got: {}",
            result.to_string_codes()
        );
    }
}
