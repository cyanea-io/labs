//! Structural biology demo datasets: protein structures, contact maps.

/// Demo protein structure — insulin chain B (30 residues).
///
/// Returns (residue, x, y, z) tuples for Cα atoms.
pub fn insulin_chain_b() -> DemoStructure {
    DemoStructure {
        name: "Insulin_B_chain",
        pdb_id: "2INS",
        chain: 'B',
        residues: vec![
            ("PHE", 1),  ("VAL", 2),  ("ASN", 3),  ("GLN", 4),  ("HIS", 5),
            ("LEU", 6),  ("CYS", 7),  ("GLY", 8),  ("SER", 9),  ("HIS", 10),
            ("LEU", 11), ("VAL", 12), ("GLU", 13), ("ALA", 14), ("LEU", 15),
            ("TYR", 16), ("LEU", 17), ("VAL", 18), ("CYS", 19), ("GLY", 20),
            ("GLU", 21), ("ARG", 22), ("GLY", 23), ("PHE", 24), ("PHE", 25),
            ("TYR", 26), ("THR", 27), ("PRO", 28), ("LYS", 29), ("THR", 30),
        ],
        ca_coords: vec![
            (13.7, 10.5, 10.1), (12.9, 11.8, 7.0), (11.5, 14.9, 5.7),
            (12.5, 18.4, 6.5),  (14.2, 18.3, 9.8), (12.0, 16.8, 12.3),
            (12.3, 13.1, 12.0), (15.2, 11.4, 13.1), (17.8, 13.7, 12.0),
            (16.5, 16.0, 9.7),  (14.4, 14.3, 7.3),  (16.6, 12.1, 5.3),
            (18.2, 14.8, 3.3),  (16.1, 17.2, 1.6),  (14.7, 15.8, -1.5),
            (17.1, 14.3, -3.8), (17.6, 17.2, -6.2), (14.8, 19.5, -5.5),
            (14.1, 19.3, -9.1), (11.8, 17.0, -11.0),(14.5, 14.9, -12.3),
            (17.1, 17.0, -13.9),(18.4, 16.3, -17.2),(15.3, 14.8, -18.5),
            (16.0, 11.1, -18.3),(18.1, 10.4, -21.4),(15.6, 8.1, -23.0),
            (13.1, 10.4, -24.5),(11.2, 8.1, -26.5), (8.5, 10.1, -28.0),
        ],
    }
}

/// A demo protein structure (Cα trace).
#[derive(Debug, Clone)]
pub struct DemoStructure {
    pub name: &'static str,
    pub pdb_id: &'static str,
    pub chain: char,
    pub residues: Vec<(&'static str, u32)>,
    pub ca_coords: Vec<(f64, f64, f64)>,
}

impl DemoStructure {
    /// Number of residues.
    pub fn len(&self) -> usize {
        self.residues.len()
    }

    pub fn is_empty(&self) -> bool {
        self.residues.is_empty()
    }

    /// Compute distance between two Cα atoms.
    pub fn ca_distance(&self, i: usize, j: usize) -> f64 {
        let (x1, y1, z1) = self.ca_coords[i];
        let (x2, y2, z2) = self.ca_coords[j];
        ((x2 - x1).powi(2) + (y2 - y1).powi(2) + (z2 - z1).powi(2)).sqrt()
    }

    /// Compute Cα contact map (distance < threshold).
    pub fn contact_map(&self, threshold: f64) -> Vec<Vec<bool>> {
        let n = self.len();
        let mut map = vec![vec![false; n]; n];
        for i in 0..n {
            for j in i + 1..n {
                if self.ca_distance(i, j) < threshold {
                    map[i][j] = true;
                    map[j][i] = true;
                }
            }
        }
        map
    }

    /// Generate simplified PDB format.
    pub fn to_pdb(&self) -> String {
        let mut lines = Vec::new();
        lines.push(format!("HEADER    DEMO STRUCTURE {}", self.pdb_id));
        for (i, ((resname, resnum), (x, y, z))) in
            self.residues.iter().zip(self.ca_coords.iter()).enumerate()
        {
            lines.push(format!(
                "ATOM  {:5}  CA  {:3} {:1}{:4}    {:8.3}{:8.3}{:8.3}  1.00  0.00           C",
                i + 1, resname, self.chain, resnum, x, y, z
            ));
        }
        lines.push("END".to_string());
        lines.join("\n")
    }
}

/// Demo Ramachandran angles (phi, psi) for secondary structure assessment.
pub fn demo_ramachandran() -> Vec<(&'static str, f64, f64)> {
    vec![
        ("alpha_helix_1", -57.0, -47.0),
        ("alpha_helix_2", -60.0, -45.0),
        ("alpha_helix_3", -55.0, -50.0),
        ("beta_sheet_1", -120.0, 130.0),
        ("beta_sheet_2", -135.0, 135.0),
        ("beta_sheet_3", -110.0, 125.0),
        ("left_alpha", 57.0, 47.0),
        ("polyproline_II", -75.0, 145.0),
        ("turn_1", -60.0, 30.0),
        ("glycine_special", 80.0, -5.0),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_insulin_chain_b() {
        let s = insulin_chain_b();
        assert_eq!(s.len(), 30);
        assert_eq!(s.pdb_id, "2INS");
        assert_eq!(s.ca_coords.len(), 30);
    }

    #[test]
    fn test_ca_distance() {
        let s = insulin_chain_b();
        // Adjacent residues should be ~3.8 Å apart
        let d = s.ca_distance(0, 1);
        assert!(d > 2.0 && d < 6.0, "Adjacent Cα distance: {:.1}", d);
    }

    #[test]
    fn test_contact_map() {
        let s = insulin_chain_b();
        let contacts = s.contact_map(8.0);
        assert_eq!(contacts.len(), 30);
        // Symmetric
        for i in 0..30 {
            for j in 0..30 {
                assert_eq!(contacts[i][j], contacts[j][i]);
            }
            assert!(!contacts[i][i]); // no self-contacts
        }
    }

    #[test]
    fn test_to_pdb() {
        let s = insulin_chain_b();
        let pdb = s.to_pdb();
        assert!(pdb.starts_with("HEADER"));
        assert!(pdb.contains("ATOM"));
        assert!(pdb.contains("END"));
        let atom_lines: Vec<&str> = pdb.lines().filter(|l| l.starts_with("ATOM")).collect();
        assert_eq!(atom_lines.len(), 30);
    }

    #[test]
    fn test_ramachandran() {
        let angles = demo_ramachandran();
        assert_eq!(angles.len(), 10);
        // Phi and psi should be in [-180, 180]
        assert!(angles.iter().all(|(_, phi, psi)| {
            *phi >= -180.0 && *phi <= 180.0 && *psi >= -180.0 && *psi <= 180.0
        }));
    }
}
