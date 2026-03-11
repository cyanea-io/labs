//! Chemistry demo datasets: small molecule library.

/// Demo small molecule library — 12 FDA-approved drugs.
///
/// Common drugs spanning different therapeutic areas with SMILES and properties.
pub fn fda_approved_drugs() -> Vec<DemoMolecule> {
    vec![
        DemoMolecule {
            name: "Aspirin",
            smiles: "CC(=O)OC1=CC=CC=C1C(=O)O",
            molecular_weight: 180.16,
            logp: 1.2,
            hbd: 1,
            hba: 4,
            therapeutic_area: "Anti-inflammatory",
        },
        DemoMolecule {
            name: "Ibuprofen",
            smiles: "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
            molecular_weight: 206.28,
            logp: 3.5,
            hbd: 1,
            hba: 2,
            therapeutic_area: "Anti-inflammatory",
        },
        DemoMolecule {
            name: "Metformin",
            smiles: "CN(C)C(=N)NC(=N)N",
            molecular_weight: 129.16,
            logp: -1.4,
            hbd: 3,
            hba: 5,
            therapeutic_area: "Antidiabetic",
        },
        DemoMolecule {
            name: "Atorvastatin",
            smiles: "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
            molecular_weight: 558.64,
            logp: 4.5,
            hbd: 4,
            hba: 6,
            therapeutic_area: "Lipid-lowering",
        },
        DemoMolecule {
            name: "Omeprazole",
            smiles: "CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=CC=CC=C3N2",
            molecular_weight: 345.42,
            logp: 2.2,
            hbd: 1,
            hba: 5,
            therapeutic_area: "Proton pump inhibitor",
        },
        DemoMolecule {
            name: "Amoxicillin",
            smiles: "CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=C(C=C3)O)N)C(=O)O)C",
            molecular_weight: 365.40,
            logp: 0.9,
            hbd: 4,
            hba: 7,
            therapeutic_area: "Antibiotic",
        },
        DemoMolecule {
            name: "Lisinopril",
            smiles: "C(CC(C(=O)O)NC(CCCCN)C(=O)O)CC1=CC=CC=C1",
            molecular_weight: 405.49,
            logp: -0.8,
            hbd: 4,
            hba: 7,
            therapeutic_area: "ACE inhibitor",
        },
        DemoMolecule {
            name: "Sertraline",
            smiles: "CNC1CCC(C2=CC=CC=C12)C3=CC(=C(C=C3)Cl)Cl",
            molecular_weight: 306.23,
            logp: 5.1,
            hbd: 1,
            hba: 1,
            therapeutic_area: "Antidepressant",
        },
        DemoMolecule {
            name: "Amlodipine",
            smiles: "CCOC(=O)C1=C(NC(=C(C1C2=CC=CC=C2Cl)C(=O)OC)C)COCCN",
            molecular_weight: 408.88,
            logp: 3.0,
            hbd: 2,
            hba: 6,
            therapeutic_area: "Calcium channel blocker",
        },
        DemoMolecule {
            name: "Doxorubicin",
            smiles: "CC1C(C(CC(O1)OC2CC(CC3=C2C(=C4C(=C3O)C(=O)C5=CC=CC(=C5C4=O)OC)O)(C(=O)CO)O)N)O",
            molecular_weight: 543.52,
            logp: 1.3,
            hbd: 6,
            hba: 12,
            therapeutic_area: "Chemotherapy",
        },
        DemoMolecule {
            name: "Imatinib",
            smiles: "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5",
            molecular_weight: 493.60,
            logp: 3.5,
            hbd: 2,
            hba: 7,
            therapeutic_area: "Targeted therapy",
        },
        DemoMolecule {
            name: "Caffeine",
            smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            molecular_weight: 194.19,
            logp: -0.1,
            hbd: 0,
            hba: 6,
            therapeutic_area: "Stimulant",
        },
    ]
}

/// A demo molecule with key properties.
#[derive(Debug, Clone)]
pub struct DemoMolecule {
    pub name: &'static str,
    pub smiles: &'static str,
    pub molecular_weight: f64,
    pub logp: f64,
    pub hbd: usize,
    pub hba: usize,
    pub therapeutic_area: &'static str,
}

impl DemoMolecule {
    /// Check Lipinski's Rule of Five compliance.
    pub fn lipinski_compliant(&self) -> bool {
        self.molecular_weight <= 500.0
            && self.logp <= 5.0
            && self.hbd <= 5
            && self.hba <= 10
    }

    /// Count Lipinski violations.
    pub fn lipinski_violations(&self) -> usize {
        let mut v = 0;
        if self.molecular_weight > 500.0 { v += 1; }
        if self.logp > 5.0 { v += 1; }
        if self.hbd > 5 { v += 1; }
        if self.hba > 10 { v += 1; }
        v
    }
}

/// Generate SDF-like text from the molecule library.
pub fn demo_sdf_summary() -> String {
    let mut lines = Vec::new();
    for mol in fda_approved_drugs() {
        lines.push(format!(">{} <NAME>", mol.name));
        lines.push(mol.name.to_string());
        lines.push(format!("> <SMILES>"));
        lines.push(mol.smiles.to_string());
        lines.push(format!("> <MW>"));
        lines.push(format!("{:.2}", mol.molecular_weight));
        lines.push("$$$$".to_string());
    }
    lines.join("\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fda_drugs() {
        let drugs = fda_approved_drugs();
        assert_eq!(drugs.len(), 12);
        assert!(drugs.iter().all(|d| !d.smiles.is_empty()));
        assert!(drugs.iter().all(|d| d.molecular_weight > 0.0));
    }

    #[test]
    fn test_lipinski() {
        let drugs = fda_approved_drugs();
        // Most FDA-approved drugs should be Lipinski-compliant or close
        let compliant = drugs.iter().filter(|d| d.lipinski_compliant()).count();
        assert!(compliant >= 8); // at least 8 of 12
    }

    #[test]
    fn test_lipinski_violations() {
        let caffeine = &fda_approved_drugs()[11];
        assert_eq!(caffeine.name, "Caffeine");
        assert_eq!(caffeine.lipinski_violations(), 0);
    }

    #[test]
    fn test_doxorubicin_violations() {
        let dox = &fda_approved_drugs()[9];
        assert_eq!(dox.name, "Doxorubicin");
        // MW > 500, HBA > 10, HBD > 5
        assert!(dox.lipinski_violations() >= 2);
    }

    #[test]
    fn test_sdf_summary() {
        let sdf = demo_sdf_summary();
        assert!(sdf.contains("$$$$"));
        let records = sdf.matches("$$$$").count();
        assert_eq!(records, 12);
    }
}
