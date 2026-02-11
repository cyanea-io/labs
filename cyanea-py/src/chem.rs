//! Python bindings for cyanea-chem: chemistry and small-molecule analysis.

use pyo3::prelude::*;

use crate::error::IntoPyResult;

// ---------------------------------------------------------------------------
// Molecule
// ---------------------------------------------------------------------------

/// A molecular graph parsed from SMILES.
#[pyclass(frozen)]
pub struct Molecule {
    inner: cyanea_chem::Molecule,
}

#[pymethods]
impl Molecule {
    /// Parse a molecule from a SMILES string.
    #[new]
    fn new(smiles: &str) -> PyResult<Self> {
        let inner = cyanea_chem::parse_smiles(smiles).into_pyresult()?;
        Ok(Self { inner })
    }

    /// Number of (heavy) atoms in the molecular graph.
    fn atom_count(&self) -> usize {
        self.inner.atom_count()
    }

    /// Number of bonds.
    fn bond_count(&self) -> usize {
        self.inner.bond_count()
    }

    /// Molecular formula in Hill system order.
    fn molecular_formula(&self) -> String {
        cyanea_chem::molecular_formula(&self.inner)
    }

    /// Molecular weight in Daltons.
    fn molecular_weight(&self) -> f64 {
        cyanea_chem::molecular_weight(&self.inner)
    }

    /// Canonical SMILES string.
    fn canonical_smiles(&self) -> String {
        cyanea_chem::canonical_smiles(&self.inner)
    }

    /// Morgan fingerprint on-bit indices.
    #[pyo3(signature = (radius=2, n_bits=2048))]
    fn morgan_fingerprint(&self, radius: usize, n_bits: usize) -> Vec<usize> {
        let fp = cyanea_chem::morgan_fingerprint(&self.inner, radius, n_bits);
        (0..n_bits).filter(|&i| fp.get_bit(i)).collect()
    }

    /// Check if this molecule contains a substructure given as SMILES.
    fn has_substructure(&self, pattern: &str) -> PyResult<bool> {
        let pattern_mol = cyanea_chem::parse_smiles(pattern).into_pyresult()?;
        Ok(cyanea_chem::has_substructure(&self.inner, &pattern_mol))
    }

    fn __repr__(&self) -> String {
        format!("Molecule('{}')", cyanea_chem::canonical_smiles(&self.inner))
    }

    fn __len__(&self) -> usize {
        self.inner.atom_count()
    }

    fn __str__(&self) -> String {
        cyanea_chem::canonical_smiles(&self.inner)
    }
}

// ---------------------------------------------------------------------------
// MolecularProperties
// ---------------------------------------------------------------------------

/// Computed molecular properties.
#[pyclass(frozen, get_all)]
pub struct MolecularProperties {
    pub formula: String,
    pub molecular_weight: f64,
    pub atom_count: usize,
    pub bond_count: usize,
    pub ring_count: usize,
    pub rotatable_bonds: usize,
    pub hbd: usize,
    pub hba: usize,
}

// ---------------------------------------------------------------------------
// Module functions
// ---------------------------------------------------------------------------

/// Parse a SMILES string into a Molecule.
#[pyfunction]
fn parse_smiles(smiles: &str) -> PyResult<Molecule> {
    let inner = cyanea_chem::parse_smiles(smiles).into_pyresult()?;
    Ok(Molecule { inner })
}

/// Compute molecular properties from a SMILES string.
#[pyfunction]
fn molecular_properties(smiles: &str) -> PyResult<MolecularProperties> {
    let mol = cyanea_chem::parse_smiles(smiles).into_pyresult()?;
    let props = cyanea_chem::compute_properties(&mol);
    Ok(MolecularProperties {
        formula: props.formula,
        molecular_weight: props.molecular_weight,
        atom_count: mol.atom_count(),
        bond_count: mol.bond_count(),
        ring_count: props.ring_count,
        rotatable_bonds: props.rotatable_bonds,
        hbd: props.hydrogen_bond_donors,
        hba: props.hydrogen_bond_acceptors,
    })
}

/// Tanimoto similarity between two molecules given as SMILES.
#[pyfunction]
#[pyo3(signature = (smiles1, smiles2, radius=2, n_bits=2048))]
fn tanimoto(smiles1: &str, smiles2: &str, radius: usize, n_bits: usize) -> PyResult<f64> {
    let mol1 = cyanea_chem::parse_smiles(smiles1).into_pyresult()?;
    let mol2 = cyanea_chem::parse_smiles(smiles2).into_pyresult()?;
    let fp1 = cyanea_chem::morgan_fingerprint(&mol1, radius, n_bits);
    let fp2 = cyanea_chem::morgan_fingerprint(&mol2, radius, n_bits);
    Ok(cyanea_chem::tanimoto_similarity(&fp1, &fp2))
}

/// Return the canonical SMILES for a given SMILES string.
#[pyfunction]
fn canonical_smiles(smiles: &str) -> PyResult<String> {
    let mol = cyanea_chem::parse_smiles(smiles).into_pyresult()?;
    Ok(cyanea_chem::canonical_smiles(&mol))
}

// ---------------------------------------------------------------------------
// Submodule registration
// ---------------------------------------------------------------------------

pub fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent.py(), "chem")?;
    m.add_class::<Molecule>()?;
    m.add_class::<MolecularProperties>()?;
    m.add_function(wrap_pyfunction!(parse_smiles, &m)?)?;
    m.add_function(wrap_pyfunction!(molecular_properties, &m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto, &m)?)?;
    m.add_function(wrap_pyfunction!(canonical_smiles, &m)?)?;
    parent.add_submodule(&m)?;
    Ok(())
}
