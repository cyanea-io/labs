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

/// Parse an SDF string (may contain multiple molecules).
#[pyfunction]
fn parse_sdf(input: &str) -> PyResult<Vec<Molecule>> {
    let results = cyanea_chem::sdf::parse_sdf(input);
    let mut molecules = Vec::new();
    for r in results {
        molecules.push(Molecule {
            inner: r.into_pyresult()?,
        });
    }
    Ok(molecules)
}

/// Parse an SDF file from disk.
#[pyfunction]
fn parse_sdf_file(path: &str) -> PyResult<Vec<Molecule>> {
    let mols = cyanea_chem::sdf::parse_sdf_file(path).into_pyresult()?;
    Ok(mols.into_iter().map(|m| Molecule { inner: m }).collect())
}

/// Compute MACCS-like 166-key structural fingerprint (on-bit indices).
#[pyfunction]
fn maccs_fingerprint(smiles: &str) -> PyResult<Vec<usize>> {
    let mol = cyanea_chem::parse_smiles(smiles).into_pyresult()?;
    let fp = cyanea_chem::maccs::maccs_fingerprint(&mol);
    Ok((0..166).filter(|&i| fp.get_bit(i)).collect())
}

// ---------------------------------------------------------------------------
// 3D Coordinate Generation
// ---------------------------------------------------------------------------

/// A 3D conformer with xyz coordinates.
#[pyclass(frozen, get_all)]
pub struct PyConformer {
    pub coords: Vec<[f64; 3]>,
    pub n_atoms: usize,
}

/// Energy components from force field calculation.
#[pyclass(frozen, get_all)]
pub struct PyEnergyComponents {
    pub bond_stretch: f64,
    pub angle_bend: f64,
    pub torsion: f64,
    pub van_der_waals: f64,
    pub electrostatic: f64,
    pub out_of_plane: f64,
    pub total: f64,
}

/// Result of energy minimization.
#[pyclass(frozen, get_all)]
pub struct PyMinimizeResult {
    pub initial_energy: f64,
    pub final_energy: f64,
    pub n_steps: usize,
    pub converged: bool,
    pub coords: Vec<[f64; 3]>,
}

/// A reaction product.
#[pyclass(frozen, get_all)]
pub struct PyReactionProduct {
    pub smiles: String,
}

/// A retrosynthetic disconnection.
#[pyclass(frozen, get_all)]
pub struct PyDisconnection {
    pub transform_name: String,
    pub smirks: String,
    pub precursors: Vec<String>,
}

/// Atom-atom mapping between two molecules.
#[pyclass(frozen, get_all)]
pub struct PyAtomMapping {
    pub mapping: Vec<(usize, usize)>,
    pub unmapped_reactant: Vec<usize>,
    pub unmapped_product: Vec<usize>,
}

/// Embed a single 3D conformer from a SMILES string.
#[pyfunction]
#[pyo3(signature = (smiles, seed=42, use_torsion_prefs=true, max_minimize_steps=200))]
fn embed_conformer(
    smiles: &str,
    seed: u64,
    use_torsion_prefs: bool,
    max_minimize_steps: usize,
) -> PyResult<PyConformer> {
    let mol = cyanea_chem::parse_smiles(smiles).into_pyresult()?;
    let config = cyanea_chem::EmbedConfig {
        random_seed: seed,
        use_torsion_prefs,
        max_minimize_steps,
        ..Default::default()
    };
    let conf = cyanea_chem::embed_molecule(&mol, &config).into_pyresult()?;
    Ok(PyConformer {
        n_atoms: conf.len(),
        coords: conf.coords,
    })
}

/// Compute UFF energy for a SMILES string (auto-embeds 3D coordinates).
#[pyfunction]
fn uff_energy(smiles: &str) -> PyResult<PyEnergyComponents> {
    let mol = cyanea_chem::parse_smiles(smiles).into_pyresult()?;
    let config = cyanea_chem::EmbedConfig {
        force_field: cyanea_chem::ForceFieldType::None,
        max_minimize_steps: 0,
        ..Default::default()
    };
    let conf = cyanea_chem::embed_molecule(&mol, &config).into_pyresult()?;
    let e = cyanea_chem::uff_energy(&mol, &conf).into_pyresult()?;
    Ok(PyEnergyComponents {
        bond_stretch: e.bond_stretch,
        angle_bend: e.angle_bend,
        torsion: e.torsion,
        van_der_waals: e.van_der_waals,
        electrostatic: e.electrostatic,
        out_of_plane: e.out_of_plane,
        total: e.total,
    })
}

/// Compute MMFF94 energy for a SMILES string (auto-embeds 3D coordinates).
#[pyfunction]
fn mmff94_energy(smiles: &str) -> PyResult<PyEnergyComponents> {
    let mol = cyanea_chem::parse_smiles(smiles).into_pyresult()?;
    let config = cyanea_chem::EmbedConfig {
        force_field: cyanea_chem::ForceFieldType::None,
        max_minimize_steps: 0,
        ..Default::default()
    };
    let conf = cyanea_chem::embed_molecule(&mol, &config).into_pyresult()?;
    let e = cyanea_chem::mmff94_energy(&mol, &conf).into_pyresult()?;
    Ok(PyEnergyComponents {
        bond_stretch: e.bond_stretch,
        angle_bend: e.angle_bend,
        torsion: e.torsion,
        van_der_waals: e.van_der_waals,
        electrostatic: e.electrostatic,
        out_of_plane: e.out_of_plane,
        total: e.total,
    })
}

/// Minimize energy using the specified force field.
#[pyfunction]
#[pyo3(signature = (smiles, force_field="uff", max_steps=500, gradient_threshold=0.1))]
fn minimize(
    smiles: &str,
    force_field: &str,
    max_steps: usize,
    gradient_threshold: f64,
) -> PyResult<PyMinimizeResult> {
    let mol = cyanea_chem::parse_smiles(smiles).into_pyresult()?;
    let embed_config = cyanea_chem::EmbedConfig {
        force_field: cyanea_chem::ForceFieldType::None,
        max_minimize_steps: 0,
        ..Default::default()
    };
    let conf = cyanea_chem::embed_molecule(&mol, &embed_config).into_pyresult()?;
    let min_config = cyanea_chem::MinimizeConfig {
        max_steps,
        gradient_threshold,
        method: cyanea_chem::MinimizeMethod::SteepestDescent,
    };
    let result = match force_field {
        "mmff94" | "MMFF94" => cyanea_chem::mmff94_minimize(&mol, &conf, &min_config).into_pyresult()?,
        _ => cyanea_chem::uff_minimize(&mol, &conf, &min_config).into_pyresult()?,
    };
    Ok(PyMinimizeResult {
        initial_energy: result.initial_energy,
        final_energy: result.final_energy,
        n_steps: result.n_steps,
        converged: result.converged,
        coords: result.conformer.coords,
    })
}

// ---------------------------------------------------------------------------
// Chemical Reactions
// ---------------------------------------------------------------------------

/// Apply a SMIRKS reaction to a molecule (SMILES).
#[pyfunction]
fn apply_reaction(smiles: &str, smirks: &str) -> PyResult<Vec<PyReactionProduct>> {
    let mol = cyanea_chem::parse_smiles(smiles).into_pyresult()?;
    let rxn = cyanea_chem::parse_reaction(smirks).into_pyresult()?;
    let products = cyanea_chem::apply_reaction(&mol, &rxn).into_pyresult()?;
    Ok(products.iter().map(|p| PyReactionProduct {
        smiles: p.smiles.clone(),
    }).collect())
}

/// Enumerate reactions from multiple reactant sets.
#[pyfunction]
fn enumerate_reactions(reactant_smiles: Vec<Vec<String>>, smirks: &str) -> PyResult<Vec<PyReactionProduct>> {
    let reactant_sets: Vec<Vec<&str>> = reactant_smiles
        .iter()
        .map(|set| set.iter().map(|s| s.as_str()).collect())
        .collect();
    let products = cyanea_chem::enumerate_reactions(&reactant_sets, smirks).into_pyresult()?;
    Ok(products.iter().map(|p| PyReactionProduct {
        smiles: p.smiles.clone(),
    }).collect())
}

/// Find retrosynthetic disconnections for a target molecule.
#[pyfunction]
fn retrosynthetic_disconnect(smiles: &str) -> PyResult<Vec<PyDisconnection>> {
    let mol = cyanea_chem::parse_smiles(smiles).into_pyresult()?;
    let disconnections = cyanea_chem::retrosynthetic_disconnections(&mol);
    Ok(disconnections.iter().map(|d| PyDisconnection {
        transform_name: d.transform_name.clone(),
        smirks: d.smirks.clone(),
        precursors: d.precursors.clone(),
    }).collect())
}

/// Compute atom-atom mapping between reactant and product SMILES.
#[pyfunction]
fn atom_atom_mapping(reactant_smiles: &str, product_smiles: &str) -> PyResult<PyAtomMapping> {
    let reactant = cyanea_chem::parse_smiles(reactant_smiles).into_pyresult()?;
    let product = cyanea_chem::parse_smiles(product_smiles).into_pyresult()?;
    let mapping = cyanea_chem::atom_atom_map(&reactant, &product);
    Ok(PyAtomMapping {
        mapping: mapping.mapping,
        unmapped_reactant: mapping.unmapped_reactant,
        unmapped_product: mapping.unmapped_product,
    })
}

/// Compute Gasteiger-Marsili partial charges.
#[pyfunction]
fn gasteiger_charges(smiles: &str) -> PyResult<Vec<f64>> {
    let mol = cyanea_chem::parse_smiles(smiles).into_pyresult()?;
    cyanea_chem::gasteiger_charges(&mol).into_pyresult()
}

// ---------------------------------------------------------------------------
// Submodule registration
// ---------------------------------------------------------------------------

pub fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent.py(), "chem")?;
    m.add_class::<Molecule>()?;
    m.add_class::<MolecularProperties>()?;
    m.add_class::<PyConformer>()?;
    m.add_class::<PyEnergyComponents>()?;
    m.add_class::<PyMinimizeResult>()?;
    m.add_class::<PyReactionProduct>()?;
    m.add_class::<PyDisconnection>()?;
    m.add_class::<PyAtomMapping>()?;
    m.add_function(wrap_pyfunction!(parse_smiles, &m)?)?;
    m.add_function(wrap_pyfunction!(molecular_properties, &m)?)?;
    m.add_function(wrap_pyfunction!(tanimoto, &m)?)?;
    m.add_function(wrap_pyfunction!(canonical_smiles, &m)?)?;
    m.add_function(wrap_pyfunction!(parse_sdf, &m)?)?;
    m.add_function(wrap_pyfunction!(parse_sdf_file, &m)?)?;
    m.add_function(wrap_pyfunction!(maccs_fingerprint, &m)?)?;
    m.add_function(wrap_pyfunction!(embed_conformer, &m)?)?;
    m.add_function(wrap_pyfunction!(uff_energy, &m)?)?;
    m.add_function(wrap_pyfunction!(mmff94_energy, &m)?)?;
    m.add_function(wrap_pyfunction!(minimize, &m)?)?;
    m.add_function(wrap_pyfunction!(apply_reaction, &m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_reactions, &m)?)?;
    m.add_function(wrap_pyfunction!(retrosynthetic_disconnect, &m)?)?;
    m.add_function(wrap_pyfunction!(atom_atom_mapping, &m)?)?;
    m.add_function(wrap_pyfunction!(gasteiger_charges, &m)?)?;
    parent.add_submodule(&m)?;
    Ok(())
}
