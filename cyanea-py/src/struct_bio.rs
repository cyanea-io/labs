//! Python bindings for cyanea-struct: protein 3D structure analysis.

use pyo3::prelude::*;

use crate::error::IntoPyResult;

// ---------------------------------------------------------------------------
// Structure
// ---------------------------------------------------------------------------

/// A macromolecular 3D structure parsed from PDB format.
#[pyclass(frozen)]
pub struct Structure {
    inner: cyanea_struct::Structure,
}

#[pymethods]
impl Structure {
    /// Number of chains.
    fn chain_count(&self) -> usize {
        self.inner.chain_count()
    }

    /// Total residues across all chains.
    fn residue_count(&self) -> usize {
        self.inner.residue_count()
    }

    /// Total atoms across all chains.
    fn atom_count(&self) -> usize {
        self.inner.atom_count()
    }

    /// PDB identifier.
    fn id(&self) -> String {
        self.inner.id.clone()
    }

    /// Assign secondary structure (simplified 4-state) for the first chain.
    fn secondary_structure(&self) -> PyResult<Vec<SecondaryStructureAssignment>> {
        let chain = self.inner.chains.first().ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err("structure has no chains")
        })?;
        let ss = cyanea_struct::assign_secondary_structure(chain).into_pyresult()?;
        Ok(ss
            .assignments
            .iter()
            .enumerate()
            .map(|(i, a)| {
                let residue = &chain.residues[i];
                SecondaryStructureAssignment {
                    residue_num: residue.seq_num as usize,
                    residue_name: residue.name.clone(),
                    structure: match a {
                        cyanea_struct::SecondaryStructure::Helix => "Helix".to_string(),
                        cyanea_struct::SecondaryStructure::Sheet => "Sheet".to_string(),
                        cyanea_struct::SecondaryStructure::Turn => "Turn".to_string(),
                        cyanea_struct::SecondaryStructure::Coil => "Coil".to_string(),
                    },
                }
            })
            .collect())
    }

    fn __repr__(&self) -> String {
        format!(
            "Structure('{}', chains={}, residues={})",
            self.inner.id,
            self.inner.chain_count(),
            self.inner.residue_count(),
        )
    }

    fn __str__(&self) -> String {
        format!(
            "Structure {} -- {} chain(s), {} residue(s), {} atom(s)",
            self.inner.id,
            self.inner.chain_count(),
            self.inner.residue_count(),
            self.inner.atom_count(),
        )
    }
}

// ---------------------------------------------------------------------------
// SecondaryStructureAssignment
// ---------------------------------------------------------------------------

/// Secondary structure assignment for a single residue.
#[pyclass(frozen, get_all)]
pub struct SecondaryStructureAssignment {
    pub residue_num: usize,
    pub residue_name: String,
    pub structure: String,
}

// ---------------------------------------------------------------------------
// Module functions
// ---------------------------------------------------------------------------

/// Parse PDB-format text into a Structure.
#[pyfunction]
fn parse_pdb(content: &str) -> PyResult<Structure> {
    let inner = cyanea_struct::parse_pdb(content).into_pyresult()?;
    Ok(Structure { inner })
}

/// Compute RMSD between two equal-length coordinate sets (no alignment).
#[pyfunction]
fn rmsd(coords1: Vec<[f64; 3]>, coords2: Vec<[f64; 3]>) -> PyResult<f64> {
    let p1: Vec<cyanea_struct::Point3D> = coords1
        .iter()
        .map(|c| cyanea_struct::Point3D::new(c[0], c[1], c[2]))
        .collect();
    let p2: Vec<cyanea_struct::Point3D> = coords2
        .iter()
        .map(|c| cyanea_struct::Point3D::new(c[0], c[1], c[2]))
        .collect();
    cyanea_struct::geometry::rmsd_points(&p1, &p2).into_pyresult()
}

/// Kabsch optimal superposition. Returns (rmsd, aligned_coords).
#[pyfunction]
fn kabsch_align(
    moving: Vec<[f64; 3]>,
    target: Vec<[f64; 3]>,
) -> PyResult<(f64, Vec<[f64; 3]>)> {
    let target_pts: Vec<cyanea_struct::Point3D> = target
        .iter()
        .map(|c| cyanea_struct::Point3D::new(c[0], c[1], c[2]))
        .collect();
    let moving_pts: Vec<cyanea_struct::Point3D> = moving
        .iter()
        .map(|c| cyanea_struct::Point3D::new(c[0], c[1], c[2]))
        .collect();
    let result =
        cyanea_struct::kabsch_points(&target_pts, &moving_pts).into_pyresult()?;
    let aligned: Vec<[f64; 3]> = result
        .transformed_coords
        .iter()
        .map(|p| [p.x, p.y, p.z])
        .collect();
    Ok((result.rmsd, aligned))
}

/// Parse mmCIF (PDBx) format text into a Structure.
#[pyfunction]
fn parse_mmcif(content: &str) -> PyResult<Structure> {
    let inner = cyanea_struct::parse_mmcif(content).into_pyresult()?;
    Ok(Structure { inner })
}

// ---------------------------------------------------------------------------
// Contact maps
// ---------------------------------------------------------------------------

/// Contact map entry: (residue_i, residue_j, distance).
#[pyclass(frozen, get_all)]
pub struct Contact {
    pub residue_i: usize,
    pub residue_j: usize,
    pub distance: f64,
}

/// Compute CA-CA contact map for a chain.
///
/// Returns contacts below the cutoff distance.
#[pyfunction]
#[pyo3(signature = (structure, chain_index=0, cutoff=8.0))]
fn contact_map(structure: &Structure, chain_index: usize, cutoff: f64) -> PyResult<Vec<Contact>> {
    let chain = structure.inner.chains.get(chain_index).ok_or_else(|| {
        pyo3::exceptions::PyIndexError::new_err(format!(
            "chain index {chain_index} out of range (structure has {} chains)",
            structure.inner.chains.len()
        ))
    })?;
    let cm = cyanea_struct::compute_contact_map(chain).into_pyresult()?;
    let pairs = cm.contacts_below(cutoff);
    Ok(pairs
        .into_iter()
        .map(|(i, j)| Contact {
            residue_i: i,
            residue_j: j,
            distance: cm.get(i, j),
        })
        .collect())
}

// ---------------------------------------------------------------------------
// Ramachandran
// ---------------------------------------------------------------------------

/// Ramachandran angles for a single residue.
#[pyclass(frozen, get_all)]
pub struct RamachandranEntry {
    pub residue_num: usize,
    pub residue_name: String,
    pub phi: f64,
    pub psi: f64,
    pub region: String,
}

/// Compute Ramachandran phi/psi angles for all residues in a structure.
#[pyfunction]
fn ramachandran(structure: &Structure) -> PyResult<Vec<RamachandranEntry>> {
    let report = cyanea_struct::ramachandran_report(&structure.inner).into_pyresult()?;
    Ok(report
        .into_iter()
        .map(|(num, name, phi, psi, region)| RamachandranEntry {
            residue_num: num,
            residue_name: name,
            phi,
            psi,
            region: format!("{:?}", region),
        })
        .collect())
}

// ---------------------------------------------------------------------------
// B-factor analysis
// ---------------------------------------------------------------------------

/// Per-residue average B-factor.
#[pyclass(frozen, get_all)]
pub struct ResidueBfactor {
    pub residue_num: usize,
    pub residue_name: String,
    pub bfactor: f64,
}

/// Compute per-residue average B-factors.
#[pyfunction]
fn residue_bfactors(structure: &Structure) -> PyResult<Vec<ResidueBfactor>> {
    let bfactors = cyanea_struct::residue_bfactors(&structure.inner).into_pyresult()?;
    Ok(bfactors
        .into_iter()
        .map(|(num, name, bf)| ResidueBfactor {
            residue_num: num,
            residue_name: name,
            bfactor: bf,
        })
        .collect())
}

/// B-factor statistics for a chain.
#[pyclass(frozen, get_all)]
pub struct ChainBfactorStats {
    pub chain_id: String,
    pub mean: f64,
    pub std_dev: f64,
    pub min: f64,
    pub max: f64,
    pub median: f64,
}

/// Compute per-chain B-factor statistics.
#[pyfunction]
fn chain_bfactors(structure: &Structure) -> PyResult<Vec<ChainBfactorStats>> {
    let stats = cyanea_struct::chain_bfactors(&structure.inner).into_pyresult()?;
    Ok(stats
        .into_iter()
        .map(|(id, s)| ChainBfactorStats {
            chain_id: id,
            mean: s.mean,
            std_dev: s.std_dev,
            min: s.min,
            max: s.max,
            median: s.median,
        })
        .collect())
}

// ---------------------------------------------------------------------------
// Submodule registration
// ---------------------------------------------------------------------------

pub fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent.py(), "struct_bio")?;
    m.add_class::<Structure>()?;
    m.add_class::<SecondaryStructureAssignment>()?;
    m.add_class::<Contact>()?;
    m.add_class::<RamachandranEntry>()?;
    m.add_class::<ResidueBfactor>()?;
    m.add_class::<ChainBfactorStats>()?;
    m.add_function(wrap_pyfunction!(parse_pdb, &m)?)?;
    m.add_function(wrap_pyfunction!(parse_mmcif, &m)?)?;
    m.add_function(wrap_pyfunction!(rmsd, &m)?)?;
    m.add_function(wrap_pyfunction!(kabsch_align, &m)?)?;
    m.add_function(wrap_pyfunction!(contact_map, &m)?)?;
    m.add_function(wrap_pyfunction!(ramachandran, &m)?)?;
    m.add_function(wrap_pyfunction!(residue_bfactors, &m)?)?;
    m.add_function(wrap_pyfunction!(chain_bfactors, &m)?)?;
    parent.add_submodule(&m)?;
    Ok(())
}
