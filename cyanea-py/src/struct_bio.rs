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

// ---------------------------------------------------------------------------
// Submodule registration
// ---------------------------------------------------------------------------

pub fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent.py(), "struct_bio")?;
    m.add_class::<Structure>()?;
    m.add_class::<SecondaryStructureAssignment>()?;
    m.add_function(wrap_pyfunction!(parse_pdb, &m)?)?;
    m.add_function(wrap_pyfunction!(rmsd, &m)?)?;
    m.add_function(wrap_pyfunction!(kabsch_align, &m)?)?;
    parent.add_submodule(&m)?;
    Ok(())
}
