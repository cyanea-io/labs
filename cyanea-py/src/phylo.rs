//! Python bindings for cyanea-phylo: phylogenetic trees and distances.

use pyo3::prelude::*;

use crate::error::IntoPyResult;

// ---------------------------------------------------------------------------
// PhyloTree
// ---------------------------------------------------------------------------

/// A rooted phylogenetic tree.
#[pyclass(frozen)]
pub struct PhyloTree {
    inner: cyanea_phylo::PhyloTree,
}

#[pymethods]
impl PhyloTree {
    /// Number of leaf nodes.
    fn leaf_count(&self) -> usize {
        self.inner.leaf_count()
    }

    /// Sorted list of leaf names.
    fn leaf_names(&self) -> Vec<String> {
        self.inner.leaf_names()
    }

    /// Serialize the tree to Newick format.
    fn to_newick(&self) -> String {
        cyanea_phylo::write_newick(&self.inner)
    }

    /// Robinson-Foulds distance to another tree.
    fn robinson_foulds(&self, other: &PhyloTree) -> PyResult<usize> {
        cyanea_phylo::robinson_foulds(&self.inner, &other.inner).into_pyresult()
    }

    fn __repr__(&self) -> String {
        format!("PhyloTree(leaves={})", self.inner.leaf_count())
    }

    fn __len__(&self) -> usize {
        self.inner.leaf_count()
    }
}

// ---------------------------------------------------------------------------
// Module functions
// ---------------------------------------------------------------------------

/// Parse a Newick format string into a PhyloTree.
#[pyfunction]
fn parse_newick(newick: &str) -> PyResult<PhyloTree> {
    let inner = cyanea_phylo::parse_newick(newick).into_pyresult()?;
    Ok(PhyloTree { inner })
}

/// Compute evolutionary distance between two aligned sequences.
///
/// Models: "p" (p-distance), "jc" (Jukes-Cantor), "k2p" (Kimura 2-parameter).
#[pyfunction]
#[pyo3(signature = (seq1, seq2, model="jc"))]
fn evolutionary_distance(seq1: &str, seq2: &str, model: &str) -> PyResult<f64> {
    let a = seq1.as_bytes();
    let b = seq2.as_bytes();
    match model {
        "p" => cyanea_phylo::p_distance(a, b).into_pyresult(),
        "jc" => {
            let p = cyanea_phylo::p_distance(a, b).into_pyresult()?;
            cyanea_phylo::jukes_cantor(p).into_pyresult()
        }
        "k2p" => {
            // For K2P we need transition and transversion proportions.
            // Count them manually from the sequences.
            let len = a.len();
            if len == 0 || a.len() != b.len() {
                return cyanea_phylo::p_distance(a, b)
                    .into_pyresult()
                    .map(|_| unreachable!());
            }
            let mut transitions = 0usize;
            let mut transversions = 0usize;
            for (&x, &y) in a.iter().zip(b.iter()) {
                let xu = x.to_ascii_uppercase();
                let yu = y.to_ascii_uppercase();
                if xu == yu {
                    continue;
                }
                match (xu, yu) {
                    (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C') => {
                        transitions += 1;
                    }
                    (b'A', b'C')
                    | (b'C', b'A')
                    | (b'A', b'T')
                    | (b'T', b'A')
                    | (b'G', b'C')
                    | (b'C', b'G')
                    | (b'G', b'T')
                    | (b'T', b'G') => {
                        transversions += 1;
                    }
                    _ => {}
                }
            }
            let total = len as f64;
            cyanea_phylo::kimura_2p(transitions as f64 / total, transversions as f64 / total)
                .into_pyresult()
        }
        _ => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "unknown distance model: {model} (expected 'p', 'jc', or 'k2p')"
        ))),
    }
}

/// Build a UPGMA tree from a distance matrix.
#[pyfunction]
fn upgma(labels: Vec<String>, matrix: Vec<Vec<f64>>) -> PyResult<PhyloTree> {
    let dm = build_distance_matrix(&labels, &matrix)?;
    let inner = cyanea_phylo::upgma(&dm, &labels).into_pyresult()?;
    Ok(PhyloTree { inner })
}

/// Build a Neighbor-Joining tree from a distance matrix.
#[pyfunction]
fn neighbor_joining(labels: Vec<String>, matrix: Vec<Vec<f64>>) -> PyResult<PhyloTree> {
    let dm = build_distance_matrix(&labels, &matrix)?;
    let inner = cyanea_phylo::neighbor_joining(&dm, &labels).into_pyresult()?;
    Ok(PhyloTree { inner })
}

/// Convert a square distance matrix (as Vec<Vec<f64>>) to a DistanceMatrix.
fn build_distance_matrix(
    labels: &[String],
    matrix: &[Vec<f64>],
) -> PyResult<cyanea_ml::DistanceMatrix> {
    let n = labels.len();
    if matrix.len() != n {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "matrix has {} rows but {} labels",
            matrix.len(),
            n
        )));
    }
    for (i, row) in matrix.iter().enumerate() {
        if row.len() != n {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "matrix row {} has {} columns, expected {}",
                i,
                row.len(),
                n
            )));
        }
    }
    // Extract upper triangle into condensed form.
    let mut condensed = Vec::with_capacity(n * (n - 1) / 2);
    for i in 0..n {
        for j in (i + 1)..n {
            condensed.push(matrix[i][j]);
        }
    }
    cyanea_ml::DistanceMatrix::from_condensed(condensed, n).into_pyresult()
}

// ---------------------------------------------------------------------------
// NEXUS I/O
// ---------------------------------------------------------------------------

/// A parsed NEXUS file containing taxa and trees.
#[pyclass(frozen, get_all)]
pub struct NexusFile {
    pub taxa: Vec<String>,
    pub tree_names: Vec<String>,
    pub tree_newicks: Vec<String>,
}

/// Parse a NEXUS format string.
#[pyfunction]
fn parse_nexus(input: &str) -> PyResult<NexusFile> {
    let nf = cyanea_phylo::nexus::parse(input).into_pyresult()?;
    let mut tree_names = Vec::new();
    let mut tree_newicks = Vec::new();
    for t in &nf.trees {
        tree_names.push(t.name.clone());
        tree_newicks.push(cyanea_phylo::write_newick(&t.tree));
    }
    Ok(NexusFile {
        taxa: nf.taxa,
        tree_names,
        tree_newicks,
    })
}

/// Write NEXUS format from taxa and Newick strings.
#[pyfunction]
fn write_nexus(taxa: Vec<String>, tree_names: Vec<String>, tree_newicks: Vec<String>) -> PyResult<String> {
    if tree_names.len() != tree_newicks.len() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "tree_names and tree_newicks must have the same length",
        ));
    }
    let mut trees = Vec::new();
    for (name, newick) in tree_names.iter().zip(tree_newicks.iter()) {
        let tree = cyanea_phylo::parse_newick(newick).into_pyresult()?;
        trees.push((name.as_str(), tree));
    }
    let refs: Vec<(&str, &cyanea_phylo::PhyloTree)> = trees.iter().map(|(n, t)| (*n, t)).collect();
    Ok(cyanea_phylo::nexus::write(&taxa, &refs))
}

// ---------------------------------------------------------------------------
// Tree comparison
// ---------------------------------------------------------------------------

/// Normalized Robinson-Foulds distance (0.0-1.0).
#[pyfunction]
fn robinson_foulds_normalized(tree1: &PhyloTree, tree2: &PhyloTree) -> PyResult<f64> {
    cyanea_phylo::robinson_foulds_normalized(&tree1.inner, &tree2.inner).into_pyresult()
}

/// Branch score distance (branch-length-aware tree distance).
#[pyfunction]
fn branch_score_distance(tree1: &PhyloTree, tree2: &PhyloTree) -> PyResult<f64> {
    cyanea_phylo::branch_score_distance(&tree1.inner, &tree2.inner).into_pyresult()
}

// ---------------------------------------------------------------------------
// Simulation classes
// ---------------------------------------------------------------------------

/// Result of simulating sequence evolution along a phylogenetic tree.
#[pyclass(frozen, get_all)]
pub struct PySimulatedAlignment {
    pub names: Vec<String>,
    pub sequences: Vec<String>,
    pub n_substitutions: usize,
}

/// A coalescent tree result (wraps a PhyloTree with metadata).
#[pyclass(frozen, get_all)]
pub struct PyCoalescentTree {
    pub newick: String,
    pub n_samples: usize,
}

/// Simulate sequence evolution along a phylogenetic tree.
///
/// Models: "jc" (Jukes-Cantor), "k2p" (Kimura 2-parameter, kappa=2.0).
#[pyfunction]
#[pyo3(signature = (tree, seq_length, model="jc", seed=42))]
fn simulate_evolution(tree: &PhyloTree, seq_length: usize, model: &str, seed: u64) -> PyResult<PySimulatedAlignment> {
    let model_obj: Box<dyn cyanea_phylo::SubstitutionModel> = match model {
        "jc" => Box::new(cyanea_phylo::Jc69Model::new()),
        "k2p" => Box::new(cyanea_phylo::Hky85Model::new(2.0, [0.25, 0.25, 0.25, 0.25]).into_pyresult()?),
        _ => return Err(pyo3::exceptions::PyValueError::new_err(format!("unknown model: {model} (expected 'jc' or 'k2p')"))),
    };
    let result = cyanea_phylo::simulate_evolution(&tree.inner, model_obj.as_ref(), seq_length, seed).into_pyresult()?;
    let sequences = result.sequences.iter().map(|s| String::from_utf8_lossy(s).to_string()).collect();
    Ok(PySimulatedAlignment { names: result.names, sequences, n_substitutions: result.n_substitutions })
}

/// Simulate a coalescent tree with constant population size.
#[pyfunction]
#[pyo3(signature = (n_samples, pop_size, seed=42))]
fn simulate_coalescent(n_samples: usize, pop_size: f64, seed: u64) -> PyResult<PyCoalescentTree> {
    let tree = cyanea_phylo::simulate_coalescent(n_samples, pop_size, seed).into_pyresult()?;
    let newick = cyanea_phylo::write_newick(&tree);
    Ok(PyCoalescentTree { newick, n_samples: tree.leaf_count() })
}

/// Simulate a coalescent tree with exponential population growth.
#[pyfunction]
#[pyo3(signature = (n_samples, pop_size, growth_rate, seed=42))]
fn simulate_coalescent_growth(n_samples: usize, pop_size: f64, growth_rate: f64, seed: u64) -> PyResult<PyCoalescentTree> {
    let tree = cyanea_phylo::simulate_coalescent_growth(n_samples, pop_size, growth_rate, seed).into_pyresult()?;
    let newick = cyanea_phylo::write_newick(&tree);
    Ok(PyCoalescentTree { newick, n_samples: tree.leaf_count() })
}

// ---------------------------------------------------------------------------
// Submodule registration
// ---------------------------------------------------------------------------

pub fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent.py(), "phylo")?;
    m.add_class::<PhyloTree>()?;
    m.add_class::<NexusFile>()?;
    m.add_class::<PySimulatedAlignment>()?;
    m.add_class::<PyCoalescentTree>()?;
    m.add_function(wrap_pyfunction!(parse_newick, &m)?)?;
    m.add_function(wrap_pyfunction!(evolutionary_distance, &m)?)?;
    m.add_function(wrap_pyfunction!(upgma, &m)?)?;
    m.add_function(wrap_pyfunction!(neighbor_joining, &m)?)?;
    m.add_function(wrap_pyfunction!(parse_nexus, &m)?)?;
    m.add_function(wrap_pyfunction!(write_nexus, &m)?)?;
    m.add_function(wrap_pyfunction!(robinson_foulds_normalized, &m)?)?;
    m.add_function(wrap_pyfunction!(branch_score_distance, &m)?)?;
    m.add_function(wrap_pyfunction!(simulate_evolution, &m)?)?;
    m.add_function(wrap_pyfunction!(simulate_coalescent, &m)?)?;
    m.add_function(wrap_pyfunction!(simulate_coalescent_growth, &m)?)?;
    parent.add_submodule(&m)?;
    Ok(())
}
