//! Python bindings for cyanea-align: pairwise sequence alignment.

use pyo3::prelude::*;

use crate::error::IntoPyResult;

// ---------------------------------------------------------------------------
// Result class
// ---------------------------------------------------------------------------

/// Result of a pairwise sequence alignment.
#[pyclass(frozen, get_all)]
pub struct AlignmentResult {
    pub score: i32,
    pub aligned_query: Vec<u8>,
    pub aligned_target: Vec<u8>,
    pub query_start: usize,
    pub query_end: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub cigar_string: String,
    pub identity: f64,
    pub matches: usize,
    pub mismatches: usize,
    pub gaps: usize,
    pub length: usize,
}

impl From<cyanea_align::AlignmentResult> for AlignmentResult {
    fn from(r: cyanea_align::AlignmentResult) -> Self {
        Self {
            score: r.score,
            aligned_query: r.aligned_query.clone(),
            aligned_target: r.aligned_target.clone(),
            query_start: r.query_start,
            query_end: r.query_end,
            target_start: r.target_start,
            target_end: r.target_end,
            cigar_string: r.cigar_string(),
            identity: r.identity(),
            matches: r.matches(),
            mismatches: r.mismatches(),
            gaps: r.gaps(),
            length: r.length(),
        }
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn parse_mode(mode: &str) -> PyResult<cyanea_align::AlignmentMode> {
    match mode {
        "local" => Ok(cyanea_align::AlignmentMode::Local),
        "global" => Ok(cyanea_align::AlignmentMode::Global),
        "semiglobal" => Ok(cyanea_align::AlignmentMode::SemiGlobal),
        _ => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "unknown alignment mode: {mode} (expected 'local', 'global', or 'semiglobal')"
        ))),
    }
}

fn parse_matrix(name: &str) -> PyResult<cyanea_align::SubstitutionMatrix> {
    match name {
        "blosum62" => Ok(cyanea_align::SubstitutionMatrix::blosum62()),
        "blosum45" => Ok(cyanea_align::SubstitutionMatrix::blosum45()),
        "blosum80" => Ok(cyanea_align::SubstitutionMatrix::blosum80()),
        "pam250" => Ok(cyanea_align::SubstitutionMatrix::pam250()),
        _ => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "unknown substitution matrix: {name} (expected 'blosum62', 'blosum45', 'blosum80', or 'pam250')"
        ))),
    }
}

// ---------------------------------------------------------------------------
// Module functions
// ---------------------------------------------------------------------------

/// Align two DNA sequences.
#[pyfunction]
#[pyo3(signature = (query, target, *, mode="local", match_score=2, mismatch_score=-1, gap_open=-5, gap_extend=-2))]
fn align_dna(
    query: &[u8],
    target: &[u8],
    mode: &str,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
) -> PyResult<AlignmentResult> {
    let mode = parse_mode(mode)?;
    let matrix =
        cyanea_align::ScoringMatrix::new(match_score, mismatch_score, gap_open, gap_extend)
            .into_pyresult()?;
    let scoring = cyanea_align::ScoringScheme::Simple(matrix);
    let result = cyanea_align::align(query, target, mode, &scoring).into_pyresult()?;
    Ok(AlignmentResult::from(result))
}

/// Align two protein sequences.
#[pyfunction]
#[pyo3(signature = (query, target, *, mode="global", matrix="blosum62"))]
fn align_protein(
    query: &[u8],
    target: &[u8],
    mode: &str,
    matrix: &str,
) -> PyResult<AlignmentResult> {
    let mode = parse_mode(mode)?;
    let sub_matrix = parse_matrix(matrix)?;
    let scoring = cyanea_align::ScoringScheme::Substitution(sub_matrix);
    let result = cyanea_align::align(query, target, mode, &scoring).into_pyresult()?;
    Ok(AlignmentResult::from(result))
}

/// Batch-align a list of (query, target) DNA pairs.
#[pyfunction]
#[pyo3(signature = (pairs, *, mode="local", match_score=2, mismatch_score=-1, gap_open=-5, gap_extend=-2))]
fn align_batch(
    pairs: Vec<(Vec<u8>, Vec<u8>)>,
    mode: &str,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
) -> PyResult<Vec<AlignmentResult>> {
    let mode = parse_mode(mode)?;
    let matrix =
        cyanea_align::ScoringMatrix::new(match_score, mismatch_score, gap_open, gap_extend)
            .into_pyresult()?;
    let scoring = cyanea_align::ScoringScheme::Simple(matrix);
    let refs: Vec<(&[u8], &[u8])> = pairs
        .iter()
        .map(|(q, t)| (q.as_slice(), t.as_slice()))
        .collect();
    let results = cyanea_align::align_batch(&refs, mode, &scoring).into_pyresult()?;
    Ok(results.into_iter().map(AlignmentResult::from).collect())
}

// ---------------------------------------------------------------------------
// Submodule registration
// ---------------------------------------------------------------------------

pub fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent.py(), "align")?;
    m.add_class::<AlignmentResult>()?;
    m.add_function(wrap_pyfunction!(align_dna, &m)?)?;
    m.add_function(wrap_pyfunction!(align_protein, &m)?)?;
    m.add_function(wrap_pyfunction!(align_batch, &m)?)?;
    parent.add_submodule(&m)?;
    Ok(())
}
