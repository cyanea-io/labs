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
// Banded alignment
// ---------------------------------------------------------------------------

/// Banded DNA alignment (restricts DP to diagonal band for speed).
#[pyfunction]
#[pyo3(signature = (query, target, *, mode="global", bandwidth=50, match_score=2, mismatch_score=-1, gap_open=-5, gap_extend=-2))]
fn align_banded(
    query: &[u8],
    target: &[u8],
    mode: &str,
    bandwidth: usize,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
) -> PyResult<AlignmentResult> {
    let matrix =
        cyanea_align::ScoringMatrix::new(match_score, mismatch_score, gap_open, gap_extend)
            .into_pyresult()?;
    let scoring = cyanea_align::ScoringScheme::Simple(matrix);
    let result = match mode {
        "global" => cyanea_align::simd::banded_nw(query, target, &scoring, bandwidth),
        "local" => cyanea_align::simd::banded_sw(query, target, &scoring, bandwidth),
        "semiglobal" => cyanea_align::simd::banded_semi_global(query, target, &scoring, bandwidth),
        _ => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "unknown alignment mode: {mode} (expected 'local', 'global', or 'semiglobal')"
            )))
        }
    };
    Ok(AlignmentResult::from(result.into_pyresult()?))
}

// ---------------------------------------------------------------------------
// MSA
// ---------------------------------------------------------------------------

/// Multiple sequence alignment result.
#[pyclass(frozen, get_all)]
pub struct MsaResult {
    pub aligned: Vec<Vec<u8>>,
    pub n_columns: usize,
}

/// Progressive multiple sequence alignment.
#[pyfunction]
#[pyo3(signature = (sequences, *, mode="dna", match_score=2, mismatch_score=-1, gap_open=-5, gap_extend=-2))]
fn progressive_msa(
    sequences: Vec<Vec<u8>>,
    mode: &str,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
) -> PyResult<MsaResult> {
    let refs: Vec<&[u8]> = sequences.iter().map(|s| s.as_slice()).collect();
    let scoring = match mode {
        "dna" => {
            let matrix = cyanea_align::ScoringMatrix::new(match_score, mismatch_score, gap_open, gap_extend)
                .into_pyresult()?;
            cyanea_align::ScoringScheme::Simple(matrix)
        }
        "protein" => cyanea_align::ScoringScheme::Substitution(cyanea_align::SubstitutionMatrix::blosum62()),
        _ => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "unknown MSA mode: {mode} (expected 'dna' or 'protein')"
            )))
        }
    };
    let result = cyanea_align::msa::progressive_msa(&refs, &scoring).into_pyresult()?;
    Ok(MsaResult {
        aligned: result.aligned,
        n_columns: result.n_columns,
    })
}

/// POA (Partial Order Alignment) consensus from multiple sequences.
#[pyfunction]
fn poa_consensus(sequences: Vec<Vec<u8>>) -> PyResult<Vec<u8>> {
    if sequences.len() < 2 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "POA requires at least 2 sequences",
        ));
    }
    let scoring = cyanea_align::poa::PoaScoring::default();
    let mut graph = cyanea_align::poa::PoaGraph::from_sequence(&sequences[0]);
    for seq in &sequences[1..] {
        graph.add_sequence(seq, &scoring).into_pyresult()?;
    }
    Ok(graph.consensus())
}

// ---------------------------------------------------------------------------
// CIGAR utilities
// ---------------------------------------------------------------------------

/// CIGAR statistics result.
#[pyclass(frozen, get_all)]
pub struct CigarStats {
    pub cigar_string: String,
    pub reference_consumed: usize,
    pub query_consumed: usize,
    pub alignment_columns: usize,
    pub identity: f64,
    pub gap_count: usize,
    pub gap_bases: usize,
    pub soft_clipped: usize,
    pub hard_clipped: usize,
}

/// Parse a SAM CIGAR string into a list of (op_char, length) tuples.
///
/// Accepts the full SAM alphabet (M, I, D, N, S, H, P, =, X) and ``*``.
///
/// >>> parse_cigar("10M3I4D")
/// [('M', 10), ('I', 3), ('D', 4)]
#[pyfunction]
fn parse_cigar(cigar: &str) -> PyResult<Vec<(char, usize)>> {
    let ops = cyanea_align::cigar::parse_cigar(cigar).into_pyresult()?;
    Ok(ops.iter().map(|op| (op.code(), op.len())).collect())
}

/// Validate a CIGAR string against SAM spec rules.
///
/// Returns ``None`` on success; raises ``ValueError`` on invalid CIGAR.
#[pyfunction]
fn validate_cigar(cigar: &str) -> PyResult<()> {
    let ops = cyanea_align::cigar::parse_cigar(cigar).into_pyresult()?;
    cyanea_align::cigar::validate_cigar(&ops).into_pyresult()
}

/// Compute statistics from a CIGAR string.
///
/// Returns a :class:`CigarStats` object with reference/query consumed,
/// identity, gap counts, and clipping totals.
#[pyfunction]
fn cigar_stats(cigar: &str) -> PyResult<CigarStats> {
    let ops = cyanea_align::cigar::parse_cigar(cigar).into_pyresult()?;
    let (soft, hard) = cyanea_align::cigar::clipped_bases(&ops);
    Ok(CigarStats {
        cigar_string: cyanea_align::cigar::cigar_string(&ops),
        reference_consumed: cyanea_align::cigar::reference_consumed(&ops),
        query_consumed: cyanea_align::cigar::query_consumed(&ops),
        alignment_columns: cyanea_align::cigar::alignment_columns(&ops),
        identity: cyanea_align::cigar::identity(&ops),
        gap_count: cyanea_align::cigar::gap_count(&ops),
        gap_bases: cyanea_align::cigar::gap_bases(&ops),
        soft_clipped: soft,
        hard_clipped: hard,
    })
}

/// Reconstruct gapped alignment from CIGAR and ungapped sequences.
///
/// Returns a tuple of ``(aligned_query, aligned_target)`` byte sequences.
#[pyfunction]
fn cigar_to_alignment(cigar: &str, query: &[u8], target: &[u8]) -> PyResult<(Vec<u8>, Vec<u8>)> {
    let ops = cyanea_align::cigar::parse_cigar(cigar).into_pyresult()?;
    cyanea_align::cigar::cigar_to_alignment(&ops, query, target).into_pyresult()
}

/// Extract a CIGAR string from a gapped alignment (using =/X distinction).
///
/// Both sequences must have the same length, with ``b'-'`` for gaps.
#[pyfunction]
fn alignment_to_cigar(query: &[u8], target: &[u8]) -> PyResult<String> {
    let ops = cyanea_align::cigar::alignment_to_cigar(query, target).into_pyresult()?;
    Ok(cyanea_align::cigar::cigar_string(&ops))
}

/// Generate a SAM MD:Z tag from CIGAR and ungapped sequences.
#[pyfunction]
fn generate_md_tag(cigar: &str, query: &[u8], reference: &[u8]) -> PyResult<String> {
    let ops = cyanea_align::cigar::parse_cigar(cigar).into_pyresult()?;
    cyanea_align::cigar::generate_md_tag(&ops, query, reference).into_pyresult()
}

/// Merge adjacent same-type CIGAR operations.
#[pyfunction]
fn merge_cigar(cigar: &str) -> PyResult<String> {
    let ops = cyanea_align::cigar::parse_cigar(cigar).into_pyresult()?;
    Ok(cyanea_align::cigar::cigar_string(
        &cyanea_align::cigar::merge_adjacent(&ops),
    ))
}

/// Reverse CIGAR operation order.
#[pyfunction]
fn reverse_cigar(cigar: &str) -> PyResult<String> {
    let ops = cyanea_align::cigar::parse_cigar(cigar).into_pyresult()?;
    Ok(cyanea_align::cigar::cigar_string(
        &cyanea_align::cigar::reverse_cigar(&ops),
    ))
}

/// Collapse =/X operations into M (alignment match).
#[pyfunction]
fn collapse_cigar(cigar: &str) -> PyResult<String> {
    let ops = cyanea_align::cigar::parse_cigar(cigar).into_pyresult()?;
    Ok(cyanea_align::cigar::cigar_string(
        &cyanea_align::cigar::collapse_matches(&ops),
    ))
}

/// Convert hard clips (H) to soft clips (S).
#[pyfunction]
fn hard_clip_to_soft(cigar: &str) -> PyResult<String> {
    let ops = cyanea_align::cigar::parse_cigar(cigar).into_pyresult()?;
    Ok(cyanea_align::cigar::cigar_string(
        &cyanea_align::cigar::hard_clip_to_soft(&ops),
    ))
}

/// Split CIGAR at a reference coordinate.
///
/// Returns a tuple of ``(left_cigar, right_cigar)`` strings.
#[pyfunction]
fn split_cigar(cigar: &str, ref_pos: usize) -> PyResult<(String, String)> {
    let ops = cyanea_align::cigar::parse_cigar(cigar).into_pyresult()?;
    let (left, right) = cyanea_align::cigar::split_at_reference(&ops, ref_pos);
    Ok((
        cyanea_align::cigar::cigar_string(&left),
        cyanea_align::cigar::cigar_string(&right),
    ))
}

// ---------------------------------------------------------------------------
// Submodule registration
// ---------------------------------------------------------------------------

pub fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent.py(), "align")?;
    m.add_class::<AlignmentResult>()?;
    m.add_class::<MsaResult>()?;
    m.add_class::<CigarStats>()?;
    m.add_function(wrap_pyfunction!(align_dna, &m)?)?;
    m.add_function(wrap_pyfunction!(align_protein, &m)?)?;
    m.add_function(wrap_pyfunction!(align_batch, &m)?)?;
    m.add_function(wrap_pyfunction!(align_banded, &m)?)?;
    m.add_function(wrap_pyfunction!(progressive_msa, &m)?)?;
    m.add_function(wrap_pyfunction!(poa_consensus, &m)?)?;
    m.add_function(wrap_pyfunction!(parse_cigar, &m)?)?;
    m.add_function(wrap_pyfunction!(validate_cigar, &m)?)?;
    m.add_function(wrap_pyfunction!(cigar_stats, &m)?)?;
    m.add_function(wrap_pyfunction!(cigar_to_alignment, &m)?)?;
    m.add_function(wrap_pyfunction!(alignment_to_cigar, &m)?)?;
    m.add_function(wrap_pyfunction!(generate_md_tag, &m)?)?;
    m.add_function(wrap_pyfunction!(merge_cigar, &m)?)?;
    m.add_function(wrap_pyfunction!(reverse_cigar, &m)?)?;
    m.add_function(wrap_pyfunction!(collapse_cigar, &m)?)?;
    m.add_function(wrap_pyfunction!(hard_clip_to_soft, &m)?)?;
    m.add_function(wrap_pyfunction!(split_cigar, &m)?)?;
    parent.add_submodule(&m)?;
    Ok(())
}
