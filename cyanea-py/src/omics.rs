//! Python bindings for cyanea-omics: genome arithmetic and liftover.

use std::collections::BTreeMap;

use pyo3::prelude::*;

use crate::error::IntoPyResult;

// ---------------------------------------------------------------------------
// Strand parsing helper
// ---------------------------------------------------------------------------

fn parse_strand(s: &str) -> PyResult<cyanea_omics::Strand> {
    match s {
        "+" | "forward" => Ok(cyanea_omics::Strand::Forward),
        "-" | "reverse" => Ok(cyanea_omics::Strand::Reverse),
        "." | "unknown" | "ignore" => Ok(cyanea_omics::Strand::Unknown),
        _ => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "invalid strand: '{s}' (expected '+', '-', '.', 'forward', 'reverse', or 'unknown')"
        ))),
    }
}

fn parse_strand_mode(s: &str) -> PyResult<cyanea_omics::StrandMode> {
    match s {
        "ignore" => Ok(cyanea_omics::StrandMode::Ignore),
        "same" => Ok(cyanea_omics::StrandMode::Same),
        "opposite" => Ok(cyanea_omics::StrandMode::Opposite),
        _ => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "invalid strand_mode: '{s}' (expected 'ignore', 'same', or 'opposite')"
        ))),
    }
}

// ---------------------------------------------------------------------------
// Result classes
// ---------------------------------------------------------------------------

/// A genomic interval.
#[pyclass(frozen, get_all)]
#[derive(Clone)]
pub struct PyGenomicInterval {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: String,
}

#[pymethods]
impl PyGenomicInterval {
    #[new]
    #[pyo3(signature = (chrom, start, end, strand="."))]
    fn new(chrom: String, start: u64, end: u64, strand: &str) -> PyResult<Self> {
        // Validate the strand string
        parse_strand(strand)?;
        if start >= end {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "start ({start}) must be less than end ({end})"
            )));
        }
        Ok(Self {
            chrom,
            start,
            end,
            strand: strand.to_string(),
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "GenomicInterval('{}', {}, {}, '{}')",
            self.chrom, self.start, self.end, self.strand
        )
    }

    fn __len__(&self) -> usize {
        (self.end - self.start) as usize
    }
}

impl PyGenomicInterval {
    fn to_rust(&self) -> PyResult<cyanea_omics::GenomicInterval> {
        let strand = parse_strand(&self.strand)?;
        cyanea_omics::GenomicInterval::with_strand(&self.chrom, self.start, self.end, strand)
            .into_pyresult()
    }

    fn from_rust(iv: &cyanea_omics::GenomicInterval) -> Self {
        let strand = match iv.strand {
            cyanea_omics::Strand::Forward => "+",
            cyanea_omics::Strand::Reverse => "-",
            cyanea_omics::Strand::Unknown => ".",
        };
        Self {
            chrom: iv.chrom.clone(),
            start: iv.start,
            end: iv.end,
            strand: strand.to_string(),
        }
    }
}

fn to_rust_intervals(intervals: Vec<PyGenomicInterval>) -> PyResult<Vec<cyanea_omics::GenomicInterval>> {
    intervals.iter().map(|iv| iv.to_rust()).collect()
}

fn from_rust_intervals(intervals: &[cyanea_omics::GenomicInterval]) -> Vec<PyGenomicInterval> {
    intervals.iter().map(PyGenomicInterval::from_rust).collect()
}

fn to_genome_info(genome: std::collections::HashMap<String, u64>) -> BTreeMap<String, u64> {
    genome.into_iter().collect()
}

/// Result of a closest-interval query.
#[pyclass(frozen, get_all)]
pub struct PyClosestResult {
    pub query: PyGenomicInterval,
    pub closest: Option<PyGenomicInterval>,
    pub distance: Option<u64>,
}

/// Jaccard statistics between two interval sets.
#[pyclass(frozen, get_all)]
pub struct PyJaccardStats {
    pub intersection_bp: u64,
    pub union_bp: u64,
    pub jaccard: f64,
    pub n_intersections: u64,
}

// ---------------------------------------------------------------------------
// Liftover
// ---------------------------------------------------------------------------

/// Parsed UCSC chain file for coordinate liftover.
#[pyclass]
pub struct PyChainFile {
    inner: cyanea_omics::ChainFile,
}

/// Result of a liftover operation.
#[pyclass(frozen, get_all)]
pub struct PyLiftoverResult {
    pub status: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub fraction_mapped: f64,
}

// ---------------------------------------------------------------------------
// Genome arithmetic functions
// ---------------------------------------------------------------------------

/// Merge overlapping intervals.
#[pyfunction]
#[pyo3(signature = (intervals, strand_mode="ignore"))]
fn merge_intervals(intervals: Vec<PyGenomicInterval>, strand_mode: &str) -> PyResult<Vec<PyGenomicInterval>> {
    let ivs = to_rust_intervals(intervals)?;
    let mode = parse_strand_mode(strand_mode)?;
    let result = cyanea_omics::merge(&ivs, mode);
    Ok(from_rust_intervals(&result))
}

/// Intersect two sets of intervals.
#[pyfunction]
#[pyo3(signature = (a, b, strand_mode="ignore"))]
fn intersect_intervals(
    a: Vec<PyGenomicInterval>,
    b: Vec<PyGenomicInterval>,
    strand_mode: &str,
) -> PyResult<Vec<PyGenomicInterval>> {
    let a = to_rust_intervals(a)?;
    let b = to_rust_intervals(b)?;
    let mode = parse_strand_mode(strand_mode)?;
    let result = cyanea_omics::intersect(&a, &b, mode).into_pyresult()?;
    Ok(from_rust_intervals(&result))
}

/// Subtract intervals in b from intervals in a.
#[pyfunction]
#[pyo3(signature = (a, b, strand_mode="ignore"))]
fn subtract_intervals(
    a: Vec<PyGenomicInterval>,
    b: Vec<PyGenomicInterval>,
    strand_mode: &str,
) -> PyResult<Vec<PyGenomicInterval>> {
    let a = to_rust_intervals(a)?;
    let b = to_rust_intervals(b)?;
    let mode = parse_strand_mode(strand_mode)?;
    let result = cyanea_omics::subtract(&a, &b, mode).into_pyresult()?;
    Ok(from_rust_intervals(&result))
}

/// Complement of intervals relative to a genome.
#[pyfunction]
fn complement_intervals(
    intervals: Vec<PyGenomicInterval>,
    genome: std::collections::HashMap<String, u64>,
) -> PyResult<Vec<PyGenomicInterval>> {
    let ivs = to_rust_intervals(intervals)?;
    let gi = to_genome_info(genome);
    let result = cyanea_omics::complement(&ivs, &gi).into_pyresult()?;
    Ok(from_rust_intervals(&result))
}

/// Union of two interval sets (merge after concatenation).
#[pyfunction]
#[pyo3(signature = (a, b, strand_mode="ignore"))]
fn union_intervals(
    a: Vec<PyGenomicInterval>,
    b: Vec<PyGenomicInterval>,
    strand_mode: &str,
) -> PyResult<Vec<PyGenomicInterval>> {
    let a = to_rust_intervals(a)?;
    let b = to_rust_intervals(b)?;
    let mode = parse_strand_mode(strand_mode)?;
    let result = cyanea_omics::union(&a, &b, mode);
    Ok(from_rust_intervals(&result))
}

/// Find closest intervals in b for each interval in a.
#[pyfunction]
#[pyo3(signature = (a, b, strand_mode="ignore"))]
fn closest_intervals(
    a: Vec<PyGenomicInterval>,
    b: Vec<PyGenomicInterval>,
    strand_mode: &str,
) -> PyResult<Vec<PyClosestResult>> {
    let a = to_rust_intervals(a)?;
    let b = to_rust_intervals(b)?;
    let mode = parse_strand_mode(strand_mode)?;
    let results = cyanea_omics::closest(&a, &b, mode);
    Ok(results
        .into_iter()
        .map(|cr| PyClosestResult {
            query: PyGenomicInterval::from_rust(&cr.query),
            closest: cr.closest.as_ref().map(PyGenomicInterval::from_rust),
            distance: cr.distance,
        })
        .collect())
}

/// Jaccard similarity between two interval sets.
#[pyfunction]
fn jaccard_index(a: Vec<PyGenomicInterval>, b: Vec<PyGenomicInterval>) -> PyResult<f64> {
    let a = to_rust_intervals(a)?;
    let b = to_rust_intervals(b)?;
    Ok(cyanea_omics::jaccard(&a, &b))
}

/// Jaccard statistics between two interval sets.
#[pyfunction]
fn jaccard_stats(
    a: Vec<PyGenomicInterval>,
    b: Vec<PyGenomicInterval>,
) -> PyResult<PyJaccardStats> {
    let a = to_rust_intervals(a)?;
    let b = to_rust_intervals(b)?;
    let s = cyanea_omics::jaccard_stats(&a, &b);
    Ok(PyJaccardStats {
        intersection_bp: s.intersection_bp,
        union_bp: s.union_bp,
        jaccard: s.jaccard,
        n_intersections: s.n_intersections,
    })
}

/// Create non-overlapping windows across a genome.
#[pyfunction]
fn make_windows(
    genome: std::collections::HashMap<String, u64>,
    window_size: u64,
) -> PyResult<Vec<PyGenomicInterval>> {
    let gi = to_genome_info(genome);
    let result = cyanea_omics::make_windows(&gi, window_size).into_pyresult()?;
    Ok(from_rust_intervals(&result))
}

/// Create sliding windows across a genome.
#[pyfunction]
fn make_sliding_windows(
    genome: std::collections::HashMap<String, u64>,
    window_size: u64,
    step: u64,
) -> PyResult<Vec<PyGenomicInterval>> {
    let gi = to_genome_info(genome);
    let result = cyanea_omics::make_sliding_windows(&gi, window_size, step).into_pyresult()?;
    Ok(from_rust_intervals(&result))
}

// ---------------------------------------------------------------------------
// Liftover functions
// ---------------------------------------------------------------------------

/// Parse a UCSC chain file from text.
#[pyfunction]
fn parse_chain(chain_text: &str) -> PyResult<PyChainFile> {
    let chain = cyanea_omics::parse_chain(chain_text).into_pyresult()?;
    Ok(PyChainFile { inner: chain })
}

/// Liftover a single genomic interval using a chain file.
#[pyfunction]
#[pyo3(signature = (chain, interval, min_match=0.95))]
fn liftover(
    chain: &PyChainFile,
    interval: PyGenomicInterval,
    min_match: f64,
) -> PyResult<PyLiftoverResult> {
    let iv = interval.to_rust()?;
    let result = cyanea_omics::liftover(&chain.inner, &iv, min_match);
    Ok(liftover_result_to_py(result))
}

/// Liftover a batch of genomic intervals.
#[pyfunction]
#[pyo3(signature = (chain, intervals, min_match=0.95))]
fn liftover_batch(
    chain: &PyChainFile,
    intervals: Vec<PyGenomicInterval>,
    min_match: f64,
) -> PyResult<Vec<PyLiftoverResult>> {
    let ivs = to_rust_intervals(intervals)?;
    let results = cyanea_omics::liftover_batch(&chain.inner, &ivs, min_match);
    Ok(results.into_iter().map(liftover_result_to_py).collect())
}

fn liftover_result_to_py(r: cyanea_omics::LiftoverResult) -> PyLiftoverResult {
    match r {
        cyanea_omics::LiftoverResult::Mapped(iv) => PyLiftoverResult {
            status: "mapped".to_string(),
            chrom: iv.chrom.clone(),
            start: iv.start,
            end: iv.end,
            fraction_mapped: 1.0,
        },
        cyanea_omics::LiftoverResult::Partial { mapped, fraction_mapped } => PyLiftoverResult {
            status: "partial".to_string(),
            chrom: mapped.chrom.clone(),
            start: mapped.start,
            end: mapped.end,
            fraction_mapped,
        },
        cyanea_omics::LiftoverResult::Unmapped => PyLiftoverResult {
            status: "unmapped".to_string(),
            chrom: String::new(),
            start: 0,
            end: 0,
            fraction_mapped: 0.0,
        },
    }
}

// ---------------------------------------------------------------------------
// Variant annotation
// ---------------------------------------------------------------------------

/// A predicted variant effect on a transcript.
#[pyclass(frozen, get_all)]
pub struct PyVariantEffect {
    pub consequence: String,
    pub gene_id: String,
    pub gene_name: String,
    pub transcript_id: String,
}

/// Annotate a variant against a list of gene definitions.
///
/// `genes` should be a list of dicts with keys: gene_id, gene_name, chrom, start, end, strand.
/// This simplified binding creates a single non-coding transcript per gene.
#[pyfunction]
#[pyo3(signature = (chrom, pos, ref_allele, alt_allele, genes))]
fn annotate_variant(
    chrom: &str,
    pos: u64,
    ref_allele: &str,
    alt_allele: &str,
    genes: Vec<PyGenomicInterval>,
) -> PyResult<Vec<PyVariantEffect>> {
    let variant = cyanea_omics::Variant::new(
        chrom,
        pos,
        ref_allele.as_bytes().to_vec(),
        vec![alt_allele.as_bytes().to_vec()],
    ).into_pyresult()?;

    let gene_defs: Vec<cyanea_omics::Gene> = genes.iter().enumerate().map(|(i, iv)| {
        let strand = match iv.strand.as_str() {
            "+" | "forward" => cyanea_omics::Strand::Forward,
            "-" | "reverse" => cyanea_omics::Strand::Reverse,
            _ => cyanea_omics::Strand::Unknown,
        };
        let gene_id = format!("gene_{}", i);
        cyanea_omics::Gene {
            gene_id: gene_id.clone(),
            gene_name: gene_id,
            chrom: iv.chrom.clone(),
            start: iv.start,
            end: iv.end,
            strand,
            gene_type: cyanea_omics::GeneType::ProteinCoding,
            transcripts: vec![cyanea_omics::Transcript {
                transcript_id: format!("tx_{}", i),
                start: iv.start,
                end: iv.end,
                exons: vec![cyanea_omics::Exon {
                    exon_number: 1,
                    start: iv.start,
                    end: iv.end,
                }],
                cds_start: Some(iv.start),
                cds_end: Some(iv.end),
            }],
        }
    }).collect();

    let config = cyanea_omics::AnnotationConfig::default();
    let effects = cyanea_omics::annotate_variant(&variant, &gene_defs, &config);
    Ok(effects.into_iter().map(|e| PyVariantEffect {
        consequence: format!("{:?}", e.consequence),
        gene_id: e.gene_id,
        gene_name: e.gene_name,
        transcript_id: e.transcript_id,
    }).collect())
}

// ---------------------------------------------------------------------------
// Methylation
// ---------------------------------------------------------------------------

/// Simulate bisulfite conversion of a DNA sequence.
///
/// Converts unmethylated C to T; positions in `methylated` remain as C.
#[pyfunction]
fn bisulfite_convert(seq: &str, methylated: Vec<u64>) -> String {
    let result = cyanea_omics::bisulfite_convert(seq.as_bytes(), &methylated);
    String::from_utf8_lossy(&result).to_string()
}

/// A CpG island detected in a sequence.
#[pyclass(frozen, get_all)]
pub struct PyCpgIsland {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub cpg_count: usize,
    pub gc_content: f64,
    pub obs_exp_ratio: f64,
}

/// Find CpG islands in a DNA sequence using a sliding window approach.
#[pyfunction]
#[pyo3(signature = (seq, chrom="chr1", min_length=200, min_gc=0.5, min_obs_exp=0.6))]
fn find_cpg_islands(seq: &str, chrom: &str, min_length: u64, min_gc: f64, min_obs_exp: f64) -> Vec<PyCpgIsland> {
    let islands = cyanea_omics::find_cpg_islands(seq.as_bytes(), chrom, min_length, min_gc, min_obs_exp);
    islands.into_iter().map(|i| PyCpgIsland {
        chrom: i.chrom,
        start: i.start,
        end: i.end,
        cpg_count: i.cpg_count,
        gc_content: i.gc_content,
        obs_exp_ratio: i.obs_exp_ratio,
    }).collect()
}

// ---------------------------------------------------------------------------
// Spatial transcriptomics
// ---------------------------------------------------------------------------

/// Moran's I spatial autocorrelation result.
#[pyclass(frozen, get_all)]
pub struct PyMoransI {
    pub morans_i: f64,
    pub expected_i: f64,
    pub variance_i: f64,
    pub z_score: f64,
    pub p_value: f64,
}

/// Compute Moran's I spatial autocorrelation statistic.
///
/// `neighbors` is a list of lists of `(neighbor_index, distance)` tuples.
#[pyfunction]
fn morans_i(values: Vec<f64>, neighbors: Vec<Vec<(usize, f64)>>) -> PyResult<PyMoransI> {
    let graph = cyanea_omics::SpatialGraph {
        n_nodes: values.len(),
        neighbors,
    };
    let result = cyanea_omics::morans_i(&values, &graph).into_pyresult()?;
    Ok(PyMoransI {
        morans_i: result.morans_i,
        expected_i: result.expected_i,
        variance_i: result.variance_i,
        z_score: result.z_score,
        p_value: result.p_value,
    })
}

/// Geary's C spatial autocorrelation result.
#[pyclass(frozen, get_all)]
pub struct PyGearysC {
    pub c: f64,
    pub expected_c: f64,
    pub z_score: f64,
    pub p_value: f64,
}

/// Compute Geary's C spatial autocorrelation statistic.
///
/// `neighbors` is a list of lists of `(neighbor_index, distance)` tuples.
#[pyfunction]
fn gearys_c(values: Vec<f64>, neighbors: Vec<Vec<(usize, f64)>>) -> PyResult<PyGearysC> {
    let graph = cyanea_omics::SpatialGraph {
        n_nodes: values.len(),
        neighbors,
    };
    let result = cyanea_omics::gearys_c(&values, &graph).into_pyresult()?;
    Ok(PyGearysC {
        c: result.c,
        expected_c: result.expected_c,
        z_score: result.z_score,
        p_value: result.p_value,
    })
}

// ---------------------------------------------------------------------------
// Submodule registration
// ---------------------------------------------------------------------------

pub fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent.py(), "omics")?;
    // Interval class
    m.add_class::<PyGenomicInterval>()?;
    m.add_class::<PyClosestResult>()?;
    m.add_class::<PyJaccardStats>()?;
    m.add_class::<PyChainFile>()?;
    m.add_class::<PyLiftoverResult>()?;
    // Genome arithmetic
    m.add_function(wrap_pyfunction!(merge_intervals, &m)?)?;
    m.add_function(wrap_pyfunction!(intersect_intervals, &m)?)?;
    m.add_function(wrap_pyfunction!(subtract_intervals, &m)?)?;
    m.add_function(wrap_pyfunction!(complement_intervals, &m)?)?;
    m.add_function(wrap_pyfunction!(union_intervals, &m)?)?;
    m.add_function(wrap_pyfunction!(closest_intervals, &m)?)?;
    m.add_function(wrap_pyfunction!(jaccard_index, &m)?)?;
    m.add_function(wrap_pyfunction!(jaccard_stats, &m)?)?;
    m.add_function(wrap_pyfunction!(make_windows, &m)?)?;
    m.add_function(wrap_pyfunction!(make_sliding_windows, &m)?)?;
    // Liftover
    m.add_function(wrap_pyfunction!(parse_chain, &m)?)?;
    m.add_function(wrap_pyfunction!(liftover, &m)?)?;
    m.add_function(wrap_pyfunction!(liftover_batch, &m)?)?;
    // Variant annotation
    m.add_class::<PyVariantEffect>()?;
    m.add_function(wrap_pyfunction!(annotate_variant, &m)?)?;
    // Methylation
    m.add_class::<PyCpgIsland>()?;
    m.add_function(wrap_pyfunction!(bisulfite_convert, &m)?)?;
    m.add_function(wrap_pyfunction!(find_cpg_islands, &m)?)?;
    // Spatial
    m.add_class::<PyMoransI>()?;
    m.add_class::<PyGearysC>()?;
    m.add_function(wrap_pyfunction!(morans_i, &m)?)?;
    m.add_function(wrap_pyfunction!(gearys_c, &m)?)?;
    parent.add_submodule(&m)?;
    Ok(())
}
