//! Python bindings for cyanea-seq: sequence types and file parsing.

use pyo3::prelude::*;

use cyanea_core::{Annotated, Sequence};

use crate::error::IntoPyResult;

// ---------------------------------------------------------------------------
// DnaSequence
// ---------------------------------------------------------------------------

/// A validated DNA sequence (IUPAC alphabet).
#[pyclass(frozen)]
pub struct DnaSequence {
    inner: cyanea_seq::DnaSequence,
}

#[pymethods]
impl DnaSequence {
    #[new]
    fn new(data: &[u8]) -> PyResult<Self> {
        let inner = cyanea_seq::DnaSequence::new(data).into_pyresult()?;
        Ok(Self { inner })
    }

    /// Return the reverse complement.
    fn reverse_complement(&self) -> Self {
        Self {
            inner: self.inner.reverse_complement(),
        }
    }

    /// Transcribe DNA → RNA (T → U).
    fn transcribe(&self) -> RnaSequence {
        RnaSequence {
            inner: self.inner.transcribe(),
        }
    }

    /// GC content as a fraction in [0.0, 1.0].
    fn gc_content(&self) -> f64 {
        self.inner.gc_content()
    }

    /// Extract k-mers as a list of bytes objects.
    fn kmers(&self, k: usize) -> PyResult<Vec<Vec<u8>>> {
        let iter = self.inner.kmers(k).into_pyresult()?;
        Ok(iter.map(|kmer| kmer.to_vec()).collect())
    }

    /// Translate DNA → protein (transcribes then translates, NCBI Table 1).
    fn translate(&self) -> PyResult<ProteinSequence> {
        let inner = self.inner.translate().into_pyresult()?;
        Ok(ProteinSequence { inner })
    }

    fn __len__(&self) -> usize {
        self.inner.as_bytes().len()
    }

    fn __bytes__(&self) -> &[u8] {
        self.inner.as_bytes()
    }

    fn __str__(&self) -> String {
        String::from_utf8_lossy(self.inner.as_bytes()).into_owned()
    }

    fn __repr__(&self) -> String {
        let seq = self.inner.as_bytes();
        let preview = if seq.len() > 20 {
            format!(
                "{}...({} bp)",
                String::from_utf8_lossy(&seq[..20]),
                seq.len()
            )
        } else {
            String::from_utf8_lossy(seq).into_owned()
        };
        format!("DnaSequence('{preview}')")
    }

    fn __eq__(&self, other: &DnaSequence) -> bool {
        self.inner == other.inner
    }
}

// ---------------------------------------------------------------------------
// RnaSequence
// ---------------------------------------------------------------------------

/// A validated RNA sequence (IUPAC alphabet).
#[pyclass(frozen)]
pub struct RnaSequence {
    inner: cyanea_seq::RnaSequence,
}

#[pymethods]
impl RnaSequence {
    #[new]
    fn new(data: &[u8]) -> PyResult<Self> {
        let inner = cyanea_seq::RnaSequence::new(data).into_pyresult()?;
        Ok(Self { inner })
    }

    /// Return the reverse complement.
    fn reverse_complement(&self) -> Self {
        Self {
            inner: self.inner.reverse_complement(),
        }
    }

    /// Translate RNA → protein (NCBI Table 1).
    fn translate(&self) -> PyResult<ProteinSequence> {
        let inner = self.inner.translate().into_pyresult()?;
        Ok(ProteinSequence { inner })
    }

    /// Reverse-transcribe RNA → DNA (U → T).
    fn reverse_transcribe(&self) -> DnaSequence {
        DnaSequence {
            inner: self.inner.reverse_transcribe(),
        }
    }

    /// Extract k-mers as a list of bytes objects.
    fn kmers(&self, k: usize) -> PyResult<Vec<Vec<u8>>> {
        let iter = self.inner.kmers(k).into_pyresult()?;
        Ok(iter.map(|kmer| kmer.to_vec()).collect())
    }

    fn __len__(&self) -> usize {
        self.inner.as_bytes().len()
    }

    fn __bytes__(&self) -> &[u8] {
        self.inner.as_bytes()
    }

    fn __str__(&self) -> String {
        String::from_utf8_lossy(self.inner.as_bytes()).into_owned()
    }

    fn __repr__(&self) -> String {
        let seq = self.inner.as_bytes();
        let preview = if seq.len() > 20 {
            format!(
                "{}...({} bp)",
                String::from_utf8_lossy(&seq[..20]),
                seq.len()
            )
        } else {
            String::from_utf8_lossy(seq).into_owned()
        };
        format!("RnaSequence('{preview}')")
    }

    fn __eq__(&self, other: &RnaSequence) -> bool {
        self.inner == other.inner
    }
}

// ---------------------------------------------------------------------------
// ProteinSequence
// ---------------------------------------------------------------------------

/// A validated protein sequence.
#[pyclass(frozen)]
pub struct ProteinSequence {
    inner: cyanea_seq::ProteinSequence,
}

#[pymethods]
impl ProteinSequence {
    #[new]
    fn new(data: &[u8]) -> PyResult<Self> {
        let inner = cyanea_seq::ProteinSequence::new(data).into_pyresult()?;
        Ok(Self { inner })
    }

    /// Molecular weight in Daltons.
    fn molecular_weight(&self) -> f64 {
        self.inner.molecular_weight()
    }

    /// Extract k-mers as a list of bytes objects.
    fn kmers(&self, k: usize) -> PyResult<Vec<Vec<u8>>> {
        let iter = self.inner.kmers(k).into_pyresult()?;
        Ok(iter.map(|kmer| kmer.to_vec()).collect())
    }

    fn __len__(&self) -> usize {
        self.inner.as_bytes().len()
    }

    fn __bytes__(&self) -> &[u8] {
        self.inner.as_bytes()
    }

    fn __str__(&self) -> String {
        String::from_utf8_lossy(self.inner.as_bytes()).into_owned()
    }

    fn __repr__(&self) -> String {
        let seq = self.inner.as_bytes();
        let preview = if seq.len() > 20 {
            format!(
                "{}...({} aa)",
                String::from_utf8_lossy(&seq[..20]),
                seq.len()
            )
        } else {
            String::from_utf8_lossy(seq).into_owned()
        };
        format!("ProteinSequence('{preview}')")
    }

    fn __eq__(&self, other: &ProteinSequence) -> bool {
        self.inner == other.inner
    }
}

// ---------------------------------------------------------------------------
// FASTA/FASTQ result classes
// ---------------------------------------------------------------------------

/// Statistics from a FASTA file.
#[pyclass(frozen, get_all)]
pub struct FastaStats {
    pub sequence_count: u64,
    pub total_bases: u64,
    pub gc_content: f64,
    pub avg_length: f64,
}

/// Statistics from a FASTQ file.
#[pyclass(frozen, get_all)]
pub struct FastqStats {
    pub sequence_count: u64,
    pub total_bases: u64,
    pub gc_content: f64,
    pub avg_length: f64,
    pub mean_quality: f64,
    pub q20_fraction: f64,
    pub q30_fraction: f64,
}

/// A single FASTQ record.
#[pyclass(frozen, get_all)]
pub struct FastqRecord {
    pub name: String,
    pub description: String,
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,
}

/// A paired-end FASTQ record (R1 + R2 mate pair).
#[pyclass(frozen, get_all)]
pub struct PairedFastqRecord {
    pub r1_name: String,
    pub r1_description: String,
    pub r1_sequence: Vec<u8>,
    pub r1_quality: Vec<u8>,
    pub r2_name: String,
    pub r2_description: String,
    pub r2_sequence: Vec<u8>,
    pub r2_quality: Vec<u8>,
}

/// Statistics from paired FASTQ files.
#[pyclass(frozen, get_all)]
pub struct PairedFastqStats {
    pub pair_count: u64,
    pub r1_total_bases: u64,
    pub r1_gc_content: f64,
    pub r1_avg_length: f64,
    pub r1_mean_quality: f64,
    pub r2_total_bases: u64,
    pub r2_gc_content: f64,
    pub r2_avg_length: f64,
    pub r2_mean_quality: f64,
}

/// Report from trimming paired FASTQ files.
#[pyclass(frozen, get_all)]
pub struct PairedTrimReport {
    pub total_input: usize,
    pub both_passed: usize,
    pub r1_only_passed: usize,
    pub r2_only_passed: usize,
    pub both_failed: usize,
    pub total_bases_input: u64,
    pub total_bases_output: u64,
    pub survival_rate: f64,
}

// ---------------------------------------------------------------------------
// Module functions
// ---------------------------------------------------------------------------

/// Compute FASTA file statistics (streaming, no full load).
#[pyfunction]
fn fasta_stats(path: &str) -> PyResult<FastaStats> {
    let stats = cyanea_seq::parse_fasta_stats(path).into_pyresult()?;
    Ok(FastaStats {
        sequence_count: stats.sequence_count,
        total_bases: stats.total_bases,
        gc_content: stats.gc_content,
        avg_length: stats.avg_length,
    })
}

/// Parse a FASTQ file into a list of records.
#[pyfunction]
fn parse_fastq(path: &str) -> PyResult<Vec<FastqRecord>> {
    let records = cyanea_seq::parse_fastq_file(path).into_pyresult()?;
    Ok(records
        .into_iter()
        .map(|r| {
            let name = r.name().to_string();
            let description = r.description().unwrap_or("").to_string();
            let quality = r.quality().as_slice().to_vec();
            let sequence = r.sequence().as_bytes().to_vec();
            FastqRecord {
                name,
                description,
                sequence,
                quality,
            }
        })
        .collect())
}

/// Compute FASTQ file statistics (streaming, no full load).
#[pyfunction]
fn fastq_stats(path: &str) -> PyResult<FastqStats> {
    let stats = cyanea_seq::parse_fastq_stats(path).into_pyresult()?;
    Ok(FastqStats {
        sequence_count: stats.sequence_count,
        total_bases: stats.total_bases,
        gc_content: stats.gc_content,
        avg_length: stats.avg_length,
        mean_quality: stats.mean_quality,
        q20_fraction: stats.q20_fraction,
        q30_fraction: stats.q30_fraction,
    })
}

// ---------------------------------------------------------------------------
// Paired-end FASTQ
// ---------------------------------------------------------------------------

/// Parse paired FASTQ files into a list of paired records.
#[pyfunction]
#[pyo3(signature = (r1_path, r2_path, validation="relaxed"))]
fn parse_paired_fastq(r1_path: &str, r2_path: &str, validation: &str) -> PyResult<Vec<PairedFastqRecord>> {
    let v = match validation {
        "strict" => cyanea_seq::MateValidation::Strict,
        "relaxed" => cyanea_seq::MateValidation::Relaxed,
        "none" => cyanea_seq::MateValidation::None,
        _ => return Err(pyo3::exceptions::PyValueError::new_err(
            format!("invalid validation mode: '{}' (expected 'strict', 'relaxed', or 'none')", validation)
        )),
    };
    let pairs = cyanea_seq::parse_paired_fastq_files(r1_path, r2_path, v).into_pyresult()?;
    Ok(pairs.into_iter().map(|p| {
        let (r1, r2) = p.into_reads();
        PairedFastqRecord {
            r1_name: r1.name().to_string(),
            r1_description: r1.description().unwrap_or("").to_string(),
            r1_sequence: r1.sequence().as_bytes().to_vec(),
            r1_quality: r1.quality().as_slice().to_vec(),
            r2_name: r2.name().to_string(),
            r2_description: r2.description().unwrap_or("").to_string(),
            r2_sequence: r2.sequence().as_bytes().to_vec(),
            r2_quality: r2.quality().as_slice().to_vec(),
        }
    }).collect())
}

/// Parse an interleaved FASTQ file into paired records.
#[pyfunction]
#[pyo3(signature = (path, validation="relaxed"))]
fn parse_interleaved_fastq(path: &str, validation: &str) -> PyResult<Vec<PairedFastqRecord>> {
    let v = match validation {
        "strict" => cyanea_seq::MateValidation::Strict,
        "relaxed" => cyanea_seq::MateValidation::Relaxed,
        "none" => cyanea_seq::MateValidation::None,
        _ => return Err(pyo3::exceptions::PyValueError::new_err(
            format!("invalid validation mode: '{}' (expected 'strict', 'relaxed', or 'none')", validation)
        )),
    };
    let pairs = cyanea_seq::parse_interleaved_fastq(path, v).into_pyresult()?;
    Ok(pairs.into_iter().map(|p| {
        let (r1, r2) = p.into_reads();
        PairedFastqRecord {
            r1_name: r1.name().to_string(),
            r1_description: r1.description().unwrap_or("").to_string(),
            r1_sequence: r1.sequence().as_bytes().to_vec(),
            r1_quality: r1.quality().as_slice().to_vec(),
            r2_name: r2.name().to_string(),
            r2_description: r2.description().unwrap_or("").to_string(),
            r2_sequence: r2.sequence().as_bytes().to_vec(),
            r2_quality: r2.quality().as_slice().to_vec(),
        }
    }).collect())
}

/// Compute statistics for paired FASTQ files.
#[pyfunction]
fn paired_fastq_stats(r1_path: &str, r2_path: &str) -> PyResult<PairedFastqStats> {
    let stats = cyanea_seq::parse_paired_fastq_stats(r1_path, r2_path).into_pyresult()?;
    Ok(PairedFastqStats {
        pair_count: stats.pair_count,
        r1_total_bases: stats.r1_stats.total_bases,
        r1_gc_content: stats.r1_stats.gc_content,
        r1_avg_length: stats.r1_stats.avg_length,
        r1_mean_quality: stats.r1_stats.mean_quality,
        r2_total_bases: stats.r2_stats.total_bases,
        r2_gc_content: stats.r2_stats.gc_content,
        r2_avg_length: stats.r2_stats.avg_length,
        r2_mean_quality: stats.r2_stats.mean_quality,
    })
}

/// Trim paired FASTQ files and return a report.
#[pyfunction]
#[pyo3(signature = (r1_path, r2_path, min_quality=20.0, window_size=4, min_length=50, max_length=None, adapters=None, orphan_policy="drop_both"))]
fn trim_paired_fastq(
    r1_path: &str,
    r2_path: &str,
    min_quality: f64,
    window_size: usize,
    min_length: usize,
    max_length: Option<usize>,
    adapters: Option<Vec<String>>,
    orphan_policy: &str,
) -> PyResult<PairedTrimReport> {
    let _ = orphan_policy; // reserved for future use

    // Build pipeline
    let mut pipeline = cyanea_seq::trim::TrimPipeline::new()
        .sliding_window(window_size, min_quality)
        .min_length(min_length);
    if let Some(max) = max_length {
        pipeline = pipeline.max_length(max);
    }
    if let Some(adapter_list) = adapters {
        for a in &adapter_list {
            pipeline = pipeline.adapter(a.as_bytes());
        }
    }

    // Parse pairs
    let pairs = cyanea_seq::parse_paired_fastq_files(r1_path, r2_path, cyanea_seq::MateValidation::None)
        .into_pyresult()?;

    // Process
    let report = pipeline.process_paired_batch_with_stats(&pairs);

    Ok(PairedTrimReport {
        total_input: report.total_input,
        both_passed: report.both_passed,
        r1_only_passed: report.r1_only_passed,
        r2_only_passed: report.r2_only_passed,
        both_failed: report.both_failed,
        total_bases_input: report.total_bases_input,
        total_bases_output: report.total_bases_output,
        survival_rate: report.survival_rate(),
    })
}

// ---------------------------------------------------------------------------
// Pattern matching
// ---------------------------------------------------------------------------

/// Boyer-Moore-Horspool exact pattern search. Returns match positions.
#[pyfunction]
fn horspool(text: &[u8], pattern: &[u8]) -> Vec<usize> {
    cyanea_seq::pattern::horspool(text, pattern)
}

/// Knuth-Morris-Pratt exact pattern search. Returns match positions.
#[pyfunction]
fn kmp(text: &[u8], pattern: &[u8]) -> Vec<usize> {
    cyanea_seq::pattern::kmp(text, pattern)
}

/// Shift-And bitparallel exact matching (patterns ≤ 64 chars).
#[pyfunction]
fn shift_and(text: &[u8], pattern: &[u8]) -> Vec<usize> {
    cyanea_seq::pattern::shift_and(text, pattern)
}

/// BNDM bitparallel exact matching (patterns ≤ 64 chars).
#[pyfunction]
fn bndm(text: &[u8], pattern: &[u8]) -> Vec<usize> {
    cyanea_seq::pattern::bndm(text, pattern)
}

/// Backward Oracle Matching exact pattern search.
#[pyfunction]
fn bom(text: &[u8], pattern: &[u8]) -> Vec<usize> {
    cyanea_seq::pattern::bom(text, pattern)
}

/// Myers bit-parallel approximate matching (patterns ≤ 64 chars).
///
/// Returns list of (end_position, edit_distance) tuples.
#[pyfunction]
fn myers_search(text: &[u8], pattern: &[u8], max_dist: usize) -> Vec<(usize, usize)> {
    cyanea_seq::pattern::myers_bitparallel(text, pattern, max_dist)
}

/// Ukkonen cut-off approximate matching.
///
/// Returns list of (end_position, edit_distance) tuples.
#[pyfunction]
fn ukkonen_search(text: &[u8], pattern: &[u8], max_dist: usize) -> Vec<(usize, usize)> {
    cyanea_seq::pattern::ukkonen(text, pattern, max_dist)
}

// ---------------------------------------------------------------------------
// ORF finder
// ---------------------------------------------------------------------------

/// An open reading frame found in a sequence.
#[pyclass(frozen, get_all)]
pub struct OrfResult {
    pub start: usize,
    pub end: usize,
    pub frame: usize,
    pub strand: String,
    pub sequence: Vec<u8>,
}

/// Find open reading frames on the forward strand.
#[pyfunction]
#[pyo3(signature = (seq, min_length=100))]
fn find_orfs(seq: &[u8], min_length: usize) -> Vec<OrfResult> {
    cyanea_seq::find_orfs(seq, min_length)
        .into_iter()
        .map(|o| OrfResult {
            start: o.start,
            end: o.end,
            frame: o.frame,
            strand: format!("{:?}", o.strand),
            sequence: o.sequence,
        })
        .collect()
}

/// Find open reading frames on both strands.
#[pyfunction]
#[pyo3(signature = (seq, min_length=100))]
fn find_orfs_both_strands(seq: &[u8], min_length: usize) -> Vec<OrfResult> {
    cyanea_seq::find_orfs_both_strands(seq, min_length)
        .into_iter()
        .map(|o| OrfResult {
            start: o.start,
            end: o.end,
            frame: o.frame,
            strand: format!("{:?}", o.strand),
            sequence: o.sequence,
        })
        .collect()
}

// ---------------------------------------------------------------------------
// FM-Index
// ---------------------------------------------------------------------------

/// An FM-Index for fast substring search on a text.
#[pyclass]
pub struct FmIndex {
    inner: cyanea_seq::FmIndex,
}

#[pymethods]
impl FmIndex {
    /// Build an FM-Index from the given text.
    #[new]
    fn new(text: &[u8]) -> Self {
        Self {
            inner: cyanea_seq::FmIndex::build(text),
        }
    }

    /// Return all positions where pattern occurs in the indexed text.
    fn search(&self, pattern: &[u8]) -> Vec<usize> {
        self.inner.search(pattern)
    }

    /// Count occurrences of pattern without returning positions.
    fn count(&self, pattern: &[u8]) -> usize {
        self.inner.count(pattern)
    }

    fn __repr__(&self) -> String {
        "FmIndex(...)".to_string()
    }
}

// ---------------------------------------------------------------------------
// MinHash sketching
// ---------------------------------------------------------------------------

/// A MinHash sketch for sequence similarity estimation.
#[pyclass]
pub struct MinHashSketch {
    inner: cyanea_seq::minhash::MinHash,
}

#[pymethods]
impl MinHashSketch {
    /// Create a MinHash sketch from a sequence.
    #[new]
    #[pyo3(signature = (seq, k=21, sketch_size=1000))]
    fn new(seq: &[u8], k: usize, sketch_size: usize) -> PyResult<Self> {
        let inner =
            cyanea_seq::minhash::MinHash::from_sequence(seq, k, sketch_size).into_pyresult()?;
        Ok(Self { inner })
    }

    /// Jaccard similarity with another sketch.
    fn jaccard(&self, other: &MinHashSketch) -> PyResult<f64> {
        self.inner.jaccard(&other.inner).into_pyresult()
    }

    /// Containment of this sketch in another.
    fn containment(&self, other: &MinHashSketch) -> PyResult<f64> {
        self.inner.containment(&other.inner).into_pyresult()
    }

    /// Average nucleotide identity estimate.
    fn ani(&self, other: &MinHashSketch) -> PyResult<f64> {
        self.inner.ani(&other.inner).into_pyresult()
    }

    /// The raw hash values in the sketch.
    fn hashes(&self) -> Vec<u64> {
        self.inner.hashes().to_vec()
    }

    /// Number of hashes in the sketch.
    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "MinHashSketch(k={}, sketch_size={})",
            self.inner.k(),
            self.inner.sketch_size()
        )
    }
}

// ---------------------------------------------------------------------------
// RNA structure prediction
// ---------------------------------------------------------------------------

/// RNA secondary structure prediction result.
#[pyclass(frozen, get_all)]
pub struct PyRnaStructure {
    pub structure: String,
    pub pairs: Vec<(usize, usize)>,
    pub max_pairs: usize,
}

/// Predict RNA secondary structure using the Nussinov algorithm (maximize base pairs).
#[pyfunction]
#[pyo3(signature = (seq, min_loop_size=3))]
fn rna_nussinov(seq: &str, min_loop_size: usize) -> PyResult<PyRnaStructure> {
    let result = cyanea_seq::nussinov(seq.as_bytes(), min_loop_size).into_pyresult()?;
    let pairs = result.structure.base_pairs();
    Ok(PyRnaStructure {
        structure: result.structure.to_dot_bracket(),
        pairs,
        max_pairs: result.max_pairs,
    })
}

/// RNA minimum free energy structure prediction result.
#[pyclass(frozen, get_all)]
pub struct PyRnaMfe {
    pub structure: String,
    pub pairs: Vec<(usize, usize)>,
    pub energy: f64,
}

/// Predict RNA secondary structure using the Zuker MFE algorithm (minimize free energy).
#[pyfunction]
fn rna_zuker(seq: &str) -> PyResult<PyRnaMfe> {
    let result = cyanea_seq::zuker_mfe(seq.as_bytes()).into_pyresult()?;
    let pairs = result.structure.base_pairs();
    Ok(PyRnaMfe {
        structure: result.structure.to_dot_bracket(),
        pairs,
        energy: result.energy,
    })
}

// ---------------------------------------------------------------------------
// Protein properties
// ---------------------------------------------------------------------------

/// Computed protein sequence properties.
#[pyclass(frozen, get_all)]
pub struct PyProteinProperties {
    pub gravy: f64,
    pub isoelectric_point: f64,
    pub extinction_reduced: f64,
    pub extinction_cystine: f64,
}

/// Compute physicochemical properties of a protein sequence.
#[pyfunction]
fn protein_properties(seq: &str) -> PyResult<PyProteinProperties> {
    let bytes = seq.as_bytes();
    let gravy = cyanea_seq::protein_properties::gravy(bytes).into_pyresult()?;
    let pi = cyanea_seq::protein_properties::isoelectric_point(bytes).into_pyresult()?;
    let ext = cyanea_seq::protein_properties::extinction_coefficient(bytes).into_pyresult()?;
    Ok(PyProteinProperties {
        gravy,
        isoelectric_point: pi,
        extinction_reduced: ext.reduced,
        extinction_cystine: ext.cystine,
    })
}

// ---------------------------------------------------------------------------
// Read simulation
// ---------------------------------------------------------------------------

/// A single simulated read.
#[pyclass(frozen, get_all)]
pub struct PySimulatedRead {
    pub name: String,
    pub sequence: String,
    pub quality: String,
    pub position: u64,
    pub is_read1: bool,
}

/// Simulate Illumina-style reads from a reference sequence.
#[pyfunction]
#[pyo3(signature = (ref_seq, read_length=150, coverage=30.0, error_rate=0.001, seed=42))]
fn simulate_reads(ref_seq: &str, read_length: usize, coverage: f64, error_rate: f64, seed: u64) -> PyResult<Vec<PySimulatedRead>> {
    let config = cyanea_seq::ReadSimConfig {
        read_length,
        coverage,
        error_rate,
        seed,
        ..Default::default()
    };
    let reads = cyanea_seq::simulate_reads(ref_seq.as_bytes(), "ref", &config);
    Ok(reads.into_iter().map(|r| PySimulatedRead {
        name: r.name,
        sequence: String::from_utf8_lossy(&r.sequence).to_string(),
        quality: String::from_utf8_lossy(&r.quality).to_string(),
        position: r.true_position,
        is_read1: r.is_read1,
    }).collect())
}

// ---------------------------------------------------------------------------
// Codon usage
// ---------------------------------------------------------------------------

/// Compute codon usage frequencies from a coding sequence.
///
/// Returns a dict mapping codon strings to counts.
#[pyfunction]
fn codon_usage(seq: &str) -> std::collections::HashMap<String, u64> {
    let usage = cyanea_seq::codon::CodonUsage::from_sequence(seq.as_bytes());
    let nucleotides = [b'A', b'C', b'G', b'T'];
    let mut result = std::collections::HashMap::new();
    for &a in &nucleotides {
        for &b in &nucleotides {
            for &c in &nucleotides {
                let codon = [a, b, c];
                let count = usage.count(&codon);
                if count > 0 {
                    result.insert(String::from_utf8_lossy(&codon).to_string(), count);
                }
            }
        }
    }
    result
}

// ---------------------------------------------------------------------------
// Assembly statistics
// ---------------------------------------------------------------------------

/// Assembly quality statistics.
#[pyclass(frozen, get_all)]
pub struct PyAssemblyStats {
    pub n_contigs: usize,
    pub total_length: usize,
    pub n50: usize,
    pub l50: usize,
    pub n90: usize,
    pub l90: usize,
    pub largest_contig: usize,
    pub smallest_contig: usize,
    pub gc_content: f64,
}

/// Compute assembly quality statistics from contig sequences.
#[pyfunction]
fn assembly_stats(contigs: Vec<Vec<u8>>) -> PyResult<PyAssemblyStats> {
    let refs: Vec<&[u8]> = contigs.iter().map(|c| c.as_slice()).collect();
    let stats = cyanea_seq::assembly_stats(&refs).into_pyresult()?;
    Ok(PyAssemblyStats {
        n_contigs: stats.n_contigs,
        total_length: stats.total_length,
        n50: stats.n50,
        l50: stats.l50,
        n90: stats.n90,
        l90: stats.l90,
        largest_contig: stats.largest_contig,
        smallest_contig: stats.smallest_contig,
        gc_content: stats.gc_content,
    })
}

// ---------------------------------------------------------------------------
// Submodule registration
// ---------------------------------------------------------------------------

pub fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent.py(), "seq")?;
    m.add_class::<DnaSequence>()?;
    m.add_class::<RnaSequence>()?;
    m.add_class::<ProteinSequence>()?;
    m.add_class::<FastaStats>()?;
    m.add_class::<FastqStats>()?;
    m.add_class::<FastqRecord>()?;
    m.add_class::<PairedFastqRecord>()?;
    m.add_class::<PairedFastqStats>()?;
    m.add_class::<PairedTrimReport>()?;
    m.add_class::<OrfResult>()?;
    m.add_class::<FmIndex>()?;
    m.add_class::<MinHashSketch>()?;
    m.add_class::<PyRnaStructure>()?;
    m.add_class::<PyRnaMfe>()?;
    m.add_class::<PyProteinProperties>()?;
    m.add_class::<PySimulatedRead>()?;
    m.add_class::<PyAssemblyStats>()?;
    m.add_function(wrap_pyfunction!(fasta_stats, &m)?)?;
    m.add_function(wrap_pyfunction!(parse_fastq, &m)?)?;
    m.add_function(wrap_pyfunction!(fastq_stats, &m)?)?;
    m.add_function(wrap_pyfunction!(parse_paired_fastq, &m)?)?;
    m.add_function(wrap_pyfunction!(parse_interleaved_fastq, &m)?)?;
    m.add_function(wrap_pyfunction!(paired_fastq_stats, &m)?)?;
    m.add_function(wrap_pyfunction!(trim_paired_fastq, &m)?)?;
    m.add_function(wrap_pyfunction!(horspool, &m)?)?;
    m.add_function(wrap_pyfunction!(kmp, &m)?)?;
    m.add_function(wrap_pyfunction!(shift_and, &m)?)?;
    m.add_function(wrap_pyfunction!(bndm, &m)?)?;
    m.add_function(wrap_pyfunction!(bom, &m)?)?;
    m.add_function(wrap_pyfunction!(myers_search, &m)?)?;
    m.add_function(wrap_pyfunction!(ukkonen_search, &m)?)?;
    m.add_function(wrap_pyfunction!(find_orfs, &m)?)?;
    m.add_function(wrap_pyfunction!(find_orfs_both_strands, &m)?)?;
    m.add_function(wrap_pyfunction!(rna_nussinov, &m)?)?;
    m.add_function(wrap_pyfunction!(rna_zuker, &m)?)?;
    m.add_function(wrap_pyfunction!(protein_properties, &m)?)?;
    m.add_function(wrap_pyfunction!(simulate_reads, &m)?)?;
    m.add_function(wrap_pyfunction!(codon_usage, &m)?)?;
    m.add_function(wrap_pyfunction!(assembly_stats, &m)?)?;
    parent.add_submodule(&m)?;
    Ok(())
}
