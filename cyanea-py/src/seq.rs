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
    m.add_function(wrap_pyfunction!(fasta_stats, &m)?)?;
    m.add_function(wrap_pyfunction!(parse_fastq, &m)?)?;
    m.add_function(wrap_pyfunction!(fastq_stats, &m)?)?;
    parent.add_submodule(&m)?;
    Ok(())
}
