//! Python bindings for cyanea-io: bioinformatics file format parsing.

use pyo3::prelude::*;

use crate::error::IntoPyResult;

// ---------------------------------------------------------------------------
// Result classes
// ---------------------------------------------------------------------------

/// Metadata about a CSV file.
#[pyclass(frozen, get_all)]
pub struct CsvInfo {
    pub row_count: u64,
    pub column_count: usize,
    pub columns: Vec<String>,
    pub has_headers: bool,
}

/// Summary statistics for a VCF file.
#[pyclass(frozen, get_all)]
pub struct VcfStats {
    pub variant_count: u64,
    pub snv_count: u64,
    pub indel_count: u64,
    pub pass_count: u64,
    pub chromosomes: Vec<String>,
}

/// Summary statistics for a BED file.
#[pyclass(frozen, get_all)]
pub struct BedStats {
    pub record_count: u64,
    pub total_bases: u64,
    pub chromosomes: Vec<String>,
}

/// Summary statistics for a GFF3 file.
#[pyclass(frozen, get_all)]
pub struct Gff3Stats {
    pub gene_count: u64,
    pub transcript_count: u64,
    pub exon_count: u64,
    pub protein_coding_count: u64,
    pub chromosomes: Vec<String>,
}

/// A single SAM/BAM alignment record.
#[pyclass(frozen, get_all)]
pub struct PySamRecord {
    pub qname: String,
    pub flag: u16,
    pub rname: String,
    pub pos: u64,
    pub mapq: u8,
    pub cigar: String,
    pub sequence: String,
    pub quality: String,
}

impl From<cyanea_io::SamRecord> for PySamRecord {
    fn from(r: cyanea_io::SamRecord) -> Self {
        Self {
            qname: r.qname,
            flag: r.flag,
            rname: r.rname,
            pos: r.pos,
            mapq: r.mapq,
            cigar: r.cigar,
            sequence: r.sequence,
            quality: r.quality,
        }
    }
}

/// Summary statistics for a SAM/BAM file.
#[pyclass(frozen, get_all)]
pub struct PySamStats {
    pub total_reads: usize,
    pub mapped: usize,
    pub unmapped: usize,
    pub avg_mapq: f64,
    pub avg_length: f64,
    pub mapq_distribution: Vec<(u8, usize)>,
}

impl From<cyanea_io::SamStats> for PySamStats {
    fn from(s: cyanea_io::SamStats) -> Self {
        Self {
            total_reads: s.total_reads,
            mapped: s.mapped,
            unmapped: s.unmapped,
            avg_mapq: s.avg_mapq,
            avg_length: s.avg_length,
            mapq_distribution: s.mapq_distribution,
        }
    }
}

// ---------------------------------------------------------------------------
// Module functions
// ---------------------------------------------------------------------------

/// Parse a CSV file and return its metadata.
#[pyfunction]
fn csv_info(path: &str) -> PyResult<CsvInfo> {
    let info = cyanea_io::parse_csv_info(path).into_pyresult()?;
    Ok(CsvInfo {
        row_count: info.row_count,
        column_count: info.column_count,
        columns: info.columns,
        has_headers: info.has_headers,
    })
}

/// Parse a VCF file and return summary statistics.
#[pyfunction]
fn vcf_stats(path: &str) -> PyResult<VcfStats> {
    let stats = cyanea_io::vcf_stats(path).into_pyresult()?;
    Ok(VcfStats {
        variant_count: stats.variant_count,
        snv_count: stats.snv_count,
        indel_count: stats.indel_count,
        pass_count: stats.pass_count,
        chromosomes: stats.chromosomes,
    })
}

/// Parse a BED file and return summary statistics.
#[pyfunction]
fn bed_stats(path: &str) -> PyResult<BedStats> {
    let stats = cyanea_io::bed_stats(path).into_pyresult()?;
    Ok(BedStats {
        record_count: stats.record_count,
        total_bases: stats.total_bases,
        chromosomes: stats.chromosomes,
    })
}

/// Parse a GFF3 file and return summary statistics.
#[pyfunction]
fn gff3_stats(path: &str) -> PyResult<Gff3Stats> {
    let stats = cyanea_io::gff3_stats(path).into_pyresult()?;
    Ok(Gff3Stats {
        gene_count: stats.gene_count,
        transcript_count: stats.transcript_count,
        exon_count: stats.exon_count,
        protein_coding_count: stats.protein_coding_count,
        chromosomes: stats.chromosomes,
    })
}

/// Parse a SAM file and return all alignment records.
#[pyfunction]
fn parse_sam(path: &str) -> PyResult<Vec<PySamRecord>> {
    let records = cyanea_io::parse_sam(path).into_pyresult()?;
    Ok(records.into_iter().map(PySamRecord::from).collect())
}

/// Compute summary statistics for a SAM file.
#[pyfunction]
fn sam_stats(path: &str) -> PyResult<PySamStats> {
    let stats = cyanea_io::sam_stats_from_path(path).into_pyresult()?;
    Ok(PySamStats::from(stats))
}

/// Parse a BAM file and return all alignment records.
#[pyfunction]
fn parse_bam(path: &str) -> PyResult<Vec<PySamRecord>> {
    let records = cyanea_io::parse_bam(path).into_pyresult()?;
    Ok(records.into_iter().map(PySamRecord::from).collect())
}

/// Compute summary statistics for a BAM file.
#[pyfunction]
fn bam_stats(path: &str) -> PyResult<PySamStats> {
    let stats = cyanea_io::bam_stats(path).into_pyresult()?;
    Ok(PySamStats::from(stats))
}

// ---------------------------------------------------------------------------
// Submodule registration
// ---------------------------------------------------------------------------

pub fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent.py(), "io")?;
    m.add_class::<CsvInfo>()?;
    m.add_class::<VcfStats>()?;
    m.add_class::<BedStats>()?;
    m.add_class::<Gff3Stats>()?;
    m.add_class::<PySamRecord>()?;
    m.add_class::<PySamStats>()?;
    m.add_function(wrap_pyfunction!(csv_info, &m)?)?;
    m.add_function(wrap_pyfunction!(vcf_stats, &m)?)?;
    m.add_function(wrap_pyfunction!(bed_stats, &m)?)?;
    m.add_function(wrap_pyfunction!(gff3_stats, &m)?)?;
    m.add_function(wrap_pyfunction!(parse_sam, &m)?)?;
    m.add_function(wrap_pyfunction!(sam_stats, &m)?)?;
    m.add_function(wrap_pyfunction!(parse_bam, &m)?)?;
    m.add_function(wrap_pyfunction!(bam_stats, &m)?)?;
    parent.add_submodule(&m)?;
    Ok(())
}
