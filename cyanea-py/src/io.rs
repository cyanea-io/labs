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
// Full record parsing
// ---------------------------------------------------------------------------

/// A VCF variant record.
#[pyclass(frozen, get_all)]
pub struct PyVcfRecord {
    pub chrom: String,
    pub position: u64,
    pub id: String,
    pub ref_allele: String,
    pub alt_alleles: Vec<String>,
    pub quality: Option<f64>,
    pub filter: String,
}

/// Parse a VCF file and return all variant records.
#[pyfunction]
fn parse_vcf(path: &str) -> PyResult<Vec<PyVcfRecord>> {
    let variants = cyanea_io::parse_vcf(path).into_pyresult()?;
    Ok(variants
        .into_iter()
        .map(|v| {
            let filter_str = format!("{:?}", v.filter);
            PyVcfRecord {
                chrom: v.chrom,
                position: v.position,
                id: v.id.unwrap_or_default(),
                ref_allele: String::from_utf8_lossy(&v.ref_allele).into_owned(),
                alt_alleles: v
                    .alt_alleles
                    .iter()
                    .map(|a| String::from_utf8_lossy(a).into_owned())
                    .collect(),
                quality: v.quality,
                filter: filter_str,
            }
        })
        .collect())
}

/// A BED record.
#[pyclass(frozen, get_all)]
pub struct PyBedRecord {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub name: Option<String>,
    pub score: Option<u32>,
}

/// Parse a BED file and return all records.
#[pyfunction]
fn parse_bed(path: &str) -> PyResult<Vec<PyBedRecord>> {
    let records = cyanea_io::parse_bed(path).into_pyresult()?;
    Ok(records
        .into_iter()
        .map(|r| PyBedRecord {
            chrom: r.interval.chrom,
            start: r.interval.start,
            end: r.interval.end,
            name: r.name,
            score: r.score,
        })
        .collect())
}

/// A GFF3 gene record.
#[pyclass(frozen, get_all)]
pub struct PyGff3Gene {
    pub id: String,
    pub name: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: String,
    pub gene_type: String,
    pub transcript_count: usize,
}

/// Parse a GFF3 file and return gene records.
#[pyfunction]
fn parse_gff3(path: &str) -> PyResult<Vec<PyGff3Gene>> {
    let genes = cyanea_io::parse_gff3(path).into_pyresult()?;
    Ok(genes
        .into_iter()
        .map(|g| PyGff3Gene {
            id: g.gene_id,
            name: g.gene_name,
            chrom: g.chrom,
            start: g.start,
            end: g.end,
            strand: format!("{:?}", g.strand),
            gene_type: format!("{:?}", g.gene_type),
            transcript_count: g.transcripts.len(),
        })
        .collect())
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
    m.add_class::<PyVcfRecord>()?;
    m.add_class::<PyBedRecord>()?;
    m.add_class::<PyGff3Gene>()?;
    m.add_function(wrap_pyfunction!(csv_info, &m)?)?;
    m.add_function(wrap_pyfunction!(vcf_stats, &m)?)?;
    m.add_function(wrap_pyfunction!(bed_stats, &m)?)?;
    m.add_function(wrap_pyfunction!(gff3_stats, &m)?)?;
    m.add_function(wrap_pyfunction!(parse_sam, &m)?)?;
    m.add_function(wrap_pyfunction!(sam_stats, &m)?)?;
    m.add_function(wrap_pyfunction!(parse_bam, &m)?)?;
    m.add_function(wrap_pyfunction!(bam_stats, &m)?)?;
    m.add_function(wrap_pyfunction!(parse_vcf, &m)?)?;
    m.add_function(wrap_pyfunction!(parse_bed, &m)?)?;
    m.add_function(wrap_pyfunction!(parse_gff3, &m)?)?;
    parent.add_submodule(&m)?;
    Ok(())
}
