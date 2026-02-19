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
    pub rnext: String,
    pub pnext: u64,
    pub tlen: i64,
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
            rnext: r.rnext,
            pnext: r.pnext,
            tlen: r.tlen,
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

/// Summary statistics for paired-end SAM/BAM data.
#[pyclass(frozen, get_all)]
pub struct PyPairedSamStats {
    pub total_reads: usize,
    pub mapped: usize,
    pub unmapped: usize,
    pub avg_mapq: f64,
    pub avg_length: f64,
    pub paired_count: usize,
    pub proper_pair_count: usize,
    pub singletons: usize,
    pub avg_insert_size: f64,
}

/// Compute paired-end statistics for a SAM file.
#[pyfunction]
fn paired_sam_stats_fn(path: &str) -> PyResult<PyPairedSamStats> {
    let records = cyanea_io::parse_sam(path).into_pyresult()?;
    let stats = cyanea_io::paired_sam_stats(&records);
    Ok(PyPairedSamStats {
        total_reads: stats.base.total_reads,
        mapped: stats.base.mapped,
        unmapped: stats.base.unmapped,
        avg_mapq: stats.base.avg_mapq,
        avg_length: stats.base.avg_length,
        paired_count: stats.paired_count,
        proper_pair_count: stats.proper_pair_count,
        singletons: stats.singletons,
        avg_insert_size: stats.avg_insert_size,
    })
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
// Pileup
// ---------------------------------------------------------------------------

/// A single pileup column.
#[pyclass(frozen, get_all)]
pub struct PyPileupColumn {
    pub rname: String,
    pub pos: u64,
    pub ref_base: String,
    pub depth: u32,
    pub base_counts: std::collections::HashMap<String, u32>,
    pub quality_sums: std::collections::HashMap<String, u64>,
}

/// Depth statistics for a reference sequence.
#[pyclass(frozen, get_all)]
pub struct PyDepthStats {
    pub rname: String,
    pub length: u64,
    pub covered: u64,
    pub breadth: f64,
    pub min_depth: u32,
    pub max_depth: u32,
    pub mean_depth: f64,
    pub median_depth: f64,
}

fn base_counts_map(counts: &[u32; 6]) -> std::collections::HashMap<String, u32> {
    let mut map = std::collections::HashMap::new();
    let labels = ["A", "C", "G", "T", "N", "del"];
    for (i, &count) in counts.iter().enumerate() {
        if count > 0 {
            map.insert(labels[i].to_string(), count);
        }
    }
    map
}

fn quality_sums_map(sums: &[u64; 6]) -> std::collections::HashMap<String, u64> {
    let mut map = std::collections::HashMap::new();
    let labels = ["A", "C", "G", "T", "N", "del"];
    for (i, &sum) in sums.iter().enumerate() {
        if sum > 0 {
            map.insert(labels[i].to_string(), sum);
        }
    }
    map
}

/// Generate pileup from a SAM file.
#[pyfunction]
fn pileup_from_sam(path: &str) -> PyResult<Vec<Vec<PyPileupColumn>>> {
    let records = cyanea_io::parse_sam(path).into_pyresult()?;
    let pileups = cyanea_io::pileup(&records, None).into_pyresult()?;
    Ok(pileups
        .iter()
        .map(|p| {
            p.columns
                .iter()
                .map(|c| PyPileupColumn {
                    rname: p.rname.clone(),
                    pos: c.pos,
                    ref_base: String::from(c.ref_base as char),
                    depth: c.depth,
                    base_counts: base_counts_map(&c.base_counts),
                    quality_sums: quality_sums_map(&c.quality_sums),
                })
                .collect()
        })
        .collect())
}

/// Compute depth statistics from a SAM file.
#[pyfunction]
fn depth_stats_from_sam(path: &str) -> PyResult<Vec<PyDepthStats>> {
    let records = cyanea_io::parse_sam(path).into_pyresult()?;
    let pileups = cyanea_io::pileup(&records, None).into_pyresult()?;
    Ok(pileups
        .iter()
        .map(|p| {
            let ds = cyanea_io::depth_stats(p);
            PyDepthStats {
                rname: ds.rname,
                length: ds.length,
                covered: ds.covered,
                breadth: ds.breadth,
                min_depth: ds.min_depth,
                max_depth: ds.max_depth,
                mean_depth: ds.mean_depth,
                median_depth: ds.median_depth,
            }
        })
        .collect())
}

/// Convert a SAM file to mpileup text format.
#[pyfunction]
fn pileup_to_mpileup_str(path: &str) -> PyResult<String> {
    let records = cyanea_io::parse_sam(path).into_pyresult()?;
    let pileups = cyanea_io::pileup(&records, None).into_pyresult()?;
    let mut result = String::new();
    for p in &pileups {
        let mpileup = cyanea_io::pileup_to_mpileup(p);
        if !result.is_empty() {
            result.push('\n');
        }
        result.push_str(&mpileup);
    }
    Ok(result)
}

// ---------------------------------------------------------------------------
// BLAST XML
// ---------------------------------------------------------------------------

/// A BLAST XML search result.
#[pyclass(frozen, get_all)]
pub struct PyBlastXmlResult {
    pub program: String,
    pub query_id: String,
    pub query_len: u64,
    pub hits: Vec<PyBlastXmlHit>,
}

/// A single hit from BLAST XML output.
#[derive(Clone)]
#[pyclass(frozen, get_all)]
pub struct PyBlastXmlHit {
    pub hit_id: String,
    pub hit_def: String,
    pub hit_accession: String,
    pub hit_len: u64,
    pub bit_score: f64,
    pub evalue: f64,
}

/// Parse BLAST XML output text.
#[pyfunction]
fn parse_blast_xml(text: &str) -> PyResult<PyBlastXmlResult> {
    let result = cyanea_io::parse_blast_xml_str(text).into_pyresult()?;
    let hits = result.iterations.into_iter().flat_map(|it| it.hits).map(|h| {
        let best_hsp = h.hsps.first();
        PyBlastXmlHit {
            hit_id: h.hit_id,
            hit_def: h.hit_def,
            hit_accession: h.hit_accession,
            hit_len: h.hit_len,
            bit_score: best_hsp.map_or(0.0, |hsp| hsp.bit_score),
            evalue: best_hsp.map_or(1.0, |hsp| hsp.evalue),
        }
    }).collect();
    Ok(PyBlastXmlResult {
        program: result.program,
        query_id: result.query_id,
        query_len: result.query_len,
        hits,
    })
}

// ---------------------------------------------------------------------------
// bedGraph
// ---------------------------------------------------------------------------

/// A bedGraph record.
#[pyclass(frozen, get_all)]
pub struct PyBedGraphRecord {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub value: f64,
}

/// Parse bedGraph format text.
#[pyfunction]
fn parse_bedgraph(text: &str) -> PyResult<Vec<PyBedGraphRecord>> {
    let records = cyanea_io::parse_bedgraph_str(text).into_pyresult()?;
    Ok(records.into_iter().map(|r| PyBedGraphRecord {
        chrom: r.chrom,
        start: r.start,
        end: r.end,
        value: r.value,
    }).collect())
}

// ---------------------------------------------------------------------------
// GFA
// ---------------------------------------------------------------------------

/// Summary of a parsed GFA graph.
#[pyclass(frozen, get_all)]
pub struct PyGfaGraph {
    pub n_segments: usize,
    pub n_links: usize,
    pub n_paths: usize,
    pub total_sequence_length: usize,
}

/// Parse GFA (Graphical Fragment Assembly) format text.
#[pyfunction]
fn parse_gfa(text: &str) -> PyResult<PyGfaGraph> {
    let graph = cyanea_io::parse_gfa_str(text).into_pyresult()?;
    let total_sequence_length = graph.segments.iter().map(|s| s.sequence.len()).sum();
    Ok(PyGfaGraph {
        n_segments: graph.segments.len(),
        n_links: graph.links.len(),
        n_paths: graph.paths.len(),
        total_sequence_length,
    })
}

// ---------------------------------------------------------------------------
// Fetch URL builders
// ---------------------------------------------------------------------------

/// Build an NCBI Entrez efetch URL.
#[pyfunction]
#[pyo3(signature = (db, ids, rettype, retmode="text"))]
fn ncbi_efetch_url(db: &str, ids: Vec<String>, rettype: &str, retmode: &str) -> String {
    let id_refs: Vec<&str> = ids.iter().map(|s| s.as_str()).collect();
    cyanea_io::EntrezUrl::efetch(db, &id_refs, rettype, retmode)
}

/// Build a UniProt entry retrieval URL.
#[pyfunction]
#[pyo3(signature = (accession, format="fasta"))]
fn uniprot_url(accession: &str, format: &str) -> String {
    cyanea_io::UniprotUrl::entry(accession, format)
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
    m.add_class::<PyPairedSamStats>()?;
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
    m.add_function(wrap_pyfunction!(paired_sam_stats_fn, &m)?)?;
    m.add_function(wrap_pyfunction!(parse_vcf, &m)?)?;
    m.add_function(wrap_pyfunction!(parse_bed, &m)?)?;
    m.add_function(wrap_pyfunction!(parse_gff3, &m)?)?;
    // Pileup
    m.add_class::<PyPileupColumn>()?;
    m.add_class::<PyDepthStats>()?;
    m.add_function(wrap_pyfunction!(pileup_from_sam, &m)?)?;
    m.add_function(wrap_pyfunction!(depth_stats_from_sam, &m)?)?;
    m.add_function(wrap_pyfunction!(pileup_to_mpileup_str, &m)?)?;
    // BLAST XML
    m.add_class::<PyBlastXmlResult>()?;
    m.add_class::<PyBlastXmlHit>()?;
    m.add_function(wrap_pyfunction!(parse_blast_xml, &m)?)?;
    // bedGraph
    m.add_class::<PyBedGraphRecord>()?;
    m.add_function(wrap_pyfunction!(parse_bedgraph, &m)?)?;
    // GFA
    m.add_class::<PyGfaGraph>()?;
    m.add_function(wrap_pyfunction!(parse_gfa, &m)?)?;
    // Fetch URL builders
    m.add_function(wrap_pyfunction!(ncbi_efetch_url, &m)?)?;
    m.add_function(wrap_pyfunction!(uniprot_url, &m)?)?;
    parent.add_submodule(&m)?;
    Ok(())
}
