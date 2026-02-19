//! In-memory SAM pileup bindings for WASM environments.
//!
//! Parses SAM text, generates pileup, and computes depth statistics.

use std::collections::HashMap;

use serde::Serialize;

use crate::error::{wasm_err, wasm_ok};

#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::*;

// ── Wrapper types ────────────────────────────────────────────────────────

/// A single pileup column for JSON serialization.
#[derive(Debug, Serialize)]
pub struct JsPileupColumn {
    pub pos: u64,
    pub ref_base: String,
    pub depth: u32,
    pub base_counts: HashMap<String, u32>,
}

/// A pileup for a single reference sequence.
#[derive(Debug, Serialize)]
pub struct JsPileup {
    pub rname: String,
    pub columns: Vec<JsPileupColumn>,
}

/// Depth statistics for a single reference sequence.
#[derive(Debug, Serialize)]
pub struct JsDepthStats {
    pub rname: String,
    pub length: u64,
    pub covered: u64,
    pub breadth: f64,
    pub min_depth: u32,
    pub max_depth: u32,
    pub mean_depth: f64,
    pub median_depth: f64,
}

// ── Helper ───────────────────────────────────────────────────────────────

/// Convert base_counts array [A, C, G, T, N, del] into a string-keyed map.
fn base_counts_map(counts: &[u32; 6]) -> HashMap<String, u32> {
    let mut map = HashMap::new();
    let labels = ["A", "C", "G", "T", "N", "del"];
    for (i, &count) in counts.iter().enumerate() {
        if count > 0 {
            map.insert(labels[i].to_string(), count);
        }
    }
    map
}

// ── JSON boundary functions ──────────────────────────────────────────────

/// Generate pileup from SAM text.
///
/// Parses SAM-formatted text and generates per-position pileup data.
/// Returns JSON array of pileups (one per reference sequence).
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn pileup_from_sam(sam_text: &str) -> String {
    let records = match cyanea_io::parse_sam_str(sam_text) {
        Ok(r) => r,
        Err(e) => return wasm_err(e),
    };
    match cyanea_io::pileup(&records, None) {
        Ok(pileups) => {
            let js: Vec<JsPileup> = pileups
                .iter()
                .map(|p| JsPileup {
                    rname: p.rname.clone(),
                    columns: p
                        .columns
                        .iter()
                        .map(|c| JsPileupColumn {
                            pos: c.pos,
                            ref_base: String::from(c.ref_base as char),
                            depth: c.depth,
                            base_counts: base_counts_map(&c.base_counts),
                        })
                        .collect(),
                })
                .collect();
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Compute depth statistics from SAM text.
///
/// Parses SAM-formatted text and returns per-reference depth statistics.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn depth_stats_from_sam(sam_text: &str) -> String {
    let records = match cyanea_io::parse_sam_str(sam_text) {
        Ok(r) => r,
        Err(e) => return wasm_err(e),
    };
    match cyanea_io::pileup(&records, None) {
        Ok(pileups) => {
            let js: Vec<JsDepthStats> = pileups
                .iter()
                .map(|p| {
                    let ds = cyanea_io::depth_stats(p);
                    JsDepthStats {
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
                .collect();
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Convert SAM text to mpileup format.
///
/// Parses SAM-formatted text and produces mpileup text output.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn pileup_to_mpileup_text(sam_text: &str) -> String {
    let records = match cyanea_io::parse_sam_str(sam_text) {
        Ok(r) => r,
        Err(e) => return wasm_err(e),
    };
    match cyanea_io::pileup(&records, None) {
        Ok(pileups) => {
            let mut result = String::new();
            for p in &pileups {
                let mpileup = cyanea_io::pileup_to_mpileup(p);
                if !result.is_empty() {
                    result.push('\n');
                }
                result.push_str(&mpileup);
            }
            wasm_ok(&result)
        }
        Err(e) => wasm_err(e),
    }
}

// ── Additional wrapper types ────────────────────────────────────────────

#[derive(Debug, Serialize)]
pub struct JsVcfVariant {
    pub chrom: String,
    pub position: u64,
    pub id: Option<String>,
    pub ref_allele: String,
    pub alt_alleles: Vec<String>,
    pub quality: Option<f64>,
    pub filter: String,
}

#[derive(Debug, Serialize)]
pub struct JsBedRecord {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub name: Option<String>,
    pub score: Option<u32>,
}

#[derive(Debug, Serialize)]
pub struct JsGff3Gene {
    pub gene_id: String,
    pub gene_name: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: String,
    pub n_transcripts: usize,
}

#[derive(Debug, Serialize)]
pub struct JsBlastXmlHit {
    pub id: String,
    pub def: String,
    pub accession: String,
    pub length: usize,
    pub n_hsps: usize,
}

#[derive(Debug, Serialize)]
pub struct JsBlastXmlResult {
    pub program: String,
    pub n_iterations: usize,
    pub hits: Vec<JsBlastXmlHit>,
}

#[derive(Debug, Serialize)]
pub struct JsBedGraphRecord {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub value: f64,
}

#[derive(Debug, Serialize)]
pub struct JsGfaGraph {
    pub n_segments: usize,
    pub n_links: usize,
    pub n_paths: usize,
    pub segments: Vec<JsGfaSegment>,
}

#[derive(Debug, Serialize)]
pub struct JsGfaSegment {
    pub name: String,
    pub sequence: String,
    pub length: usize,
}

// ── Additional JSON boundary functions ──────────────────────────────────

/// Parse VCF text and return variants as JSON.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn parse_vcf_text(text: &str) -> String {
    use cyanea_omics::variant::VariantFilter;
    match cyanea_io::parse_vcf_str(text) {
        Ok(variants) => {
            let js: Vec<JsVcfVariant> = variants
                .iter()
                .map(|v| {
                    let filter = match &v.filter {
                        VariantFilter::Pass => "PASS".to_string(),
                        VariantFilter::Missing => ".".to_string(),
                        VariantFilter::Fail(reasons) => reasons.join(";"),
                    };
                    JsVcfVariant {
                        chrom: v.chrom.clone(),
                        position: v.position,
                        id: v.id.clone(),
                        ref_allele: String::from_utf8_lossy(&v.ref_allele).to_string(),
                        alt_alleles: v
                            .alt_alleles
                            .iter()
                            .map(|a| String::from_utf8_lossy(a).to_string())
                            .collect(),
                        quality: v.quality,
                        filter,
                    }
                })
                .collect();
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Parse BED text and return records as JSON.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn parse_bed_text(text: &str) -> String {
    match cyanea_io::parse_bed_str(text) {
        Ok(records) => {
            let js: Vec<JsBedRecord> = records
                .iter()
                .map(|r| JsBedRecord {
                    chrom: r.interval.chrom.clone(),
                    start: r.interval.start,
                    end: r.interval.end,
                    name: r.name.clone(),
                    score: r.score,
                })
                .collect();
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Parse GFF3 text and return gene models as JSON.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn parse_gff3_text(text: &str) -> String {
    use cyanea_omics::genomic::Strand;
    match cyanea_io::parse_gff3_str(text) {
        Ok(genes) => {
            let js: Vec<JsGff3Gene> = genes
                .iter()
                .map(|g| {
                    let strand = match g.strand {
                        Strand::Forward => "+".to_string(),
                        Strand::Reverse => "-".to_string(),
                        Strand::Unknown => ".".to_string(),
                    };
                    JsGff3Gene {
                        gene_id: g.gene_id.clone(),
                        gene_name: g.gene_name.clone(),
                        chrom: g.chrom.clone(),
                        start: g.start,
                        end: g.end,
                        strand,
                        n_transcripts: g.transcripts.len(),
                    }
                })
                .collect();
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Parse BLAST XML output and return results as JSON.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn parse_blast_xml(xml: &str) -> String {
    match cyanea_io::parse_blast_xml_str(xml) {
        Ok(result) => {
            let hits: Vec<JsBlastXmlHit> = result
                .iterations
                .iter()
                .flat_map(|it| &it.hits)
                .map(|h| JsBlastXmlHit {
                    id: h.hit_id.clone(),
                    def: h.hit_def.clone(),
                    accession: h.hit_accession.clone(),
                    length: h.hit_len as usize,
                    n_hsps: h.hsps.len(),
                })
                .collect();
            let js = JsBlastXmlResult {
                program: result.program.clone(),
                n_iterations: result.iterations.len(),
                hits,
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Parse bedGraph text and return records as JSON.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn parse_bedgraph(text: &str) -> String {
    match cyanea_io::parse_bedgraph_str(text) {
        Ok(records) => {
            let js: Vec<JsBedGraphRecord> = records
                .iter()
                .map(|r| JsBedGraphRecord {
                    chrom: r.chrom.clone(),
                    start: r.start,
                    end: r.end,
                    value: r.value,
                })
                .collect();
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Parse GFA text and return graph summary as JSON.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn parse_gfa(text: &str) -> String {
    match cyanea_io::parse_gfa_str(text) {
        Ok(graph) => {
            let segments: Vec<JsGfaSegment> = graph
                .segments
                .iter()
                .map(|s| JsGfaSegment {
                    name: s.name.clone(),
                    sequence: String::from_utf8_lossy(&s.sequence).to_string(),
                    length: s.length.unwrap_or(s.sequence.len() as u64) as usize,
                })
                .collect();
            let js = JsGfaGraph {
                n_segments: graph.segments.len(),
                n_links: graph.links.len(),
                n_paths: graph.paths.len(),
                segments,
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Build an NCBI Entrez efetch URL for the given database and IDs.
///
/// `ids` is a comma-separated string of identifiers.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn ncbi_fetch_url(db: &str, ids: &str, rettype: &str) -> String {
    let id_vec: Vec<&str> = ids.split(',').map(|s| s.trim()).collect();
    let url = cyanea_io::EntrezUrl::efetch(db, &id_vec, rettype, "text");
    wasm_ok(&url)
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAM_TEXT: &str = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*
r2\t0\tchr1\t3\t60\t8M\t*\t0\t0\tGTACGTAC\t*
";

    #[test]
    fn pileup_from_sam_basic() {
        let json = pileup_from_sam(SAM_TEXT);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let pileups = v["ok"].as_array().unwrap();
        assert_eq!(pileups.len(), 1);
        assert_eq!(pileups[0]["rname"], "chr1");
        assert!(pileups[0]["columns"].as_array().unwrap().len() > 0);
    }

    #[test]
    fn depth_stats_from_sam_basic() {
        let json = depth_stats_from_sam(SAM_TEXT);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let stats = v["ok"].as_array().unwrap();
        assert_eq!(stats.len(), 1);
        assert_eq!(stats[0]["rname"], "chr1");
        assert!(stats[0]["max_depth"].as_u64().unwrap() >= 1);
    }

    #[test]
    fn mpileup_text_from_sam() {
        let json = pileup_to_mpileup_text(SAM_TEXT);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let text = v["ok"].as_str().unwrap();
        assert!(text.contains("chr1"));
    }

    #[test]
    fn pileup_from_invalid_sam() {
        let json = pileup_from_sam("not\tvalid\tsam");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert!(v["error"].is_string());
    }

    #[test]
    fn parse_vcf_text_basic() {
        let vcf = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\trs1\tA\tT\t30\tPASS\t.
chr1\t200\t.\tG\tC\t50\t.\t.
";
        let json = parse_vcf_text(vcf);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let variants = v["ok"].as_array().unwrap();
        assert_eq!(variants.len(), 2);
        assert_eq!(variants[0]["chrom"], "chr1");
        assert_eq!(variants[0]["position"], 100);
        assert_eq!(variants[0]["ref_allele"], "A");
        assert_eq!(variants[0]["alt_alleles"][0], "T");
        assert_eq!(variants[0]["filter"], "PASS");
        assert_eq!(variants[1]["position"], 200);
    }

    #[test]
    fn parse_vcf_text_multiallelic() {
        let vcf = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr2\t500\t.\tAC\tA,AT\t99\tPASS\t.
";
        let json = parse_vcf_text(vcf);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let variants = v["ok"].as_array().unwrap();
        assert_eq!(variants.len(), 1);
        let alts = variants[0]["alt_alleles"].as_array().unwrap();
        assert_eq!(alts.len(), 2);
        assert_eq!(alts[0], "A");
        assert_eq!(alts[1], "AT");
    }

    #[test]
    fn parse_bed_text_basic() {
        let bed = "chr1\t0\t1000\tgene1\t100\nchr2\t500\t2000\tgene2\t200\n";
        let json = parse_bed_text(bed);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let records = v["ok"].as_array().unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0]["chrom"], "chr1");
        assert_eq!(records[0]["start"], 0);
        assert_eq!(records[0]["end"], 1000);
        assert_eq!(records[0]["name"], "gene1");
        assert_eq!(records[1]["chrom"], "chr2");
    }

    #[test]
    fn parse_gff3_text_basic() {
        let gff = "\
##gff-version 3
chr1\t.\tgene\t1\t5000\t.\t+\t.\tID=gene1;Name=TP53
chr1\t.\tmRNA\t1\t5000\t.\t+\t.\tID=tx1;Parent=gene1
chr1\t.\texon\t1\t200\t.\t+\t.\tID=exon1;Parent=tx1
";
        let json = parse_gff3_text(gff);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let genes = v["ok"].as_array().unwrap();
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0]["gene_name"], "TP53");
        assert_eq!(genes[0]["strand"], "+");
        assert_eq!(genes[0]["n_transcripts"], 1);
    }

    #[test]
    fn parse_blast_xml_basic() {
        let xml = r#"<?xml version="1.0"?>
<BlastOutput><BlastOutput_program>blastn</BlastOutput_program><BlastOutput_query-ID>Query_1</BlastOutput_query-ID><BlastOutput_query-len>500</BlastOutput_query-len><BlastOutput_iterations><Iteration><Iteration_iter-num>1</Iteration_iter-num><Iteration_hits><Hit><Hit_id>hit1</Hit_id><Hit_def>test hit</Hit_def><Hit_accession>ACC1</Hit_accession><Hit_len>1000</Hit_len><Hit_hsps><Hsp><Hsp_bit-score>100</Hsp_bit-score><Hsp_score>100</Hsp_score><Hsp_evalue>1e-10</Hsp_evalue><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>100</Hsp_query-to><Hsp_hit-from>1</Hsp_hit-from><Hsp_hit-to>100</Hsp_hit-to><Hsp_identity>95</Hsp_identity><Hsp_gaps>0</Hsp_gaps><Hsp_align-len>100</Hsp_align-len><Hsp_qseq>ACGT</Hsp_qseq><Hsp_hseq>ACGT</Hsp_hseq></Hsp></Hit_hsps></Hit></Iteration_hits></Iteration></BlastOutput_iterations></BlastOutput>"#;
        let json = parse_blast_xml(xml);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let result = &v["ok"];
        assert_eq!(result["program"], "blastn");
        assert_eq!(result["n_iterations"], 1);
        let hits = result["hits"].as_array().unwrap();
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0]["id"], "hit1");
        assert_eq!(hits[0]["def"], "test hit");
        assert_eq!(hits[0]["accession"], "ACC1");
        assert_eq!(hits[0]["length"], 1000);
        assert_eq!(hits[0]["n_hsps"], 1);
    }

    #[test]
    fn parse_bedgraph_basic() {
        let bg = "chr1\t0\t100\t1.5\nchr1\t100\t200\t2.0\n";
        let json = parse_bedgraph(bg);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let records = v["ok"].as_array().unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0]["chrom"], "chr1");
        assert_eq!(records[0]["start"], 0);
        assert_eq!(records[0]["end"], 100);
        assert!((records[0]["value"].as_f64().unwrap() - 1.5).abs() < 1e-9);
        assert!((records[1]["value"].as_f64().unwrap() - 2.0).abs() < 1e-9);
    }

    #[test]
    fn parse_gfa_basic() {
        let gfa = "H\tVN:Z:1.0\nS\tseg1\tACGT\nS\tseg2\tTGCA\nL\tseg1\t+\tseg2\t+\t2M\n";
        let json = parse_gfa(gfa);
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let graph = &v["ok"];
        assert_eq!(graph["n_segments"], 2);
        assert_eq!(graph["n_links"], 1);
        assert_eq!(graph["n_paths"], 0);
        let segs = graph["segments"].as_array().unwrap();
        assert_eq!(segs[0]["name"], "seg1");
        assert_eq!(segs[0]["sequence"], "ACGT");
        assert_eq!(segs[0]["length"], 4);
    }

    #[test]
    fn ncbi_fetch_url_basic() {
        let json = ncbi_fetch_url("nucleotide", "12345,67890", "fasta");
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let url = v["ok"].as_str().unwrap();
        assert!(url.contains("db=nucleotide"));
        assert!(url.contains("12345"));
        assert!(url.contains("67890"));
        assert!(url.contains("rettype=fasta"));
    }
}
