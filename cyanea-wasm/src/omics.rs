//! Omics wrappers with JSON input/output.
//!
//! Provides genomic interval operations, liftover, variant annotation, CNV
//! segmentation, methylation analysis, and spatial statistics — all accepting
//! JSON strings and returning JSON, suitable for the WASM boundary.

use serde::Serialize;

use cyanea_omics::annotation::{Exon, Gene, GeneType, Transcript};
use cyanea_omics::cnv::{CbsConfig, circular_binary_segmentation};
use cyanea_omics::genome_arithmetic;
use cyanea_omics::genomic::{GenomicInterval, Strand};
use cyanea_omics::liftover;
use cyanea_omics::methylation;
use cyanea_omics::spatial::SpatialGraph;
use cyanea_omics::variant::Variant;
use cyanea_omics::variant_annotation::{self, AnnotationConfig, Consequence};

use crate::error::{wasm_err, wasm_ok};

#[cfg(feature = "wasm")]
use wasm_bindgen::prelude::*;

// ── Wrapper types ────────────────────────────────────────────────────────

#[derive(Debug, Serialize)]
pub struct JsGenomicInterval {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

#[derive(Debug, Serialize)]
pub struct JsVariantEffect {
    pub gene_id: String,
    pub gene_name: String,
    pub transcript_id: String,
    pub consequence: String,
    pub hgvs_c: Option<String>,
    pub hgvs_p: Option<String>,
}

#[derive(Debug, Serialize)]
pub struct JsCnvSegment {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub log2_ratio: f64,
    pub n_probes: usize,
    pub copy_number: u32,
}

#[derive(Debug, Serialize)]
pub struct JsCpgIsland {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub cpg_count: usize,
    pub obs_exp_ratio: f64,
    pub gc_content: f64,
}

#[derive(Debug, Serialize)]
pub struct JsSpatialAutocorrelation {
    pub statistic: f64,
    pub expected: f64,
    pub variance: f64,
    pub z_score: f64,
    pub p_value: f64,
}

#[derive(Debug, Serialize)]
pub struct JsLiftoverResult {
    pub status: String,
    pub chrom: Option<String>,
    pub start: Option<u64>,
    pub end: Option<u64>,
    pub fraction_mapped: Option<f64>,
}

#[derive(Debug, Serialize)]
pub struct JsJaccard {
    pub jaccard: f64,
    pub intersection: u64,
    pub union_size: u64,
}

#[derive(Debug, Serialize)]
pub struct JsClosestResult {
    pub query: JsGenomicInterval,
    pub closest: Option<JsGenomicInterval>,
    pub distance: Option<u64>,
}

#[derive(Debug, Serialize)]
pub struct JsGearysC {
    pub c: f64,
    pub expected_c: f64,
    pub z_score: f64,
    pub p_value: f64,
}

// ── Helpers ──────────────────────────────────────────────────────────────

fn parse_intervals(json: &str) -> Result<Vec<GenomicInterval>, String> {
    let arr: Vec<serde_json::Value> =
        serde_json::from_str(json).map_err(|e| format!("invalid JSON array: {e}"))?;
    let mut intervals = Vec::with_capacity(arr.len());
    for v in &arr {
        let chrom = v["chrom"]
            .as_str()
            .ok_or_else(|| "interval missing 'chrom' string".to_string())?;
        let start = v["start"]
            .as_u64()
            .ok_or_else(|| "interval missing 'start' integer".to_string())?;
        let end = v["end"]
            .as_u64()
            .ok_or_else(|| "interval missing 'end' integer".to_string())?;
        let iv = GenomicInterval::new(chrom, start, end)
            .map_err(|e| format!("invalid interval: {e}"))?;
        intervals.push(iv);
    }
    Ok(intervals)
}

fn parse_genome(json: &str) -> Result<genome_arithmetic::GenomeInfo, String> {
    let arr: Vec<serde_json::Value> =
        serde_json::from_str(json).map_err(|e| format!("invalid genome JSON: {e}"))?;
    let mut genome = genome_arithmetic::GenomeInfo::new();
    for v in &arr {
        let chrom = v["chrom"]
            .as_str()
            .ok_or_else(|| "genome entry missing 'chrom'".to_string())?;
        let length = v["length"]
            .as_u64()
            .ok_or_else(|| "genome entry missing 'length'".to_string())?;
        genome.insert(chrom.to_string(), length);
    }
    Ok(genome)
}

fn iv_to_js(iv: &GenomicInterval) -> JsGenomicInterval {
    JsGenomicInterval {
        chrom: iv.chrom.clone(),
        start: iv.start,
        end: iv.end,
    }
}

fn consequence_to_string(c: &Consequence) -> String {
    match c {
        Consequence::Missense { .. } => "missense_variant".to_string(),
        Consequence::Nonsense { .. } => "stop_gained".to_string(),
        Consequence::Synonymous { .. } => "synonymous_variant".to_string(),
        Consequence::Frameshift { .. } => "frameshift_variant".to_string(),
        Consequence::InFrame { .. } => "inframe_variant".to_string(),
        Consequence::FivePrimeUtr => "5_prime_UTR_variant".to_string(),
        Consequence::ThreePrimeUtr => "3_prime_UTR_variant".to_string(),
        Consequence::SpliceSite { donor } => {
            if *donor {
                "splice_donor_variant".to_string()
            } else {
                "splice_acceptor_variant".to_string()
            }
        }
        Consequence::Intronic => "intron_variant".to_string(),
        Consequence::Upstream => "upstream_gene_variant".to_string(),
        Consequence::Downstream => "downstream_gene_variant".to_string(),
        Consequence::NonCoding => "non_coding_transcript_variant".to_string(),
        Consequence::StartLoss { .. } => "start_lost".to_string(),
        Consequence::StopLoss { .. } => "stop_lost".to_string(),
    }
}

fn parse_f64_array(json: &str) -> Result<Vec<f64>, String> {
    serde_json::from_str::<Vec<f64>>(json).map_err(|e| format!("invalid JSON array: {e}"))
}

fn parse_u64_array(json: &str) -> Result<Vec<u64>, String> {
    serde_json::from_str::<Vec<u64>>(json).map_err(|e| format!("invalid JSON array: {e}"))
}

// ── Interval operations ──────────────────────────────────────────────────

/// Merge overlapping intervals.
///
/// Input: JSON array of `{chrom, start, end}` objects.
/// Output: JSON array of merged `JsGenomicInterval`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn merge_intervals(json: &str) -> String {
    let intervals = match parse_intervals(json) {
        Ok(v) => v,
        Err(e) => return wasm_err(e),
    };
    let merged = genome_arithmetic::merge(&intervals, genome_arithmetic::StrandMode::Ignore);
    let js: Vec<JsGenomicInterval> = merged.iter().map(iv_to_js).collect();
    wasm_ok(&js)
}

/// Intersect two interval sets, returning overlapping sub-regions.
///
/// Input: two JSON arrays of `{chrom, start, end}`.
/// Output: JSON array of `JsGenomicInterval`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn intersect_intervals(a_json: &str, b_json: &str) -> String {
    let a = match parse_intervals(a_json) {
        Ok(v) => v,
        Err(e) => return wasm_err(e),
    };
    let b = match parse_intervals(b_json) {
        Ok(v) => v,
        Err(e) => return wasm_err(e),
    };
    match genome_arithmetic::intersect(&a, &b, genome_arithmetic::StrandMode::Ignore) {
        Ok(result) => {
            let js: Vec<JsGenomicInterval> = result.iter().map(iv_to_js).collect();
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Subtract intervals in `b` from intervals in `a`.
///
/// Input: two JSON arrays of `{chrom, start, end}`.
/// Output: JSON array of `JsGenomicInterval`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn subtract_intervals(a_json: &str, b_json: &str) -> String {
    let a = match parse_intervals(a_json) {
        Ok(v) => v,
        Err(e) => return wasm_err(e),
    };
    let b = match parse_intervals(b_json) {
        Ok(v) => v,
        Err(e) => return wasm_err(e),
    };
    match genome_arithmetic::subtract(&a, &b, genome_arithmetic::StrandMode::Ignore) {
        Ok(result) => {
            let js: Vec<JsGenomicInterval> = result.iter().map(iv_to_js).collect();
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Complement of intervals with respect to a genome.
///
/// `json`: JSON array of `{chrom, start, end}`.
/// `genome_json`: JSON array of `{chrom, length}`.
/// Output: JSON array of `JsGenomicInterval`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn complement_intervals(json: &str, genome_json: &str) -> String {
    let intervals = match parse_intervals(json) {
        Ok(v) => v,
        Err(e) => return wasm_err(e),
    };
    let genome = match parse_genome(genome_json) {
        Ok(g) => g,
        Err(e) => return wasm_err(e),
    };
    match genome_arithmetic::complement(&intervals, &genome) {
        Ok(result) => {
            let js: Vec<JsGenomicInterval> = result.iter().map(iv_to_js).collect();
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Find the closest interval in `b` for each interval in `a`.
///
/// Input: two JSON arrays of `{chrom, start, end}`.
/// Output: JSON array of `JsClosestResult`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn closest_intervals(a_json: &str, b_json: &str) -> String {
    let a = match parse_intervals(a_json) {
        Ok(v) => v,
        Err(e) => return wasm_err(e),
    };
    let b = match parse_intervals(b_json) {
        Ok(v) => v,
        Err(e) => return wasm_err(e),
    };
    let results = genome_arithmetic::closest(&a, &b, genome_arithmetic::StrandMode::Ignore);
    let js: Vec<JsClosestResult> = results
        .iter()
        .map(|r| JsClosestResult {
            query: iv_to_js(&r.query),
            closest: r.closest.as_ref().map(iv_to_js),
            distance: r.distance,
        })
        .collect();
    wasm_ok(&js)
}

/// Jaccard similarity between two interval sets.
///
/// Input: two JSON arrays of `{chrom, start, end}`.
/// Output: JSON `JsJaccard`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn jaccard_intervals(a_json: &str, b_json: &str) -> String {
    let a = match parse_intervals(a_json) {
        Ok(v) => v,
        Err(e) => return wasm_err(e),
    };
    let b = match parse_intervals(b_json) {
        Ok(v) => v,
        Err(e) => return wasm_err(e),
    };
    let stats = genome_arithmetic::jaccard_stats(&a, &b);
    let js = JsJaccard {
        jaccard: stats.jaccard,
        intersection: stats.intersection_bp,
        union_size: stats.union_bp,
    };
    wasm_ok(&js)
}

/// Generate non-overlapping tiling windows across a genome.
///
/// `genome_json`: JSON array of `{chrom, length}`.
/// `window_size`: window size in bases.
/// Output: JSON array of `JsGenomicInterval`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn make_windows(genome_json: &str, window_size: u64) -> String {
    let genome = match parse_genome(genome_json) {
        Ok(g) => g,
        Err(e) => return wasm_err(e),
    };
    match genome_arithmetic::make_windows(&genome, window_size) {
        Ok(result) => {
            let js: Vec<JsGenomicInterval> = result.iter().map(iv_to_js).collect();
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

// ── Liftover ─────────────────────────────────────────────────────────────

/// Liftover a single genomic interval using a UCSC chain file.
///
/// `chain_text`: full chain file content as a string.
/// `chrom`, `start`, `end`: interval to liftover.
/// Output: JSON `JsLiftoverResult`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn liftover_interval(chain_text: &str, chrom: &str, start: u64, end: u64) -> String {
    let chain_file = match liftover::parse_chain(chain_text) {
        Ok(cf) => cf,
        Err(e) => return wasm_err(e),
    };
    let interval = match GenomicInterval::new(chrom, start, end) {
        Ok(iv) => iv,
        Err(e) => return wasm_err(e),
    };
    let result = liftover::liftover(&chain_file, &interval, 0.5);
    let js = match result {
        liftover::LiftoverResult::Mapped(mapped) => JsLiftoverResult {
            status: "mapped".to_string(),
            chrom: Some(mapped.chrom),
            start: Some(mapped.start),
            end: Some(mapped.end),
            fraction_mapped: Some(1.0),
        },
        liftover::LiftoverResult::Partial {
            mapped,
            fraction_mapped,
        } => JsLiftoverResult {
            status: "partial".to_string(),
            chrom: Some(mapped.chrom),
            start: Some(mapped.start),
            end: Some(mapped.end),
            fraction_mapped: Some(fraction_mapped),
        },
        liftover::LiftoverResult::Unmapped => JsLiftoverResult {
            status: "unmapped".to_string(),
            chrom: None,
            start: None,
            end: None,
            fraction_mapped: None,
        },
    };
    wasm_ok(&js)
}

// ── Variant annotation ───────────────────────────────────────────────────

fn parse_gene_type(s: &str) -> GeneType {
    match s {
        "protein_coding" => GeneType::ProteinCoding,
        "lncRNA" => GeneType::LncRNA,
        "miRNA" => GeneType::MiRNA,
        "rRNA" => GeneType::RRNA,
        "tRNA" => GeneType::TRNA,
        "pseudogene" => GeneType::Pseudogene,
        other => GeneType::Other(other.to_string()),
    }
}

fn parse_strand(s: &str) -> Strand {
    match s {
        "+" => Strand::Forward,
        "-" => Strand::Reverse,
        _ => Strand::Unknown,
    }
}

fn parse_genes(json: &str) -> Result<Vec<Gene>, String> {
    let arr: Vec<serde_json::Value> =
        serde_json::from_str(json).map_err(|e| format!("invalid genes JSON: {e}"))?;
    let mut genes = Vec::with_capacity(arr.len());
    for v in &arr {
        let gene_id = v["gene_id"]
            .as_str()
            .ok_or_else(|| "gene missing 'gene_id'".to_string())?
            .to_string();
        let gene_name = v["gene_name"]
            .as_str()
            .ok_or_else(|| "gene missing 'gene_name'".to_string())?
            .to_string();
        let chrom = v["chrom"]
            .as_str()
            .ok_or_else(|| "gene missing 'chrom'".to_string())?
            .to_string();
        let start = v["start"]
            .as_u64()
            .ok_or_else(|| "gene missing 'start'".to_string())?;
        let end = v["end"]
            .as_u64()
            .ok_or_else(|| "gene missing 'end'".to_string())?;
        let strand = parse_strand(v["strand"].as_str().unwrap_or("."));
        let gene_type = parse_gene_type(v["gene_type"].as_str().unwrap_or("protein_coding"));

        let transcripts_arr = v["transcripts"]
            .as_array()
            .ok_or_else(|| "gene missing 'transcripts' array".to_string())?;
        let mut transcripts = Vec::with_capacity(transcripts_arr.len());
        for tv in transcripts_arr {
            let transcript_id = tv["transcript_id"]
                .as_str()
                .ok_or_else(|| "transcript missing 'transcript_id'".to_string())?
                .to_string();
            let tx_start = tv["start"]
                .as_u64()
                .ok_or_else(|| "transcript missing 'start'".to_string())?;
            let tx_end = tv["end"]
                .as_u64()
                .ok_or_else(|| "transcript missing 'end'".to_string())?;
            let cds_start = tv["cds_start"].as_u64();
            let cds_end = tv["cds_end"].as_u64();

            let exons_arr = tv["exons"]
                .as_array()
                .ok_or_else(|| "transcript missing 'exons' array".to_string())?;
            let mut exons = Vec::with_capacity(exons_arr.len());
            for (idx, ev) in exons_arr.iter().enumerate() {
                let exon_start = ev["start"]
                    .as_u64()
                    .ok_or_else(|| format!("exon {} missing 'start'", idx))?;
                let exon_end = ev["end"]
                    .as_u64()
                    .ok_or_else(|| format!("exon {} missing 'end'", idx))?;
                exons.push(Exon {
                    exon_number: (idx + 1) as u32,
                    start: exon_start,
                    end: exon_end,
                });
            }

            transcripts.push(Transcript {
                transcript_id,
                start: tx_start,
                end: tx_end,
                exons,
                cds_start,
                cds_end,
            });
        }

        genes.push(Gene {
            gene_id,
            gene_name,
            chrom,
            start,
            end,
            strand,
            gene_type,
            transcripts,
        });
    }
    Ok(genes)
}

/// Annotate a variant against a set of genes.
///
/// `variant_json`: `{chrom, position, ref_allele, alt_alleles}`.
///   Position is 1-based (VCF convention).
///   `ref_allele` and `alt_alleles` are strings (e.g., `"A"` and `["T"]`).
/// `genes_json`: array of gene objects (see parse_genes for format).
/// Output: JSON array of `JsVariantEffect`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn annotate_variant(variant_json: &str, genes_json: &str) -> String {
    let vv: serde_json::Value = match serde_json::from_str(variant_json) {
        Ok(v) => v,
        Err(e) => return wasm_err(format!("invalid variant JSON: {e}")),
    };
    let chrom = match vv["chrom"].as_str() {
        Some(c) => c,
        None => return wasm_err("variant missing 'chrom'"),
    };
    let position = match vv["position"].as_u64() {
        Some(p) => p,
        None => return wasm_err("variant missing 'position'"),
    };
    let ref_allele_str = match vv["ref_allele"].as_str() {
        Some(r) => r,
        None => return wasm_err("variant missing 'ref_allele'"),
    };
    let alt_alleles_arr = match vv["alt_alleles"].as_array() {
        Some(a) => a,
        None => return wasm_err("variant missing 'alt_alleles' array"),
    };
    let alt_alleles: Vec<Vec<u8>> = alt_alleles_arr
        .iter()
        .filter_map(|v| v.as_str().map(|s| s.as_bytes().to_vec()))
        .collect();
    if alt_alleles.is_empty() {
        return wasm_err("alt_alleles must contain at least one string");
    }

    let variant = match Variant::new(chrom, position, ref_allele_str.as_bytes().to_vec(), alt_alleles) {
        Ok(v) => v,
        Err(e) => return wasm_err(e),
    };

    let genes = match parse_genes(genes_json) {
        Ok(g) => g,
        Err(e) => return wasm_err(e),
    };

    let config = AnnotationConfig::default();
    let effects = variant_annotation::annotate_variant(&variant, &genes, &config);
    let js: Vec<JsVariantEffect> = effects
        .iter()
        .map(|eff| JsVariantEffect {
            gene_id: eff.gene_id.clone(),
            gene_name: eff.gene_name.clone(),
            transcript_id: eff.transcript_id.clone(),
            consequence: consequence_to_string(&eff.consequence),
            hgvs_c: eff.hgvs_c.clone(),
            hgvs_p: eff.hgvs_p.clone(),
        })
        .collect();
    wasm_ok(&js)
}

// ── CNV ──────────────────────────────────────────────────────────────────

/// Segment a log2 ratio profile using Circular Binary Segmentation.
///
/// `positions_json`: JSON array of u64 positions.
/// `values_json`: JSON array of f64 log2 ratio values.
/// `chrom`: chromosome name.
/// `alpha`: significance threshold (e.g., 0.01).
/// `min_probes`: minimum probes per segment (e.g., 3).
/// Output: JSON array of `JsCnvSegment`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn cbs_segment(
    positions_json: &str,
    values_json: &str,
    chrom: &str,
    alpha: f64,
    min_probes: usize,
) -> String {
    let positions = match parse_u64_array(positions_json) {
        Ok(p) => p,
        Err(e) => return wasm_err(e),
    };
    let values = match parse_f64_array(values_json) {
        Ok(v) => v,
        Err(e) => return wasm_err(e),
    };
    let config = CbsConfig {
        alpha,
        min_probes,
        ..Default::default()
    };
    match circular_binary_segmentation(&positions, &values, chrom, &config) {
        Ok(segments) => {
            let js: Vec<JsCnvSegment> = segments
                .iter()
                .map(|s| JsCnvSegment {
                    chrom: s.chrom.clone(),
                    start: s.start,
                    end: s.end,
                    log2_ratio: s.log2_ratio,
                    n_probes: s.n_probes,
                    copy_number: s.copy_number,
                })
                .collect();
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

// ── Methylation ──────────────────────────────────────────────────────────

/// Simulate bisulfite conversion of a DNA sequence.
///
/// `seq`: DNA sequence string.
/// `methylated_json`: JSON array of 0-based positions that are methylated.
/// Output: JSON string of the converted sequence.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn bisulfite_convert(seq: &str, methylated_json: &str) -> String {
    let methylated = match parse_u64_array(methylated_json) {
        Ok(m) => m,
        Err(e) => return wasm_err(e),
    };
    let converted = methylation::bisulfite_convert(seq.as_bytes(), &methylated);
    let converted_str = String::from_utf8_lossy(&converted).to_string();
    wasm_ok(&converted_str)
}

/// Find CpG islands in a DNA sequence.
///
/// `seq`: DNA sequence string.
/// `chrom`: chromosome name.
/// Output: JSON array of `JsCpgIsland`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn find_cpg_islands(seq: &str, chrom: &str) -> String {
    let islands = methylation::find_cpg_islands(seq.as_bytes(), chrom, 200, 0.5, 0.6);
    let js: Vec<JsCpgIsland> = islands
        .iter()
        .map(|island| JsCpgIsland {
            chrom: island.chrom.clone(),
            start: island.start,
            end: island.end,
            cpg_count: island.cpg_count,
            obs_exp_ratio: island.obs_exp_ratio,
            gc_content: island.gc_content,
        })
        .collect();
    wasm_ok(&js)
}

// ── Spatial ──────────────────────────────────────────────────────────────

fn parse_spatial_graph(neighbors_json: &str, n_nodes: usize) -> Result<SpatialGraph, String> {
    let arr: Vec<Vec<Vec<serde_json::Value>>> = serde_json::from_str(neighbors_json)
        .map_err(|e| format!("invalid neighbors JSON: {e}"))?;
    if arr.len() != n_nodes {
        return Err(format!(
            "neighbors array length ({}) does not match values length ({})",
            arr.len(),
            n_nodes
        ));
    }
    let mut neighbors = Vec::with_capacity(n_nodes);
    for node_neighbors in &arr {
        let mut nb = Vec::with_capacity(node_neighbors.len());
        for pair in node_neighbors {
            if pair.len() < 2 {
                return Err("neighbor entry must have [index, distance]".to_string());
            }
            let idx = pair[0]
                .as_u64()
                .ok_or_else(|| "neighbor index must be integer".to_string())?
                as usize;
            let dist = pair[1]
                .as_f64()
                .ok_or_else(|| "neighbor distance must be number".to_string())?;
            nb.push((idx, dist));
        }
        neighbors.push(nb);
    }
    Ok(SpatialGraph {
        n_nodes,
        neighbors,
    })
}

/// Compute Moran's I spatial autocorrelation.
///
/// `values_json`: JSON array of f64 values (one per node).
/// `neighbors_json`: JSON array of arrays, where each inner array contains
///   `[neighbor_index, distance]` pairs.
/// Output: JSON `JsSpatialAutocorrelation`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn morans_i(values_json: &str, neighbors_json: &str) -> String {
    let values = match parse_f64_array(values_json) {
        Ok(v) => v,
        Err(e) => return wasm_err(e),
    };
    let graph = match parse_spatial_graph(neighbors_json, values.len()) {
        Ok(g) => g,
        Err(e) => return wasm_err(e),
    };
    match cyanea_omics::spatial::morans_i(&values, &graph) {
        Ok(result) => {
            let js = JsSpatialAutocorrelation {
                statistic: result.morans_i,
                expected: result.expected_i,
                variance: result.variance_i,
                z_score: result.z_score,
                p_value: result.p_value,
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

/// Compute Geary's C spatial autocorrelation.
///
/// `values_json`: JSON array of f64 values (one per node).
/// `neighbors_json`: JSON array of arrays of `[neighbor_index, distance]`.
/// Output: JSON `JsGearysC`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn gearys_c(values_json: &str, neighbors_json: &str) -> String {
    let values = match parse_f64_array(values_json) {
        Ok(v) => v,
        Err(e) => return wasm_err(e),
    };
    let graph = match parse_spatial_graph(neighbors_json, values.len()) {
        Ok(g) => g,
        Err(e) => return wasm_err(e),
    };
    match cyanea_omics::spatial::gearys_c(&values, &graph) {
        Ok(result) => {
            let js = JsGearysC {
                c: result.c,
                expected_c: result.expected_c,
                z_score: result.z_score,
                p_value: result.p_value,
            };
            wasm_ok(&js)
        }
        Err(e) => wasm_err(e),
    }
}

// ── Tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // Helper to parse the JSON result envelope and return the "ok" value.
    fn ok_val(json: &str) -> serde_json::Value {
        let v: serde_json::Value = serde_json::from_str(json).unwrap();
        assert!(
            v["ok"].is_array() || v["ok"].is_object() || v["ok"].is_string() || v["ok"].is_number(),
            "expected ok, got: {json}"
        );
        v["ok"].clone()
    }

    #[allow(dead_code)]
    fn err_val(json: &str) -> String {
        let v: serde_json::Value = serde_json::from_str(json).unwrap();
        v["error"].as_str().unwrap().to_string()
    }

    // 1. merge_intervals
    #[test]
    fn test_merge_intervals() {
        let input = r#"[
            {"chrom":"chr1","start":100,"end":200},
            {"chrom":"chr1","start":150,"end":300},
            {"chrom":"chr1","start":500,"end":600}
        ]"#;
        let result = merge_intervals(input);
        let ok = ok_val(&result);
        let arr = ok.as_array().unwrap();
        assert_eq!(arr.len(), 2);
        assert_eq!(arr[0]["start"], 100);
        assert_eq!(arr[0]["end"], 300);
        assert_eq!(arr[1]["start"], 500);
        assert_eq!(arr[1]["end"], 600);
    }

    // 2. intersect_intervals
    #[test]
    fn test_intersect_intervals() {
        let a = r#"[{"chrom":"chr1","start":100,"end":300}]"#;
        let b = r#"[{"chrom":"chr1","start":200,"end":400}]"#;
        let result = intersect_intervals(a, b);
        let ok = ok_val(&result);
        let arr = ok.as_array().unwrap();
        assert_eq!(arr.len(), 1);
        assert_eq!(arr[0]["start"], 200);
        assert_eq!(arr[0]["end"], 300);
    }

    // 3. subtract_intervals
    #[test]
    fn test_subtract_intervals() {
        let a = r#"[{"chrom":"chr1","start":100,"end":400}]"#;
        let b = r#"[{"chrom":"chr1","start":200,"end":300}]"#;
        let result = subtract_intervals(a, b);
        let ok = ok_val(&result);
        let arr = ok.as_array().unwrap();
        assert_eq!(arr.len(), 2);
        assert_eq!(arr[0]["start"], 100);
        assert_eq!(arr[0]["end"], 200);
        assert_eq!(arr[1]["start"], 300);
        assert_eq!(arr[1]["end"], 400);
    }

    // 4. complement_intervals
    #[test]
    fn test_complement_intervals() {
        let intervals = r#"[{"chrom":"chr1","start":100,"end":200}]"#;
        let genome = r#"[{"chrom":"chr1","length":500}]"#;
        let result = complement_intervals(intervals, genome);
        let ok = ok_val(&result);
        let arr = ok.as_array().unwrap();
        assert_eq!(arr.len(), 2);
        assert_eq!(arr[0]["start"], 0);
        assert_eq!(arr[0]["end"], 100);
        assert_eq!(arr[1]["start"], 200);
        assert_eq!(arr[1]["end"], 500);
    }

    // 5. closest_intervals
    #[test]
    fn test_closest_intervals() {
        let a = r#"[{"chrom":"chr1","start":300,"end":400}]"#;
        let b = r#"[{"chrom":"chr1","start":100,"end":200}]"#;
        let result = closest_intervals(a, b);
        let ok = ok_val(&result);
        let arr = ok.as_array().unwrap();
        assert_eq!(arr.len(), 1);
        assert_eq!(arr[0]["distance"], 100);
        assert_eq!(arr[0]["closest"]["start"], 100);
        assert_eq!(arr[0]["closest"]["end"], 200);
    }

    // 6. jaccard_intervals
    #[test]
    fn test_jaccard_intervals() {
        let a = r#"[{"chrom":"chr1","start":100,"end":300}]"#;
        let b = r#"[{"chrom":"chr1","start":200,"end":400}]"#;
        let result = jaccard_intervals(a, b);
        let ok = ok_val(&result);
        // intersection = 100bp, union = 300bp, jaccard = 100/300
        let j = ok["jaccard"].as_f64().unwrap();
        assert!((j - 100.0 / 300.0).abs() < 1e-10);
        assert_eq!(ok["intersection"], 100);
        assert_eq!(ok["union_size"], 300);
    }

    // 7. make_windows
    #[test]
    fn test_make_windows() {
        let genome = r#"[{"chrom":"chr1","length":100}]"#;
        let result = make_windows(genome, 30);
        let ok = ok_val(&result);
        let arr = ok.as_array().unwrap();
        // [0,30), [30,60), [60,90), [90,100)
        assert_eq!(arr.len(), 4);
        assert_eq!(arr[0]["start"], 0);
        assert_eq!(arr[0]["end"], 30);
        assert_eq!(arr[3]["start"], 90);
        assert_eq!(arr[3]["end"], 100);
    }

    // 8. liftover_interval - mapped
    #[test]
    fn test_liftover_mapped() {
        let chain = "chain 1000 chr1 1000 + 100 600 chrA 800 + 50 500 1\n200 100 50\n200\n";
        let result = liftover_interval(chain, "chr1", 150, 250);
        let ok = ok_val(&result);
        assert_eq!(ok["status"], "mapped");
        assert_eq!(ok["chrom"], "chrA");
        assert_eq!(ok["start"], 100);
        assert_eq!(ok["end"], 200);
        assert_eq!(ok["fraction_mapped"], 1.0);
    }

    // 9. liftover_interval - unmapped
    #[test]
    fn test_liftover_unmapped() {
        let chain = "chain 1000 chr1 1000 + 100 600 chrA 800 + 50 500 1\n200 100 50\n200\n";
        let result = liftover_interval(chain, "chrX", 100, 200);
        let ok = ok_val(&result);
        assert_eq!(ok["status"], "unmapped");
        assert!(ok["chrom"].is_null());
        assert!(ok["start"].is_null());
    }

    // 10. annotate_variant - SNV in CDS
    #[test]
    fn test_annotate_variant_snv() {
        let variant = r#"{
            "chrom": "chr1",
            "position": 105,
            "ref_allele": "T",
            "alt_alleles": ["A"]
        }"#;
        let genes = r#"[{
            "gene_id": "GENE001",
            "gene_name": "TEST1",
            "chrom": "chr1",
            "start": 100,
            "end": 400,
            "strand": "+",
            "gene_type": "protein_coding",
            "transcripts": [{
                "transcript_id": "TX001",
                "start": 100,
                "end": 400,
                "exons": [{"start": 100, "end": 400}],
                "cds_start": 100,
                "cds_end": 400
            }]
        }]"#;
        let result = annotate_variant(variant, genes);
        let ok = ok_val(&result);
        let arr = ok.as_array().unwrap();
        assert_eq!(arr.len(), 1);
        assert_eq!(arr[0]["gene_name"], "TEST1");
        assert_eq!(arr[0]["transcript_id"], "TX001");
        // Variant is in the CDS, should be a coding consequence
        let cons = arr[0]["consequence"].as_str().unwrap();
        assert!(
            cons.contains("variant") || cons.contains("lost") || cons.contains("gained"),
            "expected coding consequence, got: {cons}"
        );
    }

    // 11. annotate_variant - intronic variant
    #[test]
    fn test_annotate_variant_intronic() {
        let variant = r#"{
            "chrom": "chr1",
            "position": 1351,
            "ref_allele": "A",
            "alt_alleles": ["G"]
        }"#;
        let genes = r#"[{
            "gene_id": "GENE001",
            "gene_name": "TEST1",
            "chrom": "chr1",
            "start": 1000,
            "end": 2000,
            "strand": "+",
            "gene_type": "protein_coding",
            "transcripts": [{
                "transcript_id": "TX001",
                "start": 1000,
                "end": 2000,
                "exons": [
                    {"start": 1000, "end": 1200},
                    {"start": 1500, "end": 1800}
                ],
                "cds_start": 1050,
                "cds_end": 1750
            }]
        }]"#;
        let result = annotate_variant(variant, genes);
        let ok = ok_val(&result);
        let arr = ok.as_array().unwrap();
        assert_eq!(arr.len(), 1);
        assert_eq!(arr[0]["consequence"], "intron_variant");
    }

    // 12. cbs_segment - basic segmentation
    #[test]
    fn test_cbs_segment_basic() {
        let n = 40;
        let positions: Vec<u64> = (0..n).map(|i| i as u64 * 1000).collect();
        let mut values = vec![0.0; n as usize];
        for i in 20..n as usize {
            values[i] = 1.0;
        }
        let pos_json = serde_json::to_string(&positions).unwrap();
        let val_json = serde_json::to_string(&values).unwrap();
        let result = cbs_segment(&pos_json, &val_json, "chr1", 0.05, 3);
        let ok = ok_val(&result);
        let arr = ok.as_array().unwrap();
        assert!(arr.len() >= 2, "expected at least 2 segments, got {}", arr.len());
        // First segment should have low log2_ratio, last should have high
        let first_lr = arr[0]["log2_ratio"].as_f64().unwrap();
        let last_lr = arr.last().unwrap()["log2_ratio"].as_f64().unwrap();
        assert!(first_lr < 0.5);
        assert!(last_lr > 0.5);
    }

    // 13. cbs_segment - empty input error
    #[test]
    fn test_cbs_segment_empty() {
        let result = cbs_segment("[]", "[]", "chr1", 0.01, 3);
        let v: serde_json::Value = serde_json::from_str(&result).unwrap();
        assert!(v["error"].is_string());
    }

    // 14. bisulfite_convert
    #[test]
    fn test_bisulfite_convert() {
        let result = bisulfite_convert("ACGT", "[]");
        let ok = ok_val(&result);
        // C at position 1 should convert to T (no methylated positions)
        assert_eq!(ok.as_str().unwrap(), "ATGT");
    }

    // 15. find_cpg_islands
    #[test]
    fn test_find_cpg_islands() {
        // Build a 250 bp CG-rich sequence
        let mut seq = String::new();
        for _ in 0..125 {
            seq.push_str("CG");
        }
        let result = find_cpg_islands(&seq, "chr1");
        let ok = ok_val(&result);
        let arr = ok.as_array().unwrap();
        assert!(!arr.is_empty(), "should find at least one CpG island");
        let island = &arr[0];
        assert!(island["gc_content"].as_f64().unwrap() >= 0.5);
        assert!(island["obs_exp_ratio"].as_f64().unwrap() >= 0.6);
        assert!(island["cpg_count"].as_u64().unwrap() > 0);
    }

    // 16. morans_i - positive autocorrelation (clustered)
    #[test]
    fn test_morans_i_positive() {
        // 4x4 grid with clustered values: top-left high, bottom-right low
        let values = vec![
            10.0, 10.0, 9.0, 9.0,
            10.0, 10.0, 9.0, 9.0,
            1.0, 1.0, 2.0, 2.0,
            1.0, 1.0, 2.0, 2.0,
        ];
        // Build a 4x4 grid neighbor graph
        let neighbors = build_grid_neighbors(4, 4);
        let values_json = serde_json::to_string(&values).unwrap();
        let neighbors_json = serde_json::to_string(&neighbors).unwrap();
        let result = morans_i(&values_json, &neighbors_json);
        let ok = ok_val(&result);
        let stat = ok["statistic"].as_f64().unwrap();
        assert!(stat > 0.0, "clustered data should have positive Moran's I, got {stat}");
    }

    // 17. morans_i - random (near zero)
    #[test]
    fn test_morans_i_random() {
        let values = vec![
            12.0, 5.0, 18.0, 3.0, 15.0,
            9.0, 22.0, 1.0, 14.0, 7.0,
            20.0, 4.0, 16.0, 8.0, 11.0,
            2.0, 17.0, 6.0, 23.0, 10.0,
            13.0, 19.0, 24.0, 21.0, 25.0,
        ];
        let neighbors = build_grid_neighbors(5, 5);
        let values_json = serde_json::to_string(&values).unwrap();
        let neighbors_json = serde_json::to_string(&neighbors).unwrap();
        let result = morans_i(&values_json, &neighbors_json);
        let ok = ok_val(&result);
        let stat = ok["statistic"].as_f64().unwrap();
        assert!(
            stat.abs() < 0.6,
            "random data should have modest Moran's I, got {stat}"
        );
    }

    // 18. gearys_c - basic test
    #[test]
    fn test_gearys_c() {
        // Clustered data should give Geary's C < 1
        let values = vec![
            10.0, 10.0, 9.0, 9.0,
            10.0, 10.0, 9.0, 9.0,
            1.0, 1.0, 2.0, 2.0,
            1.0, 1.0, 2.0, 2.0,
        ];
        let neighbors = build_grid_neighbors(4, 4);
        let values_json = serde_json::to_string(&values).unwrap();
        let neighbors_json = serde_json::to_string(&neighbors).unwrap();
        let result = gearys_c(&values_json, &neighbors_json);
        let ok = ok_val(&result);
        let c = ok["c"].as_f64().unwrap();
        assert!(c < 1.0, "clustered data should have Geary's C < 1, got {c}");
        assert_eq!(ok["expected_c"], 1.0);
    }

    /// Build a grid neighbor list for testing.
    /// Returns Vec<Vec<[usize, f64]>> suitable for JSON serialization.
    fn build_grid_neighbors(rows: usize, cols: usize) -> Vec<Vec<[serde_json::Number; 2]>> {
        let n = rows * cols;
        let mut neighbors: Vec<Vec<[serde_json::Number; 2]>> = vec![Vec::new(); n];
        for r in 0..rows {
            for c in 0..cols {
                let i = r * cols + c;
                if c + 1 < cols {
                    let j = r * cols + c + 1;
                    neighbors[i].push([
                        serde_json::Number::from(j),
                        serde_json::Number::from_f64(1.0).unwrap(),
                    ]);
                    neighbors[j].push([
                        serde_json::Number::from(i),
                        serde_json::Number::from_f64(1.0).unwrap(),
                    ]);
                }
                if r + 1 < rows {
                    let j = (r + 1) * cols + c;
                    neighbors[i].push([
                        serde_json::Number::from(j),
                        serde_json::Number::from_f64(1.0).unwrap(),
                    ]);
                    neighbors[j].push([
                        serde_json::Number::from(i),
                        serde_json::Number::from_f64(1.0).unwrap(),
                    ]);
                }
            }
        }
        neighbors
    }
}
