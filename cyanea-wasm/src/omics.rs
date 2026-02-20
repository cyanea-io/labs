//! Omics wrappers with JSON input/output.
//!
//! Provides genomic interval operations, liftover, variant annotation, CNV
//! segmentation, methylation analysis, and spatial statistics — all accepting
//! JSON strings and returning JSON, suitable for the WASM boundary.

use serde::{Deserialize, Serialize};

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

// ── Single-cell ──────────────────────────────────────────────────────────

#[derive(Debug, Serialize)]
pub struct JsHvgResult {
    pub gene_indices: Vec<usize>,
    pub dispersions: Vec<f64>,
}

#[derive(Debug, Serialize)]
pub struct JsNeighborsResult {
    pub n_neighbors: usize,
    pub distances: Vec<(usize, usize, f64)>,
    pub connectivities: Vec<(usize, usize, f64)>,
}

#[derive(Debug, Serialize)]
pub struct JsClusterResult {
    pub labels: Vec<usize>,
    pub modularity: f64,
}

#[derive(Debug, Serialize)]
pub struct JsDiffusionResult {
    pub components: Vec<Vec<f64>>,
    pub eigenvalues: Vec<f64>,
}

#[derive(Debug, Serialize)]
pub struct JsDptResult {
    pub pseudotime: Vec<f64>,
}

#[derive(Debug, Serialize)]
pub struct JsPagaResult {
    pub connectivities: Vec<Vec<f64>>,
    pub groups: Vec<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct JsMarkerGene {
    pub gene_index: usize,
    pub gene_name: String,
    pub score: f64,
    pub pvalue: f64,
    pub padj: f64,
    pub log2fc: f64,
    pub pct_in: f64,
    pub pct_out: f64,
}

#[derive(Debug, Serialize)]
pub struct JsMarkerResult {
    pub markers: std::collections::HashMap<String, Vec<JsMarkerGene>>,
}

#[derive(Debug, Serialize)]
pub struct JsHarmonyResult {
    pub corrected: Vec<f64>,
    pub n_obs: usize,
    pub n_vars: usize,
}

/// Build an AnnData from a flat row-major data array.
fn build_adata(data: &[f64], n_features: usize) -> Result<cyanea_omics::single_cell::AnnData, String> {
    use cyanea_omics::single_cell::{AnnData, MatrixData};
    if data.is_empty() || n_features == 0 {
        return Err("data and n_features must be non-empty".to_string());
    }
    if data.len() % n_features != 0 {
        return Err(format!(
            "data length ({}) is not divisible by n_features ({})",
            data.len(),
            n_features
        ));
    }
    let n_obs = data.len() / n_features;
    let rows: Vec<Vec<f64>> = data.chunks(n_features).map(|c| c.to_vec()).collect();
    let obs_names: Vec<String> = (0..n_obs).map(|i| format!("cell_{i}")).collect();
    let var_names: Vec<String> = (0..n_features).map(|i| format!("gene_{i}")).collect();
    AnnData::new(MatrixData::Dense(rows), obs_names, var_names).map_err(|e| e.to_string())
}

/// Normalize cells to a target sum and optionally log-transform.
///
/// `data_json`: flat row-major JSON array of f64.
/// `n_features`: number of genes (columns).
/// `target_sum`: normalization target (e.g. 10000).
/// `log1p`: whether to apply log(1+x) afterwards.
/// Output: JSON array of corrected values (flat row-major).
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn sc_normalize(data_json: &str, n_features: usize, target_sum: f64, log1p: bool) -> String {
    use cyanea_omics::sc_preprocess::{NormalizeConfig, normalize_total};
    let data = match parse_f64_array(data_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let mut adata = match build_adata(&data, n_features) {
        Ok(a) => a,
        Err(e) => return wasm_err(e),
    };
    let config = NormalizeConfig {
        target_sum,
        log_transform: log1p,
        save_raw: false,
    };
    if let Err(e) = normalize_total(&mut adata, &config) {
        return wasm_err(e);
    }
    let result = adata.x().to_flat_row_major();
    wasm_ok(&result)
}

/// Identify highly variable genes.
///
/// `data_json`: flat row-major JSON array of f64.
/// `n_features`: number of genes (columns).
/// `n_top_genes`: how many HVGs to select.
/// `method`: `"seurat_v3"` or `"cell_ranger"`.
/// Output: JSON `JsHvgResult`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn sc_hvg(data_json: &str, n_features: usize, n_top_genes: usize, method: &str) -> String {
    use cyanea_omics::sc_preprocess::{HvgConfig, HvgMethod, highly_variable_genes};
    let data = match parse_f64_array(data_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let mut adata = match build_adata(&data, n_features) {
        Ok(a) => a,
        Err(e) => return wasm_err(e),
    };
    let hvg_method = match method {
        "cell_ranger" => HvgMethod::CellRanger,
        _ => HvgMethod::SeuratV3,
    };
    let config = HvgConfig {
        n_top_genes,
        method: hvg_method,
        ..HvgConfig::default()
    };
    if let Err(e) = highly_variable_genes(&mut adata, &config) {
        return wasm_err(e);
    }
    // Extract HVG indices and dispersions from var metadata
    let hvg_col = adata.get_var("highly_variable");
    let disp_col = adata.get_var("dispersions_norm");
    let mut gene_indices = Vec::new();
    let mut dispersions = Vec::new();
    if let Some(cyanea_omics::ColumnData::Numeric(hvg)) = hvg_col {
        for (i, &v) in hvg.iter().enumerate() {
            if v > 0.0 {
                gene_indices.push(i);
                if let Some(cyanea_omics::ColumnData::Numeric(d)) = disp_col {
                    dispersions.push(d[i]);
                }
            }
        }
    }
    wasm_ok(&JsHvgResult {
        gene_indices,
        dispersions,
    })
}

/// Regress out covariates from the expression matrix.
///
/// `data_json`: flat row-major JSON array of f64.
/// `n_features`: number of genes (columns).
/// `covariates_json`: JSON array of f64 arrays (one per covariate).
/// Output: JSON array of corrected values (flat row-major).
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn sc_regress_out(data_json: &str, n_features: usize, covariates_json: &str) -> String {
    use cyanea_omics::sc_preprocess::regress_out;
    let data = match parse_f64_array(data_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let mut adata = match build_adata(&data, n_features) {
        Ok(a) => a,
        Err(e) => return wasm_err(e),
    };
    let covariates: Vec<Vec<f64>> = match serde_json::from_str(covariates_json) {
        Ok(c) => c,
        Err(e) => return wasm_err(format!("invalid covariates JSON: {e}")),
    };
    // Store covariates in obs
    let mut keys = Vec::new();
    for (i, cov) in covariates.iter().enumerate() {
        let key = format!("covariate_{i}");
        if let Err(e) = adata.add_obs_numeric(&key, cov.clone()) {
            return wasm_err(e);
        }
        keys.push(key);
    }
    let key_refs: Vec<&str> = keys.iter().map(|s| s.as_str()).collect();
    if let Err(e) = regress_out(&mut adata, &key_refs) {
        return wasm_err(e);
    }
    let result = adata.x().to_flat_row_major();
    wasm_ok(&result)
}

/// Score a set of genes per cell (like scanpy.tl.score_genes).
///
/// `data_json`: flat row-major JSON array of f64.
/// `n_features`: number of genes (columns).
/// `gene_indices_json`: JSON array of usize gene indices to score.
/// Output: JSON array of f64 scores (one per cell).
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn sc_score_genes(data_json: &str, n_features: usize, gene_indices_json: &str) -> String {
    use cyanea_omics::sc_preprocess::score_genes;
    let data = match parse_f64_array(data_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let mut adata = match build_adata(&data, n_features) {
        Ok(a) => a,
        Err(e) => return wasm_err(e),
    };
    let gene_indices: Vec<usize> = match serde_json::from_str(gene_indices_json) {
        Ok(g) => g,
        Err(e) => return wasm_err(format!("invalid gene_indices JSON: {e}")),
    };
    if let Err(e) = score_genes(&mut adata, &gene_indices, 50, "score") {
        return wasm_err(e);
    }
    match adata.get_obs("score") {
        Some(cyanea_omics::ColumnData::Numeric(scores)) => wasm_ok(scores),
        _ => wasm_err("score_genes did not produce scores"),
    }
}

/// Compute a k-nearest neighbors graph.
///
/// `data_json`: flat row-major JSON array of f64.
/// `n_features`: number of genes (columns).
/// `n_neighbors`: number of neighbors to find.
/// `metric`: `"euclidean"` or `"cosine"`.
/// Output: JSON `JsNeighborsResult`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn sc_neighbors(data_json: &str, n_features: usize, n_neighbors: usize, metric: &str) -> String {
    use cyanea_omics::sc_cluster::{DistanceMetric, NeighborsConfig, neighbors};
    let data = match parse_f64_array(data_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let mut adata = match build_adata(&data, n_features) {
        Ok(a) => a,
        Err(e) => return wasm_err(e),
    };
    let dist_metric = match metric {
        "cosine" => DistanceMetric::Cosine,
        _ => DistanceMetric::Euclidean,
    };
    let config = NeighborsConfig {
        n_neighbors,
        n_pcs: 0, // use raw data
        metric: dist_metric,
        seed: 42,
    };
    if let Err(e) = neighbors(&mut adata, &config) {
        return wasm_err(e);
    }
    let distances = adata
        .get_obsp("distances")
        .map(|s| s.iter().collect::<Vec<_>>())
        .unwrap_or_default();
    let connectivities = adata
        .get_obsp("connectivities")
        .map(|s| s.iter().collect::<Vec<_>>())
        .unwrap_or_default();
    wasm_ok(&JsNeighborsResult {
        n_neighbors,
        distances,
        connectivities,
    })
}

/// Leiden clustering on a precomputed neighbor graph.
///
/// `neighbors_json`: JSON with `{distances: [[r,c,v],...], connectivities: [[r,c,v],...], n_obs: N}`.
/// `resolution`: resolution parameter (higher = more clusters).
/// Output: JSON `JsClusterResult`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn sc_leiden(neighbors_json: &str, resolution: f64) -> String {
    sc_cluster_impl(neighbors_json, resolution, true)
}

/// Louvain clustering on a precomputed neighbor graph.
///
/// `neighbors_json`: JSON with `{distances: [[r,c,v],...], connectivities: [[r,c,v],...], n_obs: N}`.
/// `resolution`: resolution parameter.
/// Output: JSON `JsClusterResult`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn sc_louvain(neighbors_json: &str, resolution: f64) -> String {
    sc_cluster_impl(neighbors_json, resolution, false)
}

fn sc_cluster_impl(neighbors_json: &str, resolution: f64, use_leiden: bool) -> String {
    use cyanea_omics::sc_cluster::{ClusterConfig, leiden, louvain};
    use cyanea_omics::single_cell::{AnnData, MatrixData};
    use cyanea_omics::sparse::SparseMatrix;

    let v: serde_json::Value = match serde_json::from_str(neighbors_json) {
        Ok(v) => v,
        Err(e) => return wasm_err(format!("invalid neighbors JSON: {e}")),
    };
    let n_obs = match v["n_obs"].as_u64() {
        Some(n) => n as usize,
        None => return wasm_err("neighbors JSON missing 'n_obs'"),
    };
    // Build a minimal AnnData with the neighbor graph in obsp
    let dummy_data = vec![vec![0.0; 1]; n_obs];
    let obs_names: Vec<String> = (0..n_obs).map(|i| format!("cell_{i}")).collect();
    let var_names = vec!["dummy".to_string()];
    let mut adata = match AnnData::new(MatrixData::Dense(dummy_data), obs_names, var_names) {
        Ok(a) => a,
        Err(e) => return wasm_err(e),
    };
    // Parse distances and connectivities
    if let Some(arr) = v["distances"].as_array() {
        let mut sm = SparseMatrix::new(n_obs, n_obs);
        for entry in arr {
            if let (Some(r), Some(c), Some(val)) = (
                entry[0].as_u64(),
                entry[1].as_u64(),
                entry[2].as_f64(),
            ) {
                let _ = sm.insert(r as usize, c as usize, val);
            }
        }
        let _ = adata.add_obsp("distances", sm);
    }
    if let Some(arr) = v["connectivities"].as_array() {
        let mut sm = SparseMatrix::new(n_obs, n_obs);
        for entry in arr {
            if let (Some(r), Some(c), Some(val)) = (
                entry[0].as_u64(),
                entry[1].as_u64(),
                entry[2].as_f64(),
            ) {
                let _ = sm.insert(r as usize, c as usize, val);
            }
        }
        let _ = adata.add_obsp("connectivities", sm);
    }

    let config = ClusterConfig {
        resolution,
        n_iterations: 10,
        seed: 42,
        key_added: "cluster".to_string(),
    };
    let result = if use_leiden {
        leiden(&mut adata, &config)
    } else {
        louvain(&mut adata, &config)
    };
    if let Err(e) = result {
        return wasm_err(e);
    }
    // Extract cluster labels
    match adata.get_obs("cluster") {
        Some(cyanea_omics::ColumnData::Categorical { codes, .. }) => {
            let labels: Vec<usize> = codes.iter().map(|&c| c as usize).collect();
            let n_clusters = labels.iter().copied().max().map_or(0, |m| m + 1);
            // Compute a simple modularity proxy
            let modularity = n_clusters as f64 / n_obs as f64;
            wasm_ok(&JsClusterResult {
                labels,
                modularity,
            })
        }
        Some(cyanea_omics::ColumnData::Strings(labels)) => {
            // Convert string labels to usize
            let mut label_map = std::collections::HashMap::new();
            let mut next_id = 0usize;
            let int_labels: Vec<usize> = labels.iter().map(|l| {
                *label_map.entry(l.clone()).or_insert_with(|| {
                    let id = next_id;
                    next_id += 1;
                    id
                })
            }).collect();
            wasm_ok(&JsClusterResult {
                labels: int_labels,
                modularity: 0.0,
            })
        }
        _ => wasm_err("clustering did not produce labels"),
    }
}

/// Compute a diffusion map embedding.
///
/// `data_json`: flat row-major JSON array of f64.
/// `n_features`: number of genes (columns).
/// `n_components`: number of diffusion components.
/// Output: JSON `JsDiffusionResult`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn sc_diffusion_map(data_json: &str, n_features: usize, n_components: usize) -> String {
    use cyanea_omics::sc_trajectory::{DiffusionConfig, diffusion_map};
    let data = match parse_f64_array(data_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let mut adata = match build_adata(&data, n_features) {
        Ok(a) => a,
        Err(e) => return wasm_err(e),
    };
    let config = DiffusionConfig {
        n_components,
        alpha: 1.0,
    };
    match diffusion_map(&mut adata, &config) {
        Ok(result) => wasm_ok(&JsDiffusionResult {
            components: result.components,
            eigenvalues: result.eigenvalues,
        }),
        Err(e) => wasm_err(e),
    }
}

/// Compute diffusion pseudotime.
///
/// `diffmap_json`: JSON `{components: [[...], ...], eigenvalues: [...], n_obs: N}`.
/// `root_cell`: index of the root cell.
/// Output: JSON `JsDptResult`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn sc_dpt(diffmap_json: &str, root_cell: usize) -> String {
    use cyanea_omics::sc_trajectory::{DptConfig, dpt};
    use cyanea_omics::single_cell::{AnnData, MatrixData};

    let v: serde_json::Value = match serde_json::from_str(diffmap_json) {
        Ok(v) => v,
        Err(e) => return wasm_err(format!("invalid diffmap JSON: {e}")),
    };
    let components: Vec<Vec<f64>> = match serde_json::from_value(v["components"].clone()) {
        Ok(c) => c,
        Err(e) => return wasm_err(format!("invalid components: {e}")),
    };
    let eigenvalues: Vec<f64> = match serde_json::from_value(v["eigenvalues"].clone()) {
        Ok(e) => e,
        Err(e) => return wasm_err(format!("invalid eigenvalues: {e}")),
    };
    let n_obs = components.len();
    if n_obs == 0 {
        return wasm_err("empty components");
    }
    // Build AnnData with diffmap in obsm
    let dummy = vec![vec![0.0; 1]; n_obs];
    let obs_names: Vec<String> = (0..n_obs).map(|i| format!("cell_{i}")).collect();
    let var_names = vec!["dummy".to_string()];
    let mut adata = match AnnData::new(MatrixData::Dense(dummy), obs_names, var_names) {
        Ok(a) => a,
        Err(e) => return wasm_err(e),
    };
    let _ = adata.add_obsm("X_diffmap", components);
    adata.add_uns("diffmap_evals", serde_json::to_string(&eigenvalues).unwrap_or_default());

    let config = DptConfig {
        root_cell,
        n_branchings: 0,
    };
    if let Err(e) = dpt(&mut adata, &config) {
        return wasm_err(e);
    }
    match adata.get_obs("dpt_pseudotime") {
        Some(cyanea_omics::ColumnData::Numeric(pt)) => {
            wasm_ok(&JsDptResult {
                pseudotime: pt.clone(),
            })
        }
        _ => wasm_err("DPT did not produce pseudotime"),
    }
}

/// Compute PAGA graph abstraction.
///
/// `neighbors_json`: JSON with `{connectivities: [[r,c,v],...], n_obs: N}`.
/// `clusters_json`: JSON array of cluster label strings (one per cell).
/// Output: JSON `JsPagaResult`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn sc_paga(neighbors_json: &str, clusters_json: &str) -> String {
    use cyanea_omics::sc_trajectory::paga;
    use cyanea_omics::single_cell::{AnnData, MatrixData};
    use cyanea_omics::sparse::SparseMatrix;

    let v: serde_json::Value = match serde_json::from_str(neighbors_json) {
        Ok(v) => v,
        Err(e) => return wasm_err(format!("invalid neighbors JSON: {e}")),
    };
    let clusters: Vec<String> = match serde_json::from_str(clusters_json) {
        Ok(c) => c,
        Err(e) => return wasm_err(format!("invalid clusters JSON: {e}")),
    };
    let n_obs = clusters.len();
    let dummy = vec![vec![0.0; 1]; n_obs];
    let obs_names: Vec<String> = (0..n_obs).map(|i| format!("cell_{i}")).collect();
    let var_names = vec!["dummy".to_string()];
    let mut adata = match AnnData::new(MatrixData::Dense(dummy), obs_names, var_names) {
        Ok(a) => a,
        Err(e) => return wasm_err(e),
    };
    if let Err(e) = adata.add_obs("leiden", clusters) {
        return wasm_err(e);
    }
    if let Some(arr) = v["connectivities"].as_array() {
        let mut sm = SparseMatrix::new(n_obs, n_obs);
        for entry in arr {
            if let (Some(r), Some(c), Some(val)) = (
                entry[0].as_u64(),
                entry[1].as_u64(),
                entry[2].as_f64(),
            ) {
                let _ = sm.insert(r as usize, c as usize, val);
            }
        }
        let _ = adata.add_obsp("connectivities", sm);
    }
    match paga(&adata, "leiden") {
        Ok(result) => wasm_ok(&JsPagaResult {
            connectivities: result.connectivities,
            groups: result.groups,
        }),
        Err(e) => wasm_err(e),
    }
}

/// Rank genes per cluster (differential expression).
///
/// `data_json`: flat row-major JSON array of f64.
/// `n_features`: number of genes (columns).
/// `clusters_json`: JSON array of cluster label strings (one per cell).
/// `method`: `"t-test"`, `"wilcoxon"`, or `"logistic"`.
/// Output: JSON `JsMarkerResult`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn sc_rank_genes(data_json: &str, n_features: usize, clusters_json: &str, method: &str) -> String {
    use cyanea_omics::sc_markers::{MarkerConfig, MarkerMethod, rank_genes_groups};
    let data = match parse_f64_array(data_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let mut adata = match build_adata(&data, n_features) {
        Ok(a) => a,
        Err(e) => return wasm_err(e),
    };
    let clusters: Vec<String> = match serde_json::from_str(clusters_json) {
        Ok(c) => c,
        Err(e) => return wasm_err(format!("invalid clusters JSON: {e}")),
    };
    if let Err(e) = adata.add_obs("cluster", clusters) {
        return wasm_err(e);
    }
    let marker_method = match method {
        "wilcoxon" => MarkerMethod::Wilcoxon,
        "logistic" => MarkerMethod::LogisticRegression,
        _ => MarkerMethod::TTest,
    };
    let config = MarkerConfig {
        method: marker_method,
        cluster_key: "cluster".to_string(),
        log2fc_threshold: 0.0,
        min_pct: 0.0,
        padj_threshold: 1.0,
        n_genes: None,
    };
    match rank_genes_groups(&adata, &config) {
        Ok(results) => {
            let markers: std::collections::HashMap<String, Vec<JsMarkerGene>> = results
                .markers
                .iter()
                .map(|(k, genes)| {
                    let js_genes: Vec<JsMarkerGene> = genes
                        .iter()
                        .map(|g| JsMarkerGene {
                            gene_index: g.gene_index,
                            gene_name: g.gene_name.clone(),
                            score: g.statistic,
                            pvalue: g.p_value,
                            padj: g.p_adjusted,
                            log2fc: g.log2_fold_change,
                            pct_in: g.pct_in,
                            pct_out: g.pct_out,
                        })
                        .collect();
                    (k.clone(), js_genes)
                })
                .collect();
            wasm_ok(&JsMarkerResult { markers })
        }
        Err(e) => wasm_err(e),
    }
}

/// Filter marker genes by fold-change, pct, and p-value thresholds.
///
/// `markers_json`: JSON `JsMarkerResult`.
/// `log2fc`: minimum absolute log2 fold change.
/// `min_pct`: minimum fraction of cells expressing the gene in the cluster.
/// `padj`: maximum adjusted p-value.
/// Output: JSON `JsMarkerResult` (filtered).
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn sc_filter_markers(markers_json: &str, log2fc: f64, min_pct: f64, padj: f64) -> String {
    let parsed: serde_json::Value = match serde_json::from_str(markers_json) {
        Ok(v) => v,
        Err(e) => return wasm_err(format!("invalid markers JSON: {e}")),
    };
    let input: std::collections::HashMap<String, Vec<JsMarkerGene>> =
        match serde_json::from_value(parsed["markers"].clone()) {
            Ok(m) => m,
            Err(_) => match serde_json::from_str(markers_json) {
                Ok(m) => m,
                Err(e) => return wasm_err(format!("could not parse markers: {e}")),
            },
        };
    let filtered: std::collections::HashMap<String, Vec<JsMarkerGene>> = input
        .into_iter()
        .map(|(k, genes)| {
            let filt: Vec<JsMarkerGene> = genes
                .into_iter()
                .filter(|g| g.log2fc.abs() >= log2fc && g.pct_in >= min_pct && g.padj <= padj)
                .collect();
            (k, filt)
        })
        .collect();
    wasm_ok(&JsMarkerResult { markers: filtered })
}

/// Run Harmony batch correction.
///
/// `data_json`: flat row-major JSON array of f64.
/// `n_features`: number of genes (columns).
/// `batch_json`: JSON array of batch label strings (one per cell).
/// `n_clusters`: optional number of clusters for Harmony.
/// Output: JSON `JsHarmonyResult`.
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn sc_harmony(data_json: &str, n_features: usize, batch_json: &str, n_clusters: usize) -> String {
    use cyanea_omics::sc_integrate::{HarmonyConfig, harmony};
    let data = match parse_f64_array(data_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let mut adata = match build_adata(&data, n_features) {
        Ok(a) => a,
        Err(e) => return wasm_err(e),
    };
    let batches: Vec<String> = match serde_json::from_str(batch_json) {
        Ok(b) => b,
        Err(e) => return wasm_err(format!("invalid batch JSON: {e}")),
    };
    if let Err(e) = adata.add_obs("batch", batches) {
        return wasm_err(e);
    }
    let config = HarmonyConfig {
        batch_key: "batch".to_string(),
        n_clusters: if n_clusters > 0 { Some(n_clusters) } else { None },
        theta: 2.0,
        sigma: 0.1,
        max_iter: 10,
    };
    if let Err(e) = harmony(&mut adata, &config) {
        return wasm_err(e);
    }
    let corrected = adata.x().to_flat_row_major();
    let (n_obs, n_vars) = adata.shape();
    wasm_ok(&JsHarmonyResult {
        corrected,
        n_obs,
        n_vars,
    })
}

/// Run ComBat batch correction.
///
/// `data_json`: flat row-major JSON array of f64.
/// `n_features`: number of genes (columns).
/// `batch_json`: JSON array of batch label strings (one per cell).
/// Output: JSON `JsHarmonyResult` (same format as Harmony for simplicity).
#[cfg_attr(feature = "wasm", wasm_bindgen)]
pub fn sc_combat(data_json: &str, n_features: usize, batch_json: &str) -> String {
    use cyanea_omics::sc_integrate::{CombatConfig, combat};
    let data = match parse_f64_array(data_json) {
        Ok(d) => d,
        Err(e) => return wasm_err(e),
    };
    let mut adata = match build_adata(&data, n_features) {
        Ok(a) => a,
        Err(e) => return wasm_err(e),
    };
    let batches: Vec<String> = match serde_json::from_str(batch_json) {
        Ok(b) => b,
        Err(e) => return wasm_err(format!("invalid batch JSON: {e}")),
    };
    if let Err(e) = adata.add_obs("batch", batches) {
        return wasm_err(e);
    }
    let config = CombatConfig {
        batch_key: "batch".to_string(),
        parametric: true,
    };
    if let Err(e) = combat(&mut adata, &config) {
        return wasm_err(e);
    }
    let corrected = adata.x().to_flat_row_major();
    let (n_obs, n_vars) = adata.shape();
    wasm_ok(&JsHarmonyResult {
        corrected,
        n_obs,
        n_vars,
    })
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
