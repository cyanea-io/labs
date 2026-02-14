//! UCSC chain file parsing and genomic coordinate liftover.
//!
//! Remaps coordinates between genome assemblies using UCSC chain files.
//!
//! # Chain file format
//!
//! A chain file contains pairwise alignments between a source (target) and
//! query genome. Each alignment is described by a header line followed by
//! data lines specifying ungapped alignment blocks and the gaps between them.

use std::collections::BTreeMap;

use cyanea_core::{CyaneaError, Result};

use crate::genomic::{GenomicInterval, Strand};

/// A single ungapped alignment block within a chain.
#[derive(Debug, Clone)]
#[allow(dead_code)]
struct AlignmentBlock {
    source_start: u64,
    source_end: u64,
    target_start: u64,
    target_end: u64,
}

/// A parsed chain alignment.
#[derive(Debug, Clone)]
#[allow(dead_code)]
struct Chain {
    score: u64,
    source_chrom: String,
    source_size: u64,
    source_strand: Strand,
    source_start: u64,
    source_end: u64,
    target_chrom: String,
    target_size: u64,
    target_strand: Strand,
    target_start: u64,
    target_end: u64,
    blocks: Vec<AlignmentBlock>,
}

/// A parsed chain file, indexed by source chromosome.
#[derive(Debug, Clone)]
pub struct ChainFile {
    chains: BTreeMap<String, Vec<Chain>>,
}

/// Result of a liftover operation for a single interval.
#[derive(Debug, Clone, PartialEq)]
pub enum LiftoverResult {
    /// The interval was fully mapped to the target genome.
    Mapped(GenomicInterval),
    /// The interval was partially mapped.
    Partial {
        mapped: GenomicInterval,
        fraction_mapped: f64,
    },
    /// The interval could not be mapped.
    Unmapped,
}

fn parse_strand(s: &str) -> Result<Strand> {
    match s {
        "+" => Ok(Strand::Forward),
        "-" => Ok(Strand::Reverse),
        _ => Err(CyaneaError::Parse(format!(
            "invalid strand '{s}'"
        ))),
    }
}

/// Parse a UCSC chain file from a string.
///
/// # Format
///
/// ```text
/// chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
/// size dt dq
/// size dt dq
/// size
/// ```
///
/// The header uses "target" for the reference (source) genome and "query"
/// for the destination, matching UCSC convention. In our API, "source"
/// corresponds to the chain target (the genome you're lifting FROM).
pub fn parse_chain(input: &str) -> Result<ChainFile> {
    let mut chains: BTreeMap<String, Vec<Chain>> = BTreeMap::new();
    let mut current_chain: Option<Chain> = None;
    let mut source_offset = 0u64;
    let mut target_offset = 0u64;

    for line in input.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }

        if line.starts_with("chain ") {
            // Finalize previous chain
            if let Some(chain) = current_chain.take() {
                chains
                    .entry(chain.source_chrom.clone())
                    .or_default()
                    .push(chain);
            }

            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() < 12 {
                return Err(CyaneaError::Parse(format!(
                    "chain header requires at least 12 fields, got {}: {line}",
                    fields.len()
                )));
            }

            let score: u64 = fields[1]
                .parse()
                .map_err(|_| CyaneaError::Parse(format!("invalid score: {}", fields[1])))?;

            // Target = source in our naming (genome we're lifting FROM)
            let source_chrom = fields[2].to_string();
            let source_size: u64 = fields[3].parse().map_err(|_| {
                CyaneaError::Parse(format!("invalid source size: {}", fields[3]))
            })?;
            let source_strand = parse_strand(fields[4])?;
            let source_start: u64 = fields[5].parse().map_err(|_| {
                CyaneaError::Parse(format!("invalid source start: {}", fields[5]))
            })?;
            let source_end: u64 = fields[6].parse().map_err(|_| {
                CyaneaError::Parse(format!("invalid source end: {}", fields[6]))
            })?;

            // Query = target in our naming (genome we're lifting TO)
            let target_chrom = fields[7].to_string();
            let target_size: u64 = fields[8].parse().map_err(|_| {
                CyaneaError::Parse(format!("invalid target size: {}", fields[8]))
            })?;
            let target_strand = parse_strand(fields[9])?;
            let target_start: u64 = fields[10].parse().map_err(|_| {
                CyaneaError::Parse(format!("invalid target start: {}", fields[10]))
            })?;
            let target_end: u64 = fields[11].parse().map_err(|_| {
                CyaneaError::Parse(format!("invalid target end: {}", fields[11]))
            })?;

            source_offset = source_start;
            target_offset = target_start;

            current_chain = Some(Chain {
                score,
                source_chrom,
                source_size,
                source_strand,
                source_start,
                source_end,
                target_chrom,
                target_size,
                target_strand,
                target_start,
                target_end,
                blocks: Vec::new(),
            });

            continue;
        }

        // Data line: size [dt dq]
        if let Some(ref mut chain) = current_chain {
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.is_empty() {
                continue;
            }

            let size: u64 = fields[0].parse().map_err(|_| {
                CyaneaError::Parse(format!("invalid block size: {}", fields[0]))
            })?;

            chain.blocks.push(AlignmentBlock {
                source_start: source_offset,
                source_end: source_offset + size,
                target_start: target_offset,
                target_end: target_offset + size,
            });

            if fields.len() >= 3 {
                let dt: u64 = fields[1].parse().map_err(|_| {
                    CyaneaError::Parse(format!("invalid dt: {}", fields[1]))
                })?;
                let dq: u64 = fields[2].parse().map_err(|_| {
                    CyaneaError::Parse(format!("invalid dq: {}", fields[2]))
                })?;
                source_offset += size + dt;
                target_offset += size + dq;
            }
            // Last line of chain has only `size` — no gap after
        }
    }

    // Finalize last chain
    if let Some(chain) = current_chain.take() {
        chains
            .entry(chain.source_chrom.clone())
            .or_default()
            .push(chain);
    }

    // Sort chains by score (descending) for each chromosome
    for chains_vec in chains.values_mut() {
        chains_vec.sort_by(|a, b| b.score.cmp(&a.score));
    }

    Ok(ChainFile { chains })
}

/// Convert a reverse-strand coordinate to forward-strand.
fn reverse_to_forward(pos: u64, chrom_size: u64) -> u64 {
    chrom_size - pos
}

/// Liftover a single interval using the chain file.
///
/// `min_match` is the minimum fraction of the interval that must map (0.0–1.0).
/// Use 0.95 as a reasonable default.
pub fn liftover(
    chain_file: &ChainFile,
    interval: &GenomicInterval,
    min_match: f64,
) -> LiftoverResult {
    let chains = match chain_file.chains.get(&interval.chrom) {
        Some(c) => c,
        None => return LiftoverResult::Unmapped,
    };

    // Find the best chain that overlaps this interval (chains are sorted by score desc)
    let chain = match chains.iter().find(|c| {
        c.source_start <= interval.start && c.source_end >= interval.end
    }) {
        Some(c) => c,
        None => {
            // Try partial overlap with best-scoring chain
            match chains.iter().find(|c| {
                c.source_start < interval.end && c.source_end > interval.start
            }) {
                Some(c) => c,
                None => return LiftoverResult::Unmapped,
            }
        }
    };

    // Find alignment blocks that overlap the query interval
    let mut mapped_bp = 0u64;
    let mut target_min = u64::MAX;
    let mut target_max = 0u64;
    let interval_len = interval.len();

    for block in &chain.blocks {
        if block.source_end <= interval.start || block.source_start >= interval.end {
            continue;
        }

        // Compute overlap
        let overlap_start = block.source_start.max(interval.start);
        let overlap_end = block.source_end.min(interval.end);
        let overlap_len = overlap_end - overlap_start;
        mapped_bp += overlap_len;

        // Map the overlap to target coordinates
        let offset_in_block = overlap_start - block.source_start;
        let t_start = block.target_start + offset_in_block;
        let t_end = t_start + overlap_len;

        target_min = target_min.min(t_start);
        target_max = target_max.max(t_end);
    }

    if mapped_bp == 0 {
        return LiftoverResult::Unmapped;
    }

    let fraction = mapped_bp as f64 / interval_len as f64;

    // Handle reverse strand chains: convert target coordinates
    let (final_start, final_end) = if chain.target_strand == Strand::Reverse {
        let fwd_start = reverse_to_forward(target_max, chain.target_size);
        let fwd_end = reverse_to_forward(target_min, chain.target_size);
        (fwd_start, fwd_end)
    } else {
        (target_min, target_max)
    };

    // Construct the mapped interval (use the chain's target chrom)
    let mapped = GenomicInterval {
        chrom: chain.target_chrom.clone(),
        start: final_start,
        end: final_end,
        strand: interval.strand,
    };

    if fraction >= min_match {
        if (fraction - 1.0).abs() < 1e-10 {
            LiftoverResult::Mapped(mapped)
        } else {
            LiftoverResult::Partial {
                mapped,
                fraction_mapped: fraction,
            }
        }
    } else {
        LiftoverResult::Unmapped
    }
}

/// Batch liftover of multiple intervals.
pub fn liftover_batch(
    chain_file: &ChainFile,
    intervals: &[GenomicInterval],
    min_match: f64,
) -> Vec<LiftoverResult> {
    intervals
        .iter()
        .map(|iv| liftover(chain_file, iv, min_match))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(chrom: &str, start: u64, end: u64) -> GenomicInterval {
        GenomicInterval::new(chrom, start, end).unwrap()
    }

    // A simple chain: chr1 (1000bp) -> chrA (800bp), forward strand
    // One chain with two alignment blocks:
    //   Block 1: source [100, 300) -> target [50, 250)
    //   Gap: dt=100, dq=50
    //   Block 2: source [400, 600) -> target [300, 500)
    const SIMPLE_CHAIN: &str = "\
chain 1000 chr1 1000 + 100 600 chrA 800 + 50 500 1
200 100 50
200
";

    #[test]
    fn test_parse_basic_chain() {
        let cf = parse_chain(SIMPLE_CHAIN).unwrap();
        assert!(cf.chains.contains_key("chr1"));
        let chains = &cf.chains["chr1"];
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].score, 1000);
        assert_eq!(chains[0].source_chrom, "chr1");
        assert_eq!(chains[0].target_chrom, "chrA");
        assert_eq!(chains[0].blocks.len(), 2);

        // Block 1
        assert_eq!(chains[0].blocks[0].source_start, 100);
        assert_eq!(chains[0].blocks[0].source_end, 300);
        assert_eq!(chains[0].blocks[0].target_start, 50);
        assert_eq!(chains[0].blocks[0].target_end, 250);

        // Block 2
        assert_eq!(chains[0].blocks[1].source_start, 400);
        assert_eq!(chains[0].blocks[1].source_end, 600);
        assert_eq!(chains[0].blocks[1].target_start, 300);
        assert_eq!(chains[0].blocks[1].target_end, 500);
    }

    #[test]
    fn test_parse_multiple_chains() {
        let input = "\
chain 5000 chr1 1000 + 0 500 chrA 800 + 0 400 1
200 100 50
200

chain 3000 chr2 2000 + 100 400 chrB 1500 + 200 500 2
300
";
        let cf = parse_chain(input).unwrap();
        assert!(cf.chains.contains_key("chr1"));
        assert!(cf.chains.contains_key("chr2"));
        assert_eq!(cf.chains["chr2"][0].blocks.len(), 1);
    }

    #[test]
    fn test_parse_reverse_strand() {
        let input = "\
chain 2000 chr1 1000 + 100 400 chrA 800 - 300 600 1
300
";
        let cf = parse_chain(input).unwrap();
        let chain = &cf.chains["chr1"][0];
        assert_eq!(chain.source_strand, Strand::Forward);
        assert_eq!(chain.target_strand, Strand::Reverse);
    }

    #[test]
    fn test_parse_invalid_format() {
        let input = "chain 1000 chr1";
        assert!(parse_chain(input).is_err());
    }

    #[test]
    fn test_liftover_within_single_block() {
        let cf = parse_chain(SIMPLE_CHAIN).unwrap();
        // Query [150, 250) is within block 1: source [100,300) -> target [50,250)
        // offset = 150-100 = 50, so target [100, 200)
        let result = liftover(&cf, &iv("chr1", 150, 250), 0.95);
        match result {
            LiftoverResult::Mapped(mapped) => {
                assert_eq!(mapped.chrom, "chrA");
                assert_eq!(mapped.start, 100);
                assert_eq!(mapped.end, 200);
            }
            other => panic!("expected Mapped, got {:?}", other),
        }
    }

    #[test]
    fn test_liftover_spanning_blocks() {
        let cf = parse_chain(SIMPLE_CHAIN).unwrap();
        // Query [200, 500) spans block 1 [100,300) and block 2 [400,600)
        // Block 1 overlap: [200,300) = 100bp, mapped to target [150,250)
        // Gap: source [300,400) = unmapped
        // Block 2 overlap: [400,500) = 100bp, mapped to target [300,400)
        // Total mapped = 200bp out of 300bp = 0.667
        let result = liftover(&cf, &iv("chr1", 200, 500), 0.5);
        match result {
            LiftoverResult::Partial {
                mapped,
                fraction_mapped,
            } => {
                assert_eq!(mapped.chrom, "chrA");
                assert_eq!(mapped.start, 150);
                assert_eq!(mapped.end, 400);
                assert!((fraction_mapped - 200.0 / 300.0).abs() < 1e-10);
            }
            other => panic!("expected Partial, got {:?}", other),
        }
    }

    #[test]
    fn test_liftover_in_gap() {
        let cf = parse_chain(SIMPLE_CHAIN).unwrap();
        // Query [310, 390) is in the gap between blocks
        let result = liftover(&cf, &iv("chr1", 310, 390), 0.95);
        assert_eq!(result, LiftoverResult::Unmapped);
    }

    #[test]
    fn test_liftover_reverse_strand() {
        let input = "\
chain 2000 chr1 1000 + 100 400 chrA 800 - 300 600 1
300
";
        let cf = parse_chain(input).unwrap();
        // Source [100,400) -> target reverse strand [300,600) on chrA (size=800)
        // Block: source [100,400) -> target [300,600) in reverse coords
        // Query [150,250) -> offset=50, target_start=350, target_end=450
        // Reverse conversion: fwd_start = 800-450=350, fwd_end = 800-350=450
        let result = liftover(&cf, &iv("chr1", 150, 250), 0.95);
        match result {
            LiftoverResult::Mapped(mapped) => {
                assert_eq!(mapped.chrom, "chrA");
                assert_eq!(mapped.start, 350);
                assert_eq!(mapped.end, 450);
            }
            other => panic!("expected Mapped, got {:?}", other),
        }
    }

    #[test]
    fn test_liftover_no_chain_for_chrom() {
        let cf = parse_chain(SIMPLE_CHAIN).unwrap();
        let result = liftover(&cf, &iv("chrX", 100, 200), 0.95);
        assert_eq!(result, LiftoverResult::Unmapped);
    }

    #[test]
    fn test_liftover_min_match_threshold() {
        let cf = parse_chain(SIMPLE_CHAIN).unwrap();
        // [200, 500) maps 200/300 = 0.667
        let strict = liftover(&cf, &iv("chr1", 200, 500), 0.95);
        assert_eq!(strict, LiftoverResult::Unmapped);

        let lenient = liftover(&cf, &iv("chr1", 200, 500), 0.5);
        assert!(matches!(lenient, LiftoverResult::Partial { .. }));
    }

    #[test]
    fn test_liftover_batch_consistency() {
        let cf = parse_chain(SIMPLE_CHAIN).unwrap();
        let intervals = vec![iv("chr1", 150, 250), iv("chr1", 310, 390), iv("chrX", 0, 100)];
        let results = liftover_batch(&cf, &intervals, 0.95);
        assert_eq!(results.len(), 3);
        assert!(matches!(results[0], LiftoverResult::Mapped(_)));
        assert_eq!(results[1], LiftoverResult::Unmapped);
        assert_eq!(results[2], LiftoverResult::Unmapped);
    }

    #[test]
    fn test_parse_empty_input() {
        let cf = parse_chain("").unwrap();
        assert!(cf.chains.is_empty());
    }

    #[test]
    fn test_parse_chain_score_ordering() {
        let input = "\
chain 1000 chr1 2000 + 0 500 chrA 1000 + 0 400 1
200 100 50
200

chain 5000 chr1 2000 + 0 800 chrB 1500 + 0 700 2
400 100 50
300
";
        let cf = parse_chain(input).unwrap();
        let chains = &cf.chains["chr1"];
        assert_eq!(chains.len(), 2);
        // Higher score should be first
        assert_eq!(chains[0].score, 5000);
        assert_eq!(chains[1].score, 1000);
    }
}
