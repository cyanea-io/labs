//! CRISPR analysis: guide RNA design scoring, off-target prediction,
//! CRISPR screen analysis (MAGeCK-style), and base editing outcome prediction.
//!
//! # Overview
//!
//! - **Guide scoring**: Rule Set 2 (Doench 2016) and CFD scoring for on-target activity
//! - **Off-target prediction**: Mismatch-based search with position-weighted scoring
//! - **Screen analysis**: MAGeCK-style robust rank aggregation for gene essentiality
//! - **Base editing**: Predicted editing outcomes for CBE and ABE editors

use cyanea_core::{CyaneaError, Result};
use std::collections::HashMap;

/// A CRISPR guide RNA with scoring information.
#[derive(Debug, Clone)]
pub struct GuideRna {
    /// 20-nt spacer sequence (5' to 3').
    pub spacer: String,
    /// PAM sequence (e.g., "NGG" for SpCas9).
    pub pam: String,
    /// Genomic position (chromosome).
    pub chrom: String,
    /// Genomic position (start, 0-based).
    pub start: u64,
    /// Strand ('+' or '-').
    pub strand: char,
    /// On-target activity score (0-1, higher = better).
    pub on_target_score: f64,
    /// Number of predicted off-target sites.
    pub off_target_count: u32,
    /// Specificity score (0-1, higher = fewer off-targets).
    pub specificity_score: f64,
}

/// An off-target site for a guide RNA.
#[derive(Debug, Clone)]
pub struct OffTarget {
    /// Genomic sequence at this site (20-nt).
    pub sequence: String,
    /// Chromosome.
    pub chrom: String,
    /// Position.
    pub position: u64,
    /// Strand.
    pub strand: char,
    /// Number of mismatches to the guide.
    pub mismatches: u32,
    /// CFD score for this off-target (0-1, higher = more likely cut).
    pub cfd_score: f64,
}

/// Result of a CRISPR screen analysis for one gene.
#[derive(Debug, Clone)]
pub struct ScreenGeneResult {
    /// Gene name.
    pub gene: String,
    /// Number of guides targeting this gene.
    pub n_guides: usize,
    /// Median log2 fold-change across guides.
    pub median_lfc: f64,
    /// RRA (Robust Rank Aggregation) score.
    pub rra_score: f64,
    /// Adjusted p-value (BH-corrected).
    pub fdr: f64,
    /// Classification: "essential", "enriched", or "neutral".
    pub classification: String,
}

/// A base editing outcome prediction.
#[derive(Debug, Clone)]
pub struct EditingOutcome {
    /// Position within the spacer (1-indexed from PAM-distal).
    pub position: u32,
    /// Original base.
    pub ref_base: char,
    /// Edited base.
    pub alt_base: char,
    /// Predicted editing efficiency (0-1).
    pub efficiency: f64,
    /// Whether this position is in the editing window.
    pub in_window: bool,
}

/// Editor type for base editing prediction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseEditor {
    /// Cytosine base editor (C→T).
    Cbe,
    /// Adenine base editor (A→G).
    Abe,
}

// ---------------------------------------------------------------------------
// On-target scoring
// ---------------------------------------------------------------------------

/// Score a 30-nt context (4nt upstream + 20nt spacer + 3nt PAM + 3nt downstream)
/// using a simplified Rule Set 2 model (Doench et al., 2016).
///
/// The score is in [0, 1] where higher indicates better predicted on-target activity.
/// The context must be 30 nucleotides: NNNN[20-nt spacer]NGG NNN.
pub fn score_guide_rs2(context_30: &[u8]) -> Result<f64> {
    if context_30.len() != 30 {
        return Err(CyaneaError::InvalidInput(
            format!("Expected 30-nt context, got {} nt", context_30.len()),
        ));
    }

    // Simplified scoring based on known position-dependent nucleotide preferences
    let mut score: f64 = 0.5; // baseline

    // Position 20 (spacer position 20, closest to PAM): G preferred
    if context_30[23] == b'G' { score += 0.05; }
    // Position 16 of spacer: C preferred
    if context_30[19] == b'C' { score += 0.04; }

    // GC content of spacer (positions 4-23)
    let spacer = &context_30[4..24];
    let gc_count = spacer.iter().filter(|&&b| b == b'G' || b == b'C').count();
    let gc_frac = gc_count as f64 / 20.0;

    // Optimal GC: 40-70%
    if gc_frac >= 0.4 && gc_frac <= 0.7 {
        score += 0.10;
    } else if gc_frac < 0.3 || gc_frac > 0.8 {
        score -= 0.15;
    }

    // Penalize poly-T runs (Pol III terminator signal)
    let spacer_str = std::str::from_utf8(spacer).unwrap_or("");
    if spacer_str.contains("TTTT") {
        score -= 0.20;
    }

    // Position-specific preferences (simplified from Doench 2016)
    // Seed region (positions 1-8 from PAM, i.e., spacer positions 13-20)
    // Prefer no G at position 20 of spacer (already handled above)
    // Prefer G at position 1 of spacer
    if spacer[0] == b'G' { score += 0.03; }

    // PAM-proximal nucleotide preferences
    if spacer[19] == b'C' || spacer[19] == b'G' { score += 0.03; }
    if spacer[18] == b'A' { score += 0.02; }

    // Upstream context
    if context_30[3] == b'C' { score += 0.02; }

    score = score.clamp(0.0, 1.0);
    Ok(score)
}

// ---------------------------------------------------------------------------
// CFD off-target scoring
// ---------------------------------------------------------------------------

/// Cutting Frequency Determination (CFD) score for off-target activity.
///
/// Compares a 20-nt guide with a 20-nt off-target site.
/// Returns 0-1 where 1.0 = perfect match (same as on-target),
/// and lower values indicate less likely cutting.
pub fn cfd_score(guide: &[u8], off_target: &[u8]) -> Result<f64> {
    if guide.len() != 20 || off_target.len() != 20 {
        return Err(CyaneaError::InvalidInput("Guide and off-target must be 20 nt".into()));
    }

    let mut score: f64 = 1.0;

    for i in 0..20 {
        if guide[i] != off_target[i] {
            // Position-dependent mismatch penalty
            // Seed region (positions 1-12 from PAM = indices 8-19) penalized more heavily
            let pos_from_pam = 20 - i; // 1 = PAM-proximal

            let penalty: f64 = if pos_from_pam <= 4 {
                0.0   // Critical seed: any mismatch nearly abolishes cutting
            } else if pos_from_pam <= 8 {
                0.1   // Near-seed: heavy penalty
            } else if pos_from_pam <= 12 {
                0.3   // Mid-guide
            } else if pos_from_pam <= 16 {
                0.6   // PAM-distal
            } else {
                0.8   // Very distal: mild penalty
            };

            // Specific mismatch type modifiers (rG:dT wobble is more tolerated)
            let modifier: f64 = match (guide[i], off_target[i]) {
                (b'G', b'A') | (b'A', b'G') => 1.2,  // purine-purine: slightly worse
                (b'C', b'T') | (b'T', b'C') => 1.1,  // pyrimidine-pyrimidine
                _ => 1.0,
            };

            score *= penalty / modifier.max(0.01);
        }
    }

    Ok(score.clamp(0.0, 1.0))
}

/// Count mismatches between a guide and a target sequence.
pub fn count_mismatches(guide: &[u8], target: &[u8]) -> u32 {
    guide.iter().zip(target.iter()).filter(|(a, b)| a != b).count() as u32
}

// ---------------------------------------------------------------------------
// Off-target search
// ---------------------------------------------------------------------------

/// Find off-target sites for a guide RNA in a genome sequence.
///
/// Searches for PAM (NGG or NAG) sites with up to `max_mismatches`
/// mismatches in the 20-nt upstream region.
pub fn find_off_targets(
    guide: &[u8],
    genome: &[u8],
    chrom: &str,
    max_mismatches: u32,
) -> Result<Vec<OffTarget>> {
    if guide.len() != 20 {
        return Err(CyaneaError::InvalidInput("Guide must be 20 nt".into()));
    }

    let mut results = Vec::new();

    // Search forward strand for NGG PAMs
    for i in 0..genome.len().saturating_sub(22) {
        // Check for NGG at positions i+20..i+23
        if genome[i + 21] == b'G' && genome[i + 22] == b'G' {
            let target = &genome[i..i + 20];
            let mm = count_mismatches(guide, target);
            if mm <= max_mismatches && mm > 0 {
                let cfd = cfd_score(guide, target).unwrap_or(0.0);
                results.push(OffTarget {
                    sequence: String::from_utf8_lossy(target).to_string(),
                    chrom: chrom.to_string(),
                    position: i as u64,
                    strand: '+',
                    mismatches: mm,
                    cfd_score: cfd,
                });
            }
        }
    }

    // Search reverse strand (CCN PAM on forward = NGG on reverse)
    let rc_guide = reverse_complement(guide);
    for i in 0..genome.len().saturating_sub(22) {
        if genome[i] == b'C' && genome[i + 1] == b'C' {
            let target_start = i + 3;
            if target_start + 20 <= genome.len() {
                let target = &genome[target_start..target_start + 20];
                let mm = count_mismatches(&rc_guide, target);
                if mm <= max_mismatches && mm > 0 {
                    let cfd = cfd_score(&rc_guide, target).unwrap_or(0.0);
                    results.push(OffTarget {
                        sequence: String::from_utf8_lossy(target).to_string(),
                        chrom: chrom.to_string(),
                        position: target_start as u64,
                        strand: '-',
                        mismatches: mm,
                        cfd_score: cfd,
                    });
                }
            }
        }
    }

    results.sort_by(|a, b| b.cfd_score.partial_cmp(&a.cfd_score).unwrap_or(core::cmp::Ordering::Equal));
    Ok(results)
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| match b {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        _ => b'N',
    }).collect()
}

// ---------------------------------------------------------------------------
// CRISPR screen analysis (MAGeCK-style)
// ---------------------------------------------------------------------------

/// Analyze a CRISPR screen using robust rank aggregation (MAGeCK-style).
///
/// * `guide_lfc` — log2 fold-changes per guide: `(guide_name, gene, lfc)`
/// * `negative` — if true, test for negative selection (essentiality); if false, positive selection
///
/// Returns per-gene results sorted by RRA score.
pub fn analyze_screen(
    guide_lfc: &[(String, String, f64)],
    negative: bool,
) -> Vec<ScreenGeneResult> {
    // Group guides by gene
    let mut gene_guides: HashMap<String, Vec<f64>> = HashMap::new();
    for (_, gene, lfc) in guide_lfc {
        gene_guides.entry(gene.clone()).or_default().push(*lfc);
    }

    // Rank all guides globally
    let mut all_lfcs: Vec<f64> = guide_lfc.iter().map(|(_, _, lfc)| *lfc).collect();
    let n_total = all_lfcs.len();

    // For negative selection, rank by ascending LFC (most depleted = rank 1)
    // For positive selection, rank by descending LFC (most enriched = rank 1)
    if negative {
        all_lfcs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
    } else {
        all_lfcs.sort_by(|a, b| b.partial_cmp(a).unwrap_or(core::cmp::Ordering::Equal));
    }

    // Create LFC → rank mapping
    let lfc_ranks: Vec<(f64, f64)> = all_lfcs.iter().enumerate()
        .map(|(i, &lfc)| (lfc, (i + 1) as f64 / n_total as f64))
        .collect();

    let mut results = Vec::new();

    for (gene, lfcs) in &gene_guides {
        let n_guides = lfcs.len();

        // Compute per-guide percentile ranks
        let mut guide_ranks: Vec<f64> = lfcs.iter().map(|lfc| {
            // Find closest rank
            lfc_ranks.iter()
                .min_by(|a, b| (a.0 - lfc).abs().partial_cmp(&(b.0 - lfc).abs()).unwrap_or(core::cmp::Ordering::Equal))
                .map(|r| r.1)
                .unwrap_or(0.5)
        }).collect();

        guide_ranks.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));

        // RRA: alpha-robust rank aggregation
        // For each guide rank r_k, compute p(r_k, k, n) = beta CDF approximation
        // Use the minimum of these as the RRA score
        let alpha = 0.25;
        let mut min_rra: f64 = 1.0;

        for (k, &rank) in guide_ranks.iter().enumerate() {
            if rank > alpha {
                continue;
            }
            // Simplified: uniform order statistic p-value
            // P(U_(k) <= r) ≈ Σ C(n,j) * r^j * (1-r)^(n-j) for j=k..n
            // Approximate with: r^(k+1) * C(n, k+1) for small r
            let p_val = rank.powi((k + 1) as i32) * choose(n_guides, k + 1);
            min_rra = min_rra.min(p_val);
        }

        let median_lfc = {
            let mut sorted_lfcs = lfcs.clone();
            sorted_lfcs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
            if sorted_lfcs.len() % 2 == 0 {
                (sorted_lfcs[sorted_lfcs.len() / 2 - 1] + sorted_lfcs[sorted_lfcs.len() / 2]) / 2.0
            } else {
                sorted_lfcs[sorted_lfcs.len() / 2]
            }
        };

        results.push(ScreenGeneResult {
            gene: gene.clone(),
            n_guides,
            median_lfc,
            rra_score: min_rra,
            fdr: 0.0, // placeholder, set after BH correction
            classification: String::new(),
        });
    }

    // Sort by RRA score
    results.sort_by(|a, b| a.rra_score.partial_cmp(&b.rra_score).unwrap_or(core::cmp::Ordering::Equal));

    // BH correction
    let n_genes = results.len();
    for (i, r) in results.iter_mut().enumerate() {
        r.fdr = (r.rra_score * n_genes as f64 / (i + 1) as f64).min(1.0);
    }

    // Enforce monotonicity (FDR must be non-decreasing)
    let mut min_fdr: f64 = 1.0;
    for r in results.iter_mut().rev() {
        min_fdr = min_fdr.min(r.fdr);
        r.fdr = min_fdr;
    }

    // Classify
    for r in &mut results {
        r.classification = if r.fdr < 0.05 {
            if negative {
                "essential".to_string()
            } else {
                "enriched".to_string()
            }
        } else {
            "neutral".to_string()
        };
    }

    results
}

fn choose(n: usize, k: usize) -> f64 {
    if k > n { return 0.0; }
    let mut result = 1.0;
    for i in 0..k {
        result *= (n - i) as f64 / (i + 1) as f64;
    }
    result
}

// ---------------------------------------------------------------------------
// Base editing outcome prediction
// ---------------------------------------------------------------------------

/// Predict base editing outcomes for a 20-nt spacer.
///
/// * `editor` — CBE (C→T) or ABE (A→G)
/// * `spacer` — 20-nt guide sequence
///
/// The editing window is typically positions 4-8 (CBE) or 4-7 (ABE)
/// from the PAM-distal end (1-indexed).
pub fn predict_editing(editor: BaseEditor, spacer: &[u8]) -> Result<Vec<EditingOutcome>> {
    if spacer.len() != 20 {
        return Err(CyaneaError::InvalidInput("Spacer must be 20 nt".into()));
    }

    let (target_base, edited_base, window_start, window_end) = match editor {
        BaseEditor::Cbe => (b'C', b'T', 4, 8),
        BaseEditor::Abe => (b'A', b'G', 4, 7),
    };

    let mut outcomes = Vec::new();

    for (i, &base) in spacer.iter().enumerate() {
        let pos = i + 1; // 1-indexed
        if base == target_base {
            let in_window = pos >= window_start && pos <= window_end;

            // Position-dependent efficiency within the window
            let efficiency = if in_window {
                let center = (window_start + window_end) as f64 / 2.0;
                let dist = (pos as f64 - center).abs();
                let max_dist = (window_end - window_start) as f64 / 2.0;
                // Gaussian-like falloff from center of window
                let base_eff = 0.7;
                base_eff * (-0.5 * (dist / max_dist.max(1.0)).powi(2)).exp()
            } else {
                // Bystander editing at low efficiency outside window
                let dist_to_window = if pos < window_start {
                    (window_start - pos) as f64
                } else {
                    (pos - window_end) as f64
                };
                0.05 * (-dist_to_window).exp()
            };

            outcomes.push(EditingOutcome {
                position: pos as u32,
                ref_base: base as char,
                alt_base: edited_base as char,
                efficiency,
                in_window,
            });
        }
    }

    Ok(outcomes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_score_guide_rs2_valid() {
        // Good guide: decent GC, no polyT, G at position 20 (30 nt context)
        let ctx30 = b"ACGTGCATGCTAGCTAGCTAGCAGNNGGAT";
        let score = score_guide_rs2(ctx30).unwrap();
        assert!(score > 0.0 && score <= 1.0);
    }

    #[test]
    fn test_score_guide_rs2_polyt() {
        // Guide with TTTT run should score lower (both 30 nt)
        let good = b"ACGTGCATGCTAGCTAGCGATGGNNGGATT";
        let bad  = b"ACGTTTTTGCTAGCTAGCGATGGNNGGATT";
        let s1 = score_guide_rs2(good).unwrap();
        let s2 = score_guide_rs2(bad).unwrap();
        assert!(s1 > s2);
    }

    #[test]
    fn test_score_guide_wrong_length() {
        assert!(score_guide_rs2(b"ACGT").is_err());
    }

    #[test]
    fn test_cfd_score_perfect() {
        let guide = b"ACGTACGTACGTACGTACGT";
        let score = cfd_score(guide, guide).unwrap();
        assert!((score - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_cfd_score_mismatches() {
        let guide  = b"ACGTACGTACGTACGTACGT";
        // One mismatch near PAM (position 20) — should heavily penalize
        let ot_seed = b"ACGTACGTACGTACGTACGA";
        // One mismatch far from PAM (position 1)
        let ot_distal = b"GCGTACGTACGTACGTACGT";

        let s_seed = cfd_score(guide, ot_seed).unwrap();
        let s_distal = cfd_score(guide, ot_distal).unwrap();

        assert!(s_distal > s_seed, "Distal mismatch should be more tolerated");
    }

    #[test]
    fn test_count_mismatches() {
        assert_eq!(count_mismatches(b"ACGTACGT", b"ACGTACGT"), 0);
        assert_eq!(count_mismatches(b"ACGTACGT", b"ACGTACGA"), 1);
        assert_eq!(count_mismatches(b"AAAA", b"TTTT"), 4);
    }

    #[test]
    fn test_find_off_targets() {
        let guide = b"ACGTACGTACGTACGTACGT";
        // Genome with one near-match site (1 mismatch at position 10) followed by NGG
        let mut genome = Vec::new();
        genome.extend_from_slice(b"NNNNNNNNNN"); // padding
        genome.extend_from_slice(b"ACGTACGTACATACGTACGT"); // 1 mismatch at pos 10 (G→A)
        genome.extend_from_slice(b"NGG"); // PAM
        genome.extend_from_slice(b"NNNNNNNNNN"); // padding

        let ots = find_off_targets(guide, &genome, "chr1", 3).unwrap();
        assert!(!ots.is_empty());
        assert_eq!(ots[0].mismatches, 1);
    }

    #[test]
    fn test_find_off_targets_no_match() {
        let guide = b"ACGTACGTACGTACGTACGT";
        // Genome with no NGG PAM sites near similar sequences
        let genome = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let ots = find_off_targets(guide, genome, "chr1", 3).unwrap();
        assert!(ots.is_empty());
    }

    #[test]
    fn test_analyze_screen_negative() {
        let guide_data: Vec<(String, String, f64)> = vec![
            // Essential gene — all guides depleted
            ("g1".into(), "TP53".into(), -3.5),
            ("g2".into(), "TP53".into(), -2.8),
            ("g3".into(), "TP53".into(), -4.1),
            ("g4".into(), "TP53".into(), -3.0),
            // Neutral gene
            ("g5".into(), "GFP".into(), 0.1),
            ("g6".into(), "GFP".into(), -0.2),
            ("g7".into(), "GFP".into(), 0.3),
            ("g8".into(), "GFP".into(), -0.1),
            // Another neutral
            ("g9".into(), "LUC".into(), 0.0),
            ("g10".into(), "LUC".into(), 0.1),
            ("g11".into(), "LUC".into(), -0.3),
            ("g12".into(), "LUC".into(), 0.2),
        ];

        let results = analyze_screen(&guide_data, true);
        assert_eq!(results.len(), 3);

        // TP53 should rank first (most depleted)
        assert_eq!(results[0].gene, "TP53");
        assert!(results[0].median_lfc < -2.0);
        assert!(results[0].rra_score < results[1].rra_score);
    }

    #[test]
    fn test_analyze_screen_positive() {
        let guide_data: Vec<(String, String, f64)> = vec![
            ("g1".into(), "KRAS".into(), 4.0),
            ("g2".into(), "KRAS".into(), 3.5),
            ("g3".into(), "KRAS".into(), 3.8),
            ("g4".into(), "CTRL1".into(), 0.1),
            ("g5".into(), "CTRL1".into(), -0.1),
            ("g6".into(), "CTRL1".into(), 0.2),
        ];

        let results = analyze_screen(&guide_data, false);
        assert_eq!(results[0].gene, "KRAS");
        assert!(results[0].median_lfc > 3.0);
    }

    #[test]
    fn test_predict_editing_cbe() {
        let spacer = b"AAGCACGTACGTACGTACGT";
        let outcomes = predict_editing(BaseEditor::Cbe, spacer).unwrap();

        // Should find C bases and mark those in window (positions 4-8)
        assert!(!outcomes.is_empty());
        for o in &outcomes {
            assert_eq!(o.ref_base, 'C');
            assert_eq!(o.alt_base, 'T');
        }

        // Position 4 should be in window
        let pos4 = outcomes.iter().find(|o| o.position == 4);
        if let Some(p) = pos4 {
            assert!(p.in_window);
            assert!(p.efficiency > 0.1);
        }
    }

    #[test]
    fn test_predict_editing_abe() {
        let spacer = b"ACAGACGTACGTACGTACGT";
        let outcomes = predict_editing(BaseEditor::Abe, spacer).unwrap();

        assert!(!outcomes.is_empty());
        for o in &outcomes {
            assert_eq!(o.ref_base, 'A');
            assert_eq!(o.alt_base, 'G');
        }
    }

    #[test]
    fn test_predict_editing_wrong_length() {
        assert!(predict_editing(BaseEditor::Cbe, b"ACGT").is_err());
    }

    #[test]
    fn test_predict_editing_window_efficiency() {
        // Spacer with C at positions 4, 5, 6 (in window) and 15 (out of window)
        let spacer = b"TTTCCCGGGGGGGGGCGGGG";
        let outcomes = predict_editing(BaseEditor::Cbe, spacer).unwrap();

        let in_win: Vec<_> = outcomes.iter().filter(|o| o.in_window).collect();
        let out_win: Vec<_> = outcomes.iter().filter(|o| !o.in_window).collect();

        assert!(!in_win.is_empty());
        assert!(!out_win.is_empty());

        let max_in = in_win.iter().map(|o| o.efficiency).fold(0.0f64, f64::max);
        let max_out = out_win.iter().map(|o| o.efficiency).fold(0.0f64, f64::max);
        assert!(max_in > max_out, "In-window efficiency should exceed out-of-window");
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"GCAT"), b"ATGC");
    }
}
