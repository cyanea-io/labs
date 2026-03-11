//! Protein inference from peptide-spectrum matches.
//!
//! Implements parsimony-based protein inference: given a set of PSMs,
//! determine the minimal set of proteins that explains all identified peptides.

use crate::error::Result;
use crate::search::Psm;

/// A protein group inferred from PSMs.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ProteinGroup {
    /// Protein accession(s) in this group.
    pub accessions: Vec<String>,
    /// Peptide sequences uniquely mapping to this group.
    pub unique_peptides: Vec<String>,
    /// Peptide sequences shared with other groups.
    pub shared_peptides: Vec<String>,
    /// Total number of PSMs mapping to this group.
    pub psm_count: usize,
    /// Best hyperscore among all PSMs.
    pub best_score: f64,
    /// Sequence coverage (fraction of protein covered by peptides, 0.0-1.0).
    pub coverage: f64,
    /// Whether this is a decoy protein.
    pub is_decoy: bool,
}

/// A protein entry in the database (for inference).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ProteinEntry {
    /// Accession / identifier.
    pub accession: String,
    /// Full protein sequence.
    pub sequence: Vec<u8>,
    /// Whether this is a decoy entry.
    pub is_decoy: bool,
}

impl ProteinEntry {
    pub fn new(accession: impl Into<String>, sequence: &[u8]) -> Self {
        Self {
            accession: accession.into(),
            sequence: sequence.to_vec(),
            is_decoy: false,
        }
    }

    /// Check if a peptide sequence is contained in this protein.
    pub fn contains_peptide(&self, peptide: &[u8]) -> bool {
        self.sequence.windows(peptide.len()).any(|w| w == peptide)
    }

    /// Calculate sequence coverage given a set of peptide sequences.
    pub fn coverage(&self, peptides: &[&[u8]]) -> f64 {
        if self.sequence.is_empty() {
            return 0.0;
        }
        let mut covered = vec![false; self.sequence.len()];
        for peptide in peptides {
            for i in 0..self.sequence.len().saturating_sub(peptide.len()) + 1 {
                if &self.sequence[i..i + peptide.len()] == *peptide {
                    for j in i..i + peptide.len() {
                        covered[j] = true;
                    }
                }
            }
        }
        covered.iter().filter(|&&c| c).count() as f64 / self.sequence.len() as f64
    }
}

/// Perform parsimony-based protein inference.
///
/// Given PSMs and a protein database, finds the minimal set of protein groups
/// that explains all identified peptides.
///
/// # Algorithm
///
/// 1. Map each peptide to its parent proteins
/// 2. Identify unique peptides (mapping to exactly one protein)
/// 3. Greedily select proteins that explain the most unassigned peptides
/// 4. Group proteins with identical peptide sets (indistinguishable)
pub fn infer_proteins(
    psms: &[Psm],
    proteins: &[ProteinEntry],
) -> Result<Vec<ProteinGroup>> {
    if psms.is_empty() {
        return Ok(Vec::new());
    }

    // Collect unique peptide sequences from PSMs
    let mut peptide_set: Vec<String> = psms.iter()
        .map(|p| p.peptide_sequence.clone())
        .collect();
    peptide_set.sort();
    peptide_set.dedup();

    // Map each peptide to containing proteins
    let mut peptide_to_proteins: std::collections::HashMap<String, Vec<usize>> =
        std::collections::HashMap::new();

    for (pi, protein) in proteins.iter().enumerate() {
        for pep_seq in &peptide_set {
            if protein.contains_peptide(pep_seq.as_bytes()) {
                peptide_to_proteins.entry(pep_seq.clone())
                    .or_default()
                    .push(pi);
            }
        }
    }

    // Group proteins by their peptide sets (indistinguishable proteins)
    let mut protein_peptides: std::collections::HashMap<usize, Vec<String>> =
        std::collections::HashMap::new();

    for (pep, prot_indices) in &peptide_to_proteins {
        for &pi in prot_indices {
            protein_peptides.entry(pi)
                .or_default()
                .push(pep.clone());
        }
    }

    // Sort peptide lists for comparison
    for peps in protein_peptides.values_mut() {
        peps.sort();
    }

    // Group indistinguishable proteins
    let mut groups: std::collections::HashMap<Vec<String>, Vec<usize>> =
        std::collections::HashMap::new();

    for (&pi, peps) in &protein_peptides {
        groups.entry(peps.clone())
            .or_default()
            .push(pi);
    }

    // Greedy parsimony: select groups covering the most unassigned peptides
    let mut unassigned: std::collections::HashSet<String> =
        peptide_set.iter().cloned().collect();
    let mut selected_groups = Vec::new();

    let mut group_list: Vec<(Vec<String>, Vec<usize>)> = groups.into_iter().collect();

    while !unassigned.is_empty() {
        // Find group covering the most unassigned peptides
        let best = group_list.iter().enumerate()
            .max_by_key(|(_, (peps, _))| {
                peps.iter().filter(|p| unassigned.contains(*p)).count()
            });

        match best {
            Some((idx, (peps, _))) => {
                let covered: usize = peps.iter().filter(|p| unassigned.contains(*p)).count();
                if covered == 0 {
                    break; // No more peptides can be assigned
                }
                let (peps, prot_indices) = group_list.remove(idx);
                for p in &peps {
                    unassigned.remove(p);
                }
                selected_groups.push((peps, prot_indices));
            }
            None => break,
        }
    }

    // Build ProteinGroup results
    let mut results = Vec::new();

    for (peps, prot_indices) in selected_groups {
        let accessions: Vec<String> = prot_indices.iter()
            .map(|&pi| proteins[pi].accession.clone())
            .collect();

        // Determine unique vs shared peptides
        let unique: Vec<String> = peps.iter()
            .filter(|p| {
                peptide_to_proteins.get(*p)
                    .map_or(true, |prots| prots.len() == 1)
            })
            .cloned()
            .collect();

        let shared: Vec<String> = peps.iter()
            .filter(|p| {
                peptide_to_proteins.get(*p)
                    .map_or(false, |prots| prots.len() > 1)
            })
            .cloned()
            .collect();

        // Count PSMs and best score
        let mut psm_count = 0usize;
        let mut best_score = 0.0f64;
        for psm in psms {
            if peps.contains(&psm.peptide_sequence) {
                psm_count += 1;
                if psm.hyperscore > best_score {
                    best_score = psm.hyperscore;
                }
            }
        }

        // Coverage from first protein in group
        let coverage = if let Some(&pi) = prot_indices.first() {
            let pep_bytes: Vec<&[u8]> = peps.iter().map(|s| s.as_bytes()).collect();
            proteins[pi].coverage(&pep_bytes)
        } else {
            0.0
        };

        let is_decoy = prot_indices.iter().all(|&pi| proteins[pi].is_decoy);

        results.push(ProteinGroup {
            accessions,
            unique_peptides: unique,
            shared_peptides: shared,
            psm_count,
            best_score,
            coverage,
            is_decoy,
        });
    }

    // Sort by score descending
    results.sort_by(|a, b| b.best_score.partial_cmp(&a.best_score).unwrap_or(std::cmp::Ordering::Equal));

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_psm(seq: &str, score: f64) -> Psm {
        Psm {
            spectrum_id: format!("scan_{}", seq),
            peptide_sequence: seq.to_string(),
            xcorr: 1.0,
            hyperscore: score,
            matched_b: 5,
            matched_y: 5,
            total_b: 6,
            total_y: 6,
            delta_mass: 0.0,
            charge: 2,
            is_decoy: false,
        }
    }

    #[test]
    fn test_protein_contains_peptide() {
        let prot = ProteinEntry::new("P1", b"MAAAKPEPTIDEKSEQUENCER");
        assert!(prot.contains_peptide(b"PEPTIDEK"));
        assert!(prot.contains_peptide(b"SEQUENCER"));
        assert!(!prot.contains_peptide(b"MISSING"));
    }

    #[test]
    fn test_protein_coverage() {
        let prot = ProteinEntry::new("P1", b"AAAAABBBBBCCCCC"); // 15 residues
        let coverage = prot.coverage(&[b"AAAAA".as_slice(), b"CCCCC".as_slice()]);
        assert!((coverage - 10.0 / 15.0).abs() < 1e-10);
    }

    #[test]
    fn test_infer_single_protein() {
        let psms = vec![
            make_psm("PEPTIDEK", 20.0),
            make_psm("SEQUENCER", 15.0),
        ];
        let proteins = vec![
            ProteinEntry::new("P1", b"MAAAKPEPTIDEKSEQUENCER"),
        ];

        let groups = infer_proteins(&psms, &proteins).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].accessions, vec!["P1"]);
        assert_eq!(groups[0].psm_count, 2);
    }

    #[test]
    fn test_infer_parsimony() {
        // P1 has peptide A only; P2 has peptide C only
        // Parsimony: need both proteins to explain all peptides
        let psms = vec![
            make_psm("AAAAAAK", 20.0),
            make_psm("CCCCCK", 10.0),
        ];
        let proteins = vec![
            ProteinEntry::new("P1", b"MAAAAAAK"),   // contains AAAAAAK
            ProteinEntry::new("P2", b"MCCCCCKDD"),  // contains CCCCCK
        ];

        let groups = infer_proteins(&psms, &proteins).unwrap();
        // Both proteins needed — each has a unique peptide
        let total_accessions: usize = groups.iter().map(|g| g.accessions.len()).sum();
        assert_eq!(total_accessions, 2);
        // All peptides should be explained
        let all_peps: Vec<String> = groups.iter()
            .flat_map(|g| g.unique_peptides.iter().chain(g.shared_peptides.iter()))
            .cloned()
            .collect();
        assert!(all_peps.contains(&"AAAAAAK".to_string()));
        assert!(all_peps.contains(&"CCCCCK".to_string()));
    }

    #[test]
    fn test_infer_indistinguishable() {
        // Two proteins sharing identical peptides
        let psms = vec![make_psm("PEPTIDEK", 20.0)];
        let proteins = vec![
            ProteinEntry::new("P1", b"MPEPTIDEK"),
            ProteinEntry::new("P2", b"MPEPTIDEK"),
        ];

        let groups = infer_proteins(&psms, &proteins).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].accessions.len(), 2); // Both in same group
    }

    #[test]
    fn test_infer_empty() {
        let groups = infer_proteins(&[], &[]).unwrap();
        assert!(groups.is_empty());
    }

    #[test]
    fn test_shared_peptides() {
        let psms = vec![
            make_psm("UNIQUE", 20.0),
            make_psm("SHARED", 15.0),
        ];
        let proteins = vec![
            ProteinEntry::new("P1", b"MUNIQUESHARED"),
            ProteinEntry::new("P2", b"MSHARED"),
        ];

        let groups = infer_proteins(&psms, &proteins).unwrap();
        // P1 should have UNIQUE as unique and SHARED as shared
        let p1_group = groups.iter().find(|g| g.accessions.contains(&"P1".to_string())).unwrap();
        assert!(p1_group.unique_peptides.contains(&"UNIQUE".to_string()));
        assert!(p1_group.shared_peptides.contains(&"SHARED".to_string()));
    }
}
