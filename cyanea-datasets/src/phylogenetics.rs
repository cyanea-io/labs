//! Phylogenetics demo datasets: primate tree, distance matrices.

/// Primate phylogeny in Newick format (13 species).
///
/// Based on the consensus primate phylogeny with approximate branch lengths
/// in substitutions per site.
pub fn primate_newick() -> &'static str {
    "((((((Human:0.006,Chimpanzee:0.008):0.003,Gorilla:0.009):0.008,Orangutan:0.015):0.012,(Gibbon:0.022,Siamang:0.021):0.010):0.030,((Rhesus:0.028,Baboon:0.030):0.015,(Vervet:0.025,Capuchin:0.055):0.020):0.012):0.050,(Mouse_Lemur:0.100,Ring_tailed_Lemur:0.098):0.050);"
}

/// Primate species names matching the Newick tree.
pub fn primate_species() -> Vec<&'static str> {
    vec![
        "Human", "Chimpanzee", "Gorilla", "Orangutan", "Gibbon", "Siamang",
        "Rhesus", "Baboon", "Vervet", "Capuchin",
        "Mouse_Lemur", "Ring_tailed_Lemur",
    ]
}

/// Mitochondrial cytochrome b protein sequences (partial, ~50 aa) for phylogenetics.
///
/// 6 primates — enough for demo tree building.
pub fn cytochrome_b_primates() -> Vec<(&'static str, &'static [u8])> {
    vec![
        ("Human",      b"MTNIRKSHPLFKIINHSFIDLPAPSNISAWWNFGSLLGICLIIQITTGLFLAMHYSPDAST"),
        ("Chimpanzee", b"MTNIRKSHPLFKIINHSFIDLPAPSNISAWWNFGSLLGICLIIQITTGLFLAMHYSPDAST"),
        ("Gorilla",    b"MTNIRKSHPLFKIINHSFIDLPAPSNISSWWNFGSLLGICLIIQIATGLFLAMHYSPDAST"),
        ("Orangutan",  b"MTNIRKSHPLFKIINHAFIDLPAPSNISSWWNFGSLLGICLIIQILTGLFLAMHYSPDAST"),
        ("Rhesus",     b"MTNIRKSHPLLKIINHSFIDLPAPSNISAWWNFGSLLGICLIVQILTGLFLAMHYSPDAST"),
        ("Mouse_Lemur",b"MTNIRKSHPLMKIINHSFIDLPAPSNISTWWNFGSLLGICLIAQIATGLFLAMHYTPDTST"),
    ]
}

/// Pairwise distance matrix for 6 primates (Jukes-Cantor corrected).
///
/// Returns (species_names, distance_matrix).
pub fn primate_distance_matrix() -> (Vec<&'static str>, Vec<Vec<f64>>) {
    let names = vec!["Human", "Chimpanzee", "Gorilla", "Orangutan", "Rhesus", "Mouse_Lemur"];
    let matrix = vec![
        vec![0.000, 0.012, 0.018, 0.032, 0.058, 0.165],
        vec![0.012, 0.000, 0.019, 0.033, 0.059, 0.166],
        vec![0.018, 0.019, 0.000, 0.034, 0.060, 0.168],
        vec![0.032, 0.033, 0.034, 0.000, 0.062, 0.172],
        vec![0.058, 0.059, 0.060, 0.062, 0.000, 0.175],
        vec![0.165, 0.166, 0.168, 0.172, 0.175, 0.000],
    ];
    (names, matrix)
}

/// Demo NEXUS format phylogeny block.
pub fn primate_nexus() -> String {
    format!(
        "#NEXUS\n\
         BEGIN TREES;\n\
         \tTree primate_consensus = {}\n\
         END;\n",
        primate_newick()
    )
}

/// Demo multiple sequence alignment in FASTA format (primates cytochrome b).
pub fn cytochrome_b_fasta() -> String {
    cytochrome_b_primates()
        .iter()
        .map(|(name, seq)| format!(">{}\n{}", name, std::str::from_utf8(seq).unwrap()))
        .collect::<Vec<_>>()
        .join("\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_primate_newick() {
        let tree = primate_newick();
        assert!(tree.ends_with(';'));
        // Should contain all species
        for species in primate_species() {
            assert!(tree.contains(species), "Missing: {}", species);
        }
    }

    #[test]
    fn test_primate_species() {
        let species = primate_species();
        assert_eq!(species.len(), 12);
    }

    #[test]
    fn test_cytochrome_b() {
        let seqs = cytochrome_b_primates();
        assert_eq!(seqs.len(), 6);
        // All same length (aligned)
        let len = seqs[0].1.len();
        assert!(seqs.iter().all(|s| s.1.len() == len));
    }

    #[test]
    fn test_distance_matrix() {
        let (names, matrix) = primate_distance_matrix();
        assert_eq!(names.len(), 6);
        assert_eq!(matrix.len(), 6);
        // Symmetric
        for i in 0..6 {
            for j in 0..6 {
                assert!((matrix[i][j] - matrix[j][i]).abs() < 1e-10);
            }
            assert!((matrix[i][i]).abs() < 1e-10); // zero diagonal
        }
        // Human-Chimp closest pair
        assert!(matrix[0][1] < matrix[0][4]);
    }

    #[test]
    fn test_nexus() {
        let nexus = primate_nexus();
        assert!(nexus.contains("#NEXUS"));
        assert!(nexus.contains("BEGIN TREES"));
    }

    #[test]
    fn test_cytochrome_b_fasta() {
        let fasta = cytochrome_b_fasta();
        let records: Vec<&str> = fasta.lines().filter(|l| l.starts_with('>')).collect();
        assert_eq!(records.len(), 6);
    }
}
