use cyanea_proteomics::peptide::*;
use cyanea_proteomics::search::*;
use cyanea_proteomics::spectrum::*;

fn make_spectrum_from_peptide(peptide: &Peptide, charge: i32) -> MassSpectrum {
    let ions = fragment_ions(peptide, 1);
    let peaks: Vec<Peak> = ions.iter()
        .filter(|i| i.neutral_loss.is_none() && (i.ion_type == IonType::B || i.ion_type == IonType::Y))
        .map(|i| Peak { mz: i.mz, intensity: 1000.0 })
        .collect();

    MassSpectrum::new("test", MsLevel::Ms2, 100.0, peaks)
        .unwrap()
        .with_precursor(Precursor {
            mz: peptide.mz(charge),
            charge: Some(charge),
            intensity: Some(5000.0),
            isolation_width: None,
            fragmentation: FragmentationMethod::HCD,
        })
}

#[test]
fn test_full_search_pipeline() {
    // Create protein, digest, search, score
    let protein = b"MAAAKPEPTIDEKSEQENCER";
    let config = DigestConfig {
        protease: Protease::Trypsin,
        max_missed_cleavages: 0,
        min_length: 3,
        max_length: 50,
        min_mass: 0.0,
        max_mass: 10000.0,
    };

    let peptides = digest(protein, &config).unwrap();
    assert!(peptides.len() >= 2);

    // Create spectra matching the peptides
    let spectra: Vec<MassSpectrum> = peptides.iter()
        .enumerate()
        .map(|(i, p)| {
            let mut s = make_spectrum_from_peptide(p, 2);
            s.id = format!("scan_{}", i);
            s
        })
        .collect();

    let search_config = SearchConfig {
        precursor_tolerance: 50.0,
        ..SearchConfig::default()
    };

    let psms = search_all(&spectra, &peptides, &search_config).unwrap();
    assert!(!psms.is_empty());

    // Each spectrum should match its source peptide
    for psm in &psms {
        assert!(psm.matched_b > 0 || psm.matched_y > 0);
        assert!(psm.hyperscore > 0.0);
    }
}

#[test]
fn test_decoy_generation() {
    let peptides = vec![
        Peptide::new(b"PEPTIDEK").unwrap(),
        Peptide::new(b"AACDEFGR").unwrap(),
        Peptide::new(b"SINGLE").unwrap(),
    ];

    let decoys = generate_decoys(&peptides);
    assert_eq!(decoys.len(), 3);

    // Last residue preserved
    assert_eq!(*decoys[0].sequence.last().unwrap(), b'K');
    assert_eq!(*decoys[1].sequence.last().unwrap(), b'R');
    // Single-AA peptide stays the same
    assert_eq!(*decoys[2].sequence.last().unwrap(), b'E');

    // Decoys should be different from targets (except very short ones)
    assert_ne!(decoys[0].sequence, peptides[0].sequence);
    assert_ne!(decoys[1].sequence, peptides[1].sequence);
}

#[test]
fn test_search_with_modifications() {
    let mut pep = Peptide::new(b"PEPTCIDE").unwrap();
    pep.add_modification(4, Modification::Carbamidomethyl).unwrap();

    let spec = make_spectrum_from_peptide(&pep, 2);
    let config = SearchConfig {
        precursor_tolerance: 50.0,
        ..SearchConfig::default()
    };

    let psm = score_peptide(&spec, &pep, 2, &config).unwrap();
    assert!(psm.is_some());
}

#[test]
fn test_search_multiple_candidates() {
    let target = Peptide::new(b"PEPTIDE").unwrap();
    let wrong1 = Peptide::new(b"ACDEFGH").unwrap();
    let wrong2 = Peptide::new(b"KLMNPQR").unwrap();

    let spec = make_spectrum_from_peptide(&target, 2);
    let database = vec![wrong1, target.clone(), wrong2];

    let config = SearchConfig {
        precursor_tolerance: 100.0,
        min_matched_ions: 2,
        ..SearchConfig::default()
    };

    let best = search_spectrum(&spec, &database, &config).unwrap();
    assert!(best.is_some());
    assert_eq!(best.unwrap().peptide_sequence, "PEPTIDE");
}
