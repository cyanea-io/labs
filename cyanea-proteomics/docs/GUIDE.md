# cyanea-proteomics Usage Guide

## Overview

`cyanea-proteomics` provides a complete proteomics analysis pipeline: from parsing raw MS data through peptide identification, protein inference, quantification, and FDR control.

## 1. Working with Mass Spectra

### Parse an MGF file

```rust
use cyanea_proteomics::mgf::{parse_mgf, mgf_stats};

let mgf_text = r#"
BEGIN IONS
TITLE=spectrum_1
PEPMASS=523.25 1500.0
CHARGE=2+
RTINSECONDS=120.5
100.05 1000
200.10 500
300.15 750
END IONS
"#;

let spectra = parse_mgf(mgf_text).unwrap();
let stats = mgf_stats(&spectra).unwrap();
println!("Spectra: {}, Mean peaks: {:.1}", stats.num_spectra, stats.mean_peaks);
```

### Spectrum processing

```rust
use cyanea_proteomics::spectrum::{MassSpectrum, Peak, MsLevel};

let peaks = vec![
    Peak { mz: 100.0, intensity: 1000.0 },
    Peak { mz: 200.0, intensity: 50.0 },
    Peak { mz: 300.0, intensity: 500.0 },
];

let mut spec = MassSpectrum::new("scan1", MsLevel::Ms2, 120.0, peaks).unwrap();

// Remove low-intensity noise
spec.filter_by_relative_intensity(0.05);

// Keep top 100 peaks
spec.top_n(100);

// Normalize to base peak = 100
spec.normalize();

// Remove isotope peaks
spec.deisotope(1.1);
```

## 2. Peptide Operations

### In-silico digestion

```rust
use cyanea_proteomics::peptide::{digest, DigestConfig, Protease};

let protein = b"MAAAKPEPTIDEKSEQENCER";
let config = DigestConfig {
    protease: Protease::Trypsin,
    max_missed_cleavages: 2,
    min_length: 6,
    max_length: 50,
    min_mass: 400.0,
    max_mass: 5000.0,
};

let peptides = digest(protein, &config).unwrap();
for pep in &peptides {
    println!("{}: MW={:.2} Da", pep.sequence_str(), pep.molecular_weight());
}
```

### Modifications

```rust
use cyanea_proteomics::peptide::{Peptide, Modification};

let mut pep = Peptide::new(b"PEPTCIDE").unwrap();
pep.add_modification(4, Modification::Carbamidomethyl).unwrap();
println!("Modified MW: {:.2} Da", pep.molecular_weight());
println!("m/z at charge 2: {:.4}", pep.mz(2));
```

### Fragment ion prediction

```rust
use cyanea_proteomics::peptide::{Peptide, fragment_ions, IonType};

let pep = Peptide::new(b"PEPTIDE").unwrap();
let ions = fragment_ions(&pep, 2); // up to charge 2

for ion in &ions {
    if ion.neutral_loss.is_none() {
        println!("{:?}{} (charge {}): m/z = {:.4}",
            ion.ion_type, ion.number, ion.charge, ion.mz);
    }
}
```

## 3. Database Search

```rust
use cyanea_proteomics::peptide::{digest, DigestConfig};
use cyanea_proteomics::search::{search_all, generate_decoys, SearchConfig};
use cyanea_proteomics::mgf::parse_mgf;

// 1. Parse spectra
let spectra = parse_mgf(mgf_text).unwrap();

// 2. Build peptide database from proteins
let config = DigestConfig::default();
let peptides = digest(b"MAAAKPEPTIDEKSEQENCER", &config).unwrap();

// 3. Generate decoy database
let decoys = generate_decoys(&peptides);

// 4. Search
let search_config = SearchConfig {
    fragment_tolerance: 0.02,
    precursor_tolerance: 10.0,
    max_fragment_charge: 2,
    min_matched_ions: 4,
};

let target_psms = search_all(&spectra, &peptides, &search_config).unwrap();
let decoy_psms = search_all(&spectra, &decoys, &search_config).unwrap();
```

## 4. FDR Control

```rust
use cyanea_proteomics::fdr::{estimate_fdr, filter_fdr, fdr_summary, FdrConfig};

// Estimate FDR with q-values
let config = FdrConfig { threshold: 0.01, ..Default::default() };
let results = estimate_fdr(&target_psms, &decoy_psms, &config).unwrap();

for r in results.iter().filter(|r| r.passes) {
    println!("{}: score={:.1}, q={:.4}", r.peptide_sequence, r.score, r.q_value);
}

// Or simply filter at a threshold
let filtered = filter_fdr(&target_psms, &decoy_psms, 0.01).unwrap();
println!("PSMs at 1% FDR: {}", filtered.len());

// Summary statistics
let summary = fdr_summary(&target_psms, &decoy_psms).unwrap();
println!("Passing 1%: {}, 5%: {}", summary.passing_1pct, summary.passing_5pct);
```

## 5. Protein Inference

```rust
use cyanea_proteomics::protein::{infer_proteins, ProteinEntry};

let proteins = vec![
    ProteinEntry::new("P12345", b"MAAAKPEPTIDEKSEQENCER"),
    ProteinEntry::new("Q67890", b"MDDDDEKFFFFFER"),
];

let groups = infer_proteins(&filtered, &proteins).unwrap();
for g in &groups {
    println!("Proteins: {:?}, unique peptides: {}, PSMs: {}",
        g.accessions, g.unique_peptides.len(), g.psm_count);
}
```

## 6. Quantification

### Spectral counting

```rust
use cyanea_proteomics::quantification::{spectral_counting, median_normalize};

let mut quants = spectral_counting(&groups, &filtered);
median_normalize(&mut quants);

for q in &quants {
    println!("{}: raw={:.0}, normalized={:.3}",
        q.accessions.join(","), q.raw_value, q.normalized_value);
}
```

### TMT reporter ion quantification

```rust
use cyanea_proteomics::quantification::{quantify_tmt, TmtPlex};

let tmt_quants = quantify_tmt(&spectra, TmtPlex::Tmt10, 0.01);
for q in &tmt_quants {
    println!("Spectrum {}: channels={:?}, purity={:.2}",
        q.spectrum_id, q.channel_intensities, q.purity);
}
```

## 7. Export Results

```rust
use cyanea_proteomics::mztab::{write_mztab_proteins, write_mztab_psms};

// mzTab protein section
let mztab = write_mztab_proteins(&groups, Some(&quants));

// mzTab PSM section
let psm_mztab = write_mztab_psms(&fdr_results);
```

## 8. Complete Pipeline

```rust
use cyanea_proteomics::{
    mgf::parse_mgf,
    peptide::{digest, DigestConfig},
    search::{search_all, generate_decoys, SearchConfig},
    fdr::filter_fdr,
    protein::{infer_proteins, ProteinEntry},
    quantification::{spectral_counting, median_normalize},
    mztab::write_mztab_proteins,
};

// Parse → Digest → Search → FDR → Inference → Quantify → Export
let spectra = parse_mgf(mgf_text).unwrap();
let peptides = digest(protein_seq, &DigestConfig::default()).unwrap();
let decoys = generate_decoys(&peptides);
let targets = search_all(&spectra, &peptides, &SearchConfig::default()).unwrap();
let decoy_hits = search_all(&spectra, &decoys, &SearchConfig::default()).unwrap();
let filtered = filter_fdr(&targets, &decoy_hits, 0.01).unwrap();
let db = vec![ProteinEntry::new("P1", protein_seq)];
let groups = infer_proteins(&filtered, &db).unwrap();
let mut quants = spectral_counting(&groups, &filtered);
median_normalize(&mut quants);
let output = write_mztab_proteins(&groups, Some(&quants));
```
