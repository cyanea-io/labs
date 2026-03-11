# cyanea-proteomics API Reference

## Types

### `MassSpectrum`
Mass spectrum with peaks and metadata.

| Field | Type | Description |
|-------|------|-------------|
| `id` | `String` | Spectrum identifier |
| `ms_level` | `MsLevel` | MS1, MS2, MS3, or MSN |
| `retention_time` | `f64` | Retention time in seconds |
| `peaks` | `Vec<Peak>` | Peaks sorted by m/z |
| `precursor` | `Option<Precursor>` | Precursor info (MS2+) |
| `tic` | `f64` | Total ion current |
| `base_peak_intensity` | `f64` | Base peak intensity |
| `base_peak_mz` | `f64` | Base peak m/z |

**Methods:**

| Method | Signature | Description |
|--------|-----------|-------------|
| `new` | `(id, ms_level, rt, peaks) -> Result<Self>` | Create from peaks (auto-sorts, validates) |
| `with_precursor` | `(self, Precursor) -> Self` | Attach precursor info |
| `num_peaks` | `() -> usize` | Peak count |
| `find_peak` | `(target_mz, tolerance) -> Option<&Peak>` | Find closest peak within tolerance |
| `filter_by_relative_intensity` | `(min_relative)` | Remove peaks below threshold |
| `top_n` | `(n)` | Keep only N most intense peaks |
| `deisotope` | `(tolerance)` | Merge isotope clusters |
| `normalize` | `()` | Normalize to base peak = 100 |
| `sqrt_transform` | `()` | Sqrt-transform intensities |

### `Peak`
| Field | Type | Description |
|-------|------|-------------|
| `mz` | `f64` | Mass-to-charge ratio |
| `intensity` | `f64` | Signal intensity |

### `MsLevel`
`Ms1`, `Ms2`, `Ms3`, `MsN(u8)`

### `Precursor`
| Field | Type | Description |
|-------|------|-------------|
| `mz` | `f64` | Precursor m/z |
| `charge` | `Option<i32>` | Charge state |
| `intensity` | `Option<f64>` | Precursor intensity |
| `isolation_width` | `Option<f64>` | Isolation window |
| `fragmentation` | `FragmentationMethod` | CID/HCD/ETD/EThcD/Unknown |

### `Peptide`
| Field | Type | Description |
|-------|------|-------------|
| `sequence` | `Vec<u8>` | Amino acid sequence |
| `modifications` | `Vec<ModifiedResidue>` | Post-translational modifications |
| `missed_cleavages` | `usize` | Number of missed cleavages |

**Methods:**

| Method | Signature | Description |
|--------|-----------|-------------|
| `new` | `(sequence: &[u8]) -> Result<Self>` | Create unmodified peptide |
| `add_modification` | `(position, Modification) -> Result<()>` | Add PTM |
| `molecular_weight` | `() -> f64` | Monoisotopic neutral mass |
| `mz` | `(charge: i32) -> f64` | m/z for charge state |
| `len` | `() -> usize` | Sequence length |
| `sequence_str` | `() -> String` | Sequence as string |

### `Modification`
`Carbamidomethyl` (+57.02), `Oxidation` (+15.99), `Phospho` (+79.97), `Acetyl` (+42.01), `Deamidation` (+0.98), `TMT` (+229.16), `ITRAQ4` (+144.10), `Custom(f64)`

### `Protease`
`Trypsin`, `TrypsinP`, `LysC`, `Chymotrypsin`, `AspN`, `GluC`

### `DigestConfig`
| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `protease` | `Protease` | `Trypsin` | Enzyme |
| `max_missed_cleavages` | `usize` | `2` | Max missed cleavages |
| `min_length` | `usize` | `6` | Min peptide length |
| `max_length` | `usize` | `50` | Max peptide length |
| `min_mass` | `f64` | `400.0` | Min mass (Da) |
| `max_mass` | `f64` | `5000.0` | Max mass (Da) |

### `Psm` (Peptide-Spectrum Match)
| Field | Type | Description |
|-------|------|-------------|
| `spectrum_id` | `String` | Spectrum ID |
| `peptide_sequence` | `String` | Matched peptide |
| `xcorr` | `f64` | Cross-correlation score |
| `hyperscore` | `f64` | Hyperscore |
| `matched_b` / `matched_y` | `usize` | Matched b/y ions |
| `total_b` / `total_y` | `usize` | Total possible b/y ions |
| `delta_mass` | `f64` | Mass error (Da) |
| `charge` | `i32` | Charge state |
| `is_decoy` | `bool` | Decoy hit flag |

### `ProteinGroup`
| Field | Type | Description |
|-------|------|-------------|
| `accessions` | `Vec<String>` | Protein accession(s) |
| `unique_peptides` | `Vec<String>` | Unique peptides |
| `shared_peptides` | `Vec<String>` | Shared peptides |
| `psm_count` | `usize` | Total PSMs |
| `best_score` | `f64` | Best hyperscore |
| `coverage` | `f64` | Sequence coverage (0-1) |

### `ProteinQuant`
| Field | Type | Description |
|-------|------|-------------|
| `accessions` | `Vec<String>` | Protein accession(s) |
| `raw_value` | `f64` | Raw quantification |
| `normalized_value` | `f64` | Normalized value |
| `peptide_count` | `usize` | Quantified peptides |
| `psm_count` | `usize` | PSMs used |

### `TmtPlex`
`Tmt6`, `Tmt10`, `Tmt11`, `Tmt16`

### `FdrResult`
| Field | Type | Description |
|-------|------|-------------|
| `spectrum_id` | `String` | Spectrum ID |
| `peptide_sequence` | `String` | Peptide |
| `score` | `f64` | Score used for ranking |
| `q_value` | `f64` | Estimated q-value |
| `is_decoy` | `bool` | Decoy flag |
| `passes` | `bool` | Passes FDR threshold |

## Functions

### Spectrum Operations

| Function | Signature | Description |
|----------|-----------|-------------|
| `run_stats` | `(spectra: &[MassSpectrum]) -> Result<RunStats>` | Summary statistics for a run |

### Peptide & Digestion

| Function | Signature | Description |
|----------|-----------|-------------|
| `amino_acid_mass` | `(aa: u8) -> Option<f64>` | Monoisotopic mass for amino acid |
| `fragment_ions` | `(peptide, max_charge) -> Vec<FragmentIon>` | Generate theoretical b/y/a ions |
| `digest` | `(protein, config) -> Result<Vec<Peptide>>` | In-silico protein digestion |

### File Format I/O

| Function | Signature | Description |
|----------|-----------|-------------|
| `parse_mgf` | `(text: &str) -> Result<Vec<MassSpectrum>>` | Parse MGF format |
| `write_mgf` | `(spectra: &[MassSpectrum]) -> String` | Write MGF format |
| `mgf_stats` | `(spectra: &[MassSpectrum]) -> Result<MgfStats>` | MGF file statistics |
| `parse_mzml_text` | `(text: &str) -> Result<Vec<MassSpectrum>>` | Parse simplified mzML |
| `write_mzml_text` | `(spectra: &[MassSpectrum]) -> String` | Write simplified mzML |
| `mzml_stats` | `(spectra: &[MassSpectrum]) -> Result<MzmlStats>` | mzML file statistics |

### Database Search

| Function | Signature | Description |
|----------|-----------|-------------|
| `score_peptide` | `(spectrum, peptide, charge, config) -> Result<Option<Psm>>` | Score one PSM |
| `search_spectrum` | `(spectrum, database, config) -> Result<Option<Psm>>` | Best PSM for spectrum |
| `search_all` | `(spectra, database, config) -> Result<Vec<Psm>>` | Search all MS2 spectra |
| `generate_decoys` | `(targets: &[Peptide]) -> Vec<Peptide>` | Reversed-sequence decoys |

### Protein Inference

| Function | Signature | Description |
|----------|-----------|-------------|
| `infer_proteins` | `(psms, proteins) -> Result<Vec<ProteinGroup>>` | Parsimony inference |

### Quantification

| Function | Signature | Description |
|----------|-----------|-------------|
| `spectral_counting` | `(groups, psms) -> Vec<ProteinQuant>` | Spectral count quantification |
| `intensity_quantification` | `(groups, psms, spectra, method) -> Vec<ProteinQuant>` | Intensity-based quantification |
| `quantify_tmt` | `(spectra, plex, tolerance) -> Vec<TmtQuant>` | TMT reporter ion quantification |
| `median_normalize` | `(quants: &mut [ProteinQuant])` | Median normalization |
| `log2_transform` | `(quants: &mut [ProteinQuant])` | Log2 transform |

### FDR Estimation

| Function | Signature | Description |
|----------|-----------|-------------|
| `estimate_fdr` | `(targets, decoys, config) -> Result<Vec<FdrResult>>` | Target-decoy FDR |
| `filter_fdr` | `(targets, decoys, threshold) -> Result<Vec<Psm>>` | Filter at FDR level |
| `fdr_summary` | `(targets, decoys) -> Result<FdrSummary>` | FDR summary stats |

### mzTab Output

| Function | Signature | Description |
|----------|-----------|-------------|
| `write_mztab_proteins` | `(groups, quants) -> String` | Protein-level mzTab |
| `write_mztab_psms` | `(results: &[FdrResult]) -> String` | PSM-level mzTab |
| `parse_mztab_proteins` | `(text: &str) -> Result<Vec<(String, usize)>>` | Parse PRT rows |
