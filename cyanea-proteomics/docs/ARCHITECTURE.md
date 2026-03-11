# cyanea-proteomics Architecture

## Module Map

```
cyanea-proteomics
 +-- error.rs          ProteomicsError enum (thiserror)
 +-- spectrum.rs       MassSpectrum, Peak, MsLevel, Precursor, RunStats
 +-- peptide.rs        Peptide, amino acid masses, modifications, fragment ions, digestion
 +-- mgf.rs            MGF parsing and writing, MgfStats
 +-- mzml.rs           Simplified mzML parsing and writing, MzmlStats
 +-- search.rs         PSM scoring (XCorr, hyperscore), database search, decoy generation
 +-- protein.rs        Parsimony-based protein inference, protein groups
 +-- quantification.rs Label-free (spectral counting, intensity), TMT/iTRAQ quantification
 +-- fdr.rs            Target-decoy FDR estimation, q-value calculation
 +-- mztab.rs          mzTab 1.0 output format (proteins and PSMs)
```

## Design Decisions

### Spectrum Representation

Spectra store peaks sorted by m/z for efficient binary-search-style lookups. The `find_peak` method uses linear scan with early termination (since peaks are sorted). Processing operations (filter, top-N, normalize, deisotope) mutate in place and recalculate derived statistics.

### Fragment Ion Generation

The `fragment_ions` function generates all theoretical fragment ions in a single pass:
- **Primary ions**: b, y, a series at each cleavage site
- **Neutral losses**: H2O loss from b-ions, NH3 loss from y-ions
- **Charge states**: 1 through `max_charge`

Cumulative mass arrays are precomputed for O(n) generation. Modifications at specific positions are included in the cumulative sums.

### Database Search Scoring

Two complementary scores:

1. **Hyperscore** (X! Tandem-style): `log(matched_b! * matched_y! * intensity_sum)` — rewards both ion coverage and intensity
2. **XCorr** (simplified): Sum of normalized intensities at matched positions — correlates with Sequest XCorr

Only primary b/y ions are counted for scoring (neutral losses and a-ions are generated but excluded from match counting).

### Protein Inference

The parsimony algorithm follows the standard approach:

1. Map peptides to proteins via substring search
2. Group proteins with identical peptide sets (indistinguishable proteins)
3. Greedy set cover: iteratively select the group explaining the most unassigned peptides
4. Classify peptides as unique (one protein) or shared (multiple proteins)

### Target-Decoy FDR

Implements the standard target-decoy competition:

1. Generate decoys by reversing all but the last residue (preserves tryptic C-terminus)
2. Combine target and decoy PSMs, sort by score
3. FDR at each rank = cumulative decoys / cumulative targets
4. Convert to q-values: minimum FDR at or below each rank (monotonic from bottom)

### TMT Quantification

Reporter ion extraction finds the closest peak to each theoretical reporter m/z within a tolerance window. Isolation purity is estimated as the fraction of signal from reporter ions versus total signal in the reporter mass region (126-135 Da).

### mzTab Output

Follows the mzTab 1.0 specification with:
- MTD (metadata) section with version, mode, and type
- PRH/PRT (protein header/rows) with accession, PSM count, coverage, scores
- PSH/PSM (PSM header/rows) with sequence, score, q-value, decoy flag

## Dependencies

```
cyanea-core (errors, traits)
```

No external dependencies beyond `cyanea-core` and `thiserror`. All mass calculations, scoring functions, and statistical procedures are implemented from first principles.

## Testing Strategy

- **Unit tests**: 69 tests inline in each module covering all public functions
- **Integration tests**: 17 tests in `tests/` covering full pipelines (parse→search→FDR→quantify)
- **Total**: 86 tests
- Amino acid masses validated against NIST monoisotopic values
- Digestion validated with known protein sequences
- FDR q-values validated for monotonicity
- MGF/mzML roundtrip parsing validated
