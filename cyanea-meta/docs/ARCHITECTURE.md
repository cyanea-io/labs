# cyanea-meta Architecture

## Module Map

```
cyanea-meta
 +-- error.rs          MetaError enum (thiserror)
 +-- taxonomy.rs       TaxonomyDB, TaxonNode, TaxonRank, k-mer LCA classification
 +-- profile.rs        TaxonomicProfile, abundance estimation, filtering, merging
 +-- diversity.rs      Alpha diversity (Shannon, Simpson, Chao1, ACE, Fisher's alpha),
 |                     beta diversity (Bray-Curtis), rarefaction, rarefying
 +-- composition.rs    CLR/ILR transforms, ALDEx2-style differential abundance, ANCOM
 +-- functional.rs     FunctionalProfile, taxa-to-function mapping, pathway abundance
 +-- binning.rs        Contig, Bin, tetranucleotide frequency, k-means binning, quality
 +-- assembly.rs       AssemblyStats (N50/L50/N90/L90/auN), coverage, filtering
```

## Design Decisions

### K-mer Based Classification

`TaxonomyDB.classify_sequence()` implements a Kraken-style approach:

1. Extract all k-mers (k=21 by default) from the read
2. Look up each k-mer in a HashMap to find which taxon's reference it belongs to
3. Collect all taxon hits
4. Compute the Lowest Common Ancestor (LCA) of all hits

This trades sensitivity for speed -- exact k-mer matching is fast but requires indexed reference sequences. The LCA step ensures that reads matching multiple species are assigned to their common ancestor rather than being misclassified.

### Compositional Data Analysis

Microbiome abundance data is compositional (relative abundances sum to 1.0, so changes in one taxon necessarily affect others). Raw proportions violate independence assumptions of standard statistical tests.

The composition module provides two transforms:

- **CLR (Centered Log-Ratio)**: `clr(x_i) = log(x_i / geometric_mean(x))`. Maps compositions to real space while preserving relative relationships. Used in ALDEx2-style differential abundance.
- **ILR (Isometric Log-Ratio)**: Orthonormal basis transform that maps D-part compositions to D-1 dimensional real space. More statistically rigorous but harder to interpret.

### Diversity Metrics

Alpha diversity metrics each capture different aspects of community structure:

- **Shannon** — information-theoretic entropy; sensitive to rare species
- **Simpson** — probability two random individuals are different species; emphasizes dominant species
- **Chao1** — non-parametric richness estimator using singletons and doubletons
- **ACE** — abundance-based coverage estimator; uses species with <=10 individuals
- **Fisher's alpha** — assumes log-series species abundance distribution

Rarefaction uses the hypergeometric expectation formula for computational efficiency rather than random subsampling, giving exact expected values.

### Binning Strategy

Metagenomic binning groups contigs by two signals:

1. **Tetranucleotide frequency (TNF)**: Species have characteristic 4-mer usage patterns. The 256-dimensional TNF vector (normalized to frequencies) serves as a composition signature.
2. **Coverage**: Contigs from the same genome tend to have similar sequencing depth.

The current implementation uses k-means clustering on coverage values. A more sophisticated approach would combine TNF and coverage in a joint feature space (as MetaBAT2 does), but coverage alone works reasonably well for well-separated genomes.

## Dependencies

```
cyanea-core (errors, traits)
cyanea-seq  (sequence operations for k-mer extraction)
```

No external dependencies beyond the workspace. All statistical computations are implemented from scratch to avoid bloating the dependency tree.

## Testing Strategy

- **Unit tests**: Inline in each module, testing individual functions with small synthetic data
- **Integration tests**: In `tests/` directory, testing multi-function workflows
- **Numerical validation**: Diversity metrics verified against known analytical solutions (e.g., Shannon entropy of a uniform distribution)
