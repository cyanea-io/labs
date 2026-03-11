# cyanea-epi Architecture

## Module Map

```
cyanea-epi
 +-- error.rs          EpiError enum (thiserror)
 +-- peaks.rs          Peak, PeakCallParams, PeakSet, call_peaks, call_broad_peaks
 +-- pileup.rs         TagPileup, build_pileup, normalize, smooth, correlation, fingerprint
 +-- motifs.rs         Motif, MotifMatch, PWM scanning, k-mer discovery, MEME I/O
 +-- chromatin.rs      ChromHMM-like EM state learning, Viterbi segmentation, enrichment
 +-- differential.rs   DiffResult, DESeq2-style differential binding, count reads in peaks
 +-- nucleosome.rs     NucleosomePosition, calling, NFR score, periodicity detection
 +-- accessibility.rs  ATAC-seq QC: TSS enrichment, FRiP, fragment metrics, full QC report
```

## Design Decisions

### Peak Calling (MACS2-style)

The peak caller implements the core MACS2 algorithm:

1. **Build pileup**: Extend each read to estimated fragment size, accumulate per-position coverage
2. **Estimate background**: Local lambda from control at three scales (1k, 5k, 10k windows), take the maximum
3. **Poisson test**: At each position, test enrichment over background using a Poisson p-value
4. **Merge regions**: Consecutive significant positions are merged, with a configurable max gap
5. **Find summits**: Within each merged region, find the position with highest signal
6. **BH correction**: Benjamini-Hochberg procedure on p-values to compute q-values
7. **Filter**: Keep peaks passing q-value cutoff

**Broad peak calling** uses a two-pass approach: first call narrow "seed" regions at a stringent threshold, then extend seeds into broader regions using a relaxed threshold. This captures diffuse histone marks like H3K27me3 and H3K36me3.

### Motif Discovery

Rather than implementing a full MEME-style EM algorithm (which is computationally expensive), the motif discovery uses a k-mer enrichment approach:

1. Count all k-mers of the target width in peak sequences
2. Compare to expected frequencies based on background nucleotide composition
3. Rank k-mers by enrichment score
4. Extend and refine top k-mers into PWMs by aligning surrounding sequence contexts

This is faster and suitable for finding strong motifs (like transcription factor binding sites at ChIP-seq peak summits). For weaker motifs or de novo discovery of complex patterns, a full EM-based approach would be needed.

### Chromatin State Learning

The chromatin state model is a multivariate HMM where:

- **States** represent functional chromatin types (active promoter, enhancer, repressed, etc.)
- **Emissions** are multivariate Bernoulli (independent presence/absence of each histone mark)
- **Transitions** capture the spatial structure of the genome (e.g., promoter states tend to be followed by transcribed states)

Training uses Baum-Welch EM:
1. Initialize emission probabilities randomly
2. E-step: Forward-backward algorithm to compute state responsibilities
3. M-step: Update emission and transition probabilities from responsibilities
4. Iterate until convergence (or max iterations)

Viterbi decoding then assigns the most likely state to each genomic bin.

### Differential Binding

Follows the DESeq2 approach adapted for count data at peaks:

1. **Size factors**: Median-of-ratios normalization to account for library size differences
2. **Test**: Welch's t-test on log-normalized counts (for two-condition comparisons)
3. **Correction**: Benjamini-Hochberg FDR control

This is a simplified version of the full DESeq2 negative binomial model, but works well in practice for replicated ChIP-seq experiments. For single-replicate comparisons, the module falls back to fold-change-based ranking.

### ATAC-seq QC Metrics

The `atacqc()` function computes ENCODE-standard quality metrics:

- **TSS enrichment**: Signal at TSSs vs flanking regions. Values >7 indicate good enrichment.
- **FRiP** (Fraction of Reads in Peaks): >0.3 is good for ATAC-seq.
- **NFR ratio**: Fraction of fragments in nucleosome-free range (<150 bp) vs mono-nucleosome range (150-300 bp). Higher values indicate better chromatin accessibility signal.
- **Periodicity**: Autocorrelation of fragment sizes showing ~200 bp periodicity from nucleosome arrays.

## Dependencies

```
cyanea-core (errors, traits)
```

No external dependencies. All statistical tests (Poisson, BH correction, Pearson correlation, Fisher's exact) are implemented from first principles.

## Testing Strategy

- **Unit tests**: 38 tests inline in each module
- **Integration tests**: 34 tests in `tests/` covering peak calling with synthetic data, motif round-trips, chromatin learning, and differential analysis
- **Doc tests**: 1 example in lib.rs
- Peak calling validated against synthetic enrichment signals with known properties
- Motif scanning validated with known consensus sequences
- Statistical tests validated against analytical solutions
