# cyanea-seq Architecture

## Module Map

```
cyanea-seq
 +-- alphabet.rs         DnaAlphabet, RnaAlphabet, ProteinAlphabet trait impls
 +-- types.rs            ValidatedSeq<A> generic validated sequence type
 +-- seq.rs              Sequence-specific methods (rc, transcribe, translate, gc_content)
 +-- codon.rs            7 NCBI translation tables, codon usage, CAI
 +-- kmer.rs             KmerIter (Iterator + ExactSizeIterator + DoubleEndedIterator)
 +-- quality.rs          PhredEncoding, QualityScores
 +-- fasta.rs            FASTA parsing via needletail
 +-- fastq.rs            FASTQ parsing via needletail
 +-- paired.rs           [std] Paired-end FASTQ: types, parsing, writing, interleave
 +-- trim.rs             TrimPipeline builder, adapter removal, paired-end trimming
 +-- twobit.rs           TwoBitSequence (2-bit packed DNA encoding)
 +-- suffix.rs           SuffixArray (SA-IS O(n) construction)
 +-- bwt.rs              BWT construction and inversion
 +-- fm_index.rs         FmIndex (BWT + occurrence table + C table)
 +-- fmd_index.rs        FmdIndex (bidirectional FM-index, SMEM enumeration)
 +-- pattern.rs          7 exact/approximate string matching algorithms
 +-- pssm.rs             Pssm<const N> generic PSSM, PssmDna, PssmProtein
 +-- motif.rs            Pwm, motif scanning, EM-based discovery
 +-- motif_io.rs         MEME/TRANSFAC/JASPAR I/O, motif comparison
 +-- orf.rs              ORF finding (6 frames, configurable codons)
 +-- minhash.rs          [minhash] MinHash, FracMinHash sketching
 +-- rna_structure.rs    Nussinov, Zuker MFE, McCaskill partition function
 +-- protein_properties.rs  Composition, hydrophobicity, pI, extinction, 2ary structure, disorder
 +-- debruijn.rs         De Bruijn graph, unitig extraction
 +-- assembly.rs         Assembly QC: N50/L50/N90/L90, auN
 +-- taxonomy.rs         TaxonomyTree, LCA, KmerClassifier
 +-- restriction.rs      Restriction enzymes, cut sites, in-silico digestion
 +-- masking.rs          DUST (DNA), SEG (protein), tandem repeat masking
 +-- read_sim.rs         Illumina-style read simulator
 +-- fasta_index.rs      [std] FASTA .fai index, random access reader
 +-- longread.rs         Long-read types, stats, error correction, simulation, adapter trimming
 +-- sv.rs               Structural variant detection from split alignments and CIGAR
 +-- nanopore.rs         Nanopore signal metadata, methylation calling, run QC
```

## Design Decisions

### 2-bit DNA Encoding

`TwoBitSequence` encodes unambiguous DNA bases as 2 bits each (A=00, C=01, G=10, T=11), packing 4 bases per byte. This achieves 4x compression over ASCII and enables efficient k-mer extraction as integer encodings (up to k=32 in a u64).

Complement is computed via XOR with `0xFF` (flipping all bits swaps A<->T and C<->G). The encoding only supports ACGT -- ambiguous bases (N, R, Y, etc.) must be handled before encoding.

### FM-Index: BWT + Occurrence Table + Suffix Array Sample

The FM-Index stores three data structures built from the input text:

1. **BWT** (Burrows-Wheeler Transform): a permutation of the text that groups similar suffixes together.
2. **Occurrence table** (`Occ[c][i]`): for each character `c` and position `i`, the count of `c` in `BWT[0..i]`. Sampled every 128 positions with in-between counts computed by scanning.
3. **C table** (`C[c]`): the number of text characters lexicographically smaller than `c`.

Backward search finds the suffix array interval `[lo, hi)` for a pattern in O(m) time by iterating characters from right to left. Suffix array positions are stored in full (not sampled) for O(1) position lookup from the interval.

### FMD-Index: Bidirectional FM-Index

`FmdIndex` indexes the concatenation `text # revcomp(text) $` where `#` and `$` are sentinel characters. This enables both forward and backward extension of a bidirectional interval `(lower, size, lower_rev)`:

- `extend_backward(interval, c)` prepends character `c` using standard backward search.
- `extend_forward(interval, c)` appends character `c` by mapping to the reverse complement half.

**SMEM enumeration** iterates through a query sequence, extending each position forward until it cannot be extended further, collecting super-maximal exact matches. These are used as seeds in read mapping (like BWA-MEM).

### Pattern Matching: 7 Algorithms

| Algorithm | Time | Space | Pattern limit | Best for |
|-----------|------|-------|---------------|----------|
| Horspool | O(n/m) avg | O(sigma) | None | Long patterns, large alphabets |
| KMP | O(n+m) | O(m) | None | Guaranteed linear time |
| Shift-And | O(n) | O(sigma) | 64 | Short patterns, bitwise speed |
| BNDM | O(n/m) avg | O(sigma) | 64 | Short patterns, sub-linear |
| BOM | O(n/m) avg | O(m*sigma) | None | Medium patterns |
| Myers | O(n) | O(sigma) | 64 | Approximate matching, edit distance |
| Ukkonen | O(nk) | O(m) | None | Approximate matching, longer patterns |

All functions take `&[u8]` for generality -- they work on DNA, protein, or arbitrary byte strings.

### PSSM Scanning: Log-odds with Background Model

`Pssm<const N: usize>` uses const generics for the alphabet size (N=4 for DNA, N=20 for protein). Construction from a count matrix:

1. Add pseudocounts to each cell.
2. Convert to frequencies (row normalization).
3. Divide by background frequencies.
4. Take natural logarithm to get log-odds scores.

Scanning slides the PSSM across the sequence, summing log-odds at each position. Positions scoring above the threshold are reported as hits. Information content is computed in bits using the Shannon entropy formula.

### TrimPipeline: Builder Pattern for Read Processing

`TrimPipeline` uses the builder pattern to compose trimming and filtering steps. The processing order follows the Trimmomatic convention:

1. Adapter removal (3' overlap scan with mismatch tolerance)
2. Leading quality trim (5' end)
3. Trailing quality trim (3' end)
4. Sliding window / BWA quality trim
5. Length filter (min/max)
6. Quality filter (mean Q)
7. Complexity filter (Shannon entropy)

Each step produces a `TrimRange` (half-open interval of bases to keep). Ranges are intersected to determine the final trimmed region. If the result is empty or fails a filter, the read is discarded.

For paired-end data, `OrphanPolicy` controls what happens when only one mate survives: `DropBoth` (default), `KeepFirst` (keep R1 orphans), or `KeepSecond`.

### Long-Read Analysis: Platform-Aware Types and Error Correction

`LongRead` models reads from PacBio (CLR, HiFi) and Oxford Nanopore (simplex, duplex) with platform-specific metadata (number of passes for HiFi, predicted accuracy, read group). `LongReadPlatform` is a simple enum that drives default error rates when quality scores are absent:

| Platform | Default error rate |
|----------|-------------------|
| PacBio CLR | ~13% |
| PacBio HiFi | ~0.1% |
| Nanopore simplex | ~5% |
| Nanopore duplex | ~0.5% |

**Self-correction** (`self_correct`) uses internal k-mer frequency analysis: build a k-mer count table from the read itself, identify positions where all overlapping k-mers are below a low-frequency threshold (10% of median), then try single-base substitutions to find higher-frequency alternatives. This is a lightweight heuristic; for production accuracy, use multi-read correction or POA consensus.

**Consensus** (`simple_consensus`) takes pre-aligned overlapping reads and performs position-wise majority voting across A, C, G, T. No graph or alignment is performed -- reads must already be positionally aligned. For graph-based consensus, use `cyanea_align::poa`.

**Simulation** (`simulate_long_reads`) generates reads with platform-appropriate error profiles using a simple LCG PRNG. Read lengths are sampled from a clamped normal distribution (Box-Muller). Only substitution errors are introduced (no indels in the current implementation).

**Adapter trimming** (`trim_adapters`) checks both ends of each read against a set of known adapter sequences using Hamming distance with a configurable mismatch threshold. The `common_adapters()` function provides PacBio SMRTbell 3.0 and ONT ligation/rapid adapter sequences.

### Structural Variant Detection: Split Alignment Analysis

The `sv` module detects structural variants from split/supplementary read alignments, which are the primary signal for SVs in long-read data.

**SV calling** (`call_svs`) follows this pipeline:

1. Group split alignment segments by read name.
2. Filter by minimum mapping quality.
3. Sort segments by query position within each read.
4. Classify consecutive segment pairs:
   - **Different chromosomes** -> breakend (BND)
   - **Different strands** -> inversion (INV)
   - **Same strand, reference gap > query gap** -> deletion (DEL)
   - **Same strand, query gap > reference gap** -> insertion (INS)
   - **Overlapping reference positions** -> duplication (DUP)
5. Cluster nearby SVs of the same type using positional distance.
6. Filter by minimum support (number of reads).

**CIGAR-based calling** (`svs_from_cigar`) complements split-alignment calling by extracting large insertions and deletions directly from CIGAR operations (I/D/N ops >= min_length). This catches SVs that align within a single mapping segment.

`SvCallConfig` controls sensitivity: `min_sv_length` (default 50bp), `max_sv_length` (default 100Mb), `min_support` (default 3 reads), `cluster_distance` (default 500bp), and `min_mapq` (default 20).

### Nanopore Analysis: Signal Metadata and Methylation

The `nanopore` module handles Oxford Nanopore-specific data streams.

**Signal metadata** (`parse_signal_metadata`) parses key=value formatted lines as produced by tools like `pod5 inspect` or `slow5tools stats`. Captures read/run/flow-cell IDs, channel number, timing, sampling rate, signal statistics, and basecaller information.

**Methylation calling** operates on modified base calls from basecallers like Dorado:

1. `parse_methylation_calls` parses tab-separated per-read calls (read_id, chrom, position, strand, mod_type, probability). Supports both full names (`5mC`, `5hmC`, `6mA`, `4mC`) and shorthand codes (`m`, `h`, `a`).
2. `aggregate_methylation` groups calls by (chrom, position, strand, mod_type) and computes per-site coverage, modified count (probability >= threshold), mean probability, and methylation frequency. Results are sorted by genomic coordinate.

**Run QC** (`nanopore_qc`) computes Nanopore-specific quality metrics from signal metadata and basecall statistics: total reads/bases, N50, pass/fail counts, pore occupancy (fraction of channels producing reads), reads per channel, and estimated translocation speed (bases/second).
