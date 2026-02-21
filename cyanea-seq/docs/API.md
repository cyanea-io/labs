# cyanea-seq API Reference

Sequence types, validation, manipulation, and file parsing for DNA, RNA, and protein sequences.

## Public API

### Alphabet system (`alphabet.rs`)

| Type | Description |
|------|-------------|
| `trait Alphabet` | Defines valid bytes for a sequence type (`NAME`, `VALID_BYTES`, `is_valid`) |
| `DnaAlphabet` | IUPAC DNA: `ACGTNRYSWKMBDHV` |
| `RnaAlphabet` | IUPAC RNA: `ACGUNRYSWKMBDHV` |
| `ProteinAlphabet` | 20 standard amino acids + `XBZJUO*` |

### Sequence types (`types.rs`, `seq.rs`)

All sequence types are `ValidatedSeq<A: Alphabet>` -- validated, uppercase byte sequences.

| Type | Description |
|------|-------------|
| `DnaSequence` | Validated DNA sequence |
| `RnaSequence` | Validated RNA sequence |
| `ProteinSequence` | Validated protein sequence |

**DnaSequence methods:**

| Method | Description |
|--------|-------------|
| `new(data: &[u8]) -> Result<Self>` | Validate and construct |
| `reverse_complement() -> DnaSequence` | Reverse complement |
| `transcribe() -> RnaSequence` | DNA to RNA (T -> U) |
| `translate() -> Result<ProteinSequence>` | Transcribe then translate (NCBI Table 1) |
| `gc_content() -> f64` | GC fraction in [0.0, 1.0] |
| `kmers(k) -> Result<KmerIter>` | Extract k-mers |

**RnaSequence methods:**

| Method | Description |
|--------|-------------|
| `reverse_complement() -> RnaSequence` | Reverse complement |
| `reverse_transcribe() -> DnaSequence` | RNA to DNA (U -> T) |
| `translate() -> Result<ProteinSequence>` | Translate (NCBI Table 1) |
| `translate_frames() -> [Result<ProteinSequence>; 3]` | All three reading frames |
| `kmers(k) -> Result<KmerIter>` | Extract k-mers |

**ProteinSequence methods:**

| Method | Description |
|--------|-------------|
| `molecular_weight() -> f64` | Molecular weight in Daltons |
| `kmers(k) -> Result<KmerIter>` | Extract k-mers |

### Codon translation (`codon.rs`)

| Function | Description |
|----------|-------------|
| `translate_codon(codon: &[u8]) -> Option<u8>` | Single codon to amino acid (NCBI standard code) |
| `translate_sequence(seq: &[u8]) -> Vec<u8>` | Translate full coding sequence |
| `CodonTable` | Supports 7 NCBI translation tables (1, 2, 3, 4, 5, 11, 12) |
| `codon_usage(seq: &[u8]) -> CodonUsage` | Codon usage frequency table |
| `cai(seq: &[u8], reference: &CodonUsage) -> f64` | Codon Adaptation Index |

### K-mers (`kmer.rs`)

| Type | Description |
|------|-------------|
| `KmerIter<'a>` | Iterator over k-mers; implements `Iterator`, `ExactSizeIterator`, `DoubleEndedIterator` |

### Quality scores (`quality.rs`)

| Type/Function | Description |
|---------------|-------------|
| `PhredEncoding` | Enum: `Phred33`, `Phred64` |
| `QualityScores` | Decoded Phred scores with statistics |
| `QualityScores::mean()` | Mean quality score |
| `QualityScores::fraction_above(threshold)` | Fraction of bases above a quality threshold |
| `QualityScores::error_probability(phred)` | Convert Phred score to error probability |

### FASTA parsing (`fasta.rs`)

| Type/Function | Description |
|---------------|-------------|
| `FastaStats` | `sequence_count`, `total_bases`, `gc_content`, `avg_length` |
| `parse_fasta_stats(path) -> Result<FastaStats>` | Streaming FASTA statistics |

### FASTQ parsing (`fastq.rs`)

| Type/Function | Description |
|---------------|-------------|
| `FastqRecord` | Single record: `id`, `seq`, `quality` |
| `FastqStats` | `sequence_count`, `total_bases`, `gc_content`, `avg_length`, `mean_quality`, `q20_fraction`, `q30_fraction` |
| `parse_fastq_file(path) -> Result<Vec<FastqRecord>>` | Parse all records |
| `parse_fastq_stats(path) -> Result<FastqStats>` | Streaming statistics |

### 2-bit encoding (`twobit.rs`)

Compact 2-bit-per-base representation for unambiguous DNA (A, C, G, T only). Achieves 4x compression over ASCII. Encoding: A=00, C=01, G=10, T=11.

| Type | Description |
|------|-------------|
| `TwoBitSequence` | DNA sequence stored in 2-bit packed representation (4 bases per byte) |

**TwoBitSequence methods:**

| Method | Description |
|--------|-------------|
| `encode(seq: &[u8]) -> Result<Self>` | Encode ASCII DNA into 2-bit packed form (case-insensitive, ACGT only) |
| `decode(&self) -> Vec<u8>` | Decode back to uppercase ASCII DNA bytes |
| `get(&self, index: usize) -> Option<u8>` | Get the base at a specific position |
| `len(&self) -> usize` | Number of bases in the sequence |
| `is_empty(&self) -> bool` | Whether the sequence is empty |
| `kmer(&self, pos: usize, k: usize) -> Option<u64>` | Extract a k-mer as a 2-bit integer encoding (max k=32) |
| `complement(&self) -> Self` | Bitwise complement (A<->T, C<->G) via XOR |

### Suffix array (`suffix.rs`)

O(n) suffix array construction via the SA-IS algorithm (Nong, Zhang & Chan, 2009) with O(m log n) binary search for pattern matching.

| Type | Description |
|------|-------------|
| `SuffixArray` | Suffix array built from a text, storing sorted suffix positions |

**SuffixArray methods:**

| Method | Description |
|--------|-------------|
| `build(text: &[u8]) -> Self` | Build a suffix array using SA-IS (appends sentinel, O(n) construction) |
| `search(&self, text: &[u8], pattern: &[u8]) -> Vec<usize>` | Find all occurrences of pattern in text (O(m log n), sorted positions) |
| `len(&self) -> usize` | Number of entries in the suffix array (text length + 1 for sentinel) |

### FM-Index (`fm_index.rs`)

FM-Index for O(m) exact pattern matching on DNA sequences via backward search. Built on the Burrows-Wheeler Transform (BWT), occurrence table, and C table.

| Type | Description |
|------|-------------|
| `FmIndex` | FM-Index for DNA sequences with BWT, suffix array, occurrence table, and C table |

**FmIndex methods:**

| Method | Description |
|--------|-------------|
| `build(text: &[u8]) -> Self` | Build an FM-Index from a DNA sequence (ACGT only, sentinel appended internally) |
| `search(&self, pattern: &[u8]) -> Vec<usize>` | Find all occurrence positions of pattern (O(m) search + O(k) lookup, sorted) |
| `count(&self, pattern: &[u8]) -> usize` | Count occurrences without returning positions (more efficient than `search`) |

### Pattern matching (`pattern.rs`)

7 exact and approximate string matching algorithms operating on `&[u8]` for generality. Bitparallel methods (Shift-And, BNDM, Myers) limited to patterns of length 64 or less.

**Exact matching:**

| Function | Description |
|----------|-------------|
| `horspool(text: &[u8], pattern: &[u8]) -> Vec<usize>` | Boyer-Moore-Horspool exact matching, O(n/m) average case |
| `kmp(text: &[u8], pattern: &[u8]) -> Vec<usize>` | Knuth-Morris-Pratt exact matching, O(n+m) |
| `shift_and(text: &[u8], pattern: &[u8]) -> Vec<usize>` | Shift-And bitparallel exact matching (pattern <= 64) |
| `bndm(text: &[u8], pattern: &[u8]) -> Vec<usize>` | BNDM bitparallel exact matching (pattern <= 64) |
| `bom(text: &[u8], pattern: &[u8]) -> Vec<usize>` | Backward Oracle Matching exact matching |

**Approximate matching:**

| Function | Description |
|----------|-------------|
| `myers_bitparallel(text, pattern, max_dist) -> Vec<(usize, usize)>` | Myers bit-parallel approximate matching (pattern <= 64), returns (end_position, edit_distance) |
| `ukkonen(text, pattern, max_dist) -> Vec<(usize, usize)>` | Ukkonen cut-off approximate matching, returns (end_position, edit_distance) |

### PSSM / Motif scanning (`pssm.rs`)

Position-Specific Scoring Matrix with const generic alphabet size for DNA and protein motif scanning.

| Type | Description |
|------|-------------|
| `Pssm<const N: usize>` | Position-Specific Scoring Matrix with const generic alphabet size |
| `PssmDna` | Type alias for `Pssm<4>` (DNA alphabet) |
| `PssmProtein` | Type alias for `Pssm<20>` (protein alphabet) |

**Pssm methods:**

| Method | Description |
|--------|-------------|
| `from_counts(counts, pseudocount, background, mapping) -> Result<Self>` | Build from count matrix with pseudocounts and background frequencies |
| `score(seq, mapping) -> Result<f64>` | Score a sequence window of length equal to the PSSM |
| `scan(seq, threshold, mapping) -> Result<Vec<(usize, f64)>>` | Scan for motif hits above threshold, returns (position, score) |
| `information_content() -> Vec<f64>` | Per-position information content in bits |

**Helper functions:**

| Function | Description |
|----------|-------------|
| `dna_mapping(b: u8) -> Option<usize>` | Map DNA base to column index (A=0, C=1, G=2, T=3) |
| `protein_mapping(b: u8) -> Option<usize>` | Map amino acid to column index (20 standard residues) |

### ORF finder (`orf.rs`)

Open reading frame detection across forward and reverse complement strands with configurable start/stop codons.

| Type | Description |
|------|-------------|
| `OrfResult` | ORF result: `start`, `end`, `frame`, `strand`, `sequence` |
| `Strand` | Enum: `Forward`, `Reverse` |

| Function | Description |
|----------|-------------|
| `find_orfs(seq, min_length) -> Vec<OrfResult>` | Find ORFs in all 3 forward reading frames |
| `find_orfs_both_strands(seq, min_length) -> Vec<OrfResult>` | Find ORFs in all 6 reading frames |
| `find_orfs_with_codons(seq, min_length, start_codons, stop_codons) -> Vec<OrfResult>` | Find ORFs with configurable start and stop codons |

### FASTA indexed reader (`fasta_index.rs`, `std` feature)

Random access FASTA reading using `.fai` index files.

| Type | Description |
|------|-------------|
| `FastaIndex` | Parsed `.fai` index with entries keyed by sequence name |
| `FastaIndexEntry` | Entry for one sequence: `name`, `length`, `offset`, `line_bases`, `line_width` |
| `IndexedFastaReader` | Random access FASTA reader using an index |

**FastaIndex methods:**

| Method | Description |
|--------|-------------|
| `from_file(path) -> Result<Self>` | Parse an existing `.fai` index file |
| `build(fasta_path) -> Result<Self>` | Build an index from a FASTA file |
| `write(path) -> Result<()>` | Write the index to a `.fai` file |

**IndexedFastaReader methods:**

| Method | Description |
|--------|-------------|
| `open(fasta_path, fai_path) -> Result<Self>` | Open a FASTA file with its index |
| `fetch(name, start, end) -> Result<Vec<u8>>` | Fetch a region by sequence name and coordinates |
| `fetch_all(name) -> Result<Vec<u8>>` | Fetch the full sequence by name |

### FMD-Index (`fmd_index.rs`)

Bidirectional FM-Index for strand-aware search on DNA sequences. Indexes `text#reverse_complement(text)$` to support both forward and backward extension of search intervals, enabling super-maximal exact match (SMEM) enumeration.

| Type | Description |
|------|-------------|
| `FmdIndex` | Bidirectional FM-Index over `text#revcomp(text)$` |
| `BiInterval` | Bidirectional interval: `lower`, `size`, `lower_rev` |

**FmdIndex methods:**

| Method | Description |
|--------|-------------|
| `new(seq) -> Self` | Build from a DNA sequence |
| `init_interval(c) -> BiInterval` | Initialize a bidirectional interval with a single character |
| `extend_backward(interval, c) -> BiInterval` | Extend search backward (prepend character) |
| `extend_forward(interval, c) -> BiInterval` | Extend search forward (append character) |
| `locate(interval) -> Vec<usize>` | Get text positions from a bidirectional interval |
| `backward_search(pattern) -> BiInterval` | Full backward search for a pattern |
| `smems(query, min_len) -> Vec<(usize, usize, BiInterval)>` | Enumerate super-maximal exact matches (query_start, query_end, interval) |

### BWT utilities (`bwt.rs`)

Standalone Burrows-Wheeler Transform construction and inversion.

| Type | Description |
|------|-------------|
| `Bwt` | Burrows-Wheeler Transform with primary index |

**Bwt methods:**

| Method | Description |
|--------|-------------|
| `build(text) -> Self` | Build BWT from text (SA-via-sort, sentinel appended internally) |
| `as_bytes() -> &[u8]` | The BWT string |
| `primary_index() -> usize` | Position of sentinel in the BWT |
| `len() -> usize` | Length of BWT (text length + 1) |
| `invert() -> Vec<u8>` | Reconstruct original text via LF-mapping |

### Paired-end FASTQ (`paired.rs`, `std` feature)

Paired-end read support for Illumina short-read workflows.

| Type | Description |
|------|-------------|
| `PairedFastqRecord` | Paired-end record containing R1 (forward) and R2 (reverse) `FastqRecord` |
| `MateValidation` | Enum: `Strict`, `Relaxed`, `None` |
| `PairedFastqStats` | Per-file `FastqStats` for R1 and R2, plus `pair_count` |

**Parsing:**

| Function | Description |
|----------|-------------|
| `parse_paired_fastq_files(r1_path, r2_path, validation) -> Result<Vec<PairedFastqRecord>>` | Parse separate R1/R2 files |
| `parse_interleaved_fastq(path, validation) -> Result<Vec<PairedFastqRecord>>` | Parse interleaved file |
| `parse_paired_fastq_stats(r1_path, r2_path) -> Result<PairedFastqStats>` | Streaming statistics |

**Writing:**

| Function | Description |
|----------|-------------|
| `write_paired_fastq(pairs, r1_path, r2_path, encoding) -> Result<()>` | Write to separate R1/R2 files |
| `write_interleaved_fastq(pairs, path, encoding) -> Result<()>` | Write to interleaved file |
| `interleave_fastq_files(r1_path, r2_path, output_path, validation) -> Result<u64>` | Interleave two files |
| `deinterleave_fastq_file(input_path, r1_path, r2_path, validation) -> Result<u64>` | Split interleaved file |

### Quality trimming & filtering (`trim.rs`)

Quality trimming, adapter removal, and read filtering with `TrimPipeline` builder.

| Type | Description |
|------|-------------|
| `TrimRange` | Half-open range `[start, end)` describing which bases to keep |
| `TrimPipeline` | Configurable read-processing pipeline builder |
| `TrimReport` | Batch processing statistics |

**Low-level functions** (return `TrimRange`):

| Function | Description |
|----------|-------------|
| `trim_sliding_window(quality, window_size, threshold) -> TrimRange` | Trimmomatic SLIDINGWINDOW |
| `trim_leading(quality, threshold) -> TrimRange` | Remove 5' bases below threshold |
| `trim_trailing(quality, threshold) -> TrimRange` | Remove 3' bases below threshold |
| `trim_quality_3prime(quality, threshold) -> TrimRange` | BWA-style 3' trim |
| `find_adapter_3prime(seq, adapter, max_mismatches) -> usize` | Adapter detection |
| `shannon_entropy(seq) -> f64` | Base composition entropy |

**TrimPipeline builder methods:**

| Method | Description |
|--------|-------------|
| `new() -> Self` | Create empty pipeline |
| `adapter(self, adapter) -> Self` | Add adapter sequence |
| `illumina_adapters(self) -> Self` | Add all Illumina adapters |
| `leading(self, threshold) -> Self` | Trim 5' low-quality bases |
| `trailing(self, threshold) -> Self` | Trim 3' low-quality bases |
| `sliding_window(self, window_size, threshold) -> Self` | Sliding window quality trim |
| `bwa_quality(self, threshold) -> Self` | BWA-style 3' trim |
| `min_length(self, len) -> Self` | Minimum read length filter |
| `max_length(self, len) -> Self` | Maximum read length filter |
| `min_mean_quality(self, quality) -> Self` | Minimum mean quality filter |
| `min_entropy(self, entropy) -> Self` | Low-complexity filter |
| `process(record) -> Option<FastqRecord>` | Process single record |
| `process_batch(records) -> Vec<FastqRecord>` | Process batch |
| `process_batch_with_stats(records) -> TrimReport` | Process batch with stats |

**Paired-end trimming** (`std` feature):

| Type | Description |
|------|-------------|
| `OrphanPolicy` | Enum: `DropBoth`, `KeepFirst`, `KeepSecond` |
| `PairedTrimResult` | Enum: `BothPassed`, `OnlyFirst`, `OnlySecond`, `Dropped` |
| `PairedTrimReport` | Batch stats with survival rate and orphan counts |

### MinHash sketching (`minhash.rs`, `minhash` feature)

Bottom-k MinHash and scaled FracMinHash sketching for rapid genome comparison.

| Type | Description |
|------|-------------|
| `MinHash` | Bottom-k MinHash sketch |
| `FracMinHash` | Scaled (fractional) MinHash sketch |

**Shared methods** (both types):

| Method | Description |
|--------|-------------|
| `from_sequence(seq, k, size_or_scale) -> Result<Self>` | Build from DNA sequence |
| `add_sequence(seq)` | Add k-mers to sketch |
| `jaccard(other) -> Result<f64>` | Estimate Jaccard similarity |
| `containment(other) -> Result<f64>` | Estimate containment |
| `ani(other) -> Result<f64>` | Estimate average nucleotide identity |

### RNA secondary structure prediction (`rna_structure.rs`)

| Type | Description |
|------|-------------|
| `RnaSecondaryStructure` | Pair table: `pairs[i] = Some(j)` if paired, `None` if unpaired |
| `NussinovResult` | Nussinov prediction: `structure` + `max_pairs` |
| `MfeResult` | Zuker MFE prediction: `structure` + `energy` (kcal/mol) |
| `PartitionResult` | McCaskill partition function: `pair_probabilities` + `ensemble_energy` |

| Function | Description |
|----------|-------------|
| `nussinov(seq, min_loop_size) -> Result<NussinovResult>` | Maximize base pair count (O(n^3) DP) |
| `zuker_mfe(seq) -> Result<MfeResult>` | Minimum free energy with Turner parameters |
| `mccaskill(seq, temperature) -> Result<PartitionResult>` | Base pair probabilities via inside-outside |
| `base_pair_distance(a, b) -> Result<usize>` | Symmetric difference of base pair sets |
| `mountain_distance(a, b) -> Result<f64>` | L1 distance of mountain representations |

### Protein sequence properties (`protein_properties.rs`)

| Type | Description |
|------|-------------|
| `HydrophobicityScale` | Enum: `KyteDoolittle`, `HoppWoods` |
| `AminoAcidComposition` | Residue counts, fractions, and length |
| `ExtinctionCoefficient` | Molar absorptivity at 280 nm |
| `SecondaryStructurePrediction` | Per-residue states, scores, and fraction summaries |
| `DisorderPrediction` | Per-residue disorder scores and calls |

| Function | Description |
|----------|-------------|
| `amino_acid_composition(seq) -> Result<AminoAcidComposition>` | Residue counts and fractions |
| `hydrophobicity_profile(seq, window, scale) -> Result<Vec<f64>>` | Sliding window hydrophobicity |
| `gravy(seq) -> Result<f64>` | Grand average of hydropathicity |
| `isoelectric_point(seq) -> Result<f64>` | pI via bisection |
| `extinction_coefficient(seq) -> Result<ExtinctionCoefficient>` | Molar extinction at 280 nm |
| `chou_fasman(seq) -> Result<SecondaryStructurePrediction>` | Chou-Fasman prediction |
| `gor(seq) -> Result<SecondaryStructurePrediction>` | GOR prediction |
| `predict_disorder(seq, window) -> Result<DisorderPrediction>` | Disorder propensity |

### De Bruijn graph (`debruijn.rs`)

| Type | Description |
|------|-------------|
| `DeBruijnGraph` | Node-centric De Bruijn graph built from k-mers |
| `Unitig` | A maximal non-branching path with coverage |

| Method | Description |
|--------|-------------|
| `from_sequences(sequences, k) -> Result<Self>` | Build from input sequences |
| `from_kmers(kmers) -> Result<Self>` | Build from pre-extracted k-mers |
| `node_count() -> usize` | Number of (k-1)-mer nodes |
| `edge_count() -> usize` | Number of k-mer edges |
| `unitigs() -> Vec<Unitig>` | Extract all maximal non-branching paths |

### Assembly QC metrics (`assembly.rs`)

| Type/Function | Description |
|---------------|-------------|
| `AssemblyStats` | N50, L50, N90, L90, GC content, auN, contig counts |
| `assembly_stats(contigs) -> Result<AssemblyStats>` | Compute all assembly statistics |
| `nx_values(contigs, x) -> Result<(usize, usize)>` | Compute Nx and Lx for arbitrary x |

### Taxonomy (`taxonomy.rs`)

| Type | Description |
|------|-------------|
| `TaxonRank` | Enum: Domain, Phylum, Class, Order, Family, Genus, Species, Unranked |
| `TaxonomyNode` | Node with id, name, rank, parent link |
| `TaxonomyTree` | Rooted taxonomy tree with LCA queries |
| `KmerClassifier` | Kraken-style k-mer taxonomic classifier |

### Restriction enzymes (`restriction.rs`)

| Type | Description |
|------|-------------|
| `RestrictionEnzyme` | Enzyme with recognition site and cut offsets |
| `CutSite` | A located cut site with overhang type |
| `Fragment` | A restriction fragment with start, end, length |

| Function | Description |
|----------|-------------|
| `common_enzymes() -> Vec<RestrictionEnzyme>` | 20 common enzymes (EcoRI, BamHI, etc.) |
| `find_cut_sites(seq, enzyme) -> Vec<CutSite>` | Find all cut sites |
| `digest(seq, enzymes) -> Vec<Fragment>` | Digest with one or more enzymes |

### Motif discovery (`motif.rs`)

| Type | Description |
|------|-------------|
| `Pwm` | Position Weight Matrix for DNA |
| `MotifMatch` | A motif match: position, score, strand |
| `DiscoveredMotif` | Discovered motif with PWM, sites, and score |

| Function | Description |
|----------|-------------|
| `discover_motifs(sequences, motif_length, n_motifs, max_iter) -> Result<Vec<DiscoveredMotif>>` | MEME-style EM motif discovery |

### Motif format I/O (`motif_io.rs`)

| Type | Description |
|------|-------------|
| `Motif` | Sequence motif with name, alphabet, PWM matrix, and metadata |
| `MotifAlphabet` | Enum: `Dna`, `Rna`, `Protein` |

| Function | Description |
|----------|-------------|
| `parse_meme(input) -> Result<Vec<Motif>>` | Parse MEME minimal format |
| `write_meme(motifs) -> String` | Write MEME minimal format |
| `parse_transfac(input) -> Result<Vec<Motif>>` | Parse TRANSFAC format |
| `write_transfac(motifs) -> String` | Write TRANSFAC format |
| `parse_jaspar(input) -> Result<Vec<Motif>>` | Parse JASPAR format |
| `write_jaspar(motifs) -> String` | Write JASPAR format |
| `motif_similarity(a, b) -> f64` | Pearson correlation-based PWM similarity |

### Sequence masking (`masking.rs`)

| Type | Description |
|------|-------------|
| `MaskMode` | Enum: `Soft` (lowercase), `Hard` (N/X replacement) |
| `MaskSource` | Enum: `Dust`, `Seg`, `TandemRepeat` |
| `MaskedRegion` | Region with start, end, score, source |
| `MaskResult` | Masked sequence, regions, masked fraction |

| Function | Description |
|----------|-------------|
| `dust_mask(seq, threshold, mode) -> Result<MaskResult>` | DUST low-complexity DNA masking |
| `seg_mask(seq, window, threshold, mode) -> Result<MaskResult>` | SEG low-complexity protein masking |
| `tandem_repeat_mask(seq, min_period, min_copies, mode) -> Result<MaskResult>` | Tandem repeat detection |

### Read simulation (`read_sim.rs`)

| Type | Description |
|------|-------------|
| `ReadSimConfig` | Configuration: read_length, coverage, fragment_mean/std, error_rate, paired, seed |
| `SimulatedRead` | Simulated read with name, sequence, quality, and true origin |

| Function | Description |
|----------|-------------|
| `simulate_reads(reference, config) -> Result<Vec<SimulatedRead>>` | Generate synthetic Illumina-style reads |

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support |
| `wasm` | No | WASM target |
| `serde` | No | Serialization support |
| `minhash` | No | MinHash/FracMinHash sketching |

## Dependencies

- `cyanea-core` -- traits, error types
- `needletail` -- high-performance FASTA/FASTQ parsing

## Tests

474 unit + 26 doc tests across 31 source files.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 140 | Module declarations, re-exports |
| `alphabet.rs` | 94 | Alphabet trait and implementations |
| `types.rs` | 283 | ValidatedSeq generic type |
| `seq.rs` | 197 | Sequence-specific methods |
| `codon.rs` | 159 | Codon translation tables, usage, CAI |
| `kmer.rs` | 113 | K-mer iterator |
| `quality.rs` | 164 | Phred quality scores |
| `fasta.rs` | 98 | FASTA parsing |
| `fastq.rs` | 256 | FASTQ parsing |
| `twobit.rs` | 373 | 2-bit packed DNA encoding |
| `suffix.rs` | 635 | Suffix array (SA-IS algorithm) |
| `fm_index.rs` | 355 | FM-Index (BWT backward search) |
| `fmd_index.rs` | 732 | Bidirectional FM-Index (FMD-Index) |
| `bwt.rs` | 225 | Standalone BWT construction and inversion |
| `pattern.rs` | 631 | Exact and approximate pattern matching (7 algorithms) |
| `pssm.rs` | 370 | Position-Specific Scoring Matrix / motif scanning |
| `motif.rs` | ~400 | DNA motif PWM, scanning, EM discovery |
| `motif_io.rs` | ~350 | MEME/TRANSFAC/JASPAR I/O, motif comparison |
| `orf.rs` | 288 | Open reading frame finder |
| `fasta_index.rs` | 421 | FASTA indexed reader (.fai) |
| `paired.rs` | 1084 | Paired-end FASTQ: types, parsing, writing, interleave/deinterleave |
| `trim.rs` | 1535 | Quality trimming, adapter removal, filtering, TrimPipeline |
| `minhash.rs` | 833 | MinHash and FracMinHash sketching |
| `rna_structure.rs` | 1110 | RNA secondary structure: Nussinov, Zuker MFE, McCaskill |
| `protein_properties.rs` | 1223 | Protein properties: composition, hydrophobicity, pI, etc. |
| `debruijn.rs` | ~300 | De Bruijn graph and unitig extraction |
| `assembly.rs` | ~200 | Assembly QC metrics |
| `taxonomy.rs` | ~400 | Taxonomic trees and k-mer classification |
| `restriction.rs` | ~350 | Restriction enzyme digestion |
| `masking.rs` | ~400 | DUST/SEG/tandem repeat masking |
| `read_sim.rs` | ~300 | Illumina-style read simulator |
