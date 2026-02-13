# cyanea-seq

Sequence types, validation, manipulation, and file parsing for DNA, RNA, and protein sequences.

## Status: Complete

All planned functionality is implemented and tested. FASTA/FASTQ parsing uses `needletail` for streaming performance.

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

FM-Index for O(m) exact pattern matching on DNA sequences via backward search. Built on the Burrows-Wheeler Transform (BWT), occurrence table, and C table. Supports the DNA alphabet (A, C, G, T) plus sentinel '$'.

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
| `myers_bitparallel(text: &[u8], pattern: &[u8], max_dist: usize) -> Vec<(usize, usize)>` | Myers bit-parallel approximate matching (pattern <= 64), returns (end_position, edit_distance) |
| `ukkonen(text: &[u8], pattern: &[u8], max_dist: usize) -> Vec<(usize, usize)>` | Ukkonen cut-off approximate matching, returns (end_position, edit_distance) |

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
| `from_counts(counts: &[Vec<f64>], pseudocount: f64, background: &[f64], mapping: fn(u8) -> Option<usize>) -> Result<Self>` | Build from count matrix with pseudocounts and background frequencies |
| `score(seq: &[u8], mapping: fn(u8) -> Option<usize>) -> Result<f64>` | Score a sequence window of length equal to the PSSM |
| `scan(seq: &[u8], threshold: f64, mapping: fn(u8) -> Option<usize>) -> Result<Vec<(usize, f64)>>` | Scan for motif hits above threshold, returns (position, score) |
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
| `find_orfs(seq: &[u8], min_length: usize) -> Vec<OrfResult>` | Find ORFs in all 3 forward reading frames |
| `find_orfs_both_strands(seq: &[u8], min_length: usize) -> Vec<OrfResult>` | Find ORFs in all 6 reading frames (forward + reverse complement) |
| `find_orfs_with_codons(seq: &[u8], min_length: usize, start_codons: &[&[u8]], stop_codons: &[&[u8]]) -> Vec<OrfResult>` | Find ORFs with configurable start and stop codons |

### FASTA indexed reader (`fasta_index.rs`)

Random access FASTA reading using `.fai` index files. Feature-gated behind `std`.

| Type | Description |
|------|-------------|
| `FastaIndex` | Parsed `.fai` index with entries keyed by sequence name |
| `FastaIndexEntry` | Entry for one sequence: `name`, `length`, `offset`, `line_bases`, `line_width` |
| `IndexedFastaReader` | Random access FASTA reader using an index |

**FastaIndex methods:**

| Method | Description |
|--------|-------------|
| `from_file(path: &Path) -> Result<Self>` | Parse an existing `.fai` index file |
| `build(fasta_path: &Path) -> Result<Self>` | Build an index from a FASTA file |
| `write(path: &Path) -> Result<()>` | Write the index to a `.fai` file |

**IndexedFastaReader methods:**

| Method | Description |
|--------|-------------|
| `open(fasta_path: &Path, fai_path: &Path) -> Result<Self>` | Open a FASTA file with its index |
| `fetch(name: &str, start: usize, end: usize) -> Result<Vec<u8>>` | Fetch a region by sequence name and coordinates |
| `fetch_all(name: &str) -> Result<Vec<u8>>` | Fetch the full sequence by name |

### FMD-Index (`fmd_index.rs`)

Bidirectional FM-Index for strand-aware search on DNA sequences. Indexes `text#reverse_complement(text)$` to support both forward and backward extension of search intervals, enabling super-maximal exact match (SMEM) enumeration.

| Type | Description |
|------|-------------|
| `FmdIndex` | Bidirectional FM-Index over `text#revcomp(text)$` |
| `BiInterval` | Bidirectional interval: `lower`, `size`, `lower_rev` |

**FmdIndex methods:**

| Method | Description |
|--------|-------------|
| `new(seq: &[u8]) -> Self` | Build from a DNA sequence |
| `init_interval(c: u8) -> BiInterval` | Initialize a bidirectional interval with a single character |
| `extend_backward(&self, interval: &BiInterval, c: u8) -> BiInterval` | Extend search backward (prepend character) |
| `extend_forward(&self, interval: &BiInterval, c: u8) -> BiInterval` | Extend search forward (append character) |
| `locate(&self, interval: &BiInterval) -> Vec<usize>` | Get text positions from a bidirectional interval |
| `backward_search(&self, pattern: &[u8]) -> BiInterval` | Full backward search for a pattern |
| `smems(&self, query: &[u8], min_len: usize) -> Vec<(usize, usize, BiInterval)>` | Enumerate super-maximal exact matches (query_start, query_end, interval) |

### MinHash sketching (`minhash.rs`)

Bottom-k MinHash and scaled FracMinHash sketching for rapid genome comparison. Estimates Jaccard similarity, containment, and average nucleotide identity (ANI) between DNA sequences without full alignment. Uses canonical k-mers (min of forward and reverse complement hash) for strand-agnostic sketching.

| Type | Description |
|------|-------------|
| `MinHash` | Bottom-k MinHash sketch: keeps the `sketch_size` smallest canonical k-mer hash values |
| `FracMinHash` | Scaled (fractional) MinHash sketch: keeps all hash values below `u64::MAX / scale` |

**MinHash methods:**

| Method | Description |
|--------|-------------|
| `new(k: usize, sketch_size: usize) -> Result<Self>` | Create an empty bottom-k sketch |
| `from_sequence(seq: &[u8], k: usize, sketch_size: usize) -> Result<Self>` | Build a sketch from a DNA sequence |
| `add_sequence(&mut self, seq: &[u8])` | Add k-mers from a DNA sequence to the sketch |
| `jaccard(&self, other: &MinHash) -> Result<f64>` | Estimate Jaccard similarity (merge-based estimator) |
| `containment(&self, other: &MinHash) -> Result<f64>` | Estimate containment C(A,B) = \|A intersect B\| / \|A\| |
| `ani(&self, other: &MinHash) -> Result<f64>` | Estimate ANI from Jaccard (Mash formula) |
| `len(&self) -> usize` | Number of hash values in the sketch |
| `is_empty(&self) -> bool` | Whether the sketch is empty |
| `k(&self) -> usize` | The k-mer size |
| `sketch_size(&self) -> usize` | The target sketch size (bottom-k parameter) |
| `hashes(&self) -> &[u64]` | Reference to the sorted hash values |

**FracMinHash methods:**

| Method | Description |
|--------|-------------|
| `new(k: usize, scale: u64) -> Result<Self>` | Create an empty scaled sketch |
| `from_sequence(seq: &[u8], k: usize, scale: u64) -> Result<Self>` | Build a scaled sketch from a DNA sequence |
| `add_sequence(&mut self, seq: &[u8])` | Add k-mers from a DNA sequence to the sketch |
| `jaccard(&self, other: &FracMinHash) -> Result<f64>` | Estimate Jaccard similarity |
| `containment(&self, other: &FracMinHash) -> Result<f64>` | Estimate containment C(A,B) = \|A intersect B\| / \|A\| |
| `ani(&self, other: &FracMinHash) -> Result<f64>` | Estimate ANI from Jaccard (Mash formula) |
| `len(&self) -> usize` | Number of hash values in the sketch |
| `is_empty(&self) -> bool` | Whether the sketch is empty |
| `k(&self) -> usize` | The k-mer size |
| `scale(&self) -> u64` | The scale factor |
| `hashes(&self) -> &[u64]` | Reference to the sorted hash values |

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

205 tests across 18 source files.

## Source Files

| File | Lines | Purpose |
|------|-------|---------|
| `lib.rs` | 81 | Module declarations, re-exports |
| `alphabet.rs` | 94 | Alphabet trait and implementations |
| `types.rs` | 283 | ValidatedSeq generic type |
| `seq.rs` | 197 | Sequence-specific methods |
| `codon.rs` | 159 | Codon translation table |
| `kmer.rs` | 113 | K-mer iterator |
| `quality.rs` | 164 | Phred quality scores |
| `fasta.rs` | 98 | FASTA parsing |
| `fastq.rs` | 256 | FASTQ parsing |
| `twobit.rs` | 373 | 2-bit packed DNA encoding |
| `suffix.rs` | 635 | Suffix array (SA-IS algorithm) |
| `fm_index.rs` | 355 | FM-Index (BWT backward search) |
| `minhash.rs` | 833 | MinHash and FracMinHash sketching |
| `pattern.rs` | 631 | Exact and approximate pattern matching (7 algorithms) |
| `pssm.rs` | 370 | Position-Specific Scoring Matrix / motif scanning |
| `orf.rs` | 288 | Open reading frame finder |
| `fasta_index.rs` | 421 | FASTA indexed reader (.fai) |
| `fmd_index.rs` | 732 | Bidirectional FM-Index (FMD-Index) |
