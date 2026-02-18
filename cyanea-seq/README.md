# cyanea-seq

> Sequence types, I/O, and analysis for DNA, RNA, and protein sequences.

## What's Inside

- **Sequence types** -- validated `DnaSequence`, `RnaSequence`, `ProteinSequence` with reverse complement, transcription, translation
- **FASTA/FASTQ parsing** -- streaming statistics via needletail, paired-end support, interleave/deinterleave
- **K-mers** -- iterator-based extraction, 2-bit packed encoding (4x compression)
- **Indexing** -- suffix array (SA-IS), FM-index, FMD-index (bidirectional, SMEM enumeration)
- **Pattern matching** -- 7 algorithms: Horspool, KMP, Shift-And, BNDM, BOM, Myers bitparallel, Ukkonen
- **MinHash** -- bottom-k and scaled FracMinHash sketching for rapid genome comparison
- **Quality trimming** -- `TrimPipeline` builder with adapter removal, sliding window, BWA-style, paired-end
- **Motif scanning** -- PSSM with const-generic alphabet, PWM, EM-based motif discovery
- **ORF finding** -- all 6 reading frames, configurable start/stop codons
- **Codon tables** -- 7 NCBI translation tables, codon usage analysis, CAI
- **Sequence masking** -- DUST (low-complexity), SEG (protein), tandem repeat detection
- **RNA structure** -- Nussinov, Zuker MFE, McCaskill partition function, dot-bracket notation
- **Protein properties** -- composition, hydrophobicity, pI, extinction coefficient, Chou-Fasman/GOR secondary structure, disorder prediction
- **De Bruijn graphs** -- k-mer graph construction, unitig extraction
- **Assembly QC** -- N50/L50/N90/L90, auN statistics
- **Taxonomy** -- taxonomic trees, LCA queries, Kraken-style k-mer classification
- **Restriction enzymes** -- 20 common enzymes, cut-site finding, in-silico digestion

## Quick Start

```toml
[dependencies]
cyanea-seq = { version = "0.1", features = ["minhash"] }
```

```rust
use cyanea_seq::{DnaSequence, MinHash};

let seq = DnaSequence::new(b"ACGTACGTACGT").unwrap();
println!("GC: {:.1}%", seq.gc_content() * 100.0);
println!("RevComp: {:?}", seq.reverse_complement());

let sketch = MinHash::from_sequence(b"ACGTACGT", 4, 100).unwrap();
```

## Feature Flags

| Flag | Default | Description |
|------|---------|-------------|
| `std` | Yes | Standard library support |
| `wasm` | No | WASM target marker |
| `serde` | No | Serialize/Deserialize derives |
| `minhash` | No | MinHash/FracMinHash sketching |

## Modules

| Module | Description |
|--------|-------------|
| `alphabet` | `DnaAlphabet`, `RnaAlphabet`, `ProteinAlphabet` |
| `types` / `seq` | `DnaSequence`, `RnaSequence`, `ProteinSequence` |
| `codon` | Codon translation tables |
| `kmer` | K-mer iterator |
| `quality` | Phred quality scores |
| `fasta` / `fastq` | FASTA/FASTQ parsing and statistics |
| `paired` | Paired-end FASTQ support |
| `trim` | Quality trimming, adapter removal, `TrimPipeline` |
| `twobit` | 2-bit packed DNA encoding |
| `suffix` | Suffix array (SA-IS algorithm) |
| `fm_index` | FM-Index (BWT backward search) |
| `fmd_index` | Bidirectional FM-Index (SMEM enumeration) |
| `bwt` | Burrows-Wheeler Transform |
| `pattern` | 7 exact/approximate string matching algorithms |
| `pssm` | Position-Specific Scoring Matrix |
| `motif` | DNA motif PWM, scanning, EM discovery |
| `orf` | Open reading frame finder |
| `minhash` | MinHash/FracMinHash sketching (feature-gated) |
| `rna_structure` | RNA secondary structure prediction |
| `protein_properties` | Protein physicochemical analysis |
| `debruijn` | De Bruijn graph and unitig extraction |
| `assembly` | Assembly QC metrics |
| `taxonomy` | Taxonomic trees and k-mer classification |
| `restriction` | Restriction enzyme digestion |
| `fasta_index` | FASTA indexed reader (.fai) |

## See Also

- [API Reference (STATUS.md)](docs/STATUS.md)
- [Architecture](../ARCHITECTURE.md)
- [Build Guide](../docs/BUILDING.md)
