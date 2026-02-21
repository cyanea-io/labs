# cyanea-wasm Usage Guide

## Installation

```bash
npm install @cyanea/bio
```

## Quick Start

```javascript
import init, { gc_content_json, align_dna, describe } from '@cyanea/bio';

await init();

// Sequence analysis
const gc = JSON.parse(gc_content_json("ATCGATCG"));
console.log(gc.ok); // 0.5

// Alignment
const result = JSON.parse(align_dna("ACGT", "ACGTT", "global"));
console.log(result.ok.score);

// Statistics
const stats = JSON.parse(describe("[1.0, 2.0, 3.0, 4.0, 5.0]"));
console.log(stats.ok.mean); // 3.0
```

## JSON Envelope

All functions return JSON strings with a consistent envelope:

```javascript
// Success
{ "ok": <value> }

// Error
{ "error": "<message>" }
```

Always parse the result and check for the `ok` or `error` key:

```javascript
const result = JSON.parse(some_function(args));
if (result.error) {
  console.error(result.error);
} else {
  console.log(result.ok);
}
```

## TypeScript Types

Import from `@cyanea/bio/types` for type-safe JSON parsing:

```typescript
import type { AlignmentResult, DescriptiveStats, MolecularProperties } from '@cyanea/bio/types';
```

## Module Examples

### Sequence Analysis

```javascript
import { parse_fasta, reverse_complement, translate, rna_fold_nussinov } from '@cyanea/bio';

const fasta = JSON.parse(parse_fasta(">seq1\nATCGATCG\n>seq2\nGGCCTTAA"));
const rc = reverse_complement("ATCG"); // "CGAT"
const protein = JSON.parse(translate("ATGAAAGGG"));
const fold = JSON.parse(rna_fold_nussinov("GCGCAAUAGCGC"));
```

### Alignment

```javascript
import { align_dna_custom, progressive_msa, parse_cigar } from '@cyanea/bio';

const aln = JSON.parse(align_dna_custom("ACGTACGT", "ACGTGCGT", "local", 2, -1, -5, -2));
const msa = JSON.parse(progressive_msa('["ACGT","ACGTT","ACGT"]', 1, -1, -2, -1));
const cigar = JSON.parse(parse_cigar("4M1I3M"));
```

### Statistics

```javascript
import { describe, t_test_two_sample, kaplan_meier, fst_hudson } from '@cyanea/bio';

const stats = JSON.parse(describe("[1,2,3,4,5]"));
const ttest = JSON.parse(t_test_two_sample("[1,2,3]", "[4,5,6]", false));
const km = JSON.parse(kaplan_meier("[1,2,3,5,7]", "[1,1,0,1,0]"));
const fst = JSON.parse(fst_hudson("[0,1,2,0,1]", "[2,2,1,2,1]"));
```

### Machine Learning

```javascript
import { pca, kmeans, random_forest_classify } from '@cyanea/bio';

const pca_result = JSON.parse(pca("[1,2,3,4,5,6,7,8]", 4, 2));
const clusters = JSON.parse(kmeans("[1,2,3,4,5,6,7,8]", 4, 2, 300));
```

### Chemistry

```javascript
import { smiles_properties, canonical, tanimoto, maccs_fingerprint } from '@cyanea/bio';

const props = JSON.parse(smiles_properties("c1ccccc1O"));
const canon = JSON.parse(canonical("OC1=CC=CC=C1"));
const sim = JSON.parse(tanimoto("c1ccccc1", "c1ccccc1O"));
const maccs = JSON.parse(maccs_fingerprint("c1ccccc1O"));
```

### Structural Biology

```javascript
import { pdb_info, pdb_secondary_structure, contact_map, ramachandran_analysis } from '@cyanea/bio';

const info = JSON.parse(pdb_info(pdb_text));
const ss = JSON.parse(pdb_secondary_structure(pdb_text));
const contacts = JSON.parse(contact_map(pdb_text, 8.0, "ca"));
const rama = JSON.parse(ramachandran_analysis(pdb_text));
```

### Phylogenetics

```javascript
import { newick_info, build_nj, simulate_coalescent } from '@cyanea/bio';

const tree = JSON.parse(newick_info("((A:0.1,B:0.2):0.3,C:0.4);"));
const nj = JSON.parse(build_nj('["A","B","C"]', '[[0,1,2],[1,0,1.5],[2,1.5,0]]'));
const coal = JSON.parse(simulate_coalescent(10, 1000, 42));
```

### I/O

```javascript
import { parse_vcf_text, parse_bed_text, parse_blast_xml } from '@cyanea/bio';

const vcf = JSON.parse(parse_vcf_text(vcf_string));
const bed = JSON.parse(parse_bed_text(bed_string));
const blast = JSON.parse(parse_blast_xml(xml_string));
```

### Omics

```javascript
import { merge_intervals, annotate_variant, morans_i } from '@cyanea/bio';

const merged = JSON.parse(merge_intervals('[{"chrom":"chr1","start":100,"end":200},{"chrom":"chr1","start":150,"end":300}]'));
const vep = JSON.parse(annotate_variant(variant_json, transcripts_json));
const moran = JSON.parse(morans_i(values_json, weights_json));
```

## Web Worker Usage

For heavy computation (alignment of long sequences, large PCA), run in a Web Worker to avoid blocking the main thread:

```javascript
// worker.js
import init, { align_dna_custom } from '@cyanea/bio';
await init();
self.onmessage = (e) => {
  const result = align_dna_custom(e.data.query, e.data.target, "local", 2, -1, -5, -2);
  self.postMessage(JSON.parse(result));
};
```

## Feature Flags

Build with the `wasm` feature to enable wasm-bindgen annotations:

```bash
wasm-pack build cyanea-wasm --features wasm
```
