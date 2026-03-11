# cyanea-datasets API Reference

## Modules

### `genomics`
Genomics demo datasets: E. coli genome subset, SARS-CoV-2 spike, demo VCF.

### `alignment`
Alignment demo datasets: pre-aligned sequences, reference pairs.

### `epigenomics`
Epigenomics demo datasets: ChIP-seq peaks, methylation sites.

### `single_cell`
Single-cell demo dataset: small PBMC-like expression matrix.

### `chemistry`
Chemistry demo datasets: small molecule library.

### `phylogenetics`
Phylogenetics demo datasets: primate tree, distance matrices.

### `metagenomics`
Metagenomics demo datasets: OTU tables, taxonomy profiles.

### `structural`
Structural biology demo datasets: protein structures, contact maps.

### `protocols`
Structured protocol templates: 10 wet-lab and 6 dry-lab workflows with steps, timing, tips, and markdown rendering.

### `protocols::Protocol`
A protocol template with structured steps and metadata.

| Field | Type | Description |
|-------|------|-------------|
| `title` | `&'static str` | Protocol title |
| `slug` | `&'static str` | URL-friendly identifier |
| `category` | `ProtocolCategory` | Wet lab or dry lab |
| `description` | `&'static str` | Overview of the protocol |
| `estimated_time` | `&'static str` | Estimated total time |
| `difficulty` | `Difficulty` | Beginner, Intermediate, or Advanced |
| `requirements` | `Vec<&'static str>` | Materials or software needed |
| `steps` | `Vec<ProtocolStep>` | Ordered protocol steps |
| `expected_outputs` | `Vec<&'static str>` | What the protocol produces |
| `references` | `Vec<&'static str>` | Literature references |

**Methods:**

| Method | Signature | Description |
|--------|-----------|-------------|
| `to_markdown` | `(&self) -> String` | Render the full protocol as a markdown document |

### `protocols::ProtocolCategory`
Protocol category enum.

| Variant | Description |
|---------|-------------|
| `WetLab` | Bench/laboratory protocol |
| `DryLab` | Computational/bioinformatics pipeline |

### `protocols::Difficulty`
Protocol difficulty level enum.

| Variant | Description |
|---------|-------------|
| `Beginner` | Suitable for new researchers |
| `Intermediate` | Requires some domain experience |
| `Advanced` | Requires significant expertise |

### `protocols::ProtocolStep`
A single step in a protocol.

| Field | Type | Description |
|-------|------|-------------|
| `number` | `u32` | Step number (1-based) |
| `title` | `&'static str` | Step title |
| `description` | `&'static str` | Detailed instructions |
| `duration` | `Option<&'static str>` | Time estimate for this step |
| `tips` | `Vec<&'static str>` | Practical tips |
| `caution` | `Option<&'static str>` | Safety or critical warnings |

## Types

### `genomics::DemoVariant`
A simple variant record.

| Field | Type | Description |
|-------|------|-------------|
| `chrom` | `&'static str` | Chromosome |
| `pos` | `u64` | 1-based position |
| `id` | `&'static str` | Variant ID (e.g., rsID) |
| `ref_a` | `&'static str` | Reference allele |
| `alt` | `&'static str` | Alternate allele |
| `qual` | `f64` | Quality score |
| `af` | `f64` | Allele frequency |

**Methods:**

| Method | Signature | Description |
|--------|-----------|-------------|
| `to_vcf_line` | `(&self) -> String` | Format as a VCF data line |

### `genomics::DemoGene`
A simplified gene annotation.

| Field | Type | Description |
|-------|------|-------------|
| `name` | `&'static str` | Gene symbol |
| `chrom` | `&'static str` | Chromosome |
| `start` | `u64` | Start position |
| `end` | `u64` | End position |
| `strand` | `char` | Strand ('+' or '-') |
| `exon_count` | `usize` | Number of exons |
| `transcript` | `&'static str` | Transcript accession |
| `description` | `&'static str` | Gene description |

### `alignment::DemoAlignment`
A pre-computed pairwise alignment result.

| Field | Type | Description |
|-------|------|-------------|
| `seq_a_name` | `&'static str` | Name of first sequence |
| `seq_b_name` | `&'static str` | Name of second sequence |
| `aligned_a` | `&'static [u8]` | Aligned first sequence (with gaps) |
| `aligned_b` | `&'static [u8]` | Aligned second sequence (with gaps) |
| `score` | `i64` | Alignment score |
| `identity` | `f64` | Fraction identity (0.0-1.0) |
| `gaps` | `usize` | Number of gap positions |

### `epigenomics::DemoPeak`
A ChIP-seq narrow peak.

| Field | Type | Description |
|-------|------|-------------|
| `chrom` | `&'static str` | Chromosome |
| `start` | `u64` | Peak start |
| `end` | `u64` | Peak end |
| `name` | `&'static str` | Peak name |
| `score` | `f64` | Peak score |
| `signal` | `f64` | Signal value |
| `pvalue` | `f64` | P-value |
| `qvalue` | `f64` | Q-value (FDR) |
| `summit` | `u64` | Summit position |

**Methods:**

| Method | Signature | Description |
|--------|-----------|-------------|
| `to_narrowpeak` | `(&self) -> String` | Format as narrowPeak (BED6+4) line |
| `to_bed` | `(&self) -> String` | Format as BED line |

### `epigenomics::DemoCpgSite`
A CpG methylation site.

| Field | Type | Description |
|-------|------|-------------|
| `chrom` | `&'static str` | Chromosome |
| `pos` | `u64` | Position |
| `methylated` | `u32` | Methylated read count |
| `total` | `u32` | Total read count |
| `beta` | `f64` | Beta value (methylated / total) |

### `single_cell::DemoSingleCell`
A demo single-cell dataset.

| Field | Type | Description |
|-------|------|-------------|
| `genes` | `Vec<String>` | Gene names (rows) |
| `cells` | `Vec<String>` | Cell barcodes (columns) |
| `matrix` | `Vec<Vec<f64>>` | Count matrix (genes x cells) |
| `cell_types` | `Vec<String>` | Ground-truth cell type labels |

**Methods:**

| Method | Signature | Description |
|--------|-----------|-------------|
| `num_genes` | `(&self) -> usize` | Number of genes |
| `num_cells` | `(&self) -> usize` | Number of cells |
| `gene_expression` | `(&self, gene_name: &str) -> Option<&Vec<f64>>` | Get expression vector for a gene by name |
| `total_counts_per_cell` | `(&self) -> Vec<f64>` | Total UMI counts per cell |
| `genes_per_cell` | `(&self) -> Vec<usize>` | Number of detected genes (count > 0) per cell |

### `chemistry::DemoMolecule`
A demo molecule with key properties.

| Field | Type | Description |
|-------|------|-------------|
| `name` | `&'static str` | Drug name |
| `smiles` | `&'static str` | SMILES string |
| `molecular_weight` | `f64` | Molecular weight (g/mol) |
| `logp` | `f64` | Partition coefficient |
| `hbd` | `usize` | Hydrogen bond donors |
| `hba` | `usize` | Hydrogen bond acceptors |
| `therapeutic_area` | `&'static str` | Therapeutic classification |

**Methods:**

| Method | Signature | Description |
|--------|-----------|-------------|
| `lipinski_compliant` | `(&self) -> bool` | Check Lipinski's Rule of Five compliance |
| `lipinski_violations` | `(&self) -> usize` | Count number of Lipinski violations (0-4) |

### `metagenomics::DemoOtuTable`
A demo OTU table.

| Field | Type | Description |
|-------|------|-------------|
| `taxa` | `Vec<&'static str>` | Taxon names (rows) |
| `samples` | `Vec<&'static str>` | Sample names (columns) |
| `counts` | `Vec<Vec<u64>>` | Counts matrix (taxa x samples) |
| `metadata` | `Vec<(&'static str, &'static str)>` | Sample metadata: (sample_name, group) |

**Methods:**

| Method | Signature | Description |
|--------|-----------|-------------|
| `sample_totals` | `(&self) -> Vec<u64>` | Total counts per sample |
| `relative_abundance` | `(&self) -> Vec<Vec<f64>>` | Relative abundance matrix (columns sum to 1.0) |
| `shannon_diversity` | `(&self) -> Vec<f64>` | Shannon diversity index per sample |

### `metagenomics::DemoTaxonomy`
Taxonomy lineage for a species.

| Field | Type | Description |
|-------|------|-------------|
| `species` | `&'static str` | Species name |
| `genus` | `&'static str` | Genus |
| `family` | `&'static str` | Family |
| `order` | `&'static str` | Order |
| `class` | `&'static str` | Class |
| `phylum` | `&'static str` | Phylum |

### `structural::DemoStructure`
A demo protein structure (Cα trace).

| Field | Type | Description |
|-------|------|-------------|
| `name` | `&'static str` | Structure name |
| `pdb_id` | `&'static str` | PDB identifier |
| `chain` | `char` | Chain ID |
| `residues` | `Vec<(&'static str, u32)>` | Residues as (name, number) pairs |
| `ca_coords` | `Vec<(f64, f64, f64)>` | Cα coordinates (x, y, z) in angstroms |

**Methods:**

| Method | Signature | Description |
|--------|-----------|-------------|
| `len` | `(&self) -> usize` | Number of residues |
| `is_empty` | `(&self) -> bool` | Check if structure has no residues |
| `ca_distance` | `(&self, i: usize, j: usize) -> f64` | Euclidean distance between two Cα atoms |
| `contact_map` | `(&self, threshold: f64) -> Vec<Vec<bool>>` | Compute Cα contact map (distance < threshold) |
| `to_pdb` | `(&self) -> String` | Generate simplified PDB format output |

## Functions

### Genomics

| Function | Signature | Description |
|----------|-----------|-------------|
| `genomics::ecoli_16s_rrna` | `() -> (&'static str, &'static [u8])` | E. coli K-12 MG1655 rrsA (16S rRNA) gene, 1,542 bp |
| `genomics::sars_cov2_spike_rbd` | `() -> (&'static str, &'static [u8])` | SARS-CoV-2 Spike RBD nucleotide sequence, ~670 bp |
| `genomics::sars_cov2_spike_rbd_protein` | `() -> (&'static str, &'static [u8])` | SARS-CoV-2 Spike RBD protein sequence, ~222 aa |
| `genomics::demo_variants_chr22` | `() -> Vec<DemoVariant>` | 10 demo variants on human chr22 (SNVs, indels) |
| `genomics::demo_vcf_chr22` | `() -> String` | Full VCF v4.2 string from demo variants |
| `genomics::tp53_gene` | `() -> DemoGene` | Demo gene annotation for TP53 |
| `genomics::brca1_gene` | `() -> DemoGene` | Demo gene annotation for BRCA1 |

### Alignment

| Function | Signature | Description |
|----------|-----------|-------------|
| `alignment::spike_alignment_seqs` | `() -> Vec<(&'static str, &'static [u8])>` | 4 SARS-CoV-2 spike RBD sequences (Wuhan, Alpha, Delta, Omicron), ~60 aa |
| `alignment::hemoglobin_alpha` | `() -> Vec<(&'static str, &'static [u8])>` | Hemoglobin alpha from 5 species (Human, Chimpanzee, Dog, Chicken, Zebrafish), ~141 aa |
| `alignment::demo_pairwise_alignment` | `() -> DemoAlignment` | Pre-computed Human vs. Chicken hemoglobin alignment |
| `alignment::demo_cigar_strings` | `() -> Vec<(&'static str, &'static str)>` | 6 demo CIGAR strings for SAM/BAM examples |

### Epigenomics

| Function | Signature | Description |
|----------|-----------|-------------|
| `epigenomics::demo_chipseq_peaks` | `() -> Vec<DemoPeak>` | 10 H3K27ac ChIP-seq narrow peaks on chr1 |
| `epigenomics::demo_narrowpeak` | `() -> String` | narrowPeak format output from demo peaks |
| `epigenomics::demo_cpg_methylation` | `() -> Vec<DemoCpgSite>` | 12 CpG methylation sites in a promoter region |

### Single-Cell

| Function | Signature | Description |
|----------|-----------|-------------|
| `single_cell::demo_pbmc_50` | `() -> DemoSingleCell` | 50-cell PBMC dataset: 20 genes, 3 cell types (T cell, B cell, Monocyte) |

### Chemistry

| Function | Signature | Description |
|----------|-----------|-------------|
| `chemistry::fda_approved_drugs` | `() -> Vec<DemoMolecule>` | 12 FDA-approved drugs with SMILES and properties |
| `chemistry::demo_sdf_summary` | `() -> String` | SDF-like text summary of the molecule library |

### Phylogenetics

| Function | Signature | Description |
|----------|-----------|-------------|
| `phylogenetics::primate_newick` | `() -> &'static str` | 12-species primate phylogeny in Newick format |
| `phylogenetics::primate_species` | `() -> Vec<&'static str>` | Species names matching the Newick tree |
| `phylogenetics::cytochrome_b_primates` | `() -> Vec<(&'static str, &'static [u8])>` | Cytochrome b protein sequences for 6 primates, ~60 aa |
| `phylogenetics::primate_distance_matrix` | `() -> (Vec<&'static str>, Vec<Vec<f64>>)` | Jukes-Cantor distance matrix for 6 primates |
| `phylogenetics::primate_nexus` | `() -> String` | Primate phylogeny in NEXUS format |
| `phylogenetics::cytochrome_b_fasta` | `() -> String` | Primate cytochrome b in FASTA format |

### Metagenomics

| Function | Signature | Description |
|----------|-----------|-------------|
| `metagenomics::demo_otu_table` | `() -> DemoOtuTable` | Gut microbiome OTU table: 8 taxa x 6 samples (3 healthy, 3 IBD) |
| `metagenomics::demo_taxonomy` | `() -> Vec<DemoTaxonomy>` | Full taxonomy lineage for the 8 OTU table taxa |

### Structural Biology

| Function | Signature | Description |
|----------|-----------|-------------|
| `structural::insulin_chain_b` | `() -> DemoStructure` | Insulin chain B Cα trace (PDB 2INS), 30 residues |
| `structural::demo_ramachandran` | `() -> Vec<(&'static str, f64, f64)>` | 10 demo (label, phi, psi) angles for Ramachandran plot |

### Protocols

| Function | Signature | Description |
|----------|-----------|-------------|
| `protocols::all_protocols` | `() -> Vec<Protocol>` | All 16 protocol templates |
| `protocols::wet_lab_protocols` | `() -> Vec<Protocol>` | 10 wet-lab protocols |
| `protocols::dry_lab_protocols` | `() -> Vec<Protocol>` | 6 dry-lab (computational) protocols |
| `protocols::elisa` | `() -> Protocol` | ELISA (Enzyme-Linked Immunosorbent Assay) |
| `protocols::qpcr` | `() -> Protocol` | qPCR / RT-qPCR |
| `protocols::immunofluorescence` | `() -> Protocol` | Immunofluorescence staining |
| `protocols::cell_viability_mtt` | `() -> Protocol` | MTT Cell Viability Assay |
| `protocols::site_directed_mutagenesis` | `() -> Protocol` | Site-Directed Mutagenesis |
| `protocols::bacterial_transformation` | `() -> Protocol` | Bacterial Transformation |
| `protocols::coimmunoprecipitation` | `() -> Protocol` | Co-Immunoprecipitation (Co-IP) |
| `protocols::chipseq_library_prep` | `() -> Protocol` | ChIP-seq Library Preparation |
| `protocols::atacseq_library_prep` | `() -> Protocol` | ATAC-seq Library Preparation |
| `protocols::hic_library_prep` | `() -> Protocol` | Hi-C Library Preparation |
| `protocols::gwas_pipeline` | `() -> Protocol` | GWAS Pipeline |
| `protocols::longread_genome_assembly` | `() -> Protocol` | Long-Read Genome Assembly |
| `protocols::proteomics_dia_search` | `() -> Protocol` | DIA Proteomics Analysis |
| `protocols::spatial_transcriptomics_analysis` | `() -> Protocol` | Spatial Transcriptomics Analysis |
| `protocols::hic_analysis` | `() -> Protocol` | Hi-C Analysis Pipeline |
| `protocols::alphafold_structure_prediction` | `() -> Protocol` | AlphaFold Structure Prediction |
