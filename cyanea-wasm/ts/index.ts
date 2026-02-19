// ---------------------------------------------------------------------------
// @cyanea/bio -- Typed wrapper around cyanea-wasm WASM bindings
// ---------------------------------------------------------------------------
// Imports the raw WASM functions (which return JSON strings), parses the
// JSON envelope, and exposes a clean, fully typed TypeScript API.
//
// Usage:
//   import init from "../pkg/cyanea_wasm.js";
//   import { Seq, Align, Stats, ML, Chem, StructBio, Phylo, IO, Omics, Core } from "./index.js";
//
//   await init();
//   const stats = Seq.parseFasta(">seq1\nACGT\n");
//   console.log(stats.gc_content);
// ---------------------------------------------------------------------------

import type {
  FastaStats,
  FastqRecord,
  CigarOp,
  CigarStats,
  AlignmentResult,
  AlignmentMode,
  SubstitutionMatrix,
  SeqPair,
  DescriptiveStats,
  TestResult,
  KmerCounts,
  UmapResult,
  PcaResult,
  TsneResult,
  KmeansResult,
  DistanceMetric,
  MinHashSketch,
  MinHashComparison,
  Pileup,
  DepthStats,
  MolecularProperties,
  Fingerprint,
  SubstructureResult,
  StructureInfo,
  SecondaryStructure,
  TreeInfo,
  RFDistance,
  DistanceModel,
  Alphabet,
  WasmResult,
  // New T15 types
  RnaStructure,
  ProteinProperties,
  SimulatedRead,
  CodonUsage,
  AssemblyStats,
  MsaResult,
  PoaConsensus,
  KmResult,
  LogRankResult,
  CoxPhResult,
  WrightFisherResult,
  TajimaD,
  FstResult,
  RandomForestResult,
  GbdtRegressionResult,
  GbdtClassifyResult,
  HmmViterbiResult,
  ConfusionMatrix,
  RocCurve,
  PrCurve,
  CvResult,
  FeatureSelection,
  SdfMolecule,
  MaccsFingerprint,
  ContactMap,
  RamachandranEntry,
  MmcifInfo,
  KabschResult,
  NexusFile,
  SimulatedAlignment,
  CoalescentTree,
  VcfVariant,
  BedRecord,
  Gff3Gene,
  BlastXmlResult,
  BedGraphRecord,
  GfaGraph,
  GenomicInterval,
  VariantEffect,
  CnvSegment,
  CpgIsland,
  SpatialAutocorrelation,
  GearysC,
  LiftoverResult,
  JaccardResult,
  ClosestResult,
} from "./types.js";

export type {
  FastaStats,
  FastqRecord,
  CigarOp,
  CigarStats,
  AlignmentResult,
  AlignmentMode,
  SubstitutionMatrix,
  SeqPair,
  DescriptiveStats,
  TestResult,
  KmerCounts,
  UmapResult,
  PcaResult,
  TsneResult,
  KmeansResult,
  DistanceMetric,
  MinHashSketch,
  MinHashComparison,
  Pileup,
  DepthStats,
  MolecularProperties,
  Fingerprint,
  SubstructureResult,
  StructureInfo,
  SecondaryStructure,
  TreeInfo,
  RFDistance,
  DistanceModel,
  Alphabet,
  // New T15 types
  RnaStructure,
  ProteinProperties,
  SimulatedRead,
  CodonUsage,
  AssemblyStats,
  MsaResult,
  PoaConsensus,
  KmResult,
  LogRankResult,
  CoxPhResult,
  WrightFisherResult,
  TajimaD,
  FstResult,
  RandomForestResult,
  GbdtRegressionResult,
  GbdtClassifyResult,
  HmmViterbiResult,
  ConfusionMatrix,
  RocCurve,
  PrCurve,
  CvResult,
  FeatureSelection,
  SdfMolecule,
  MaccsFingerprint,
  ContactMap,
  RamachandranEntry,
  MmcifInfo,
  KabschResult,
  NexusFile,
  SimulatedAlignment,
  CoalescentTree,
  VcfVariant,
  BedRecord,
  Gff3Gene,
  BlastXmlResult,
  BedGraphRecord,
  GfaGraph,
  GenomicInterval,
  VariantEffect,
  CnvSegment,
  CpgIsland,
  SpatialAutocorrelation,
  GearysC,
  LiftoverResult,
  JaccardResult,
  ClosestResult,
} from "./types.js";

// Re-export the envelope types for advanced users who want raw access.
export type { WasmOk, WasmErr, WasmResult } from "./types.js";

// ── Raw WASM imports ───────────────────────────────────────────────────────
// These are the raw wasm-bindgen functions that return JSON strings.
// The import path assumes wasm-pack output in ../pkg/.
// All imports are aliased with a _raw prefix to avoid shadowing by the
// namespace wrapper functions that share similar names.

import {
  // seq
  parse_fasta as _raw_parse_fasta,
  parse_fastq as _raw_parse_fastq,
  gc_content_json as _raw_gc_content_json,
  reverse_complement as _raw_reverse_complement,
  transcribe as _raw_transcribe,
  translate as _raw_translate,
  validate as _raw_validate,
  minhash_sketch as _raw_minhash_sketch,
  minhash_compare as _raw_minhash_compare,
  rna_fold_nussinov as _raw_rna_fold_nussinov,
  rna_fold_zuker as _raw_rna_fold_zuker,
  protein_props as _raw_protein_props,
  simulate_reads as _raw_simulate_reads,
  codon_usage as _raw_codon_usage,
  assembly_stats_json as _raw_assembly_stats_json,
  // align
  align_dna as _raw_align_dna,
  align_dna_custom as _raw_align_dna_custom,
  align_protein as _raw_align_protein,
  align_batch as _raw_align_batch,
  parse_cigar as _raw_parse_cigar,
  validate_cigar as _raw_validate_cigar,
  cigar_stats as _raw_cigar_stats,
  cigar_to_alignment as _raw_cigar_to_alignment,
  alignment_to_cigar as _raw_alignment_to_cigar,
  generate_md_tag as _raw_generate_md_tag,
  merge_cigar as _raw_merge_cigar,
  reverse_cigar as _raw_reverse_cigar,
  collapse_cigar as _raw_collapse_cigar,
  hard_clip_to_soft as _raw_hard_clip_to_soft,
  split_cigar as _raw_split_cigar,
  progressive_msa as _raw_progressive_msa,
  poa_consensus as _raw_poa_consensus,
  align_banded as _raw_align_banded,
  // stats
  describe as _raw_describe,
  pearson as _raw_pearson,
  spearman as _raw_spearman,
  t_test as _raw_t_test,
  t_test_two_sample as _raw_t_test_two_sample,
  mann_whitney_u as _raw_mann_whitney_u,
  bonferroni as _raw_bonferroni,
  benjamini_hochberg as _raw_benjamini_hochberg,
  kaplan_meier as _raw_kaplan_meier,
  log_rank_test as _raw_log_rank_test,
  cox_ph as _raw_cox_ph,
  wright_fisher as _raw_wright_fisher,
  permutation_test as _raw_permutation_test,
  bootstrap_ci as _raw_bootstrap_ci,
  shannon_index as _raw_shannon_index,
  simpson_index as _raw_simpson_index,
  bray_curtis as _raw_bray_curtis,
  fst_hudson as _raw_fst_hudson,
  tajimas_d as _raw_tajimas_d,
  // ml
  kmer_count as _raw_kmer_count,
  euclidean_distance as _raw_euclidean_distance,
  manhattan_distance as _raw_manhattan_distance,
  hamming_distance as _raw_hamming_distance,
  cosine_similarity as _raw_cosine_similarity,
  umap as _raw_umap,
  pca as _raw_pca,
  tsne as _raw_tsne,
  kmeans as _raw_kmeans,
  random_forest_classify as _raw_random_forest_classify,
  gbdt_regression as _raw_gbdt_regression,
  gbdt_classify as _raw_gbdt_classify,
  hmm_viterbi as _raw_hmm_viterbi,
  hmm_likelihood as _raw_hmm_likelihood,
  confusion_matrix as _raw_confusion_matrix,
  roc_curve as _raw_roc_curve,
  pr_curve as _raw_pr_curve,
  cross_validate_rf as _raw_cross_validate_rf,
  feature_importance_variance as _raw_feature_importance_variance,
  // chem
  smiles_properties as _raw_smiles_properties,
  canonical as _raw_canonical,
  smiles_fingerprint as _raw_smiles_fingerprint,
  tanimoto as _raw_tanimoto,
  smiles_substructure as _raw_smiles_substructure,
  parse_sdf as _raw_parse_sdf,
  maccs_fingerprint as _raw_maccs_fingerprint,
  tanimoto_maccs as _raw_tanimoto_maccs,
  // struct_bio
  pdb_info as _raw_pdb_info,
  pdb_secondary_structure as _raw_pdb_secondary_structure,
  rmsd as _raw_rmsd,
  contact_map as _raw_contact_map,
  ramachandran_analysis as _raw_ramachandran_analysis,
  parse_mmcif as _raw_parse_mmcif,
  kabsch_align as _raw_kabsch_align,
  // phylo
  newick_info as _raw_newick_info,
  evolutionary_distance as _raw_evolutionary_distance,
  build_upgma as _raw_build_upgma,
  build_nj as _raw_build_nj,
  rf_distance as _raw_rf_distance,
  parse_nexus as _raw_parse_nexus,
  write_nexus as _raw_write_nexus,
  simulate_evolution as _raw_simulate_evolution,
  simulate_coalescent as _raw_simulate_coalescent,
  simulate_coalescent_growth as _raw_simulate_coalescent_growth,
  // io
  pileup_from_sam as _raw_pileup_from_sam,
  depth_stats_from_sam as _raw_depth_stats_from_sam,
  pileup_to_mpileup_text as _raw_pileup_to_mpileup_text,
  parse_vcf_text as _raw_parse_vcf_text,
  parse_bed_text as _raw_parse_bed_text,
  parse_gff3_text as _raw_parse_gff3_text,
  parse_blast_xml as _raw_parse_blast_xml,
  parse_bedgraph as _raw_parse_bedgraph,
  parse_gfa as _raw_parse_gfa,
  ncbi_fetch_url as _raw_ncbi_fetch_url,
  // omics
  merge_intervals as _raw_merge_intervals,
  intersect_intervals as _raw_intersect_intervals,
  subtract_intervals as _raw_subtract_intervals,
  complement_intervals as _raw_complement_intervals,
  closest_intervals as _raw_closest_intervals,
  jaccard_intervals as _raw_jaccard_intervals,
  make_windows as _raw_make_windows,
  liftover_interval as _raw_liftover_interval,
  annotate_variant as _raw_annotate_variant,
  cbs_segment as _raw_cbs_segment,
  bisulfite_convert as _raw_bisulfite_convert,
  find_cpg_islands as _raw_find_cpg_islands,
  morans_i as _raw_morans_i,
  gearys_c as _raw_gearys_c,
  // core
  sha256 as _raw_sha256,
  zstd_compress as _raw_zstd_compress,
  zstd_decompress as _raw_zstd_decompress,
} from "../pkg/cyanea_wasm.js";

// ── Error type ─────────────────────────────────────────────────────────────

/**
 * Error thrown when a WASM function returns an error envelope.
 * The `message` property contains the error string from the Rust side.
 */
export class CyaneaError extends Error {
  override readonly name = "CyaneaError";

  constructor(message: string) {
    super(message);
  }
}

// ── JSON envelope unwrapper ────────────────────────────────────────────────

/**
 * Parse a JSON envelope string and return the success value, or throw
 * a `CyaneaError` if the envelope contains an error.
 */
function unwrap<T>(json: string): T {
  const parsed: WasmResult<T> = JSON.parse(json);
  if ("error" in parsed) {
    throw new CyaneaError(parsed.error);
  }
  return parsed.ok;
}

// ── Seq ────────────────────────────────────────────────────────────────────

export namespace Seq {
  export function parseFasta(data: string): FastaStats {
    return unwrap<FastaStats>(_raw_parse_fasta(data));
  }

  export function parseFastq(data: string): FastqRecord[] {
    return unwrap<FastqRecord[]>(_raw_parse_fastq(data));
  }

  export function gcContent(seq: string): number {
    return unwrap<number>(_raw_gc_content_json(seq));
  }

  export function reverseComplement(seq: string): string {
    return unwrap<string>(_raw_reverse_complement(seq));
  }

  export function transcribe(seq: string): string {
    return unwrap<string>(_raw_transcribe(seq));
  }

  export function translate(seq: string): string {
    return unwrap<string>(_raw_translate(seq));
  }

  export function validate(seq: string, alphabet: Alphabet): boolean {
    return unwrap<boolean>(_raw_validate(seq, alphabet));
  }

  export function minhashSketch(seq: string, k: number, sketchSize: number): MinHashSketch {
    return unwrap<MinHashSketch>(_raw_minhash_sketch(seq, k, sketchSize));
  }

  export function minhashCompare(
    seqA: string, seqB: string, k: number, sketchSize: number,
  ): MinHashComparison {
    return unwrap<MinHashComparison>(_raw_minhash_compare(seqA, seqB, k, sketchSize));
  }

  /** Predict RNA secondary structure using Nussinov algorithm. */
  export function rnaFoldNussinov(seq: string): RnaStructure {
    return unwrap<RnaStructure>(_raw_rna_fold_nussinov(seq));
  }

  /** Predict RNA secondary structure using Zuker MFE algorithm. */
  export function rnaFoldZuker(seq: string): RnaStructure {
    return unwrap<RnaStructure>(_raw_rna_fold_zuker(seq));
  }

  /** Compute protein sequence properties (MW, pI, GRAVY, etc.). */
  export function proteinProperties(seq: string): ProteinProperties {
    return unwrap<ProteinProperties>(_raw_protein_props(seq));
  }

  /** Simulate sequencing reads from a reference sequence. */
  export function simulateReads(refSeq: string, configJson: string): SimulatedRead[] {
    return unwrap<SimulatedRead[]>(_raw_simulate_reads(refSeq, configJson));
  }

  /** Compute codon usage from a coding DNA sequence. */
  export function codonUsage(seq: string): CodonUsage {
    return unwrap<CodonUsage>(_raw_codon_usage(seq));
  }

  /** Compute assembly statistics for a set of contigs. */
  export function assemblyStats(contigsJson: string): AssemblyStats {
    return unwrap<AssemblyStats>(_raw_assembly_stats_json(contigsJson));
  }
}

// ── Align ──────────────────────────────────────────────────────────────────

export namespace Align {
  export function alignDna(
    query: string, target: string, mode: AlignmentMode,
  ): AlignmentResult {
    return unwrap<AlignmentResult>(_raw_align_dna(query, target, mode));
  }

  export function alignDnaCustom(
    query: string, target: string, mode: AlignmentMode,
    matchScore: number, mismatchScore: number, gapOpen: number, gapExtend: number,
  ): AlignmentResult {
    return unwrap<AlignmentResult>(
      _raw_align_dna_custom(query, target, mode, matchScore, mismatchScore, gapOpen, gapExtend),
    );
  }

  export function alignProtein(
    query: string, target: string, mode: AlignmentMode, matrix: SubstitutionMatrix,
  ): AlignmentResult {
    return unwrap<AlignmentResult>(_raw_align_protein(query, target, mode, matrix));
  }

  export function alignBatch(
    pairs: SeqPair[], mode: AlignmentMode,
    matchScore: number, mismatchScore: number, gapOpen: number, gapExtend: number,
  ): AlignmentResult[] {
    return unwrap<AlignmentResult[]>(
      _raw_align_batch(JSON.stringify(pairs), mode, matchScore, mismatchScore, gapOpen, gapExtend),
    );
  }

  export function parseCigar(cigar: string): CigarOp[] {
    return unwrap<CigarOp[]>(_raw_parse_cigar(cigar));
  }

  export function validateCigar(cigar: string): boolean {
    return unwrap<boolean>(_raw_validate_cigar(cigar));
  }

  export function cigarStats(cigar: string): CigarStats {
    return unwrap<CigarStats>(_raw_cigar_stats(cigar));
  }

  export function cigarToAlignment(
    cigar: string, query: string, target: string,
  ): { aligned_query: number[]; aligned_target: number[] } {
    return unwrap<{ aligned_query: number[]; aligned_target: number[] }>(
      _raw_cigar_to_alignment(cigar, query, target),
    );
  }

  export function alignmentToCigar(query: string, target: string): string {
    return unwrap<string>(_raw_alignment_to_cigar(query, target));
  }

  export function generateMdTag(cigar: string, query: string, reference: string): string {
    return unwrap<string>(_raw_generate_md_tag(cigar, query, reference));
  }

  export function mergeCigar(cigar: string): string {
    return unwrap<string>(_raw_merge_cigar(cigar));
  }

  export function reverseCigar(cigar: string): string {
    return unwrap<string>(_raw_reverse_cigar(cigar));
  }

  export function collapseCigar(cigar: string): string {
    return unwrap<string>(_raw_collapse_cigar(cigar));
  }

  export function hardClipToSoft(cigar: string): string {
    return unwrap<string>(_raw_hard_clip_to_soft(cigar));
  }

  export function splitCigar(cigar: string, refPos: number): { left: string; right: string } {
    return unwrap<{ left: string; right: string }>(_raw_split_cigar(cigar, refPos));
  }

  /** Progressive multiple sequence alignment. */
  export function progressiveMsa(
    seqsJson: string, matchScore: number, mismatchScore: number,
    gapOpen: number, gapExtend: number,
  ): MsaResult {
    return unwrap<MsaResult>(
      _raw_progressive_msa(seqsJson, matchScore, mismatchScore, gapOpen, gapExtend),
    );
  }

  /** Partial-order alignment consensus. */
  export function poaConsensus(seqsJson: string): PoaConsensus {
    return unwrap<PoaConsensus>(_raw_poa_consensus(seqsJson));
  }

  /** Banded pairwise alignment. */
  export function alignBanded(
    query: string, target: string, mode: AlignmentMode, bandwidth: number,
    matchScore: number, mismatchScore: number, gapOpen: number, gapExtend: number,
  ): AlignmentResult {
    return unwrap<AlignmentResult>(
      _raw_align_banded(query, target, mode, bandwidth, matchScore, mismatchScore, gapOpen, gapExtend),
    );
  }
}

// ── Stats ──────────────────────────────────────────────────────────────────

export namespace Stats {
  export function describe(data: number[]): DescriptiveStats {
    return unwrap<DescriptiveStats>(_raw_describe(JSON.stringify(data)));
  }

  export function pearson(x: number[], y: number[]): number {
    return unwrap<number>(_raw_pearson(JSON.stringify(x), JSON.stringify(y)));
  }

  export function spearman(x: number[], y: number[]): number {
    return unwrap<number>(_raw_spearman(JSON.stringify(x), JSON.stringify(y)));
  }

  export function tTest(data: number[], mu: number): TestResult {
    return unwrap<TestResult>(_raw_t_test(JSON.stringify(data), mu));
  }

  export function tTestTwoSample(x: number[], y: number[], equalVar: boolean): TestResult {
    return unwrap<TestResult>(
      _raw_t_test_two_sample(JSON.stringify(x), JSON.stringify(y), equalVar),
    );
  }

  export function mannWhitneyU(x: number[], y: number[]): TestResult {
    return unwrap<TestResult>(
      _raw_mann_whitney_u(JSON.stringify(x), JSON.stringify(y)),
    );
  }

  export function bonferroni(pValues: number[]): number[] {
    return unwrap<number[]>(_raw_bonferroni(JSON.stringify(pValues)));
  }

  export function benjaminiHochberg(pValues: number[]): number[] {
    return unwrap<number[]>(_raw_benjamini_hochberg(JSON.stringify(pValues)));
  }

  /** Kaplan-Meier survival analysis. */
  export function kaplanMeier(timesJson: string, statusJson: string): KmResult {
    return unwrap<KmResult>(_raw_kaplan_meier(timesJson, statusJson));
  }

  /** Log-rank test comparing two survival curves. */
  export function logRankTest(
    t1Json: string, s1Json: string, t2Json: string, s2Json: string,
  ): LogRankResult {
    return unwrap<LogRankResult>(_raw_log_rank_test(t1Json, s1Json, t2Json, s2Json));
  }

  /** Cox proportional hazards model. */
  export function coxPh(
    timesJson: string, statusJson: string, covariatesJson: string, nCovariates: number,
  ): CoxPhResult {
    return unwrap<CoxPhResult>(_raw_cox_ph(timesJson, statusJson, covariatesJson, nCovariates));
  }

  /** Wright-Fisher population simulation. */
  export function wrightFisher(popSize: number, initFreq: number, nGens: number, seed: number): WrightFisherResult {
    return unwrap<WrightFisherResult>(_raw_wright_fisher(popSize, initFreq, nGens, seed));
  }

  /** Permutation test for group differences. */
  export function permutationTest(
    valuesJson: string, groupSizesJson: string, nPerms: number, seed: number,
  ): TestResult {
    return unwrap<TestResult>(_raw_permutation_test(valuesJson, groupSizesJson, nPerms, seed));
  }

  /** Bootstrap confidence interval. */
  export function bootstrapCi(dataJson: string, nBootstrap: number, seed: number): number[] {
    return unwrap<number[]>(_raw_bootstrap_ci(dataJson, nBootstrap, seed));
  }

  /** Shannon diversity index. */
  export function shannonIndex(countsJson: string): number {
    return unwrap<number>(_raw_shannon_index(countsJson));
  }

  /** Simpson diversity index. */
  export function simpsonIndex(countsJson: string): number {
    return unwrap<number>(_raw_simpson_index(countsJson));
  }

  /** Bray-Curtis dissimilarity. */
  export function brayCurtis(aJson: string, bJson: string): number {
    return unwrap<number>(_raw_bray_curtis(aJson, bJson));
  }

  /** Hudson's Fst estimator. */
  export function fstHudson(pop1Json: string, pop2Json: string): FstResult {
    return unwrap<FstResult>(_raw_fst_hudson(pop1Json, pop2Json));
  }

  /** Tajima's D neutrality test. */
  export function tajimasD(genotypesJson: string): TajimaD {
    return unwrap<TajimaD>(_raw_tajimas_d(genotypesJson));
  }
}

// ── ML ─────────────────────────────────────────────────────────────────────

export namespace ML {
  export function kmerCount(seq: string, k: number): KmerCounts {
    return unwrap<KmerCounts>(_raw_kmer_count(seq, k));
  }

  export function euclideanDistance(a: number[], b: number[]): number {
    return unwrap<number>(_raw_euclidean_distance(JSON.stringify(a), JSON.stringify(b)));
  }

  export function manhattanDistance(a: number[], b: number[]): number {
    return unwrap<number>(_raw_manhattan_distance(JSON.stringify(a), JSON.stringify(b)));
  }

  export function hammingDistance(a: string, b: string): number {
    return unwrap<number>(_raw_hamming_distance(a, b));
  }

  export function cosineSimilarity(a: number[], b: number[]): number {
    return unwrap<number>(_raw_cosine_similarity(JSON.stringify(a), JSON.stringify(b)));
  }

  export function umap(
    data: number[], nFeatures: number, nComponents: number,
    nNeighbors: number, minDist: number, nEpochs: number, metric: DistanceMetric,
  ): UmapResult {
    return unwrap<UmapResult>(
      _raw_umap(JSON.stringify(data), nFeatures, nComponents, nNeighbors, minDist, nEpochs, metric),
    );
  }

  export function pca(data: number[], nFeatures: number, nComponents: number): PcaResult {
    return unwrap<PcaResult>(_raw_pca(JSON.stringify(data), nFeatures, nComponents));
  }

  export function tsne(
    data: number[], nFeatures: number, nComponents: number,
    perplexity: number, learningRate: number, nIter: number, seed: number,
  ): TsneResult {
    return unwrap<TsneResult>(
      _raw_tsne(JSON.stringify(data), nFeatures, nComponents, perplexity, learningRate, nIter, seed),
    );
  }

  export function kmeans(
    data: number[], nFeatures: number, nClusters: number, maxIter: number, seed: number,
  ): KmeansResult {
    return unwrap<KmeansResult>(
      _raw_kmeans(JSON.stringify(data), nFeatures, nClusters, maxIter, seed),
    );
  }

  /** Random forest classification. */
  export function randomForestClassify(
    dataJson: string, nFeatures: number, labelsJson: string, configJson: string,
  ): RandomForestResult {
    return unwrap<RandomForestResult>(
      _raw_random_forest_classify(dataJson, nFeatures, labelsJson, configJson),
    );
  }

  /** Gradient-boosted tree regression. */
  export function gbdtRegression(
    dataJson: string, nFeatures: number, targetsJson: string, configJson: string,
  ): GbdtRegressionResult {
    return unwrap<GbdtRegressionResult>(
      _raw_gbdt_regression(dataJson, nFeatures, targetsJson, configJson),
    );
  }

  /** Gradient-boosted tree classification. */
  export function gbdtClassify(
    dataJson: string, nFeatures: number, labelsJson: string, configJson: string,
  ): GbdtClassifyResult {
    return unwrap<GbdtClassifyResult>(
      _raw_gbdt_classify(dataJson, nFeatures, labelsJson, configJson),
    );
  }

  /** HMM Viterbi decoding. */
  export function hmmViterbi(
    nStates: number, nSymbols: number,
    initJson: string, transJson: string, emissJson: string, obsJson: string,
  ): HmmViterbiResult {
    return unwrap<HmmViterbiResult>(
      _raw_hmm_viterbi(nStates, nSymbols, initJson, transJson, emissJson, obsJson),
    );
  }

  /** HMM log-likelihood. */
  export function hmmLikelihood(
    nStates: number, nSymbols: number,
    initJson: string, transJson: string, emissJson: string, obsJson: string,
  ): number {
    return unwrap<number>(
      _raw_hmm_likelihood(nStates, nSymbols, initJson, transJson, emissJson, obsJson),
    );
  }

  /** Confusion matrix from actual and predicted labels. */
  export function confusionMatrix(actualJson: string, predictedJson: string): ConfusionMatrix {
    return unwrap<ConfusionMatrix>(_raw_confusion_matrix(actualJson, predictedJson));
  }

  /** ROC curve from scores and binary labels. */
  export function rocCurve(scoresJson: string, labelsJson: string): RocCurve {
    return unwrap<RocCurve>(_raw_roc_curve(scoresJson, labelsJson));
  }

  /** Precision-recall curve from scores and binary labels. */
  export function prCurve(scoresJson: string, labelsJson: string): PrCurve {
    return unwrap<PrCurve>(_raw_pr_curve(scoresJson, labelsJson));
  }

  /** K-fold cross-validation with random forest. */
  export function crossValidateRf(
    dataJson: string, nFeatures: number, labelsJson: string, k: number, seed: number,
  ): CvResult {
    return unwrap<CvResult>(_raw_cross_validate_rf(dataJson, nFeatures, labelsJson, k, seed));
  }

  /** Variance-threshold feature selection. */
  export function featureImportanceVariance(
    dataJson: string, nFeatures: number, threshold: number,
  ): FeatureSelection {
    return unwrap<FeatureSelection>(
      _raw_feature_importance_variance(dataJson, nFeatures, threshold),
    );
  }
}

// ── Chem ───────────────────────────────────────────────────────────────────

export namespace Chem {
  export function properties(smiles: string): MolecularProperties {
    return unwrap<MolecularProperties>(_raw_smiles_properties(smiles));
  }

  export function canonical(smiles: string): string {
    return unwrap<string>(_raw_canonical(smiles));
  }

  export function fingerprint(smiles: string, radius: number, nBits: number): Fingerprint {
    return unwrap<Fingerprint>(_raw_smiles_fingerprint(smiles, radius, nBits));
  }

  export function tanimoto(smiles1: string, smiles2: string): number {
    return unwrap<number>(_raw_tanimoto(smiles1, smiles2));
  }

  export function substructure(molecule: string, pattern: string): SubstructureResult {
    return unwrap<SubstructureResult>(_raw_smiles_substructure(molecule, pattern));
  }

  /** Parse SDF V2000/V3000 text and return molecules. */
  export function parseSdf(sdfText: string): SdfMolecule[] {
    return unwrap<SdfMolecule[]>(_raw_parse_sdf(sdfText));
  }

  /** Compute MACCS 166-key fingerprint. */
  export function maccsFingerprint(smiles: string): MaccsFingerprint {
    return unwrap<MaccsFingerprint>(_raw_maccs_fingerprint(smiles));
  }

  /** Tanimoto similarity using MACCS fingerprints. */
  export function tanimotoMaccs(smiles1: string, smiles2: string): number {
    return unwrap<number>(_raw_tanimoto_maccs(smiles1, smiles2));
  }
}

// ── StructBio ──────────────────────────────────────────────────────────────

export namespace StructBio {
  export function pdbInfo(pdbText: string): StructureInfo {
    return unwrap<StructureInfo>(_raw_pdb_info(pdbText));
  }

  export function secondaryStructure(pdbText: string): SecondaryStructure {
    return unwrap<SecondaryStructure>(_raw_pdb_secondary_structure(pdbText));
  }

  export function rmsd(
    coords1: [number, number, number][], coords2: [number, number, number][],
  ): number {
    return unwrap<number>(_raw_rmsd(JSON.stringify(coords1), JSON.stringify(coords2)));
  }

  /** Compute CA-CA contact map from PDB text. */
  export function contactMap(pdbText: string, cutoff: number): ContactMap {
    return unwrap<ContactMap>(_raw_contact_map(pdbText, cutoff));
  }

  /** Ramachandran analysis from PDB text. */
  export function ramachandran(pdbText: string): RamachandranEntry[] {
    return unwrap<RamachandranEntry[]>(_raw_ramachandran_analysis(pdbText));
  }

  /** Parse mmCIF text and return structure info. */
  export function parseMmcif(text: string): MmcifInfo {
    return unwrap<MmcifInfo>(_raw_parse_mmcif(text));
  }

  /** Kabsch superposition on two coordinate sets. */
  export function kabschAlign(
    coords1Json: string, coords2Json: string,
  ): KabschResult {
    return unwrap<KabschResult>(_raw_kabsch_align(coords1Json, coords2Json));
  }
}

// ── Phylo ──────────────────────────────────────────────────────────────────

export namespace Phylo {
  export function newickInfo(newick: string): TreeInfo {
    return unwrap<TreeInfo>(_raw_newick_info(newick));
  }

  export function evolutionaryDistance(seq1: string, seq2: string, model: DistanceModel): number {
    return unwrap<number>(_raw_evolutionary_distance(seq1, seq2, model));
  }

  export function buildUpgma(labels: string[], matrix: number[][]): string {
    return unwrap<string>(_raw_build_upgma(JSON.stringify(labels), JSON.stringify(matrix)));
  }

  export function buildNj(labels: string[], matrix: number[][]): string {
    return unwrap<string>(_raw_build_nj(JSON.stringify(labels), JSON.stringify(matrix)));
  }

  export function rfDistance(newick1: string, newick2: string): RFDistance {
    return unwrap<RFDistance>(_raw_rf_distance(newick1, newick2));
  }

  /** Parse NEXUS format text. */
  export function parseNexus(text: string): NexusFile {
    return unwrap<NexusFile>(_raw_parse_nexus(text));
  }

  /** Write NEXUS format from taxa and trees. */
  export function writeNexus(taxaJson: string, treesJson: string): string {
    return unwrap<string>(_raw_write_nexus(taxaJson, treesJson));
  }

  /** Simulate sequence evolution along a phylogenetic tree. */
  export function simulateEvolution(
    newick: string, seqLength: number, model: string, seed: number,
  ): SimulatedAlignment {
    return unwrap<SimulatedAlignment>(_raw_simulate_evolution(newick, seqLength, model, seed));
  }

  /** Simulate a coalescent tree. */
  export function simulateCoalescent(nSamples: number, popSize: number, seed: number): CoalescentTree {
    return unwrap<CoalescentTree>(_raw_simulate_coalescent(nSamples, popSize, seed));
  }

  /** Simulate a coalescent tree with exponential growth. */
  export function simulateCoalescentGrowth(
    nSamples: number, popSize: number, growthRate: number, seed: number,
  ): CoalescentTree {
    return unwrap<CoalescentTree>(
      _raw_simulate_coalescent_growth(nSamples, popSize, growthRate, seed),
    );
  }
}

// ── IO ─────────────────────────────────────────────────────────────────

export namespace IO {
  export function pileup(samText: string): Pileup[] {
    return unwrap<Pileup[]>(_raw_pileup_from_sam(samText));
  }

  export function depthStats(samText: string): DepthStats[] {
    return unwrap<DepthStats[]>(_raw_depth_stats_from_sam(samText));
  }

  export function mpileup(samText: string): string {
    return unwrap<string>(_raw_pileup_to_mpileup_text(samText));
  }

  /** Parse VCF text. */
  export function parseVcf(text: string): VcfVariant[] {
    return unwrap<VcfVariant[]>(_raw_parse_vcf_text(text));
  }

  /** Parse BED text. */
  export function parseBed(text: string): BedRecord[] {
    return unwrap<BedRecord[]>(_raw_parse_bed_text(text));
  }

  /** Parse GFF3 text. */
  export function parseGff3(text: string): Gff3Gene[] {
    return unwrap<Gff3Gene[]>(_raw_parse_gff3_text(text));
  }

  /** Parse BLAST XML text. */
  export function parseBlastXml(xml: string): BlastXmlResult {
    return unwrap<BlastXmlResult>(_raw_parse_blast_xml(xml));
  }

  /** Parse bedGraph text. */
  export function parseBedgraph(text: string): BedGraphRecord[] {
    return unwrap<BedGraphRecord[]>(_raw_parse_bedgraph(text));
  }

  /** Parse GFA assembly graph text. */
  export function parseGfa(text: string): GfaGraph {
    return unwrap<GfaGraph>(_raw_parse_gfa(text));
  }

  /** Build an NCBI E-utilities fetch URL. */
  export function ncbiFetchUrl(db: string, ids: string, rettype: string): string {
    return unwrap<string>(_raw_ncbi_fetch_url(db, ids, rettype));
  }
}

// ── Omics ──────────────────────────────────────────────────────────────────

export namespace Omics {
  /** Merge overlapping/adjacent intervals. */
  export function mergeIntervals(intervalsJson: string): GenomicInterval[] {
    return unwrap<GenomicInterval[]>(_raw_merge_intervals(intervalsJson));
  }

  /** Intersect two interval sets. */
  export function intersectIntervals(aJson: string, bJson: string): GenomicInterval[] {
    return unwrap<GenomicInterval[]>(_raw_intersect_intervals(aJson, bJson));
  }

  /** Subtract interval set B from A. */
  export function subtractIntervals(aJson: string, bJson: string): GenomicInterval[] {
    return unwrap<GenomicInterval[]>(_raw_subtract_intervals(aJson, bJson));
  }

  /** Complement intervals relative to genome. */
  export function complementIntervals(intervalsJson: string, genomeJson: string): GenomicInterval[] {
    return unwrap<GenomicInterval[]>(_raw_complement_intervals(intervalsJson, genomeJson));
  }

  /** Find closest intervals. */
  export function closestIntervals(aJson: string, bJson: string): ClosestResult[] {
    return unwrap<ClosestResult[]>(_raw_closest_intervals(aJson, bJson));
  }

  /** Jaccard similarity between interval sets. */
  export function jaccardIntervals(aJson: string, bJson: string): JaccardResult {
    return unwrap<JaccardResult>(_raw_jaccard_intervals(aJson, bJson));
  }

  /** Generate genomic windows. */
  export function makeWindows(genomeJson: string, windowSize: number): GenomicInterval[] {
    return unwrap<GenomicInterval[]>(_raw_make_windows(genomeJson, windowSize));
  }

  /** Liftover a genomic interval using chain file. */
  export function liftoverInterval(
    chainText: string, chrom: string, start: number, end: number,
  ): LiftoverResult {
    return unwrap<LiftoverResult>(_raw_liftover_interval(chainText, chrom, start, end));
  }

  /** Annotate a variant against gene definitions. */
  export function annotateVariant(variantJson: string, genesJson: string): VariantEffect[] {
    return unwrap<VariantEffect[]>(_raw_annotate_variant(variantJson, genesJson));
  }

  /** Circular binary segmentation for CNV detection. */
  export function cbsSegment(
    positionsJson: string, valuesJson: string, chrom: string, configJson: string,
  ): CnvSegment[] {
    return unwrap<CnvSegment[]>(_raw_cbs_segment(positionsJson, valuesJson, chrom, configJson));
  }

  /** Bisulfite convert a DNA sequence. */
  export function bisulfiteConvert(seq: string, methylatedJson: string): string {
    return unwrap<string>(_raw_bisulfite_convert(seq, methylatedJson));
  }

  /** Find CpG islands in a DNA sequence. */
  export function findCpgIslands(seq: string, chrom: string): CpgIsland[] {
    return unwrap<CpgIsland[]>(_raw_find_cpg_islands(seq, chrom));
  }

  /** Moran's I spatial autocorrelation. */
  export function moransI(valuesJson: string, neighborsJson: string): SpatialAutocorrelation {
    return unwrap<SpatialAutocorrelation>(_raw_morans_i(valuesJson, neighborsJson));
  }

  /** Geary's C spatial autocorrelation. */
  export function gearysC(valuesJson: string, neighborsJson: string): GearysC {
    return unwrap<GearysC>(_raw_gearys_c(valuesJson, neighborsJson));
  }
}

// ── Core Utilities ─────────────────────────────────────────────────────────

export namespace Core {
  export function sha256(data: string): string {
    return unwrap<string>(_raw_sha256(data));
  }

  export function zstdCompress(data: string, level: number): number[] {
    return unwrap<number[]>(_raw_zstd_compress(data, level));
  }

  export function zstdDecompress(data: number[]): string {
    return unwrap<string>(_raw_zstd_decompress(JSON.stringify(data)));
  }
}
