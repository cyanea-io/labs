// ---------------------------------------------------------------------------
// @cyanea/bio -- Worker dispatch registry
// ---------------------------------------------------------------------------
// Maps "Namespace.method" strings to the corresponding sync functions from
// index.ts.  The worker uses this to dispatch incoming call messages.
// ---------------------------------------------------------------------------

import { Seq, Align, Stats, ML, Chem, StructBio, Phylo, Core } from "./index.js";

/** A single entry in the dispatch registry. */
export interface RegistryEntry {
  /** The function to call. */
  fn: (...args: unknown[]) => unknown;
  /** If true, the worker will decompose the call for progress reporting. */
  progressive: boolean;
}

/**
 * Flat dispatch table: "Namespace.method" → { fn, progressive }.
 *
 * Progressive functions (alignBatch, umap) are decomposed by the worker
 * into smaller units of work with progress reporting between them.
 */
export const registry = new Map<string, RegistryEntry>([
  // ── Seq ──────────────────────────────────────────────────────────────────
  ["Seq.parseFasta",       { fn: (d: unknown) => Seq.parseFasta(d as string), progressive: false }],
  ["Seq.parseFastq",       { fn: (d: unknown) => Seq.parseFastq(d as string), progressive: false }],
  ["Seq.gcContent",        { fn: (s: unknown) => Seq.gcContent(s as string), progressive: false }],
  ["Seq.reverseComplement",{ fn: (s: unknown) => Seq.reverseComplement(s as string), progressive: false }],
  ["Seq.transcribe",       { fn: (s: unknown) => Seq.transcribe(s as string), progressive: false }],
  ["Seq.translate",        { fn: (s: unknown) => Seq.translate(s as string), progressive: false }],
  ["Seq.validate",         { fn: (s: unknown, a: unknown) => Seq.validate(s as string, a as "dna" | "rna" | "protein"), progressive: false }],

  // ── Align ────────────────────────────────────────────────────────────────
  ["Align.alignDna",       { fn: (q: unknown, t: unknown, m: unknown) => Align.alignDna(q as string, t as string, m as "local" | "global" | "semiglobal"), progressive: false }],
  ["Align.alignDnaCustom", { fn: (q: unknown, t: unknown, m: unknown, ms: unknown, mm: unknown, go: unknown, ge: unknown) => Align.alignDnaCustom(q as string, t as string, m as "local" | "global" | "semiglobal", ms as number, mm as number, go as number, ge as number), progressive: false }],
  ["Align.alignProtein",   { fn: (q: unknown, t: unknown, m: unknown, mx: unknown) => Align.alignProtein(q as string, t as string, m as "local" | "global" | "semiglobal", mx as "blosum62" | "blosum45" | "blosum80" | "pam250"), progressive: false }],
  ["Align.alignBatch",     { fn: (...args: unknown[]) => Align.alignBatch(args[0] as {query: string; target: string}[], args[1] as "local" | "global" | "semiglobal", args[2] as number, args[3] as number, args[4] as number, args[5] as number), progressive: true }],

  // ── Stats ────────────────────────────────────────────────────────────────
  ["Stats.describe",           { fn: (d: unknown) => Stats.describe(d as number[]), progressive: false }],
  ["Stats.pearson",            { fn: (x: unknown, y: unknown) => Stats.pearson(x as number[], y as number[]), progressive: false }],
  ["Stats.spearman",           { fn: (x: unknown, y: unknown) => Stats.spearman(x as number[], y as number[]), progressive: false }],
  ["Stats.tTest",              { fn: (d: unknown, m: unknown) => Stats.tTest(d as number[], m as number), progressive: false }],
  ["Stats.tTestTwoSample",     { fn: (x: unknown, y: unknown, e: unknown) => Stats.tTestTwoSample(x as number[], y as number[], e as boolean), progressive: false }],
  ["Stats.mannWhitneyU",       { fn: (x: unknown, y: unknown) => Stats.mannWhitneyU(x as number[], y as number[]), progressive: false }],
  ["Stats.bonferroni",         { fn: (p: unknown) => Stats.bonferroni(p as number[]), progressive: false }],
  ["Stats.benjaminiHochberg",  { fn: (p: unknown) => Stats.benjaminiHochberg(p as number[]), progressive: false }],

  // ── ML ───────────────────────────────────────────────────────────────────
  ["ML.kmerCount",         { fn: (s: unknown, k: unknown) => ML.kmerCount(s as string, k as number), progressive: false }],
  ["ML.euclideanDistance",  { fn: (a: unknown, b: unknown) => ML.euclideanDistance(a as number[], b as number[]), progressive: false }],
  ["ML.manhattanDistance",  { fn: (a: unknown, b: unknown) => ML.manhattanDistance(a as number[], b as number[]), progressive: false }],
  ["ML.hammingDistance",    { fn: (a: unknown, b: unknown) => ML.hammingDistance(a as string, b as string), progressive: false }],
  ["ML.cosineSimilarity",   { fn: (a: unknown, b: unknown) => ML.cosineSimilarity(a as number[], b as number[]), progressive: false }],
  ["ML.umap",              { fn: (...args: unknown[]) => ML.umap(args[0] as number[], args[1] as number, args[2] as number, args[3] as number, args[4] as number, args[5] as number, args[6] as "euclidean" | "manhattan" | "cosine"), progressive: true }],

  // ── Chem ─────────────────────────────────────────────────────────────────
  ["Chem.properties",    { fn: (s: unknown) => Chem.properties(s as string), progressive: false }],
  ["Chem.canonical",     { fn: (s: unknown) => Chem.canonical(s as string), progressive: false }],
  ["Chem.fingerprint",   { fn: (s: unknown, r: unknown, n: unknown) => Chem.fingerprint(s as string, r as number, n as number), progressive: false }],
  ["Chem.tanimoto",      { fn: (a: unknown, b: unknown) => Chem.tanimoto(a as string, b as string), progressive: false }],
  ["Chem.substructure",  { fn: (m: unknown, p: unknown) => Chem.substructure(m as string, p as string), progressive: false }],

  // ── StructBio ────────────────────────────────────────────────────────────
  ["StructBio.pdbInfo",             { fn: (p: unknown) => StructBio.pdbInfo(p as string), progressive: false }],
  ["StructBio.secondaryStructure",  { fn: (p: unknown) => StructBio.secondaryStructure(p as string), progressive: false }],
  ["StructBio.rmsd",                { fn: (a: unknown, b: unknown) => StructBio.rmsd(a as [number, number, number][], b as [number, number, number][]), progressive: false }],

  // ── Phylo ────────────────────────────────────────────────────────────────
  ["Phylo.newickInfo",           { fn: (n: unknown) => Phylo.newickInfo(n as string), progressive: false }],
  ["Phylo.evolutionaryDistance", { fn: (a: unknown, b: unknown, m: unknown) => Phylo.evolutionaryDistance(a as string, b as string, m as "p" | "jc" | "k2p"), progressive: false }],
  ["Phylo.buildUpgma",          { fn: (l: unknown, m: unknown) => Phylo.buildUpgma(l as string[], m as number[][]), progressive: false }],
  ["Phylo.buildNj",             { fn: (l: unknown, m: unknown) => Phylo.buildNj(l as string[], m as number[][]), progressive: false }],
  ["Phylo.rfDistance",           { fn: (a: unknown, b: unknown) => Phylo.rfDistance(a as string, b as string), progressive: false }],

  // ── Core ─────────────────────────────────────────────────────────────────
  ["Core.sha256",          { fn: (d: unknown) => Core.sha256(d as string), progressive: false }],
  ["Core.zstdCompress",    { fn: (d: unknown, l: unknown) => Core.zstdCompress(d as string, l as number), progressive: false }],
  ["Core.zstdDecompress",  { fn: (d: unknown) => Core.zstdDecompress(d as number[]), progressive: false }],
]);
