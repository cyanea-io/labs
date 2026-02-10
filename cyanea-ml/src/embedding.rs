//! Sequence embedding via composition vectors and k-mer frequency profiles.
//!
//! Provides methods for representing biological sequences as fixed-length
//! numerical vectors suitable for machine learning:
//!
//! - **K-mer frequency embedding** — represent a sequence by its normalized
//!   k-mer frequency distribution
//! - **Composition vector** — nucleotide or amino acid composition as a
//!   fixed-length vector
//! - **Batch embedding** — embed multiple sequences for downstream distance
//!   computation or clustering

use cyanea_core::{CyaneaError, Result};

use crate::encoding::Alphabet;
use crate::kmer::KmerCounter;

/// Configuration for sequence embedding.
#[derive(Debug, Clone)]
pub struct EmbeddingConfig {
    /// K-mer size for frequency embedding.
    pub k: usize,
    /// Alphabet to use (determines vector dimensionality).
    pub alphabet: Alphabet,
    /// Whether to L2-normalize the output vectors.
    pub normalize: bool,
}

impl Default for EmbeddingConfig {
    fn default() -> Self {
        Self {
            k: 3,
            alphabet: Alphabet::Dna,
            normalize: true,
        }
    }
}

/// A fixed-length embedding vector for a sequence.
#[derive(Debug, Clone)]
pub struct SequenceEmbedding {
    /// The embedding vector.
    pub vector: Vec<f64>,
    /// Dimensionality of the embedding.
    pub dim: usize,
}

/// Embed a single sequence as a normalized k-mer frequency vector.
///
/// The output dimension is `alphabet.size().pow(k)`.
///
/// # Errors
///
/// Returns an error if the sequence is empty or k is invalid.
pub fn kmer_embedding(seq: &[u8], config: &EmbeddingConfig) -> Result<SequenceEmbedding> {
    if seq.is_empty() {
        return Err(CyaneaError::InvalidInput("empty sequence".into()));
    }

    let counter = KmerCounter::new(config.k)?;
    let counts = counter.count_sequence(seq);
    let mut vector = counts.to_frequency_vector(config.alphabet);

    // Normalize to frequencies
    let total: f64 = vector.iter().sum();
    if total > 0.0 {
        for v in vector.iter_mut() {
            *v /= total;
        }
    }

    if config.normalize {
        l2_normalize(&mut vector);
    }

    let dim = vector.len();
    Ok(SequenceEmbedding { vector, dim })
}

/// Embed a sequence as a simple composition vector.
///
/// Each position in the vector corresponds to a symbol in the alphabet,
/// giving the frequency of that symbol in the sequence.
///
/// # Errors
///
/// Returns an error if the sequence is empty.
pub fn composition_vector(seq: &[u8], alphabet: Alphabet) -> Result<SequenceEmbedding> {
    if seq.is_empty() {
        return Err(CyaneaError::InvalidInput("empty sequence".into()));
    }

    let syms = alphabet.symbols();
    let mut vector = vec![0.0; syms.len()];
    let mut total = 0usize;

    for &b in seq {
        let upper = b.to_ascii_uppercase();
        if let Some(pos) = syms.iter().position(|&s| s == upper) {
            vector[pos] += 1.0;
            total += 1;
        }
    }

    if total > 0 {
        let t = total as f64;
        for v in vector.iter_mut() {
            *v /= t;
        }
    }

    let dim = vector.len();
    Ok(SequenceEmbedding { vector, dim })
}

/// Batch-embed multiple sequences using k-mer frequency vectors.
///
/// Returns a matrix where each row is a sequence embedding.
///
/// # Errors
///
/// Returns an error if any sequence is empty or config is invalid.
pub fn batch_embed(sequences: &[&[u8]], config: &EmbeddingConfig) -> Result<Vec<SequenceEmbedding>> {
    if sequences.is_empty() {
        return Err(CyaneaError::InvalidInput("empty sequence list".into()));
    }
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        sequences
            .par_iter()
            .enumerate()
            .map(|(i, seq)| {
                kmer_embedding(seq, config).map_err(|e| {
                    CyaneaError::InvalidInput(format!("sequence {}: {}", i, e))
                })
            })
            .collect()
    }
    #[cfg(not(feature = "parallel"))]
    sequences
        .iter()
        .enumerate()
        .map(|(i, seq)| {
            kmer_embedding(seq, config).map_err(|e| {
                CyaneaError::InvalidInput(format!("sequence {}: {}", i, e))
            })
        })
        .collect()
}

/// Compute pairwise cosine distances between embedding vectors.
///
/// Returns a condensed distance vector of length `n*(n-1)/2`.
pub fn pairwise_cosine_distances(embeddings: &[SequenceEmbedding]) -> Result<Vec<f64>> {
    let n = embeddings.len();
    if n < 2 {
        return Err(CyaneaError::InvalidInput(
            "need at least 2 embeddings".into(),
        ));
    }

    #[cfg(feature = "parallel")]
    let distances = {
        use rayon::prelude::*;
        (0..n)
            .into_par_iter()
            .map(|i| {
                ((i + 1)..n)
                    .map(|j| {
                        let sim = cosine_sim(&embeddings[i].vector, &embeddings[j].vector);
                        1.0 - sim
                    })
                    .collect::<Vec<f64>>()
            })
            .collect::<Vec<_>>()
            .into_iter()
            .flatten()
            .collect::<Vec<f64>>()
    };
    #[cfg(not(feature = "parallel"))]
    let distances = {
        let mut distances = Vec::with_capacity(n * (n - 1) / 2);
        for i in 0..n {
            for j in (i + 1)..n {
                let sim = cosine_sim(&embeddings[i].vector, &embeddings[j].vector);
                distances.push(1.0 - sim);
            }
        }
        distances
    };
    Ok(distances)
}

fn cosine_sim(a: &[f64], b: &[f64]) -> f64 {
    let mut dot = 0.0;
    let mut na = 0.0;
    let mut nb = 0.0;
    for (x, y) in a.iter().zip(b) {
        dot += x * y;
        na += x * x;
        nb += y * y;
    }
    let denom = na.sqrt() * nb.sqrt();
    if denom == 0.0 {
        0.0
    } else {
        dot / denom
    }
}

fn l2_normalize(v: &mut [f64]) {
    let norm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    if norm > 0.0 {
        for val in v.iter_mut() {
            *val /= norm;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kmer_embedding_basic() {
        let config = EmbeddingConfig {
            k: 2,
            alphabet: Alphabet::Dna,
            normalize: false,
        };
        let emb = kmer_embedding(b"ACGTACGT", &config).unwrap();
        assert_eq!(emb.dim, 16); // 4^2
        let sum: f64 = emb.vector.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10); // frequencies sum to 1
    }

    #[test]
    fn kmer_embedding_normalized() {
        let config = EmbeddingConfig::default();
        let emb = kmer_embedding(b"ACGTACGT", &config).unwrap();
        let norm: f64 = emb.vector.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!((norm - 1.0).abs() < 1e-10);
    }

    #[test]
    fn kmer_embedding_empty_error() {
        let config = EmbeddingConfig::default();
        assert!(kmer_embedding(b"", &config).is_err());
    }

    #[test]
    fn composition_vector_dna() {
        let emb = composition_vector(b"AACG", Alphabet::Dna).unwrap();
        assert_eq!(emb.dim, 4);
        assert!((emb.vector[0] - 0.5).abs() < 1e-10); // A = 2/4
        assert!((emb.vector[1] - 0.25).abs() < 1e-10); // C = 1/4
    }

    #[test]
    fn composition_vector_protein() {
        let emb = composition_vector(b"AAW", Alphabet::Protein).unwrap();
        assert_eq!(emb.dim, 20);
    }

    #[test]
    fn batch_embed_basic() {
        let config = EmbeddingConfig::default();
        let seqs: Vec<&[u8]> = vec![b"ACGTACGT", b"TTTTAAAA", b"GCGCGCGC"];
        let embeddings = batch_embed(&seqs, &config).unwrap();
        assert_eq!(embeddings.len(), 3);
        assert_eq!(embeddings[0].dim, embeddings[1].dim);
    }

    #[test]
    fn pairwise_distances_basic() {
        let config = EmbeddingConfig {
            k: 2,
            alphabet: Alphabet::Dna,
            normalize: true,
        };
        let seqs: Vec<&[u8]> = vec![b"ACGTACGT", b"ACGTACGT", b"TTTTTTTT"];
        let embeddings = batch_embed(&seqs, &config).unwrap();
        let dists = pairwise_cosine_distances(&embeddings).unwrap();
        assert_eq!(dists.len(), 3); // 3*(3-1)/2
        assert!(dists[0] < 1e-10); // identical sequences → distance 0
        assert!(dists[1] > 0.0); // different sequences
    }

    #[test]
    fn pairwise_distances_too_few() {
        let emb = vec![SequenceEmbedding {
            vector: vec![1.0],
            dim: 1,
        }];
        assert!(pairwise_cosine_distances(&emb).is_err());
    }
}
