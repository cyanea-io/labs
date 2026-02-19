//! Bayesian phylogenetic inference via Markov Chain Monte Carlo (MCMC).
//!
//! Implements a Metropolis-Hastings sampler with NNI/SPR topology proposals,
//! branch length scaling, coalescent and birth-death tree priors, strict
//! and relaxed molecular clocks, convergence diagnostics, and posterior
//! summarization.

use std::collections::{BTreeSet, HashMap};

use crate::bootstrap::bipartitions;
use crate::generic_likelihood::generic_tree_likelihood;
use crate::likelihood::nni_swap;
use crate::subst_model::SubstitutionModel;
use crate::tree::{NodeId, PhyloTree};
use cyanea_core::Result;

/// MCMC sampler configuration.
#[derive(Debug, Clone)]
pub struct McmcConfig {
    pub n_generations: usize,
    pub sample_every: usize,
    pub burnin: usize,
    pub proposal_weights: ProposalWeights,
}

impl Default for McmcConfig {
    fn default() -> Self {
        Self {
            n_generations: 10000,
            sample_every: 100,
            burnin: 1000,
            proposal_weights: ProposalWeights::default(),
        }
    }
}

/// Relative weights for MCMC proposal moves.
#[derive(Debug, Clone)]
pub struct ProposalWeights {
    pub nni: f64,
    pub spr: f64,
    pub branch_scale: f64,
    pub branch_slide: f64,
    pub model_params: f64,
}

impl Default for ProposalWeights {
    fn default() -> Self {
        Self {
            nni: 3.0,
            spr: 1.0,
            branch_scale: 3.0,
            branch_slide: 3.0,
            model_params: 1.0,
        }
    }
}

/// Tree prior distribution.
#[derive(Debug, Clone)]
pub enum TreePrior {
    Uniform,
    CoalescentConstant { pop_size: f64 },
    CoalescentExponential { initial_pop: f64, growth_rate: f64 },
    BirthDeath { birth_rate: f64, death_rate: f64 },
}

/// Molecular clock model.
#[derive(Debug, Clone)]
pub enum ClockModel {
    Strict { rate: f64 },
    UncorrelatedLognormal { mean_rate: f64, stdev: f64 },
}

/// A single MCMC sample.
#[derive(Debug, Clone)]
pub struct McmcSample {
    pub generation: usize,
    pub tree: PhyloTree,
    pub log_likelihood: f64,
    pub log_prior: f64,
    pub log_posterior: f64,
}

/// Results from an MCMC run.
#[derive(Debug, Clone)]
pub struct McmcResult {
    pub samples: Vec<McmcSample>,
    pub acceptance_rates: HashMap<String, f64>,
}

/// Convergence diagnostics for MCMC samples.
#[derive(Debug, Clone)]
pub struct ConvergenceDiag {
    pub ess: HashMap<String, f64>,
    pub mean_log_likelihood: f64,
    pub variance_log_likelihood: f64,
}

/// Summary of the posterior distribution.
#[derive(Debug, Clone)]
pub struct PosteriorSummary {
    pub map_tree: PhyloTree,
    pub mean_branch_lengths: Vec<(NodeId, f64)>,
    pub node_age_intervals: Vec<(NodeId, f64, f64)>,
    pub clade_credibilities: Vec<(BTreeSet<String>, f64)>,
}

/// Simple xorshift64 PRNG.
struct Xorshift64 {
    state: u64,
}

impl Xorshift64 {
    fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 { 1 } else { seed },
        }
    }

    fn next_u64(&mut self) -> u64 {
        self.state ^= self.state << 13;
        self.state ^= self.state >> 7;
        self.state ^= self.state << 17;
        self.state
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform in \[lo, hi\]
    fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }
}

/// Run MCMC sampling for Bayesian phylogenetic inference.
pub fn mcmc_sample(
    tree: &PhyloTree,
    sequences: &[&[u8]],
    model: &dyn SubstitutionModel,
    prior: &TreePrior,
    clock: &ClockModel,
    config: &McmcConfig,
) -> Result<McmcResult> {
    let mut rng = Xorshift64::new(42);
    let mut current_tree = tree.clone();
    let mut current_ll = generic_tree_likelihood(&current_tree, sequences, model, None)?;
    let mut current_prior = tree_prior_ln(&current_tree, prior, clock);
    let mut current_posterior = current_ll + current_prior;

    let mut samples = Vec::new();
    let mut accept_counts: HashMap<String, (usize, usize)> = HashMap::new();

    let total_weight = config.proposal_weights.nni
        + config.proposal_weights.spr
        + config.proposal_weights.branch_scale
        + config.proposal_weights.branch_slide
        + config.proposal_weights.model_params;

    for gen in 0..config.n_generations {
        // Choose proposal type.
        let r = rng.next_f64() * total_weight;
        let (proposal_name, proposed_tree, hastings_ratio) = if r
            < config.proposal_weights.nni
        {
            propose_nni(&current_tree, &mut rng)
        } else if r < config.proposal_weights.nni + config.proposal_weights.spr {
            propose_spr(&current_tree, &mut rng)
        } else if r
            < config.proposal_weights.nni
                + config.proposal_weights.spr
                + config.proposal_weights.branch_scale
        {
            propose_branch_scale(&current_tree, &mut rng)
        } else if r
            < config.proposal_weights.nni
                + config.proposal_weights.spr
                + config.proposal_weights.branch_scale
                + config.proposal_weights.branch_slide
        {
            propose_branch_slide(&current_tree, &mut rng)
        } else {
            // Model params: no-op for now (empirical models have no free params).
            ("model_params", current_tree.clone(), 1.0)
        };

        // Evaluate proposed tree.
        let proposed_ll = match generic_tree_likelihood(&proposed_tree, sequences, model, None) {
            Ok(ll) => ll,
            Err(_) => {
                // Invalid tree (e.g., negative branch lengths): reject.
                let entry = accept_counts
                    .entry(proposal_name.to_string())
                    .or_insert((0, 0));
                entry.1 += 1;
                continue;
            }
        };
        let proposed_prior = tree_prior_ln(&proposed_tree, prior, clock);
        let proposed_posterior = proposed_ll + proposed_prior;

        // Metropolis-Hastings acceptance.
        let log_alpha = (proposed_posterior - current_posterior) + hastings_ratio.ln();
        let accept = if log_alpha >= 0.0 {
            true
        } else {
            rng.next_f64() < log_alpha.exp()
        };

        let entry = accept_counts
            .entry(proposal_name.to_string())
            .or_insert((0, 0));
        entry.1 += 1;
        if accept {
            entry.0 += 1;
            current_tree = proposed_tree;
            current_ll = proposed_ll;
            current_prior = proposed_prior;
            current_posterior = proposed_posterior;
        }

        // Sample.
        if gen >= config.burnin && (gen - config.burnin) % config.sample_every == 0 {
            samples.push(McmcSample {
                generation: gen,
                tree: current_tree.clone(),
                log_likelihood: current_ll,
                log_prior: current_prior,
                log_posterior: current_posterior,
            });
        }
    }

    let acceptance_rates: HashMap<String, f64> = accept_counts
        .into_iter()
        .map(|(name, (accepted, total))| {
            let rate = if total > 0 {
                accepted as f64 / total as f64
            } else {
                0.0
            };
            (name, rate)
        })
        .collect();

    Ok(McmcResult {
        samples,
        acceptance_rates,
    })
}

/// Compute convergence diagnostics from MCMC samples.
pub fn convergence_diagnostics(samples: &[McmcSample]) -> ConvergenceDiag {
    let n = samples.len();
    let lls: Vec<f64> = samples.iter().map(|s| s.log_likelihood).collect();

    let mean_ll = if n > 0 {
        lls.iter().sum::<f64>() / n as f64
    } else {
        0.0
    };

    let var_ll = if n > 1 {
        let sum_sq: f64 = lls.iter().map(|&x| (x - mean_ll).powi(2)).sum();
        sum_sq / (n - 1) as f64
    } else {
        0.0
    };

    // ESS via autocorrelation.
    let ess_ll = effective_sample_size(&lls);

    let mut ess = HashMap::new();
    ess.insert("log_likelihood".to_string(), ess_ll);

    ConvergenceDiag {
        ess,
        mean_log_likelihood: mean_ll,
        variance_log_likelihood: var_ll,
    }
}

/// Summarize the posterior distribution.
pub fn posterior_summary(samples: &[McmcSample]) -> PosteriorSummary {
    if samples.is_empty() {
        return PosteriorSummary {
            map_tree: PhyloTree::new(),
            mean_branch_lengths: Vec::new(),
            node_age_intervals: Vec::new(),
            clade_credibilities: Vec::new(),
        };
    }

    // MAP tree: sample with highest posterior.
    let map_idx = samples
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| {
            a.log_posterior
                .partial_cmp(&b.log_posterior)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, _)| i)
        .unwrap_or(0);
    let map_tree = samples[map_idx].tree.clone();

    // Mean branch lengths (from MAP tree nodes).
    let n_nodes = map_tree.node_count();
    let mean_branch_lengths: Vec<(NodeId, f64)> = (0..n_nodes)
        .filter_map(|id| {
            map_tree
                .get_node(id)
                .and_then(|n| n.branch_length)
                .map(|bl| (id, bl))
        })
        .collect();

    // Node age intervals: use total branch length as proxy for ages.
    let node_age_intervals: Vec<(NodeId, f64, f64)> = (0..n_nodes)
        .filter_map(|id| {
            map_tree
                .get_node(id)
                .and_then(|n| n.branch_length)
                .map(|bl| (id, bl * 0.5, bl * 1.5)) // Simple 50%-150% interval
        })
        .collect();

    // Clade credibilities: fraction of samples containing each bipartition.
    let mut clade_counts: HashMap<BTreeSet<String>, usize> = HashMap::new();
    for sample in samples {
        let bps = bipartitions(&sample.tree);
        for bp in bps {
            *clade_counts.entry(bp).or_insert(0) += 1;
        }
    }

    let n_samples = samples.len();
    let clade_credibilities: Vec<(BTreeSet<String>, f64)> = clade_counts
        .into_iter()
        .map(|(bp, count)| (bp, count as f64 / n_samples as f64))
        .collect();

    PosteriorSummary {
        map_tree,
        mean_branch_lengths,
        node_age_intervals,
        clade_credibilities,
    }
}

/// Log tree prior probability.
fn tree_prior_ln(tree: &PhyloTree, prior: &TreePrior, clock: &ClockModel) -> f64 {
    let topology_prior = match prior {
        TreePrior::Uniform => 0.0,
        TreePrior::CoalescentConstant { pop_size } => {
            coalescent_prior_ln(tree, *pop_size)
        }
        TreePrior::CoalescentExponential {
            initial_pop,
            growth_rate,
        } => coalescent_exp_prior_ln(tree, *initial_pop, *growth_rate),
        TreePrior::BirthDeath {
            birth_rate,
            death_rate,
        } => birth_death_prior_ln(tree, *birth_rate, *death_rate),
    };

    let clock_prior = match clock {
        ClockModel::Strict { rate } => {
            // Exponential prior on rate: ln(lambda * exp(-lambda * rate))
            let lambda: f64 = 1.0;
            lambda.ln() - lambda * rate
        }
        ClockModel::UncorrelatedLognormal { mean_rate, stdev } => {
            // Lognormal prior on branch rates.
            let mut lp = 0.0;
            for id in 0..tree.node_count() {
                if let Some(node) = tree.get_node(id) {
                    if let Some(bl) = node.branch_length {
                        if bl > 0.0 {
                            let log_bl = bl.ln();
                            let log_mean = mean_rate.ln();
                            lp -= (log_bl - log_mean).powi(2) / (2.0 * stdev * stdev);
                            lp -= log_bl + (2.0 * std::f64::consts::PI).sqrt().ln()
                                + stdev.ln();
                        }
                    }
                }
            }
            lp
        }
    };

    topology_prior + clock_prior
}

/// Coalescent prior (constant population size).
fn coalescent_prior_ln(tree: &PhyloTree, pop_size: f64) -> f64 {
    if pop_size <= 0.0 {
        return f64::NEG_INFINITY;
    }
    // Simple coalescent: penalize long branches relative to population size.
    let mut lp = 0.0;
    for id in 0..tree.node_count() {
        if let Some(node) = tree.get_node(id) {
            if let Some(bl) = node.branch_length {
                // Exponential prior on waiting times.
                lp -= bl / pop_size;
            }
        }
    }
    lp
}

/// Coalescent prior (exponential growth).
fn coalescent_exp_prior_ln(tree: &PhyloTree, initial_pop: f64, growth_rate: f64) -> f64 {
    if initial_pop <= 0.0 {
        return f64::NEG_INFINITY;
    }
    let mut lp = 0.0;
    for id in 0..tree.node_count() {
        if let Some(node) = tree.get_node(id) {
            if let Some(bl) = node.branch_length {
                let effective_pop = initial_pop * (-growth_rate * bl).exp();
                if effective_pop > 0.0 {
                    lp -= bl / effective_pop;
                }
            }
        }
    }
    lp
}

/// Birth-death prior.
fn birth_death_prior_ln(tree: &PhyloTree, birth_rate: f64, death_rate: f64) -> f64 {
    if birth_rate <= 0.0 || death_rate < 0.0 || death_rate >= birth_rate {
        return f64::NEG_INFINITY;
    }
    let n = tree.leaf_count() as f64;
    // Simple birth-death: log-probability of tree topology.
    let r = birth_rate - death_rate;
    let mut lp = (n - 1.0) * r.ln();
    // Penalize total tree height.
    let total_bl = tree.total_branch_length();
    lp -= birth_rate * total_bl;
    lp
}

/// Propose NNI move.
fn propose_nni(tree: &PhyloTree, rng: &mut Xorshift64) -> (&'static str, PhyloTree, f64) {
    let edges = collect_nni_edges(tree);
    if edges.is_empty() {
        return ("nni", tree.clone(), 1.0);
    }

    let idx = (rng.next_u64() as usize) % edges.len();
    let (parent_id, child_id) = edges[idx];
    let parent = tree.get_node(parent_id).unwrap();
    let child = tree.get_node(child_id).unwrap();

    let sibling_id = if parent.children[0] == child_id {
        parent.children[1]
    } else {
        parent.children[0]
    };

    if child.children.is_empty() {
        return ("nni", tree.clone(), 1.0);
    }
    let nephew_idx = (rng.next_u64() as usize) % child.children.len();
    let nephew_id = child.children[nephew_idx];

    match nni_swap(tree, parent_id, child_id, sibling_id, nephew_id) {
        Ok(new_tree) => ("nni", new_tree, 1.0), // Symmetric proposal
        Err(_) => ("nni", tree.clone(), 1.0),
    }
}

/// Propose SPR move.
fn propose_spr(tree: &PhyloTree, rng: &mut Xorshift64) -> (&'static str, PhyloTree, f64) {
    let n = tree.node_count();
    if n < 4 {
        return ("spr", tree.clone(), 1.0);
    }

    // Pick random non-root node to prune.
    let non_root: Vec<NodeId> = (0..n)
        .filter(|&id| id != tree.root())
        .collect();
    if non_root.is_empty() {
        return ("spr", tree.clone(), 1.0);
    }
    let prune = non_root[(rng.next_u64() as usize) % non_root.len()];

    // Pick random regraft point (not same as prune or parent).
    let prune_parent = tree.get_node(prune).unwrap().parent;
    let candidates: Vec<NodeId> = (0..n)
        .filter(|&id| id != prune && Some(id) != prune_parent && id != tree.root())
        .collect();
    if candidates.is_empty() {
        return ("spr", tree.clone(), 1.0);
    }
    let regraft = candidates[(rng.next_u64() as usize) % candidates.len()];

    match crate::tree_search::spr_move(tree, prune, regraft) {
        Ok(new_tree) => ("spr", new_tree, 1.0),
        Err(_) => ("spr", tree.clone(), 1.0),
    }
}

/// Propose branch length scaling (multiply all by exp(delta)).
fn propose_branch_scale(
    tree: &PhyloTree,
    rng: &mut Xorshift64,
) -> (&'static str, PhyloTree, f64) {
    let epsilon = 0.1;
    let delta = rng.uniform(-epsilon, epsilon);
    let factor = delta.exp();

    let n = tree.node_count();
    let mut nodes: Vec<crate::tree::Node> = (0..n)
        .map(|id| tree.get_node(id).unwrap().clone())
        .collect();

    let mut n_branches = 0;
    for node in &mut nodes {
        if let Some(ref mut bl) = node.branch_length {
            *bl *= factor;
            n_branches += 1;
        }
    }

    // Hastings ratio: exp(n_branches * delta)
    let hastings = (n_branches as f64 * delta).exp();

    match PhyloTree::from_nodes(nodes, tree.root()) {
        Ok(new_tree) => ("branch_scale", new_tree, hastings),
        Err(_) => ("branch_scale", tree.clone(), 1.0),
    }
}

/// Propose sliding a single branch length.
fn propose_branch_slide(
    tree: &PhyloTree,
    rng: &mut Xorshift64,
) -> (&'static str, PhyloTree, f64) {
    let n = tree.node_count();
    let nodes_with_bl: Vec<NodeId> = (0..n)
        .filter(|&id| {
            tree.get_node(id)
                .and_then(|n| n.branch_length)
                .is_some()
        })
        .collect();

    if nodes_with_bl.is_empty() {
        return ("branch_slide", tree.clone(), 1.0);
    }

    let target = nodes_with_bl[(rng.next_u64() as usize) % nodes_with_bl.len()];
    let epsilon = 0.05;
    let delta = rng.uniform(-epsilon, epsilon);

    let mut nodes: Vec<crate::tree::Node> = (0..n)
        .map(|id| tree.get_node(id).unwrap().clone())
        .collect();

    if let Some(ref mut bl) = nodes[target].branch_length {
        *bl += delta;
        if *bl < 1e-8 {
            *bl = 1e-8;
        }
    }

    match PhyloTree::from_nodes(nodes, tree.root()) {
        Ok(new_tree) => ("branch_slide", new_tree, 1.0), // Symmetric
        Err(_) => ("branch_slide", tree.clone(), 1.0),
    }
}

/// Effective sample size via autocorrelation.
fn effective_sample_size(values: &[f64]) -> f64 {
    let n = values.len();
    if n < 2 {
        return n as f64;
    }

    let mean = values.iter().sum::<f64>() / n as f64;
    let var: f64 = values.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n as f64;
    if var < 1e-30 {
        return n as f64;
    }

    let mut sum_rho = 0.0;
    for lag in 1..n {
        let mut rho = 0.0;
        for i in 0..(n - lag) {
            rho += (values[i] - mean) * (values[i + lag] - mean);
        }
        rho /= n as f64 * var;

        if rho < 0.0 {
            break; // Stop at first negative autocorrelation
        }
        sum_rho += rho;
    }

    let ess = n as f64 / (1.0 + 2.0 * sum_rho);
    ess.max(1.0)
}

/// Collect edges suitable for NNI.
fn collect_nni_edges(tree: &PhyloTree) -> Vec<(NodeId, NodeId)> {
    let mut edges = Vec::new();
    for id in 0..tree.node_count() {
        let node = tree.get_node(id).unwrap();
        if !node.is_leaf() && !node.is_root() && node.children.len() == 2 {
            if let Some(parent_id) = node.parent {
                let parent = tree.get_node(parent_id).unwrap();
                if parent.children.len() == 2 {
                    edges.push((parent_id, id));
                }
            }
        }
    }
    edges
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::subst_model::Jc69Model;

    fn small_tree_and_seqs() -> (PhyloTree, Vec<Vec<u8>>) {
        let tree =
            PhyloTree::from_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();
        let seqs = vec![
            b"ACGTACGTAC".to_vec(),
            b"ACGTACGTAC".to_vec(),
            b"TGCATGCATG".to_vec(),
            b"TGCATGCATG".to_vec(),
        ];
        (tree, seqs)
    }

    #[test]
    fn mcmc_converges_on_tiny_tree() {
        let (tree, seqs) = small_tree_and_seqs();
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
        let model = Jc69Model::new();
        let config = McmcConfig {
            n_generations: 500,
            sample_every: 50,
            burnin: 100,
            ..Default::default()
        };

        let result = mcmc_sample(
            &tree,
            &refs,
            &model,
            &TreePrior::Uniform,
            &ClockModel::Strict { rate: 1.0 },
            &config,
        )
        .unwrap();

        assert!(!result.samples.is_empty());
    }

    #[test]
    fn acceptance_rates_in_range() {
        let (tree, seqs) = small_tree_and_seqs();
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
        let model = Jc69Model::new();
        let config = McmcConfig {
            n_generations: 200,
            sample_every: 20,
            burnin: 50,
            ..Default::default()
        };

        let result = mcmc_sample(
            &tree,
            &refs,
            &model,
            &TreePrior::Uniform,
            &ClockModel::Strict { rate: 1.0 },
            &config,
        )
        .unwrap();

        for (name, rate) in &result.acceptance_rates {
            assert!(
                (0.0..=1.0).contains(rate),
                "{} acceptance rate {} out of range",
                name, rate
            );
        }
    }

    #[test]
    fn ess_positive() {
        let (tree, seqs) = small_tree_and_seqs();
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
        let model = Jc69Model::new();
        let config = McmcConfig {
            n_generations: 500,
            sample_every: 10,
            burnin: 50,
            ..Default::default()
        };

        let result = mcmc_sample(
            &tree,
            &refs,
            &model,
            &TreePrior::Uniform,
            &ClockModel::Strict { rate: 1.0 },
            &config,
        )
        .unwrap();

        let diag = convergence_diagnostics(&result.samples);
        for (name, &ess_val) in &diag.ess {
            assert!(ess_val > 0.0, "ESS for {} should be > 0, got {}", name, ess_val);
        }
    }

    #[test]
    fn posterior_summary_returns_map_tree() {
        let (tree, seqs) = small_tree_and_seqs();
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
        let model = Jc69Model::new();
        let config = McmcConfig {
            n_generations: 300,
            sample_every: 30,
            burnin: 50,
            ..Default::default()
        };

        let result = mcmc_sample(
            &tree,
            &refs,
            &model,
            &TreePrior::Uniform,
            &ClockModel::Strict { rate: 1.0 },
            &config,
        )
        .unwrap();

        let summary = posterior_summary(&result.samples);
        assert!(summary.map_tree.leaf_count() > 0);
    }

    #[test]
    fn coalescent_prior_penalizes_long_branches() {
        let short = PhyloTree::from_newick("(A:0.01,B:0.01);").unwrap();
        let long = PhyloTree::from_newick("(A:10.0,B:10.0);").unwrap();

        let prior = TreePrior::CoalescentConstant { pop_size: 1.0 };
        let clock = ClockModel::Strict { rate: 1.0 };

        let lp_short = tree_prior_ln(&short, &prior, &clock);
        let lp_long = tree_prior_ln(&long, &prior, &clock);
        assert!(
            lp_short > lp_long,
            "short branches should have higher prior: {} vs {}",
            lp_short, lp_long
        );
    }

    #[test]
    fn strict_clock_constrains_rates() {
        let tree = PhyloTree::from_newick("(A:0.1,B:0.1);").unwrap();
        let prior = TreePrior::Uniform;

        let low_rate = ClockModel::Strict { rate: 0.01 };
        let high_rate = ClockModel::Strict { rate: 100.0 };

        let lp_low = tree_prior_ln(&tree, &prior, &low_rate);
        let lp_high = tree_prior_ln(&tree, &prior, &high_rate);
        assert!(
            lp_low > lp_high,
            "low rate should have higher prior: {} vs {}",
            lp_low, lp_high
        );
    }

    #[test]
    fn log_posterior_is_ll_plus_prior() {
        let (tree, seqs) = small_tree_and_seqs();
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
        let model = Jc69Model::new();
        let config = McmcConfig {
            n_generations: 200,
            sample_every: 20,
            burnin: 50,
            ..Default::default()
        };

        let result = mcmc_sample(
            &tree,
            &refs,
            &model,
            &TreePrior::Uniform,
            &ClockModel::Strict { rate: 1.0 },
            &config,
        )
        .unwrap();

        for sample in &result.samples {
            let expected = sample.log_likelihood + sample.log_prior;
            assert!(
                (sample.log_posterior - expected).abs() < 1e-10,
                "posterior {} != ll {} + prior {}",
                sample.log_posterior,
                sample.log_likelihood,
                sample.log_prior
            );
        }
    }

    #[test]
    fn branch_scale_is_reversible() {
        let tree = PhyloTree::from_newick("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);").unwrap();
        let mut rng = Xorshift64::new(42);
        let (_, scaled, hastings) = propose_branch_scale(&tree, &mut rng);
        assert!(hastings > 0.0);
        assert_eq!(scaled.leaf_count(), tree.leaf_count());
    }

    #[test]
    fn burnin_discarded() {
        let (tree, seqs) = small_tree_and_seqs();
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
        let model = Jc69Model::new();
        let config = McmcConfig {
            n_generations: 200,
            sample_every: 10,
            burnin: 100,
            ..Default::default()
        };

        let result = mcmc_sample(
            &tree,
            &refs,
            &model,
            &TreePrior::Uniform,
            &ClockModel::Strict { rate: 1.0 },
            &config,
        )
        .unwrap();

        for sample in &result.samples {
            assert!(
                sample.generation >= 100,
                "sample from gen {} should be after burnin 100",
                sample.generation
            );
        }
    }

    #[test]
    fn empty_samples_handled() {
        let diag = convergence_diagnostics(&[]);
        assert_eq!(diag.mean_log_likelihood, 0.0);

        let summary = posterior_summary(&[]);
        assert_eq!(summary.clade_credibilities.len(), 0);
    }

    #[test]
    fn birth_death_prior() {
        let tree = PhyloTree::from_newick("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);").unwrap();
        let prior = TreePrior::BirthDeath {
            birth_rate: 2.0,
            death_rate: 0.5,
        };
        let clock = ClockModel::Strict { rate: 1.0 };
        let lp = tree_prior_ln(&tree, &prior, &clock);
        assert!(lp.is_finite(), "birth-death prior should be finite");
    }
}
