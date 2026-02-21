# cyanea-ml Architecture

Internal design of the machine learning modules.

## Data Layout

All data is passed as flat `&[f64]` arrays in row-major order with a separate `n_features` parameter. This avoids allocation overhead from nested `Vec<Vec<f64>>` and is compatible with both native and WASM targets. A row (sample) of length `n_features` starts at index `i * n_features`.

## Clustering

### K-Means (Lloyd's Algorithm)

1. Initialize centroids by selecting k random data points (deterministic seed via xorshift64).
2. Assignment step: assign each point to its nearest centroid (Euclidean distance).
3. Update step: recompute centroids as the mean of assigned points.
4. Repeat until convergence (centroid shift < tolerance) or max iterations reached.

Complexity: O(n * k * d * iterations) where n = samples, k = clusters, d = features.

### DBSCAN

Density-based clustering without specifying k:
1. Compute pairwise distances (O(n^2 * d)).
2. For each unvisited point, find its eps-neighborhood.
3. If the neighborhood has >= min_samples points, expand the cluster by recursively visiting all density-reachable points.
4. Points not reachable from any core point are labeled as noise (-1).

### Hierarchical (Agglomerative)

Bottom-up clustering with configurable linkage:
- **Single**: distance = min pairwise distance between clusters
- **Complete**: distance = max pairwise distance
- **Average**: distance = mean pairwise distance
- **Ward**: distance = increase in total within-cluster variance

Uses a O(n^3) naive approach: at each step, find the two closest clusters, merge them, update the distance matrix. The merge history is stored as `MergeStep` records for dendrogram reconstruction.

## Distance Metrics

Ten metrics available via the `DistanceMetric` enum dispatch:
- Euclidean, Manhattan, Cosine, Hamming (core four)
- Additional metrics available through `compute_distance()`

`pairwise_distances()` computes the full n x n distance matrix. With the `parallel` feature, rows are computed in parallel using Rayon.

## Decision Trees and Random Forests

### CART Decision Tree

Each internal node stores a `(feature_index, threshold)` split chosen by minimizing Gini impurity:
- Gini = 1 - sum(p_k^2) where p_k is the proportion of class k
- At each node, iterate over features and candidate thresholds (midpoints between sorted unique values)
- Select the split with lowest weighted Gini of the two children
- Recurse until max_depth or pure leaves

### Random Forest

Ensemble of CART trees with two sources of randomness:
1. **Bootstrap sampling**: each tree is trained on a random sample with replacement (same size as original)
2. **Feature bagging**: at each split, only consider a random subset of sqrt(n_features) features

Prediction: majority vote across all trees. Feature importance: normalized split frequency.

## Gradient Boosted Decision Trees (GBDT)

Additive ensemble of regression trees, fit by gradient descent in function space:

1. Initialize with the mean (regression) or log-odds (classification).
2. For each boosting round:
   - Compute negative gradient of the loss function (residuals for MSE, pseudo-residuals for log-loss)
   - Fit a regression tree to the pseudo-residuals
   - Line search for optimal leaf values
   - Add the tree scaled by the learning rate
3. Early stopping: hold out a validation fraction, stop if validation loss does not improve for `early_stopping_rounds`.

Loss functions:
- **Regression**: MSE (squared error)
- **Binary classification**: log-loss (logistic)
- **Multiclass**: cross-entropy with softmax (one tree per class per round)

Internal regression trees use a simplified CART that minimizes MSE at each split.

## t-SNE

Implementation of t-distributed Stochastic Neighbor Embedding:

1. Compute pairwise affinities P_ij using Gaussian kernels with per-point perplexity calibration (binary search for sigma).
2. Symmetrize: P_ij = (P_i|j + P_j|i) / 2n.
3. Initialize embedding randomly or from PCA.
4. Gradient descent on the KL divergence between P (high-D) and Q (low-D, Student-t with 1 df):
   - Early exaggeration: multiply P by 4 for the first 250 iterations
   - Momentum: use 0.5 initially, switch to 0.8 after 250 iterations
   - Adaptive learning rate per point

Complexity: O(n^2) per iteration (pairwise Q computation). No Barnes-Hut approximation in the current implementation.

## UMAP

Uniform Manifold Approximation and Projection:

1. **Simplicial set construction**: For each point, find k nearest neighbors. Compute fuzzy set membership strengths using the local distance scale (rho = distance to nearest neighbor, sigma calibrated via binary search).
2. **Symmetrization**: Combine the directed kNN graph into an undirected fuzzy simplicial set via: w_sym = w_ij + w_ji - w_ij * w_ji.
3. **Initialization**: Spectral embedding of the fuzzy graph Laplacian (or PCA fallback).
4. **SGD optimization**: Minimize the cross-entropy between the high-dimensional fuzzy set and the low-dimensional embedding. Attractive forces for connected pairs, repulsive forces for negative samples (default 5 per positive edge).

The `min_dist` and `spread` parameters control the a/b curve parameters for the low-dimensional similarity kernel.

## HMM

All computations use log-space arithmetic to prevent underflow on long observation sequences.

- **Forward**: alpha[t][s] = log P(o_1..o_t, q_t = s). Computed left-to-right using log-sum-exp.
- **Backward**: beta[t][s] = log P(o_{t+1}..o_T | q_t = s). Computed right-to-left.
- **Viterbi**: delta[t][s] = max log P(o_1..o_t, q_1..q_t | q_t = s). Uses max instead of log-sum-exp, with backtracking for the optimal path.
- **Baum-Welch**: E-step computes posterior state/transition probabilities using forward and backward matrices. M-step updates initial, transition, and emission parameters. Iterates until log-likelihood converges.

## Cross-Validation

All CV methods use a closure `eval_fn(train_indices, test_indices) -> Result<f64>` that the user provides. This avoids coupling to any specific model type.

- **K-fold**: Shuffle indices (Fisher-Yates with seed), split into k contiguous folds.
- **Stratified K-fold**: Group indices by class label, then interleave groups into folds to preserve class proportions.
- **Leave-one-out**: Special case of k-fold with k = n. No shuffling needed.
