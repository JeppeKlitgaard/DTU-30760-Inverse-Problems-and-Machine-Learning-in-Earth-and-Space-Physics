#import "common/common.typ": *

#show: setup.with(4, "Unsupervised Machine Learning")
#show: show-page1

#let questions = [
= Questions
- What is meant by *machine learning*?
  - _Algorithms that can incorporate experience to improve predictive performance._
- What are *supervised and unsupervised machine learning*?
  - Supervised: _Learning from labelled data, e.g. regression and classification._
  - Unsupervised: _Learning from unlabelled data, e.g. clustering and dimensionality reduction._
    - Examples: K-means, GMM, PCA
- Choose an *unsupervised machine learning algorithm* and present this in detail.
- Briefly *describe a possible application* in Earth and Space physics
]
#questions

= Gaussian Mixture Model (GMM)
- $i∈{1, ..., N}$ datapoints
- $j∈{1, ..., M}$ mixture components, $M$ is hyperparameter
- Let likelihood be a mixture of multivariate Gaussians:
  - $p(vv(d)_i|vv(m)) = ∑_(j=1)^M α_j cal(N) (vv(d)_i|mm(μ)_j, mm(Σ)_j)$
  - $α_j ∈ [0, 1]$ are mixture weights, $∑ α_j = 1$
- Assume IID:
  - #inline-block-math(<eq:gmm-loglik>)[$ln p(vv(d)|vv(α), vv(μ), vv(Σ)) = ∑_(i=1)^N ln [∑_(j=1)^M α_j cal(N) (vv(d)_i|mm(μ)_j, mm(Σ)_j)]$]
- Find parameters $vv(m) = (vv(α), vv(μ), vv(Σ))$ through Expectation-Maximisation algorithm:
  + Initialise $vv(α), vv(μ), vv(Σ)$
  + Expectation: evaluate posterior probability: $γ_(i k) = (α_k cal(N) (vv(d)_i|mm(μ)_k, mm(Σ)_k)) / (∑_(j=1)^M α_j cal(N) (vv(d)_i|mm(μ)_j, mm(Σ)_j))$
    - Comes from Bayes theorem: $p(z_i=k|vv(d)_i) = (p(vv(d)_i|z_i=k) p(z_i=k)) / p(vv(d)_i)$
    - Marginal+product: $p(vv(d)_i) = ∑_(j=1)^M p(z_i = j) p(vv(d)_i|z_i = j)$
    - $p(z_i = k) = α_k$
    - $p(vv(d)_i|z_i = k) = cal(N) (vv(d)_i|mm(μ)_k, mm(Σ)_k)$
    - Plug in to get result
  + Maximisation: update parameters
    - $N_k = ∑_(i=1)^N γ_(i k)$
    - $vv(μ)'_k = 1/N_k ∑_(i=1)^N γ_(i k) vv(d)_i$
    - $vv(Σ)'_k = 1/N_k ∑_(i=1)^N γ_(i k) (vv(d)_i - vv(μ)'_k)(vv(d)_i - vv(μ)'_k)^T$
    - $α'_k = N_k / N$
  + Evaluate log-likelihood @eq:gmm-loglik, repeat from (2.) until convergence
- Could pick correct $M$ based on AIC, BIC, or other criterion.

= Application
- Here 2D stellar metallicity $[α\/"Fe"]$ vs. $["Fe"\/"H"]$
  - $α$ is abundance of alpha process elements: O, Ne, Mg, Si, S, Ar, Ca, Ti
- Used for cluster and characterisation of stellar populations in Milky Way
  - Stars can either belong to high-$α$ or low-$α$ disks

#colbreak(weak: true)
= K-means
- Only clustering
- Parittion dataset into $K$ subsets using _centroids_:
  + Initialise $K$ centroids
  + Assign each datapoint to closest centroid
  + Update centroids to be mean of assigned datapoints
  + Repeat from (2.) until convergence
- Solution is not unique or global minimum
  - repeat with different initialisation and pick best solution

= PCA
- Transform data into coordinate system of orthogonal principal components
- Order by variance
- Truncate to reduce dimensionality while retaining most variance
  - Can calculate "explained variance" by truncated basis.
- Dimensionality reduction/compression technique

#v(1fr)
#rotate(
  -90deg,
  box(
    figure(
      image("gmm.png"),
      caption: [
        Left: Data; Middle: GMM, $M=2$; Right: GMM, $M=5$.
      ],
    ),
    width: 14cm,
  ),
  reflow: true,
)
