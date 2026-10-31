#import "common/common.typ": *

#show: setup.with(1, "Robust Estimation")
#show: show-page1

#let questions = [
= Questions
- Describe *how one can make robust estimates* of model parameters from observed data in the presence of non-Gaussian noise.
- What are the *assumptions* underlying this approach?
- Under what *circumstances* would you recommend to use robust estimation?
]

#questions


= Approach
- Draw distributions (Gaussian, heavy-tailed)
- Keep track on assumptions on board
- Summarise with example (1st assignment?)

= Setup
- Data/observations: $vv(d)$
- Model parameters: $vv(m)$
- Likelihood: $p(vv(d)|vv(m))$
- Forward model $g(vv(m)) = hat(vv(d))$

= Maximum Likelihood Estimation (MLE)
- Choose MLE as estimator: $hat(vv(m)) = arg max_vv(m) p(vv(d)|vv(m))$
- Equivalently: $hat(vv(m)) = arg min_vv(m) -log p(vv(d)|vv(m))$
  - (log-transform preserves ranking since log is monotonic)
- For IID data: $hat(vv(m)) = arg min_vv(m) - sum_i log p(d_i|vv(m))$
  - Since $p(vv(d)|vv(m)) = ∏_i p(d_i|vv(m))$
- Let loss function be NLL: $ϕ(vv(m)) ≔ -log p(vv(d)|vv(m)) = - sum_i log p(d_i|vv(m))$
- Assume $ϕ$ is differentiable, convex, (lower-bounded, coercive)
  - Then _unique, global_ minimum is at $∇_vv(m) ϕ(vv(hat(m))) = 0$
  - $vv(hat(m))$ is _MLE_
- Estimators, $E{hat(vv(m))}$, desirable properties are:
  - Unbiased: $E{hat(vv(m))} = vv(m)$
  - Consistent: $hat(vv(m)) → vv(m)$ as number of data points $N → ∞$
  - Efficient: Has lowest possible variance among all unbiased estimators
- Summary of assumptions:
  - IID data
  - Convex, differentiable loss function (e.g. negative log-likelihood)

= Robust Estimation
- Motivation: We need loss function $ϕ(vv(m))$ faithfully represent underlying data distribution
- Gaussian likelihood assumption heavily penalises outliers ⇒ wrong model if outliers present
- We restrict variational family of likelihoods to be those characterised by:
  - Location parameter, $vv(μ)$
  - Scale parameter, $σ$
  - These are most common, continuous distributions (Guassian, Laplace, ...)
- Suppose _a priori_ known scale parameter, $σ$
- We can express location as a _(normalised) residual_, $vv(ε) = (vv(d) - μ)\/σ$
- Then we can express loss as a function of residuals: $ϕ(vv(ε))$
- Under prev. assumptions: $∑_i ϕ'(vv(ε_i)) = 0$ gives optimal $hat(μ)$
  - Influence function $ψ(ε_i) ≔ dv(ϕ(ε_i), ε_i)$
  - Optimality condition: $∑_i ψ(ε_i) = 0$

== Linearisation
- $ψ(ε_i)$ are generally _non-linear_.
- Linearise: $ψ(ε_i) = w_i ε_i$
- $w_i$ is _robust weight_ for data point
- Determined using _Iteratively Reweighted Least Squares_
- Optimality:
$
  ∑_i ψ(ε_i) &= ∑_i w_i ε_i = ∑_i w_i (d_i - hat(μ))cancel(\/σ) = 0\
  hat(μ) &= (∑_i w_i d_i)/(∑_i w_i)
$

== IRLS
-  Initialise $vv(w) = vv(1)$
+ Estimate $hat(vv(m)) = g^(-1)(hat(vv(d)))$
+ Calculate residuals $vv(ε) = (vv(d) - g(hat(vv(m))))\/σ$
+ Update weights $vv(w) = ψ(vv(ε))\/vv(ε)$
+ Repeat until convergence

= Assumption summary
- IID data (factorisable likelihood)
- Likelihood characterised by location and scale parameters
- (A priori scale parameter)
- Differentiable, convex loss function (also coercive, lower-bounded)
- (Linearisation of influence function)

= Perspectives
- Suitable for data with outliers or non-Gaussian noise
- Full Bayesian approach instead of point estimation with MLE

= Distributions
#table(
  columns: 4,
  [*Dist.*], [*Loss*, $φ(vv(ε))$], [*Influence*, $ψ(vv(ε))$], [*Weight*, $w(ε)$],

  [*Gauss.*],
  $ε^2$,
  $ε$,
  $1$,

  [*Lapla.*],
  $|ε|$,
  $"sgn"(ε)$,
  $1\/|ε|$,

  [*Huber*],
  $∝ cases(
    c abs(ε) &quad ε > c,
    ε^2 &quad abs(ε) ≤ c,
    c abs(ε) &quad ε < c,
  )$,
  $∝ cases(
    c &quad ε > c,
    ε &quad abs(ε) ≤ c,
    -c &quad ε < -c,
  )$,
  $∝ cases(
    c\/|ε| &quad ε > c,
    1 &quad abs(ε) ≤ c,
    c\/|ε| &quad ε < -c,
  )$,
)

= Linear Algebra
$
  hat(vv(m)) = (mm(G)^T mm(W) mm(G))^(-1) mm(G)^T mm(W) vv(d)
$
Where:
- $mm(G)$ is forward model
- $mm(W)$ is diagonal matrix of robust weights
  - $mm(W) = mm(I)$ for non-robust case

// ---
#show: show-page2
#questions
#LINE

#unwritten[
Robust Estimation:
- Estimator: given data, give estimate of parameter
- Robust: changes in individual data points have limited influence on estimate
- Important when outliers are present or (equiv) non-Gaussian noise
]

#unwritten[Setup: ][$
vv(d), quad vv(m), quad p(vv(d)|vv(m))
$]

#unwritten[
  #box(width: 100%)[
  #fill-in-drawing(16em, 10em, caption: [Gaussian + Heavy-Tailed])
  ]
]

#unwritten[MLE: ][$
hat(vv(m)) = arg max_vv(m) p(vv(d)|vv(m))
$]

#unwritten[Introduce Loss: ][$
ϕ(vv(m)) ≔ -log p(vv(d)|vv(m)) =^"IID" -∑_i log p(d_i|vv(m))
$]

#unwritten[Optimality: ][$
∇_vv(m) ϕ(vv(hat(m))) = 0
$]

#unwritten[
Restrict to location-scale family:\
Note: Heteroscedasticity makes $vv(σ)$ scalar
][$
p(vv(d)|vv(μ), vv(σ))
$]

#unwritten[Define _normalized_ residuals: ][$
vv(ε) = (vv(d) - vv(μ))\/vv(σ)
$]

#unwritten[
Define influence function:\
Note: non-linear!][$
ψ(ε_i) ≔ ϕ'(ε_i)
$]
#unwritten[Optimality: ][$
∑_i ψ(ε_i) = 0
$]

#unwritten[= IRLS]

#unwritten[Assume linear influence: ][$
ψ(ε_i) = w_i ε_i
$]
#unwritten[Optimality: ][$
hat(μ) = (∑_i w_i d_i)/(∑_i w_i)
$]

#unwritten[IRLS: ][
Initialise $vv(w) = vv(1)$\
Estimate $hat(vv(m)) = g^(-1)(hat(vv(d)))$\
Calculate residuals $vv(ε) = (vv(d) - g(hat(vv(m))))\/σ$\
Update weights $vv(w) = ψ(vv(ε))\/vv(ε)$\
]

#unwritten[Huber: ][$
ψ(ε) ∝ cases(c &quad ε > c, ε &quad abs(ε) ≤ c, -c &quad ε < -c,)
$]

#colbreak(weak: true)
#unwritten[
= Assumptions
- IID data (factorisable likelihood)
- Likelihood characterised by location and scale parameters
- (A priori scale parameter)
- Differentiable, convex loss function (also coercive, lower-bounded)
- (Linearisation of influence function)
]

#unwritten[
= Perspectives
- Suitable for data with outliers or non-Gaussian noise
- Perhps better yet: Full Bayesian approach instead of point estimation with MLE
]
