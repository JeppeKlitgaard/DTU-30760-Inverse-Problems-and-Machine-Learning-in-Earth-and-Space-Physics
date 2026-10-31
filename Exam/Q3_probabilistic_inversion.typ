#import "common/common.typ": *

#show: setup.with(3, "Probabilistic Inversion")
#show: show-page1

#let questions = [
= Questions
- What is *Bayes theorem*?
- *How can it be used* to solve an inverse problem?
- What are *MCMC methods*?
- Briefly *give an example* of their use in finding the probabilistic solution to an inverse problem in Earth and Space Physics.
]
#questions

= Approach
- Bayes' theorem: $σ_m (vv(m)) = k L(vv(m)) ρ_m (vv(m))$
  - $k = 1\/p(vv(d))$ is scaling constant (evidence in terms of standard notation)
- Intuition around MCMC and Bayes
  - Consider MCMC determination of $π$, similar here:
  - Assume how often we end up in some volume of parameter space is proportional to the probability of that volume (Lebesgue _/le-beg/_ measure)
  - Thus $σ_(m, i) = (L_i ρ_(m, i))/(∑ L_i ρ_(m, i))$
- Example: Assignment 2
- Benefit: Easily handles non-linearity
  - Can do so efficiently with right prior

= Probability Basics
- Expectation: $𝔼[f(x)] = ∫ f(x) p(x) dif x$
- Product rule: $p(x, y) = p(x|y) p(y)$
- Sum rule (marginalisation): $p(x) = ∫ p(x, y) dif y$
- Independence: $p(x, y) = p(x) p(y)$
- Bayes' theorem: $p(x|y) = (p(y|x) p(x)) / p(y)$

= Setup
- Forward model: $vv(d) = g(vv(m)) + vv(e)$
- Joint prior: $ρ(vv(d), vv(m)) = ρ_d (vv(d)) ρ_m (vv(m))$
- Distribution linking data and model: $θ(vv(d), vv(m)) = θ(vv(d)|vv(m)) μ_m (vv(m))$
  - $μ_m (vv(m))$ is _homogenous probability density function_, ensures coordinate transforms preserve probability volume
- Joint posterior: $σ(vv(d), vv(m)) = k (ρ(vv(d), vv(m)) θ(vv(d), vv(m))) / μ(vv(d), vv(m)) = k L(vv(m)) ρ_m (vv(m))$
  - $L(vv(m)) = ∫ (ρ_d (vv(d)) θ(vv(d)|vv(m))) / (μ_d (vv(d))) dif vv(d)$ is _likelihood_, $k$ is normalisation constant

#colbreak(weak: true)
= Monte Carlo Sampling
- Markov property: $p(x_n|x_(n-1)) = p(x_n|x_(n-1), x_(n-2), ...)$
- $p_i = p(vv(m)_i) Δ m^(1) Δ m^(2) ...$ where $Δ m^(j)$ are discretisation lengths in parameter space

== MCMC Methods
- Converge to target distribution as $N→∞$ (ergodicity)
- Algorithms must satisfy _detailed balance_: $p_i P(i→j) = p_j P(j→i)$
- Does not seek to produce _independent_ samples, instead seeks to do _random walk_ through parameter space

== Extended Metropolis Hastings
+ Initialise: $vv(m)_0 ~ ρ_m (vv(m))$
+ Propose: $vv(m)' ~ q(vv(m)'|vv(m)_n)$
  - $q$ is _proposal distribution_; easy to sample from
+ Accept with probability: $A = min(1, L(vv(m)')\/L(vv(m)_0))$
+ Repeat from (2.)

== Sampling Prior
- For proposal distribution, we can:
  - Gibbs Sequential: Update random subset of parameters with conditional distribution.
  - Pre-Conditioned Crank-Nicolson: Propose new sample as linear combination of current sample and random draw from prior.

= Encoding Prior Information
- Can encode prior information through choice of prior distribution
- 1-point statistics: Gaussian with diagonal covariance
- 2-point statistics: Gaussian with non-diagonal covariance
  - Gaussian Process and covariance functions can be used to encode assumptions about smoothness and correlation structure of model parameters
  - Can encode some spatial correlation
- Higher-order statistics: could encode more complex structures

#v(1fr)
#figure(image("inversion_test_synthetic_100.png"))

#show: show-page2
#questions
#LINE

#unwritten[
  Bayes Theorem
  - TM version
  - $k$: normalization, inverse of evidence
  - Likelihood: $L(vv(m)) = ∫ (ρ_d (vv(d)) θ(vv(d)|vv(m))) / (μ_d (vv(d))) dif vv(d)$
][$
  σ_m (vv(m)) = k L(vv(m)) ρ_m (vv(m))
$]

#unwritten[
Compare to prev:
- Regularisation through prior distribution
- Data fit through likelihood
- Posterior is solution

Benefit: Easily handles non-linearity, gives full distribution
]

#unwritten[
Forward model:
][$vv(d) = g(vv(m)) + vv(e)$]

#unwritten[
= MCMC
- Lebesgue: Probability measure must $π$ must be _absolutely continuous_ w.r.t Lebesgue measure $λ$.
- Or: how often we end up in a volume of parameter space is proportional to the probability of that volume
- This is also what implies the _homogenous probability density functions_.
][Monte Carlo intuition: estimation of $π$]

#unwritten[
- Markov property: $p(x_n|x_(n-1)) = p(x_n|x_(n-1), x_(n-2), ...)$
  - Next state depends only on current state, not history
- Detailed balance: $p_i P(i→j) = p_j P(j→i)$
- Ergodicity: Converge to target distribution as $N→∞$
][Markov chains: do away with independent sampling, do random walk at detailed balance]


#unwritten[
== Extended Metropolis Hastings
]
+ Initialise: $vv(m)_0 ~ ρ_m (vv(m))$
+ Propose: $vv(m)' ~ q(vv(m)'|vv(m)_n)$
  - $q$ is _proposal distribution_; easy to sample from
+ Accept with probability: $A = min(1, L(vv(m)')\/L(vv(m)_0))$
+ Repeat from (2.)

Extended vs. Standard:
- Only likelihoods enter acceptance probability, not priors
- Requires that proposal is symmetric
- Requires that proposal is _ergodic_ w.r.t prior (samples all of prior)

#unwritten[== Sampling Prior]
- Gibbs: Update random subset of parameters with conditional distribution.
- Pre-Conditioned Crank-Nicolson: Propose new sample as linear combination of current sample and random draw from prior.
  - $vv(m)' = sqrt(1-β^2) vv(m)_n + β vv(z)$
  - $vv(z) ~ ρ_m (vv(m))$ is random draw from prior
  - $β ∈(0, 1]$ is step size

#unwritten[= Encoding Prior Information]
- Encode prior info in prior dist.
- 1-point stats: Gaussian with diagonal covariance
- 2-point stats: Gaussian with non-diagonal covariance
  - Spatial correlation

#v(1fr)
#figure(image("inversion_test_synthetic_100.png"))