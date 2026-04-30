#import "@preview/codly:1.3.0": *
// #import "@preview/cetz:0.4.2"
#import "@preview/unify:0.7.1": num, qty
#import "@preview/algorithmic:0.1.0"
#import "@preview/physica:0.9.8": curl, grad, tensor, pdv, dv, TT

#import "@preview/meander:0.4.1"
#import "@preview/algorithmic:1.0.7"
#import "@preview/subpar:0.2.2"
#import algorithmic: style-algorithm, algorithm-figure, algorithm

#import "preamble_funcs.typ": mref, mathformatter
#import "style.typ": style, frontpage-1

#show: style-algorithm
#show: style
#set page(numbering: "i")
#show figure.caption: set text(size: 9.5pt)

// Fix smallcaps in math
#show math.equation: it => {
  show smallcaps: set text(font: "Libertinus Serif")
  it
}

#let num = num.with(multiplier: "×")

#let vv = mathformatter(underbars: 0, bold: true, upright: false)
#let mm = mathformatter(underbars: 0, bold: false, upright: true)
#let vvh = x => math.hat(vv(x))

#set math.mat(delim:"[")

#let Cov = math.op("Cov")
#let diag = math.op("diag")
#let median = math.op("median")
#let TWIP = math.op("TWIP")

#let wider = h(3em)

#let HH = math.sans(math.upright($H$))

#let nobreak(body) = block(breakable: false)[#body]
#let steal-margin(body) = block[
  #show: pad.with(right: -10mm)
  #body
]
#let right-figs(inner) = steal-margin[#meander.reflow(inner)]

#set document(title: "Report on Probablistic Inversion and Cross Hole GPR Tomography")
#let authors = (
  [Jeppe Klitgaard #link("emailto:s250250@dtu.dk", raw("<s250250@dtu.dk>"))],
)

#let abstract = [
  #set heading(numbering: none)
  #show heading: set text(size: 16pt)
  #set par(justify: true)
  #set text(size: 12pt)
  = Abstract

  The existence of the South Atlantic Anomaly (SAA) in the geomagnetic field remains a topic of debate within the literature. We present reconstructions of the radial geomagnetic field at the Core-Mantle Boundary (CMB) which to high certainty support findings of a region of reversed flux. The employed model solves the inverse problem of determining the Gauss coefficients of a spherical harmonic expansion using data from ESA's Swarm mission. Two inversions are carried out: one with unweighted data and a penalty on the forward projection onto the CMB in the $L_2$-norm, and the other with robustly reweighted data employing the $L_1$-norm instead. The latter model is determined using the Iteratively Reweighted Least Squares algorithm with Huber weights and a perturbed $L_p$-norm. Both of the produced models reveal flux reversal beneath the South Atlantic Ocean. Subsequent analysis of the prediction errors and the spatial resolution of the model leads to the conclusion that the null-hypothesis of purely dipolar geomagnetism must be rejected.
]

#frontpage-1(front-content: abstract, authors: authors, date: datetime(year: 2026, month: 3, day: 19))
#counter(page).update(1)
#counter(heading).update(0)

#set heading(numbering: none)
#outline()

#v(10mm)
= Artificial Intelligence Declaration
The body of work presented below is solely authored by the author of this document, Jeppe Klitgaard, and none of the presented material has been produced by means of _Generative AI_ in the manner requiring citation as described by @ref:dtu-genai outside the entries listed below.

Generative AI has, however, been used selectively in the generation of limited parts of the code. GitHub Copilot tab-completion has been used during the development of the attached code and in the generation of an initial version of the following entries:
- Some parts of `nanopinv` (which falls outside the examined material of this report), in particular:
  - A number of implementations of the Fast Sweeping Method
  - Boiler-plate code around benchmarking and testing these implementations
  - Initial implementations for the plotting methods supported by the `History` class

The generated code output has in all instances been checked and often revised by the authors and is _not_ considered to be a _significant_ contribution to total intellectual work undertaken and presented.

#set heading(numbering: "1.1.a")
#pagebreak()
#counter(page).update(1)
#set page(numbering: "1/1")

// // Measure paper width
// #let get-page-width(id) = layout(size => [
//   #let obj = line(length: 100%)
//   #let (width, height) = measure(..size, obj)
//   #obj
//   #id: Available page width: #width.mm()mm
// ])

// #get-page-width([`normal`])

// #steal-margin[
//   #get-page-width([`steal-margin`])
// ]

= Introduction <sec:intro>

This report seeks to demonstrate the use of probabilistic inversion techniques to solve the non-linear inverse problem of Ground Penetrating Radar (GPR) Tomography in which the spatial distribution of the electromagnetic phase velocities through the subsurface are reconstructed from observed travel times of radar waves between sources and receivers located in two nearby boreholes.

In order to produce the results presented in this report, we have developed a performant Python library, `nanopinv`, which is able to achieve accurate reconstructions of the modelled field by sampling the posterior distribution of the modelled field using Markov Chain Monte Carlo (MCMC) methods. The library is implemented in JAX @bib:jax and can efficiently leverage multiple CPU cores or GPU hardware.

The performance of the library is further enhanced by implementations of the Fast Sweeping Method (FSM) @bib:fsm-zhao for solving the eikonal equation and enables to user to rapidly obtain independent samples by leveraging the Parallel Tempering method.

The structure of the report deviates from the typical layout of an experimental report and instead interleaves theory and results. This choice is made in the interest of sparing the reader the rather hefty theory section featured in the first draft of the report. Further, the exploratory and interdisciplinary nature of the project incentivises a less rigid structure in which theory and design choices are introduced as they become relevant to the discussion of the results. Alternatively, the reader may simply regard @sec:theory as a lengthy theory section with accompanying examples, which is followed by more focused results and discussions sections in @sec:results and @sec:discussion respectively.

= Data and Experimental Design <sec:experiment-data>
#let fig-data-rays = [
  #box(width: 50%)[
    #figure(
      image("export/data_rays_50.png"),
      caption: "Source and receiver locations with observed travel times overlaid on the model space."
    ) <fig:data-rays>
  ]
]
#let fig-data-observations = [
  #box(width: 50%)[
    #figure(
      image("export/data_observations_50.png"),
      caption: "Observed travel times for all source-receiver pairs."
    ) <fig:data-observations>
  ]
]

#right-figs({
  import meander: *

  placed(top + right, stack(fig-data-rays, v(1em), fig-data-observations))
  container()
  content([

The dataset used in this report consists of $N_d = 758$ travel time measurements taken using pairs of $N_"source" = 52$ source locations and $N_"receiver" = 235$ receiver locations in two boreholes at the Boise Hydrogeophysical Research Site (BHRS) near Boise, Idaho in the United States. The boreholes are separated by approximately 10 metres and measurements are obtained for elevations in the range $y∈[#qty("831.51", "m"), #qty("846.79", "m")]$ as shown on @fig:data-rays.

During the experiment, measurements were obtained are using the technique of Cross Hole Ground Penetrating Radar (GPR), in which the delays in the observed travel times of electromagnetic radar waves emitted by the source and received by its corresponding network of receivers were recorded for each source-receiver pair.

Additionally, we are given 3 samples of subsurface velocity fields shown in @fig:empirical-prior, which server as a source of prior information about the spatial distribution of velocities that we may expect to find using the dataset. Lastly, we are informed that the water table at the site is known to be located at an elevation of $y_* = #qty("844", "m")$.
])})

#figure(
  image("export/empirical_prior_samples.png"),
  caption: "The three samples representing the prior knowledge of the velocity field."
) <fig:empirical-prior>



= Theory and Exploration <sec:theory>
In order to arrive at a suitable model and framework for obtaining it, we will first lay out a brief summary of the relevant theory from the fields of physics and probabilistic inversion.



TODO Much more introduction


== Probabilistic Inversion <sec:theory-probinv>
In order to solve the inverse problem of estimating the velocity distribution from the observed travel times, we reach for the framework of _probabilistic inversion_, in which we adopt the notation due to Tarantola @bib:tarantola2005. This choice of framework is convenient for travel-time tomography problems as it handles non-linearities in the forward model naturally, while also offering a principled way to incorporate prior information and uncertainty quantification into the inversion.

We may tersely summarise the framework as presented in @ref:course_book[4.2] by letting $ρ_m (vv(m)) : cal(M) → ℝ$ and $ρ_d (vv(d)) : cal(D) → ℝ$ denote the _prior model probability distribution_ and _prior data probability distribution_ respectively, with their arguments denoting the model parameters and data. For the problem covered by this report, the model space is 2-dimensional, $cal(M) = ℝ^(N_x × N_y)$, while the data space is 1-dimensional, $cal(D) = ℝ^(N_d)$. By assuming the data to be independent of the model parameters, we construct a joint prior distribution:
$
  ρ(vv(d), vv(m)) = ρ_d (vv(d)) ρ_m (vv(m)).
$ <eq:joint-prior>

The forward model is then introduced through $vv(d) = g(vv(m)) + vv(e)$ where the error $vv(e)$ captures any numerical errors in the computation of the forward model $g(vv(m))$. Additionally, any potential stochasticity in the forward model and systematic biases such as those introduced by the simplifying assumptions made in the physical model as discussed in @sec:theory-physics enter through this error term. We may capture these effects by another joint probability distribution, $θ(vv(d), vv(m))$ that links the model and data spaces through the forward model.

We introduce the notion of _homogenous probability density functions_, $μ(⋅)$, which have probability densities proportionate to the volume element of the space on which they are defined. As such, they may be used to ensure correct normalization of probability distributions under coordinate transformations and notably do not carry any information on their own.

With this in mind, we relate the joint linking distribution $θ(vv(d), vv(m))$ to its conditional through the product rule to obtain
$
  θ(vv(d), vv(m))
  &= θ(vv(d)|vv(m)) θ(vv(m))\
  &= θ(vv(d)|vv(m)) μ_m (vv(m))\
$ <eq:linking-distribution>
where $θ(vv(m))$ is identified with $μ_m (vv(m))$ on the grounds that neither carry any information about the model parameters. This can be understood as $θ(vv(d), vv(m))$ not describing the model parameters themselves, only how the data is linked to them.

As shown in @bib:tarantola2005[Eq. 1.83], these distributions may be combined to obtain the _joint posterior model distribution_:
$
  σ(vv(m), vv(d)) = k (ρ(vv(m), vv(d)) θ(vv(d), vv(m)))/(μ(vv(m), vv(d))),
$
where $k$ is a normalisation constant and $μ(vv(m), vv(d))$ is the _joint homogenous probability density function_. This further simplies to the _posterior model distribution_ under the following assumptions
$
  σ_m (vv(m))
  &= k ∫_cal(D) (ρ(vv(m), vv(d)) θ(vv(d), vv(m)))/(μ(vv(m), vv(d))) dif vv(d)
  &wider #[Margin. over data space, $cal(D)$]
  \
  &= k ρ_m (vv(m)) ∫_cal(D) (ρ_d (vv(d)) θ(vv(d), vv(m)))/(μ(vv(m), vv(d))) dif vv(d)
  &wider #[@eq:joint-prior]
  \
  &= k ρ_m (vv(m)) ∫_cal(D) (ρ_d (vv(d)) θ(vv(d), vv(m)))/(μ_m (vv(m)) μ_d (vv(d))) dif vv(d)
  &wider #[Assume $μ(vv(m)), μ(vv(d))$ independent]
  \
  &= k ρ_m (vv(m)) ∫_cal(D) (ρ_d (vv(d)) θ(vv(d)|vv(m)) cancel(μ_m (vv(m))))/(cancel(μ_m (vv(m))) μ_d (vv(d))) dif vv(d)
  &wider #[@eq:linking-distribution]
  \
  &= k L(vv(m)) ρ_m (vv(m))
  &wider L(vv(m)) ≔ ∫_cal(D) (ρ_d (vv(d)) θ(vv(d)|vv(m)))/(μ_d (vv(d))) dif vv(d)
$ <eq:posterior>
where $L(vv(m)) $ is the _likelihood probability density_ describing how likely the observed data $d$ is given the model parameters $vv(m)$.

=== Likelihood

If we assume the data to be independently and normally distributed and additionally assume the same of the errors introduced by the forward model, both $ρ_d (vv(d))$ and $θ(vv(d)|vv(m))$ can be expressed as independent normal distributions:
$
  ρ_d (vv(d)) &= cal(N)(vv(d)|vv(μ)_d, Σ_d)\
  θ(vv(d)|vv(m)) &= cal(N)(vv(d)|g(vv(m)), Σ_e)\
$
which convolve to give the likelihood as another normal distribution:
$
  L(vv(m)) ∝ cal(N)(g(vv(m))|vv(μ)_d, Σ_d + Σ_e).
$

We further assume the errors of the forward model to be negligible compared to the variance of the data distribution, which is known a priori to be #qty("4", "ns"). That is, given $Σ_e ≪ Σ_d$, the likelihood can be approximated as
$
  L(vv(m)) ∝ cal(N)(g(vv(m))|vv(μ)_d, Σ_d).
$ <eq:likelihood>

== Prior Information and 2-Point Statistics <sec:theory-prior>

In order to turn the given velocity field samples shown in @fig:empirical-prior into a usable prior distribution, realise that 1-point statistics will not suffice due to the strong spatial correlation in the field. Instead, we turn to 2-point statistics wherein the covariance between pairs of points in the field can be modelled. By assuming that the covariance of the field is isotropic and the process it is drawn from to be stationary, we propose that the covariance between any two points in the field depends only their separation, $h$. This leads to the notion of _covariance functions_ of the form $C(h)$, with a notable example being the _spherical covariance function_:
$
  C_"sph" (h) = cases(
    σ^2 (1 - (3h)/(2s) + (h^3)/(2s^3)) &wide h ≤ s,
    0 &wide h > s,
  )
$
where the _variance_, $σ^2$, and _range_, $s$, are hyperparameters of the covariance function.

In order to assess the suitability of a covariance function for a given set of samples, we may estimate the empirical _semi-variogram_, $hat(γ)(h)$ of the samples by binning the separations into $N$ bins of width $Δ h$:
$
  hat(γ)(h_k ± Δ h) = 1/(2 N_k) ∑_((i, j) ∈ S_k) (vv(v)_i - vv(v)_j)^2
$
where $h_k$ is the bin centre, $S_k = {(i, j) : h_k - Δ h < norm(vv(v)_i - vv(v)_j) < h_k + Δ h}$ is the set of point pairs in the bin and $N_k$ is the number of elements in $S_k$. The empirical semi-variogram can then be obtained by computing $hat(γ)(h_k ± Δ h)$ for each bin and averaging across the samples.

Such a semi-variogram is related to the covariance function by
$
  γ(h) = 1/2 (σ^2 - C(h))
$
where notably $γ(h)$ here represents an analytical, exact semi-variogram.

=== Obtaining a prior distribution <sec:prior-dist>

In order to obtain a covariance function that suitably captures the spatial correlations in the prior samples, we may perform a fit of candidate correlation functions $C(h; θ)$ in which their hyperparameters $θ$ are determined by minimising the mean squared error between the empirical semi-variogram and the analytical semi-variogram implied by the covariance function:
$
  θ^* = arg min_θ norm(hat(γ)(h) - 1/2(σ^2 - C(h; θ)))_2^2
$

By first computing the empirical semi-variogram of the prior samples from @fig:empirical-prior and subsequently fitting a variety of candidate covariance functions to it, we are able to pick a suitable covarinace function. The empirical semi-variogram is estimated using `gstools` @bib:gstools, which is also used to perform the fitting. Notably, the global tolerance parameter of the underlying fitting routine is decreased to $10^(-14)$ as the fits were otherwise observed to be poor.

#figure(
  image("export/empirical_variogram_fits.png"),
  caption: "Empirical semi-variogram of the prior samples and fits of candidate covariance functions."
) <fig:empirical-variogram-fits>

While the circular covariance function slightly outperforms the spherical covariance function in terms of the $R^2$ metric, we select the spherical covariance function for our prior distribution as it appears to be a more standard choice in the literature. The parameters of the fitted spherical covariance function are found to be $σ^2 = num("8.06e-6"), s=num("5.16"), τ^2 = num("1.83e-6")$ where $τ^2$ is the _nugget_, which represents the variance at zero separation and is simply added to the covariance function as $C(h) + τ^2$.

In order to obtain a prior distribution using two-point statistics, we also need to obtain an estimate of the mean of the distribution. As we have already assumed isotropy, we let the mean be constant across the field and estimate it by averaging across the prior samples to obtain $μ = num("0.0849")$.

With estimates of both the mean and covariance function, we are able to construct a Gaussian prior distribution over the model space. In practice, this involves constructing the full covariance matrix, which for a model discretisation of $Δ x = Δ y = qty("0.25", "m")$ gives a model grid of $N_x = 41, N_y=63$ and covariance of size $2583×2583$. The efficiency of the prior sampling in `nanopinv` is enhanced by performing a Cholesky decomposition of the covariance matrix to obtain the lower triangular matrix $mm(L)$ during initialisation. Samples are then produced by sampling a standard normal distribution $z∼cal(N)(0, 1)$ and applying the affine transformation $x = μ + L z$.

To gauge whether the obtained prior distribution is a plausibly represents the prior samples, we draw 3 samples from the distribution and compare them to the original samples in @fig:empirical-vs-estimated-prior-samples. The drawn samples appear to be similar in the characteristic length scale and magnitude of the variations, though notably does not reproduce the small-scale, directional patterns that appear to be present in the original samples. These violate the assumption of isotropy and may further be speculated to be systematic errors from the experiment.

#figure(
  image("export/empirical_vs_estimated_prior_samples.png"),
  caption: "Comparison of samples from the obtained Gaussian prior distribution and those of the empirical prior."
) <fig:empirical-vs-estimated-prior-samples>

== Physics and Forward Model <sec:theory-physics>

Having produced a prior distribution over the model space, we turn to the physics involved in the ground penetrating radar experiment in order to obtain a suitable forward model, $g : cal(M) -> cal(D)$, that maps the model parameters to the space of observations.

=== Eikonal Equation and Wave Travel Time <sec:physics-eikonal>
In order to relate the observed travel times to the spatial distribution of the phase velocity, we refer to @born-optics[3.1.1] in which the _eikonal equation_, @eq:eikonal, is derived from Maxwell's equations under a number of relevant assumptions.
$
  norm(∇ S(vv(r)))_2^2 = n^2 (vv(r))
$ <eq:eikonal>
where $S : ℝ^N -> ℝ$ is the _optical path_ as a function of position $vv(r)$ and $n$ is the _refractive index_.
We may relate this to the travel time $T(vv(r))$ and phase velocity $v(vv(r))$ by noting that the optical path is defined as the distance travelled in vacuum, $S(vv(r)) = c_0 T(vv(r))$ and recalling the definition of the refractive index, $n(vv(r)) = c_0 / v(vv(r))$, which yields
$
  norm(∇ T(vv(r)))_2 = 1 / (v (vv(r))).
$ <eq:eikonal-travel-time>

The derivation due to Born and Wolf relies on the assumption that the characteristic length scale of the variations in the refractive index — or equivalently, the phase velocities — are much larger than the wavelengths of the electromagnetic waves. If we assume the data to be collected with similar equipment to that used as part of the study described in @boise-gpr, we may obtain estimates of the wavelengths of the radar waves with which we can reason about the validity of the high-frequency assumption, $λ ≪ L$, where $L$ is the characteristic length scale of variations in the modelled field through the medium.

The study cites an antenna frequency of $f = #qty("100", "MHz")$, with relative permittivities $ε\/ε_0$ of $≈ {1, 4, 81} med ["F/m"]$ for air, soil, and water respectively. We may find the wavelengths associated with these parameters using the relation
$
  λ = v/f = c_0 / (f sqrt(ε\/ε_0)),
$
which gives $λ_"soil" = #qty("1.5", "m")$ and $λ_"water" = #qty("0.3", "m")$. These effectively set a lower bound on the size of the features that can be resolved using the acquired data.
This is a useful metric to keep in mind, as it also provides a useful scale of discretisation for the model space. In particular, no additional resolution can be obtained by further refining the model discretisation beyond this limit without also increasing the frequency of the radar waves or resorting to a more complex physical model that relaxes the high-frequency assumption. As such, the previously proposed discretisation of $Δ x = Δ y = #qty("0.25", "m")$ may already be below the limit of resolution set by the GPR frequency and the properties of the medium.

=== Numerical Solution of the Forward Problem <sec:theory-eikonal-solution>
In order to compute the forward projection of a given model, here understood to be a discretised spatial distribution of phase velocities, we must solve @eq:eikonal-travel-time, and ideally do so as efficiently as possible.

A popular choice of algorithm dedicated to solving the family of equations to which the eikonal equation belongs is the Fast Marching Method (FMM) due to Sethian @bib:fmm, which are implemented as first and second order methods in @bib:skfmm.
The finer details of the FFM algorithm is beyond the scope of this report, but a notable consequence of the formulation of the algorithm is that it is not amenable to parallelisation. While multiple source-receiver pairs may be computed in parallel using CPU cores, the method is not able utilise modern, massively-parallel hardware such as GPUs efficiently.

For this reason, the alternative Fast Sweeping Method (FSM) due to Zhao @bib:fsm-zhao @bib:fsm-zhao-parallel is considered a superior choice for solving the eikonal equation on GPUs, as it more naturally lends itself to parallelisation. Taking the paper @bib:fsm-2nd-order as a reference, `nanopinv` implements both a first and second order version of the parallel FSM algorithms suitable for use with CPU or GPU hardware.

Common to all the implementations is that they produce an approximation of the travel time $T$ from a source $vv(r)_s$ to all other points in the model space when given a propagation velocity field $v(vv(r))$.
The estimated travel time to the receivers, $vv(r)_r$, may then be obtained by interpolation. @fig:forward-model-prior-samples shows the forward projections obtained by the 2nd order FSM solver for three samples drawn from the prior distribution obtained in @sec:prior-dist.

#figure(
  image("export/forward_model_prior_samples.png"),
  caption: "Comparison of the travel time fields obtained by the first and second order FSM implementations for a given velocity field."
) <fig:forward-model-prior-samples>

== Markov Chain Monte Carlo Methods <sec:theory-mcmc>

In order to sample the posterior, $σ_m (vv(m))$, we turn to the theory of Markov Chains and Monte Carlo sampling, in which a chain of samples is generated by a Markov process — that is, a process in which the next state depends only on the current state.

While many algorithms for simulating such processes exist, we will focus on the _Extended Metropolis algorithm_ due to Mosegaard and Tarantola @bib:extended-metropolis whose transition procedure is described in @alg:extended-metropolis. In order to efficiently sample the prior model distribution even for large model spaces, we leverage the Preconditioned Crank-Nicolson proposal algorithm given in @alg:preconditioned-crank-nicolson.

It should be noted that implementations of @alg:extended-metropolis typically work with the log-transformed likelihood and ensure the forward model and likelihood are only evaluated once per model by to improve performance.

#nobreak[
  #algorithm-figure(
    "Preconditioned Crank-Nicolson Proposal",
    vstroke: .5pt + luma(200),
    inset: 0.35em,
    {
      import algorithmic: *
      let Solve = Call.with("Solve")
      let Input(doc) = Line([*input*~#context box(doc, baseline: 100% - par.leading)])
      let Output(doc) = Line([*output*~#context box(doc, baseline: 100% - par.leading)])
      Procedure(
        "PreconditionedCrankNicolsonProposal",
        ([$vv(m)_0$], [$δ$], [$ρ_m$]),
        {
          Input([
            Initial model $vv(m)_0 ∈ cal(M)$,\
            Step size $δ ∈ (0, 1)$,\
            Prior distribution $ρ_m : cal(M) → ℝ$
          ])
          Output([
            Proposed model $vv(m)_p ∈ cal(M)$
          ])
          LineComment(
            Assign[$vv(μ)$][$𝔼[ρ_m]$],
            [Mean of prior model distribution]
          )
          LineComment(
            Assign[$vv(m)_*$][$vv(x) ∼ ρ_m$],
            [Sample from prior model distribution]
          )
          LineComment(
            Assign[$vv(m)_p$][$vv(μ) + sqrt(1 - δ^2) (vv(m)_0 - vv(μ)) + δ (vv(m)_* - vv(μ))$],
            [Preconditioned Crank-Nicolson]
          )
          Return[$vv(m)_p$]
        }
      )
    },
  ) <alg:preconditioned-crank-nicolson>
]

#nobreak[
  #algorithm-figure(
    "Extended Metropolis Algorithm",
    vstroke: .5pt + luma(200),
    inset: 0.35em,
    {
      import algorithmic: *
      let Solve = Call.with("Solve")
      let Input(doc) = Line([*input*~#context box(doc, baseline: 100% - par.leading)])
      let Output(doc) = Line([*output*~#context box(doc, baseline: 100% - par.leading)])
      let True = [#smallcaps("true")]
      let False = [#smallcaps("false")]
      Procedure(
        "ExtendedMetropolisStep",
        ([$vv(m)_0$], [$β$], [$q$], [$g$], [$L$]),
        {
          Input([
            Initial model $vv(m)_0 ∈ cal(M)$, Inverse temperature $β ∈ (0, 1]$,\
            Proposal algorithm $q : cal(M) → cal(M)$,\
            Forward model $g : cal(M) → cal(D)$, Likelihood $L : cal(D) → ℝ$,\
          ])
          Output([
            Final model $vv(m) ∈ cal(M)$, Acceptance $∈ {True, False}$
          ])
          LineComment(
            Assign[$vv(m)_p$][$q(vv(m)_0)$],
            [Propose new model]
          )
          LineComment(
            Assign[$P_"accept"$][$min(1, (L(g(vv(m)_p)) / L(g(vv(m)_0)))^β)$],
            [Compute acceptance probability]
          )
          LineComment(
            Assign[$α$][$x ~ cal(U)(0, 1)$],
            [Sample uniform distribution]
          )

          IfElseChain(
            $α < P_"accept"$,
            {
              LineComment(
                Return[$vv(m)_p, True$],
                [Accept proposal]
              )
            },
            {
            LineComment(
              Return[$vv(m)_0, False$],
              [Reject proposal]
            )
            }
          )
        }
      )
    }
  ) <alg:extended-metropolis>
]



=== Sampling the Posterior <sec:sampling-posterior>

By using these algorithms with sufficient care, we will be able to sample from the posterior distribution of the model parameters, $σ_m (vv(m))$. Firstly, the choice of the step size, $δ$, will have significant influence on the efficiency with which we can sample. By inspecting the proposal algorithm in @alg:preconditioned-crank-nicolson, it is clear that a larger step size will retain less of the current model and instead be more heavily influenced by the sample drawn from the prior distribution. As such, a larger step size will more aggressively explore the model space. This is of course desirable, but comes at the cost of a lower acceptance rate as many of the proposed models will inevitably lie in regions of low likelihood inside the model space. Conversely, if the chosen step size is too small, the acceptance rate will be high, but the chain will be slow to explore the model space and will be more likely to get stuck in local valleys within the modal space, thus failing to accurately sample the posterior. The literature generally favours a step size corresponding to an acceptance rate of approximately $25%$.

The second important hyperparameter is the inverse temperature, $β ∈ (0, 1]$, which directly affects the acceptance probability of the proposed models. As $β$ goes towards zero, the acceptance probability approaches unity, and as such the inverse temperature may be understood to effectively flatten the likelihood landscape thus making it easier to traverse. This comes with the unfortunate consequence that samples drawn during simulations away from unit inverse temperature are not distributed according to the true posterior distrubution. As such, we concern ourselves only with _cold sampling_ at $β=1$ for now.

Further, we must only consider _accepted_ samples as being drawn from the posterior distribution, and must additionally take care that the drawn samples are sufficiently decorrelated from one another as they would otherwise bias our estimate of the posterior distribution.

#let fig-emc-1-burn-in = [
  #box(width: 50%)[
    #figure(
      image("export/emc_1_burn-in_50.png"),
      caption: [
        Trace log likelihood during the burn-in phase of the Extended Metropolis algorithm.
        ]
    ) <fig:emc-1-burn-in>
  ]
]

#let fig-emc-1-samples-likelihood = [
  #box(width: 50%)[
    #figure(
      image("export/emc_1_samples-likelihood_50.png"),
      caption: [
        Trace of log likelihood values for 10 chains during sampling after burn-in.
      ]
    ) <fig:emc-1-samples-likelihood>
  ]
]

#let fig-emc-1-samples-acceptance = [
  #box(width: 50%)[
    #figure(
      image("export/emc_1_samples-acceptance_50.png"),
      caption: [
        Trace of acceptance rate for samples during sampling. Red dashed line indicates $25%$ and black line denotes mean acceptance rate across all chains.
      ]
    ) <fig:emc-1-samples-acceptance>
  ]
]

#let fig-emc-1-samples-autocorrelation = [
  #box(width: 50%)[
    #figure(
      image("export/emc_1_samples-autocorrelation_50.png"),
      caption: [
        Autocorrelation of the samples during sampling. Black line denotes mean across all chains.
      ]
    ) <fig:emc-1-samples-autocorrelation>
  ]
]

#pagebreak(weak: true)
#right-figs({
  import meander: *

  placed(top + right, stack(fig-emc-1-burn-in, v(2em), fig-emc-1-samples-likelihood, v(2em), fig-emc-1-samples-acceptance, v(2em), fig-emc-1-samples-autocorrelation))
  container()
  content([
We set up the `ExtendedMetropolisChain` sampler from `nanopinv` with the proposal algorithm given in @alg:preconditioned-crank-nicolson using the previously obtained prior distribution $ρ_m (vv(m))$ and a step size of $δ = 0.05$. An initial model $vv(m)_0$ is drawn directly from the prior distribution and the $N_"iter" = num("1000")$ steps are taken while recording intermediate samples and acceptance decisions. It is highly likely that the sample drawn from the prior distribution will lie in a region of low likelihood, which leads to the burn-in phenomenon shown in @fig:emc-1-burn-in, where the log#{sym.hyph.nobreak}likelihood is observed to decrease rapidly during the initial steps of the chain before stabilising as the chain diffuses towards a higher likelihood region of the model space.

Using JAX it is particularly easy to run multiple chains in parallel using the `vmap` function. This is useful as it allows us to obtain multiple independent samples, and additionally is a useful way to gauge whether the chains get stuck in local minima by comparing the results of multiple chains initialised with different random seeds. By performing $N_"burn-in" = num("3000")$ steps on $N_"chains" = 10$ chains, and subsequently recording another $N_"samples" = num("10 000")$ samples, we are able to produce a large number of samples. Running the cumulative #num("100 000") steps takes approximately 6 minutes and 30 seconds on an NVIDIA RTX 4090 GPU.

Inspecting @fig:emc-1-samples-likelihood reveals that all chains appear to have converged to a similar level set of the likelihood landscape, which is an indication that the chains are not getting stuck in local minima.
@fig:emc-1-samples-acceptance shows that the acceptance rate is observed to be around $45%$ during the sampling, which is higher than the $25%$ target. This can negatively affect the correlation times and heightens the risk of poor mixing.

As observed in @fig:emc-1-samples-autocorrelation, the autocorrelation of the likelihood roughly follows a decaying exponential with a characteristic decay length of roughly #num("700") steps. We note with caution that the autocorrelation of the _likelihood_ is not perfectly representative of the autocorrelation of the _model parameters_ and may be artificially longer in cases where the proposal along contours of the likelihood landscape. In order to obtain uncorrelated samples we must thin the accepted samples such that the number of steps between recorded samples is greater than some small multiple of the decay length.
])})

By requiring a minimum interval of #num("2000") steps between recorded samples we obtain a total of #num("50") uncorrelated samples from the posterior distribution, 3 of which are shown in @fig:emc-1-posterior-samples.
#figure(
  image("export/emc_1_posterior_samples_100.png"),
  caption: "Posterior samples obtained by the Extended Metropolis algorithm."
) <fig:emc-1-posterior-samples>

=== Step size tuning <sec:step-size-tuning>

#let fig-emc-1-tuning = [
  #box(width: 50%)[
    #figure(
      image("export/emc_1_tuning.png"),
      caption: [
        Trace of step size during tuning procedure. Red dashed line indicates the initial step size and black dashed line indicates the final step size after tuning.
      ]
    ) <fig:emc-1-tuning>
  ]
]

As mentioned previously, the choice of step size has significant implications for the efficiency of the sampling. For this reason, `nanopinv` additionally implements a simple step size tuning procedure in which the step size is adjusted iteratively by using batches of $N_"batch"$ samples to compute the update according to the Robbins-Monro algorithm:
$
  δ_(k+1) = δ_k + ξ_k (dash(P)_k - P_"target") wide "where" wide ξ_k = ξ_0 / (k + 1)^(-κ)
$
#right-figs({
  import meander: *

  placed(top + right, fig-emc-1-tuning)
  container()
  content([

where $ξ_k ∈ (0, 1]$ is a _learning rate_ which we decay from a base value $ξ_0$ according to the batch number $k$ and a decay parameter $κ ∈ (1\/2, 1]$. $P_"target"$ denotes the target acceptance rate and $dash(P)_k$ is the observed acceptance rate for the $k$th batch of samples.

This is implemented in `nanopinv` via the `tune` method of the `ExtendedMetropolisChain` class, which for $κ=0.5, ξ_0=0.1, N_"batch"=num("1000")$ produces the tuning history shown in @fig:emc-1-tuning, which converges to a step size of $δ = num("0.808")$ and achieves a acceptance rate of approximately $25%$.

])})

=== Parallel Tempering

The primary bottle-neck of our current simulations are the long correlation times between the samples. While running multiple chains in parallel is a useful workaround, the performance can be further improved by implementing a _parallel tempering_ scheme. In @sec:sampling-posterior we discussed how the inverse temperature, $β$, can flatten the likelihood landscape, thus allowing much larger steps to be taken without compromising the acceptance rate. With this in mind, we can arrange a ladder of chains with increasing temperatures and individual step sizes. Hotter chains will be able to easily traverse any hills in the likelihood landscape and reach all parts of the model space. Unfortunately, the samples drawn from these hotter chains still cannot be used to garner information about the posterior distribution.

Instead, we propose to swap the states of neighbouring chains on the ladder according to a criterion that preserves _detailed balance_ such that the chains remain in equilibrium at all times. This allows the hot chains to lead the charge in exploring the model space while the coldest chain is still able to draw samples from the true posterior distribution. While the full theory behind the parallel tempering scheme is beyond the scope of this report, an accessible description by Paweł Czyż can be found in @bib:czyz.

In particular, the variant of parallel tempering that is implemented in `nanopinv` is the _Deterministic Even Odd_ (DEO) _non-reversible parallel tempering_ (PT) scheme due to Syed et al. @bib:syed, which is outlined in @alg:non-rev-pt where the DEO variant arises when the parity offset $Ξ$ is alternated between $0$ and $1$ for each swap. This approach has the limitation that chain states can only diffuse by at most one rung on the ladder per swap, which puts a lower bound on the number of swaps required to traverse the ladder and thus the correlation time of the samples.

With this in mind, the authors of @bib:syed introduce the notion of the _communication barrier_, which is used to produce a tuning scheme for the temperature ladder, which is optimally tuned when the swapping rate is uniform across all rungs of the ladder. An implementation of this tuning scheme alongside a Robbins-Monro-based tuning scheme akin to that discussed in @sec:step-size-tuning is implemented in `nanopinv`, closely following that presented in @bib:czyz.

#nobreak[
  #algorithm-figure(
    "Non-Reversible Parallel Tempering Swap Procedure",
    vstroke: .5pt + luma(200),
    inset: 0.35em,
    {
      import algorithmic: *
      let Input(doc) = Line([*input*~#context box(doc, baseline: 100% - par.leading)])
      let Output(doc) = Line([*output*~#context box(doc, baseline: 100% - par.leading)])
      let True = smallcaps("true")
      let False = smallcaps("false")

      Procedure(
        "AdjacentSwap",
        ([$N$], [$Ξ$], [$X$], [$L$], [$β$],),
        {
          Input([
            Number of chains $N$, Parity offset $Ξ ∈ {0, 1}$,\
            Chain states $X = (x_0, ..., x_(N-1))$, Likelihoods $L = (L_0, ..., L_(N-1))$,\
            Inverse temperatures $β = (β_0, ..., β_(N-1))$
          ])
          Output([
            Updated states $X'$, Updated likelihoods $L'$, Acceptances $A ∈ 𝟙^N$
          ])

          LineComment(
            Assign[$A_i$][False $∀ i ∈ {0, ..., N-2}$],
            [Initialize acceptances]
          )
          LineComment(
            Assign[$X', L'$][$X, L$],
            [Initialize outputs]
          )

          For(
            [$i = 0, ..., N-2$],
            {
              If(
                [$(i mod 2) = Ξ$],
                {
                  LineComment(
                    Assign[$α_i$][$min(1, (L_(i+1) / L_i)^(β_i - β_(i+1)))$],
                    [Compute acceptance probability]
                  )
                  LineComment(
                    Assign[$u_i$][$u ~ cal(U)(0, 1)$],
                    [Sample uniform distribution]
                  )
                  If(
                    [$u_i < α_i$],
                    {
                      LineComment(
                        Assign[$A_i$][True],
                        [Accept swap proposal]
                      )
                      LineComment(
                        Assign[$X'_i, X'_(i+1)$][$X'_(i+1), X'_i$],
                        [Apply state swap]
                      )
                      LineComment(
                        Assign[$L'_i, L'_(i+1)$][$L'_(i+1), L'_i$],
                        [Apply likelihood swap]
                      )
                    }
                  )
                }
              )
            }
          )
          LineComment(
            Return[$X', L', A$],
            [Return updated chains and acceptance outcomes]
          )
        }
      )
    }
  ) <alg:non-rev-pt>
]

The parallel tempering scheme combines with @alg:preconditioned-crank-nicolson and @alg:extended-metropolis by first performing a step of the Extended Metropolis algorithm on each chain in parallel and then performing the non-reversible parallel tempering swap procedure on the ladder of chains.

Importantly, a chain state is considered an _acceptable_ sample of the posterior distribution if it resides on the _cold_ chain after the step and swap and has either been locally updated by the proposal algorithm _or_ has been swapped with its neighbour. Even if the state has just arrived from a hotter chain, it is still considered an acceptable sample as the swap procedure obeys _detailed balance_ and thus the state is still distributed according to the posterior distribution.

While a more elegant approach to jointly tune the step sizes and temperatures of the ladder may exist, we resort to iteratively tuning the step sizes and temperatures separately, which is found to converge well for the data used in this report.

=== Obtaining a Tuned Parallel Tempering Sampler

In order to obtain a tuned parallel tempering sampler we instantiate a ladder of chains and evolve them according to the following procedure:
+ Construct $N_"chains" = 16$ chains:
  + Initial inverse temperatures: $β_k = ω^(-k)$ where $ω = 2.0 ∈ [1, ∞)$ is a hyperparameter controlling the initial temperature spacing.
    - As a heuristic, $ω$ should be chosen such that the hottest initial temperature is just hot enough that a chain acceptance rate of $25%$ is reached for a step size of $δ = 1$.
    - This ensures that the hottest chain is able to explore the model space entirely unhindered and effectively corresponds to drawing proposals directly from the prior distribution at each step.
    - There is no goldilocks acceptance rate for non-reversible parallel tempering — the acceptance rate will generally be larger for more densely packed temperature ladders, which further lowers the correlation time by improving mixing between the chains.
  + Initial step sizes: linearly spaced over the range $[0.1, 1.0]$.
  + Initial states: Draw initial states $vv(m)_0$ from prior $ρ_m (vv(m))$.
+ Burn in: Perform $N_"burn" = num("500")$ steps of the PT scheme.
+ Tune step sizes: Perform 3 batches of $N_"batch" = 300$ steps of the PT scheme and tune the step sizes according to the Robbins-Monro algorithm with $κ=1\/2, ξ_0=1$ and target acceptance rate of $25%$.
+ Tune temperatures: Perform a single batch of $N_"batch" = num("1000")$ steps of the PT scheme and tune the temperatures according to the method outlined in @bib:syed.
+ Tune step sizes: Perform another 5 batches of tuning.
+ Tune temperatures: Perform another 2 batches of tuning.
+ Tune step sizes: Perform another 5 batches of tuning.
+ Tune temperatures: Perform another 2 batches of tuning.

After tuning the sampler is run for another #num("1000") steps to ensure the chains are in equilibrium at the final parameters after which further steps can be recorded as potential samples from the posterior distribution.

=== Validating the `nanopinv` Implementation <sec:validation>

In order to validate the implementation of these algorithms in `nanopinv`, we first tune a parallel tempering sampler as described above and fit an empirical variogram to it in a similar the manner to that of @sec:prior-dist. This variogram is then fitted to a spherical covariance function and found to have a variance $σ^2=num("1.18e-5")$, range $s=5.72$, and a nugget $τ^2=num("1.72e-6")$. The mean of the posterior samples is found to be $μ=num("0.086")$. These are then used to produce a new distribution similar to how the prior distribution was obtained and samples are drawn from it. If the implementation is correct, these samples should have similar isotropic spatial characteristics to true posterior. As such, we may pick a sample that appears physically plausible and where the velocity field is distributed in a manner that can be resolved by the rays traced by the forward model.




== Porosity and Total Water in Place <sec:physics-porosity>
In order to obtain estimates of the Total Water In Place (TWIP) over a region of interest, we first relate the phase velocity of the radar waves to the porosity of the medium, $Φ$, through the following empirical relation:
$
  Φ = exp(-41.7 v + 2.03).
$

The origin of the relation is unknown to the author and as such is simply assumed sufficiently accurate for the purposes of this report.

The TWIP may be obtained by integrating the porosity over the area of interest:
$
  TWIP
  &= ∫ Φ(vv(r)) dif vv(r)\
  &= ∑_(i=1)^(N_x) ∑_(j=1)^(N_y) Φ_(i j) Δ x Δ y\
$
where $N_x, N_y$ refer to the number of discretisation points in the $x$ and $y$ directions respectively, and $Δ x, Δ y$ are the corresponding grid spacings.

For the purposes of this report, we define the forward model to be only the solution to the eikonal equation and thus relegate the porosity and TWIP calculations to the post-inversion analysis of the obtained model.
= Results <sec:results>

= Discussion <sec:discussion>


= Conclusion <sec:conclusion>
We have employed an Iteratively Reweighted Least Squares (IRLS) algorithm with robust reweighting and Ekblom relaxation of the $L_1$-norm to solve the regularised inverse problem of obtaining the Gauss coefficients of the radial magnetic field expansion in a truncated spherical harmonics basis of degree $N=20$. The forward model relies on the assumption of interior sources, and the regularisation is applied to the projection of the model onto the Core-Mantle Boundary (CMB). This model was constructed using data gathered by the ESA Swarm mission in addition to a similar model using $L_2$-regularisation without data reweighting.

While the inversion is found to be predominantly determined by the a priori choice of regularisation for components of high spatial frequency, both models confidently produce a region of reversed flux beneath the southern Atlantic when evaluated at the CMB. The prediction error of the $L_1$-regularised model was studied under the assumption of univariate and uncorrelated data errors, and the estimated field was determined to be above a $5σ$ threshold over the entire surface with the exception of higher-frequency edges in the reconstruction.

Having considered the influence of the unmodelled ionosphere and magnetosphere, we reject the null-hypothesis of a purely dipolar geomagnetic model and find evidence supporting the existence of the South Atlantic Anomaly (SAA). We do so with reference to the obtained confidence of prediction and the ability of the proposed model to resolve features of spatial frequencies similar to that of the flux reversal region.

#pagebreak(weak: true)
= Bibliography <bib>
#bibliography("references.bib", title: none)

#pagebreak(weak: true)
#counter(heading).update(0)
#set heading(numbering: "A.1", supplement: "Appendix")
= Appendix <app>

== Student Feedback on Assignment
The assignment description is relatively sparse with details about the experiment and the experimental design. It is difficult to come up with this after the fact without knowing where the data came from.

In particular:
- The assignment incorrectly states that there are 758 sources and receivers. The correct numbers are given in @sec:experiment-data.
- The GPR carrier frequencies are not given, which makes it difficult to reason about the validity of the high-frequency assumption that arises in the derivation of the forward model through the eikonal equation. See @born-optics[3.1.1].
- The dates and location of the data collection are not available.
- Maybe assumed/obivous, but is elevation relative to sea level?
- I suspect the truncation error of the solvers may be significant for $Δ x = Δ y = #qty("0.25", "m")$

// == Large Plot of Solution to $L_2$-regularised Problem <app:l2-plot>
// #figure(
//   image("export/L2_field_map_v2.png"),
//   caption: [
//     // Large version of the plot shown in @fig:l2-solutions (b)
//   ]
// )
