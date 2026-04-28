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

#set document(title: "Report on Geomagnetic Field Reconstruction using Iterative Regularisation Techniques")
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

TODO:

The body of work presented below is solely authored by the author of this document, Jeppe Klitgaard,, and none of the presented material has been produced by means of _Generative AI_ in the manner requiring citation as described by @ref:dtu-genai outside the entries listed below.

Generative AI has, however, been used selectively in the generation of limited parts of the code. GitHub Copilot tab-completion has been used during the development of the attached code and in the generation of an initial version of the following, which is also clearly sign-posted by comments in the attached code:
- Some parts of `nanopinv` (which falls outside the examined material of this report), in particular:
  - A number of implementations of the Fast Sweeping Method
  - Boiler-plate code around benchmarking and testing these implementations

- Initial version of the plotting code for the residual plot
- Alternative L-curve corner finding algorithm: Normalised Menger Curvature
- Alternative L-curve corner finding algorithm: Chord Distance

The generated code output has in all instances been checked and often revised by the authors and is _not_ considered to be a _significant_ contribution to total intellectual work undertaken and presented.

#set heading(numbering: "1.1.a")
#pagebreak()
#counter(page).update(1)
#set page(numbering: "1/1")

= Introduction <sec:intro>
// This report seeks to determine the radial magnetic field, $B_r$, at the Core-Mantle Boundary (CMB) of the Earth using data collected as part of the three-satellite ESA Swarm mission. This may be considered as an _inverse problem_ in which the model space is spanned by the Gauss coefficients of the expansion of the projection onto the spherical harmonics.

// A derivation of the forward model is presented using the Ampère-Maxwell law and the assumption of interior sources after which the relevant theory from regularisation theory and the study of inverse problems is laid out. The inversion, which produces a model from a set of observations, is proposed to be carried out using either $L_2$ or $L_1$ regularisation applied to the projection of the model onto the CMB. The $L_1$-regularised case is solved using the proposed Iteratively Reweighted Least Squares (IRLS) algorithm with robust weights.

// Suitable candidates for the regularisation parameters, $α^2$, are determined by grid searches and the familiar L-curve heuristic after which the resulting models are used to reconstruct predictions of the radial field at the CMB. Estimates of the model uncertainties and their resolution are carried out alongside an estimated prediction confidence found using the plug-in principle and the assumption of univariate, uncorrelated noise.

// The results in the form of model projections onto the CMB and associated power spectra are then evaluated against the dipole model of the geomagnetic field, which serves as the null-hypothesis of this report. The outcomes of both the $L_2$ and robust $L_1$-regularised models are compared with particular consideration of the bias introduced by the regularisation terms.

// Lastly, the findings of the report are summarised and the conclusion is presented.

= Data and Experimental Design <sec:experiment-data>

The dataset used in this report consists of $N_d = 758$ travel time measurements taken using pairs of $N_"source" = 52$ source locations and $N_"receiver" = 235$ receiver locations in two boreholes at the Boise Hydrogeophysical Research Site (BHRS) near Boise, Idaho in the United States. The boreholes are separated by approximately 10 metres and measurements are obtained for elevations in the range $y∈[#qty("831.51", "m"), #qty("846.79", "m")]$.

During the experiment, measurements were obtained are using the technique of Cross Hole Ground Penetrating Radar (GPR), in which the delays in the observed travel times of electromagnetic radar waves emitted by the source and received by its corresponding network of receivers were recorded for each source-receiver pair.

The water table at the site is known to be located at an elevation of $y_* = #qty("844", "m")$.

= Theory <sec:theory>
In order to arrive at a suitable model and framework for obtaining it, we will first lay out a brief summary of the relevant theory from the fields of physics and probabilistic inversion.

== Physics and Forward Model <sec:theory-physics>
The data may then be used to formulate a classic tomographic inverse problem in which the spatial distribution of the electromagnetic phase velocity is reconstructed from the observed travel times through the medium between the boreholes.

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
This is a useful metric to keep in mind, as it also provides a useful scale of discretisation for the model space. In particular, no additional resolution can be obtained by further refining the model discretisation beyond this limit without also increasing the frequency of the radar waves or resorting to a more complex physical model that relaxes the high-frequency assumption.

=== Porosity and Total Water in Place <sec:physics-porosity>

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

== Numerical Solution of the Forward Problem <sec:theory-eikonal-solution>
In order to compute the forward projection of a given model, here understood to be a discretised spatial distribution of phase velocities, we must solve @eq:eikonal-travel-time, and ideally do so as efficiently as possible.

A popular choice of algorithm dedicated to solving the family of equations to which the eikonal equation belongs is the Fast Marching Method (FMM) due to Sethian @bib:fmm, which are implemented as first and second order methods in @bib:skfmm.
The finer details of the FFM algorithm is beyond the scope of this report, but a notable consequence of the formulation of the algorithm is that it is not amenable to parallelisation. While multiple source-receiver pairs may be computed in parallel using CPU cores, the method is not able utilise modern, massively-parallel hardware such as GPUs efficiently.

For this reason, the alternative Fast Sweeping Method (FSM) due to Zhao @bib:fsm-zhao @bib:fsm-zhao-parallel is considered a superior choice for solving the eikonal equation on GPUs, as it more naturally lends itself to parallelisation. In particular, have used the paper @bib:fsm-2nd-order as a reference for the implementation of both first and second order parallel FSM algorithms, which are now available as part of the `nanopinv` library developed concurrently with the writing of this report.

Common to all the implementations is that they produce an approximation of the travel time $T$ from a source $vv(r)_s$ to all other points in the model space when given a propagation velocity field $v(vv(r))$.
The estimated travel time to the receivers, $vv(r)_r$, may then be obtained by interpolation.
== Probabilistic Inversion <sec:theory-probinv>
In order to solve the inverse problem of estimating the velocity distribution from the observed travel times, we reach for the framework of _probabilistic inversion_, in which we adopt the notation due to Tarantola @bib:tarantola2005. This choice of framework is a convenient for travel-time tomography problems, as it handles non-linearities in the forward model naturally, while also offering a principled way to incorporate prior information and uncertainty quantification into the inversion.

We may tersely summarise the framework as presented in @ref:course_book[4.2] by letting $ρ_m (vv(m)) : cal(M) → ℝ$ and $ρ_d (vv(d)) : ℝ^(N_d) → ℝ$ denote the _prior model probability distribution_ and _prior data probability distribution_ respectively, with their arguments denoting the model parameters and data. For a 2-dimensional model space, $cal(M) = ℝ^(N_x × N_y)$. By assuming the data to be independent of the model parameters, we construct a joint prior distribution:
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

=== Probablistic Model

The Bayes-like formula in @eq:posterior may be re-expressed for the inverse problem at hand by identifying terms with those described in @sec:theory-physics:
$
  σ_v (vv(v)) = L(vv(v)) ρ_v (vv(v))
$
$
  L(vv(v)) = ∫_cal(D) (ρ_d (vv(d)) θ(vv(d)|vv(v)))/(μ_d (vv(d))) dif vv(d)
$

TODO THIS SUCKS
TODO: Map symbols onto physics symbols

== Markov Chain Monte Carlo Methods <sec:theory-mcmc>

In order to sample the posterior, $σ_m (vv(m))$, we turn to the theory of Markov Chains and Monte Carlo sampling, in which chain of samples are generated by a Markov process — that is, a process in which the next state depends only on the current state.

While many algorithms for simulating such processes exist, we will focus on the _Extended Metropolis algorithm_ due to Mosegaard and Tarantola @bib:extended-metropolis, which is combined with a Preconditioned Crank-Nicolson proposal distribution in order to efficiently sample the

TODO Proposal algorithm

which is summarised in @alg:extended-metropolis


#pagebreak(weak: true)
#algorithm-figure(
  "Extended Metropolis Algorithm",
  vstroke: .5pt + luma(200),
  inset: 0.35em,
  {
    import algorithmic: *
    let Solve = Call.with("Solve")
    let Input(doc) = Line([*input*~#context box(doc, baseline: 100% - par.leading)])
    let Output(doc) = Line([*output*~#context box(doc, baseline: 100% - par.leading)])
    Procedure(
      "ExtendedMetropolis",
      ([$mm(G)$], [$vv(d)$], [$mm(H)_k$], [$α_k$], [$δ_k$], [$p_k$], [$cal(V)$], [$cal(S)$], [$N_"iter"$], [$τ$]),
      {
        Input([
          Model $mm(G) ∈ ℝ^(N_m × N_d)$, data $vv(d) ∈ ℝ^(N_d)$,\
          Regularisers $mm(H)_(k∈cal(K)) ∈ ℝ^(N_k × N_m)$, reg. parameters $α_(k ∈ cal(K))$, norm orders $p_(k∈cal(K))$,\
          Ekblom perturbations $δ_(k∈cal(K))$,\
          M-estimator $cal(S)(vv(ε); vv(θ)_cal(S))$, scale estimator $cal(V)(vv(ε); vv(θ)_cal(V))$,\
          Maximal iterations $N_"iter"$, tolerance $τ$
        ])
        Comment[$cal(K) = {1, …, K} ⊂ ℕ$: regularisation term indices]
        Output([
          Model parameters $vv(m) ∈ ℝ^(N_m)$,\
          Robust weights $vv(w) ∈ ℝ^(N_d)$, Ekblom measures $vv(γ)_(k∈cal(K))$
        ])
        Comment[Initialisation]
        Assign[$vv(m)$][$bb(0)_(N_m)$]
        Assign[$vv(w)$][$bb(1)_(N_d)$]
        Assign[$vv(γ)_k$][$bb(0)_(N_k) quad ∀ k ∈ cal(K)$]

        For(
          $i = 1 " to " N_"iter"$,
          {
            Comment[Construct system matrix and right-hand side, $mm(A) vv(m)_"new" = vv(b)$]
            LineComment(
              Assign[$vv(w)_(m k)$][$p_k (vv(γ)_k^2 + δ_k^2)^(p_k\/2 -1) quad ∀ k$], [Compute Ekblom norm weights]
            )
            LineComment(
              Assign[$mm(A)$][$mm(G)^TT diag(vv(w)) mm(G) + ∑_(k=1)^K α_k^2 mm(H)_k^TT diag(vv(w)_(m k)) mm(H)_k$],
              [System matrix, $mm(A)$]
            )
            LineComment(
              Assign[$vv(b)$][$mm(G)^TT (vv(w) dot.o vv(d))$],
              [RHS, $vv(b)$; $dot.o$ is element-wise mult.]
            )
            LineComment(
              Assign($vv(m)_"new"$, Solve[$mm(A)$, $vv(b)$]),
              [$mm(A) succ 0, mm(A)^HH = mm(A), $ ∴ use Cholesky]
            )
            LineComment(
              Assign[$Δ_2$][$norm(vv(m) - vv(m)_"new")_2 \/ norm(vv(m))_2$],
              [Convergence criterion]
            )

            // Update model
            Comment[Update model]
            Assign[$vv(m)$][$vv(m)_"new"$]

            Comment[Update weights]
            LineComment(
              Assign[$vv(ε)$][$vv(vv(d) - mm(G) vv(m))$],
              [Compute unweighted residual]
            )

            LineComment(
              Assign[$hat(σ)$][$cal(V)(vv(ε))$],
              [Estimate data dispersion scale]
            )

            LineComment(
              Assign[$vv(w)$][$cal(S)(vv(ε)\/hat(σ))$],
              [Compute robust weights]
            )

            LineComment(
              Assign[$vv(γ)_k$][$mm(H)_k vv(m) quad ∀ k$],
              [Compute Ekblom measures]
            )

            // Convergence check
            LineComment(
                If($Δ_2 < τ$, {
                Break
              }),
              [Stop if converged]
            )
          }
        )
        Return[$vv(m), vv(w), vv(γ)_(k∈cal(K))$]
      }
    )
  }
) <alg:extended-metropolis>

== Parallel Tempering

== Hyperparameter Tuning

---

// Before an inverse problem can be solved, we must present the corresponding forward model of the system, which may be formulated in the linear, discretised case as
// $
//   vvh(d) = vv(d) + vv(e) = mm(G) vv(m) + vv(e),
// $ <eq:forward-model>

// where $mm(G) : cal(M) → cal(D)$ is the _system matrix_ which maps model parameters $vv(m)$ from the discrete model space, $cal(M) = ℝ^(N_m)$, to a vector $vv(d) = mm(G) vv(m)$ in the observation space, $cal(D) = ℝ^(N_d)$. Here, $vv(d)$ is understood to be the _deterministic_ component of the forward projection, and is colloquially — and somewhat inaccurately — referred to as the _noise-free data_. Naturally, we do not observe the forward projection of some idealised physical model, but rather observations $vvh(d)$ which also contain contributions from both the realisations of influencing stochastic processes and any systematic biases, both of which are captured by $vv(e)$. As we will later discover, $vv(e)$ cannot reasonably be assumed to be normally distributed.

// In order to determine the system matrix $mm(G)$, we must consider carefully the physical model of the system. To this end, we start from the Ampère-Maxwell law and the assumption that the region of measurement — the orbits upon which the radial component $B_r$ is measured using the _Swarm_ satellites — is devoid of free charges and currents. That is, we assert that the magnetic field is curl-free:
// $
//   ∇ × vv(B) = μ_0 (vv(J) + ε_0 pdv(vv(E), t)) = vv(0).
// $

// Thus, we may define a _scalar magnetic potential_, $V$, such that
// $
//   vv(B) = -∇V.
// $ <eq:scalar_potential>

// Deeming it unlikely that any magnetic monopoles will present themselves, we may conclude that the field is also divergence-free, $∇ ⋅ vv(B) = 0$, which yields the Laplace equation for the scalar potential:
// $
//   ∇ ⋅ vv(B) = -∇^2 V = 0.
// $

// If we make the simplifying assumption that a spherical boundary surface of radius $a$ exists which delineates the electromagnetic sources from empty space, we find that the spatial domain of our model is the sphere, $S^2$. While such a spherical boundary is reasonable for the case of the Earth, assuming that no sources — here understood to be magnetic moments and electrical currents — exist outside of it is more difficult to justify. In particular, this neglects the charge carriers in the ionosphere and the influence of the magnetosphere, a discussion which we will return to later.

// With this in mind, we may choose the spherical harmonics $Y_n^m (ϑ, λ)$ as a convenient, orthogonal basis for the model space and continue our analysis in spherical coordinates $(r, ϑ, λ)$, in which we recall Laplace's equation to be of the form:
// $
//   ∇^2 V = 1/r^2 pdv(,r)(r^2 pdv(V, r)) + 1/(r^2 sin ϑ) pdv(, ϑ)(sin ϑ pdv(V, ϑ)) + 1/(r^2 sin^2 ϑ) pdv(V, λ, 2) = 0.
// $

// Separation of variables and the assumption of interior sources yield the solution @ref:geomath[p. 507, eq. 13]
// $
//   V(r, ϑ, λ)
//   &=a ∑_(n=1)^∞ (a/r)^(n+1) ∑_(m=0)^n (g_n^m cos m λ + h_n^m sin m λ ) P_n^m (cos ϑ)\
// $ <eq:lapl_sol>

// in which $P_n^m$ are the associated Legendre polynomials and we have used the spherical harmonics as given by:
// $
// Y_(n,"even")^m = cos (m λ) P_n^m (cos ϑ)
// wider
// Y_(n,"odd")^m = sin (m λ) P_n^m (cos ϑ).
// $

// The parameters $g_n^m, h_n^m$ of @eq:lapl_sol are referred to as the _Gauss coefficients_ and the inverse problem is thus "simply" to determine these given the field measurements.

// We do not, however, have access to the potential $V$ directly, which is introduced more as a mathematical convenience than a physical quantity. Instead, we measure the radial component of its gradient as given by @eq:scalar_potential,
// $
//   B_r (r, ϑ, λ)
//   = - pdv(V, r)
//   = ∑_(n=1)^∞ (n+1) (a/r)^(n+2) ∑_(m=0)^n (g_n^m cos m λ + h_n^m sin m λ ) P_n^m (cos ϑ).
// $ <eq:B_r>

// By identifying this as forward projection and truncating the spherical harmonic expansion to degree $N$, we obtain the semi-discrete model in terms of coordinates:
// $
//   B_r (r, ϑ, λ) = ∑_(n=1)^N (n+1) (a/r)^(n+2) ∑_(m=0)^n (g_n^m cos m λ + h_n^m sin m λ ) P_n^m (cos ϑ),
// $ <eq:semi-discrete-model>
// while the finite-dimensional model parameters are consequently given by
// $
// vv(m) = (g_1^0, g_1^1, h_1^1, g_2^0, ..., g_N^N, h_N^N)^T  ∈ cal(M)
// $
// where notably the $h_n^0$ coefficients are omitted as their terms vanish due to the $sin(m λ)$ factor in @eq:lapl_sol. This truncation is sensible under the assumption that the variations of the magnetic fields are smooth and that the features of interest are of sufficiently low frequency to be captured by the remaining basis functions.

// By constructing a system of equations using @eq:semi-discrete-model and the associated discrete model parameters $vv(m)$ over an idealised, noise-free set $D = {d_i = B_r (r_i, ϑ_i, λ_i)}_(i=1)^(N_d)$, we produce the discrete forward operator as the system matrix, $mm(G)$, which is built as follows:
// $
//   vv(d) = mm(G) vv(m)\
//   vec(
//     B_r (vv(x)_1),
//     B_r (vv(x)_2),
//     ⋮,
//     B_r (vv(x)_(N_d)),
//   )
//   =
//   mat(
//     // g=cos, n=1, m=0
//     Ξ_1^0 (vv(x)_1),
//     // h=sin, n=1, m=0, this vanishes
//     // g=cos, n=1, m=1
//     Ξ_1^1 (vv(x)_1),
//     // h=sin, n=1, m=1
//     ξ_1^1 (vv(x)_1),
//     …,
//     // g=cos, n=N, m=N
//     Ξ_N^N (vv(x)_1),
//     // h=sin, n=N, m=N
//     ξ_N^N (vv(x)_1)
//     ;
//     Ξ_1^0 (vv(x)_2),
//     Ξ_1^1 (vv(x)_2),
//     ξ_1^1 (vv(x)_2),
//     …,
//     Ξ_N^N (vv(x)_2),
//     ξ_N^N (vv(x)_2)
//     ;
//     ⋮,⋮,⋮,⋮,⋮,⋮
//     ;
//     Ξ_1^0 (vv(x)_(N_d)),
//     Ξ_1^1 (vv(x)_(N_d)),
//     ξ_1^1 (vv(x)_(N_d)),
//     …,
//     Ξ_N^N (vv(x)_(N_d)),
//     ξ_N^N (vv(x)_(N_d))
//     ;

//   )
//   vec(
//     g_1^0,
//     g_1^1,
//     h_1^1,
//     ⋮,
//     g_N^N,
//     h_N^N,
//   )
// $ <eq:forward-model-big>

// where we have introduced some shorthand notation
// $
// vv(x)_i &= (r_i, ϑ_i, λ_i)\
// Ξ_n^m (r, ϑ, λ) &= (n+1) (a/r)^(n+2) cos (m λ) P_n^m (cos ϑ)\
// ξ_n^m (r, ϑ, λ) &= (n+1) (a/r)^(n+2) sin (m λ) P_n^m (cos ϑ).
// $

// In conclusion, we can build the discrete forward operator $mm(G) ∈ ℝ^(N_d × N_m), N_m = N(N+2)$ simply by evaluating the basis functions $Ξ_n^m$ and $ξ_n^m$ at the measurement locations, ${vv(x)_i}_(i=1)^(N_d)$. We note that a final assumption made is that the temporal variations in the field are locally much slower than the period over which the data has been collected, or equivalently that $B_r (r, ϑ, λ)$ is constant with respect to time.

// = Theory — Mathematical Treatment of the Inverse Problem <sec:theory-math>

// Having now produced a discrete forward model that may be used with our experimental data, we turn our attention to the associated _inverse problem_, that is, the estimation of model parameters $vv(m)$ given the observations $hat(D) = {hat(d)_i = B_r (r_i, ϑ_i, λ_i) + e_i}_(i=1)^(N_d)$.

// == Least Squares Solution <sec:math-ls>
// Naïvely, we may obtain a solution to the inverse problem by the computing the inverse of $mm(G)$ with which we would left-multiply to obtain an estimate
// $
// vvh(m) = mm(G)^(-1) vvh(d),
// $
// except of course such an inversion does not exist since $mm(G)$ is not generally square. Instead, we seek the model parameters which minimise the misfit, $ϕ$, between the observed data $vvh(d)$ and the forward projection in a norm, here taken to be the $L_2$-norm:
// $
// vvh(m)
// = arg min_vv(m) ϕ =
// arg min_vv(m) 1/2 norm(mm(G) vv(m) - vvh(d))_2^2
// $

// By inspection, the misfit $ϕ$ is convex, differentiable with respect to $vv(m)$, and coercive provided $mm(G)$ is full-rank, and thus the (global) minimum may be found by solving
// $
// ∇_vv(m) ϕ = 0 &= ∇_vv(m) [1/2 (mm(G) vv(m) - vvh(d))^TT (mm(G) vv(m) - vvh(d))]\
// &= 1/2 ∇_vv(m) [-vv(m)^TT mm(G)^TT vvh(d) + vv(m)^TT mm(G)^TT mm(G) vv(m) - vvh(d)^TT mm(G) vv(m) + vvh(d)^TT vvh(d)]\
// &= 1/2 ∇_vv(m) [-2 vvh(d)^TT mm(G) vv(m) + vv(m)^TT mm(G)^TT mm(G) vv(m)] = - vvh(d)^TT mm(G) + vv(m)^TT mm(G)^TT mm(G) ⇔ mm(G)^TT mm(G) vvh(m) = mm(G)^TT vvh(d)\
// vvh(m) &= (mm(G)^TT mm(G))^(-1) mm(G)^TT vvh(d) = (mm(G)^TT mm(G))^(-1) mm(G)^TT (vv(d) + vv(e)) = (mm(G)^TT mm(G))^(-1) mm(G)^TT (mm(G) vv(m) + vv(e)) = vv(m) + (mm(G)^TT mm(G))^(-1) mm(G)^TT vv(e)\
// &= vv(m) + mm(G)^† vv(e)
// $ <eq:m-hat-ls>

// Clearly, our estimate will be rather poor unless $norm(mm(G)^† vv(e)) ≪ norm(vv(m))$. We can more clearly see how errors may enter our estimate if we consider using the constructed model to compute the radial magnetic field at the Earth's Core-Mantle Boundary (CMB).

// == Stability Analysis <sec:math-illposedness>
// Since the model is fully determined by the expansion coefficients in the angular coordinates, the projection onto the CMB using the obtained model may be carried out simply by constructing a new system matrix $mm(G)_"CMB"$ as before, but replacing the ${r_i}_(i=1)^(N_d)$ radial coordinates with the radius of the CMB, $r_"CMB" = qty("3480.0", "km")$. By inspection of @eq:forward-model-big, this will yield an estimate of the field at the boundary, $vvh(d)_"CMB"$:
// $
// vvh(d)_"CMB" = mm(G)_"CMB" vvh(m) = mm(G)_"CMB" vv(m) + mm(G)_"CMB" mm(G)^† vv(e) = vv(d)_"CMB" + vvh(e)
// $

// Where $vvh(e) = mm(G)_"CMB" mm(G)^† vv(e)$ is the error of our prediction $vvh(d)_"CMB"$ given the realisation $vv(e)$. If we consider a specific spherical degree $n$ and @eq:semi-discrete-model, we find that the radial dependence factors out as $(n+1) (a\/r)^(n+2)$ which by identity $mm(G)^† mm(G) = mm(I)$ implies that $mm(G)^†$ must scale by the reciprocal, $(n+1)^(-1) (a\/r)^(-n-2)$. If we model all the measurements as having been taken at a single, mean radius, $r_"sat"$, we thus find
// $
// vvh(e)_n ∝ (n+1)^(-1) (a/r_"sat")^(-n-2) (n+1) (a/r_"CMB")^(n+2) vv(e) = (r_"sat" / r_"CMB")^(n+2) vv(e)
// $

// In our data, the average radius is $r_"sat" = qty("6891.51", "km")$, giving a spectral noise amplification factor $η_n=1.9803^(n+2), η_20 ≈ num("3 370 000")$. Clearly, even minute errors in the high-frequency components of the measurements may dominate our estimate at the core.

// This sensitivity to small changes causes the inverse problem to be severely _ill-posed_ in the sense of Hadamard and requires us to regularise our inversion.

// == Regularisation Theory <sec:math-regularisation>
// A common regularisation scheme is to introduce a penalty on the norm of model parameters, or some linear transformation thereof. In a general case, such an $L_p$-regularised solution may be found using the variational formulation
// $
// arg min_vv(m) 1/2 norm(mm(G) vv(m) - vvh(d))_mm(W)^2 + α^2/2 norm(mm(H) vv(m))_p^p,
// $ <eq:robust-variational>

// where $p$ is the _norm order_, $α^2$ is the _regularisation parameter_, and $mm(H) ∈ ℝ^(N_k × N_m)$ is a linear regularisation operator, taken here to be the CMB projection, $mm(G)_"CMB"$. The norm $norm(vv(v))_mm(W)^2 ≔ vv(v)^TT mm(W) vv(v)$ indicates the _weighted_ $L_2$ norm where $mm(W)$ is a diagonal matrix of _robust weights_ that must be determined jointly in the absence of an _a priori_ estimate.

// To ensure differentiability at zero, the $L_p$ norm may be approximated using the relaxation scheme due to Ekblom:
// $
// norm(mm(H) vv(m))_p^p ≈ sum_(i=1)^(N_k) [((mm(H) vv(m))_i)^2 + δ^2]^(p\/2),
// $
// which yields an approximate solution to @eq:robust-variational of the form @ref:course_book[eq. 51]
// $
// vvh(m)_(r α p) = (mm(G)^TT mm(W) mm(G) + α^2 mm(H)^TT mm(W)_(m p) mm(H))^(-1) mm(G)^TT mm(W) vvh(d) = tilde(mm(G))_(α p)^(-1) mm(G)^TT mm(W) vvh(d) = mm(G)_(α p)^(-1) vvh(d).
// $

// Here, both the robust data weights $mm(W)$ and the (diagonal) regularisation weights may be obtained through the Iteratively Reweighted Least Squares (IRLS) algorithm, as described by @alg:irls-r. Specifically, the diagonal of $mm(W)_(m p)$ is populated by the vector $vv(w)_m$, which is updated element-wise at each iteration by
// $
// vv(w)_m = p (vv(γ)^2 + δ^2)^(p\/2 - 1), quad "where" quad vv(γ) = mm(H) vv(m),
// $
// with $δ$ denoting the _Ekblom perturbation_ which must satisfy $δ^2 ≪ norm(vv(γ))_2^2 \/ N_k.$

== Robust Iteratively Reweighted Least Squares Algorithm <sec:theory-irls>
We propose a version of the standard Iteratively Reweighted Least Squares algorithm which features a number of additional considerations. While we will only use a single regularisation term in this report, the proposal in @alg:irls-r allows an arbitrary number $k$ of regularisation terms to be employed. It also allows for different M-estimators and scale estimators to be used for the computation of the robust weights and the dispersion scale $hat(σ)$ respectively.

For the purposes of this report, we choose the Huber weights @ref:course_book[p.9] and normalised Median Absolute Deviation (MAD) @ref:robust-stats[p.36] as the M-estimator and scale estimator respectively, which are given by $cal(S)(vv(ε);vv(θ)_cal(S))$ and $cal(V)(vv(ε);vv(θ)_cal(V))$ in @alg:irls-r:
$
w_"Huber" (ε; c=1.345) = cases(
  c \/ abs(ε) &wide abs(ε) > c,
  1 &wide abs(ε) ≤ c,
)
$ <eq:huber-weights>

$
  hat(σ)_"MAD" (vv(ε)) = median(sum_(i=1)^(N_d) abs(ε_i - median(vv(ε))) ) / (Φ^(-1)(0.75))
$
where $Φ(p)$ is the CDF of a standard Gaussian distribution such that $Φ^(-1)(0.75) ≈ 0.675$ @ref:robust-stats[p.36].

As a termination criterion for the algorithm, the relative size of the model parameter update in the $L_2$-norm is compared against the given tolerance $τ$ at each time step to enable stopping before the iteration limit $N_"iter"$.

// The implementation of @alg:irls-r is provided in the attached code and may also be found in @app:irls. This features some additional parameters pertaining to debugging and observability into the performance of the algorithm.

#pagebreak(weak: true)
#algorithm-figure(
  "Robust Iteratively Reweighted Least Squares (IRLS-R)",
  vstroke: .5pt + luma(200),
  inset: 0.35em,
  {
    import algorithmic: *
    let Solve = Call.with("Solve")
    let Input(doc) = Line([*input*~#context box(doc, baseline: 100% - par.leading)])
    let Output(doc) = Line([*output*~#context box(doc, baseline: 100% - par.leading)])
    Procedure(
      "IRLS-R",
      ([$mm(G)$], [$vv(d)$], [$mm(H)_k$], [$α_k$], [$δ_k$], [$p_k$], [$cal(V)$], [$cal(S)$], [$N_"iter"$], [$τ$]),
      {
        Input([
          Model $mm(G) ∈ ℝ^(N_m × N_d)$, data $vv(d) ∈ ℝ^(N_d)$,\
          Regularisers $mm(H)_(k∈cal(K)) ∈ ℝ^(N_k × N_m)$, reg. parameters $α_(k ∈ cal(K))$, norm orders $p_(k∈cal(K))$,\
          Ekblom perturbations $δ_(k∈cal(K))$,\
          M-estimator $cal(S)(vv(ε); vv(θ)_cal(S))$, scale estimator $cal(V)(vv(ε); vv(θ)_cal(V))$,\
          Maximal iterations $N_"iter"$, tolerance $τ$
        ])
        Comment[$cal(K) = {1, …, K} ⊂ ℕ$: regularisation term indices]
        Output([
          Model parameters $vv(m) ∈ ℝ^(N_m)$,\
          Robust weights $vv(w) ∈ ℝ^(N_d)$, Ekblom measures $vv(γ)_(k∈cal(K))$
        ])
        Comment[Initialisation]
        Assign[$vv(m)$][$bb(0)_(N_m)$]
        Assign[$vv(w)$][$bb(1)_(N_d)$]
        Assign[$vv(γ)_k$][$bb(0)_(N_k) quad ∀ k ∈ cal(K)$]

        For(
          $i = 1 " to " N_"iter"$,
          {
            Comment[Construct system matrix and right-hand side, $mm(A) vv(m)_"new" = vv(b)$]
            LineComment(
              Assign[$vv(w)_(m k)$][$p_k (vv(γ)_k^2 + δ_k^2)^(p_k\/2 -1) quad ∀ k$], [Compute Ekblom norm weights]
            )
            LineComment(
              Assign[$mm(A)$][$mm(G)^TT diag(vv(w)) mm(G) + ∑_(k=1)^K α_k^2 mm(H)_k^TT diag(vv(w)_(m k)) mm(H)_k$],
              [System matrix, $mm(A)$]
            )
            LineComment(
              Assign[$vv(b)$][$mm(G)^TT (vv(w) dot.o vv(d))$],
              [RHS, $vv(b)$; $dot.o$ is element-wise mult.]
            )
            LineComment(
              Assign($vv(m)_"new"$, Solve[$mm(A)$, $vv(b)$]),
              [$mm(A) succ 0, mm(A)^HH = mm(A), $ ∴ use Cholesky]
            )
            LineComment(
              Assign[$Δ_2$][$norm(vv(m) - vv(m)_"new")_2 \/ norm(vv(m))_2$],
              [Convergence criterion]
            )

            // Update model
            Comment[Update model]
            Assign[$vv(m)$][$vv(m)_"new"$]

            Comment[Update weights]
            LineComment(
              Assign[$vv(ε)$][$vv(vv(d) - mm(G) vv(m))$],
              [Compute unweighted residual]
            )

            LineComment(
              Assign[$hat(σ)$][$cal(V)(vv(ε))$],
              [Estimate data dispersion scale]
            )

            LineComment(
              Assign[$vv(w)$][$cal(S)(vv(ε)\/hat(σ))$],
              [Compute robust weights]
            )

            LineComment(
              Assign[$vv(γ)_k$][$mm(H)_k vv(m) quad ∀ k$],
              [Compute Ekblom measures]
            )

            // Convergence check
            LineComment(
                If($Δ_2 < τ$, {
                Break
              }),
              [Stop if converged]
            )
          }
        )
        Return[$vv(m), vv(w), vv(γ)_(k∈cal(K))$]
      }
    )
  }
) <alg:irls-r>

== Covariance and Resolution Matrices <sec:math-cov-res-matrix>
We wish to be able to compute estimates of the covariance and resolution matrices of our model, which will give insight into respectively the uncertainty and bias of our obtained model, $vvh(m)$. Additionally, we may propagate the model covariance onto the prediction to produce estimates of the uncertainty of the obtained radial magnetic field.

We first note that $tilde(mm(G))_(α p)^(-1) = (mm(G)^TT mm(W) mm(G) + α^2 mm(H)^TT mm(W)_(m p) mm(H))^(-1)$ is symmetric by the diagonality of $mm(W), mm(W)_(m p)$.
Thus, by identity $Cov(mm(A) vv(b)) = mm(A) Cov(vv(b)) mm(A)^TT$, we obtain an estimated model covariance matrix
$
mm(Σ)_vvh(m) = Cov(vvh(m)_(r α p))
&≈ mm(G)_(α p)^(-1) Cov(vvh(d)) (mm(G)_(α p)^(-1))^TT\
&= tilde(mm(G))_(α p)^(-1) mm(G)^TT mm(W) Cov(vvh(d)) (tilde(mm(G))_(α p)^(-1) mm(G)^TT mm(W))^TT\
&= tilde(mm(G))_(α p)^(-1) mm(G)^TT mm(W) Cov(vvh(d)) mm(W) mm(G) tilde(mm(G))_(α p)^(-1).
$ <eq:cov-matrix>

Here $mm(W)$ denotes the _frozen_ robust weights obtained at the final iteration of @alg:irls-r, which provide a reasonable, linear approximation in this case.

The variance of the measurements are not known a priori and we resort to a _plug-in_ estimate in which we assume the data to be uncorrelated and homoscedastic with a variance approximated by the Median Absolute Deviation of the residuals, which has been jointly estimated as part of @alg:irls-r:
$
Cov(vvh(d)) ≈ [hat(σ)_"MAD" (vvh(d) - mm(G) vvh(m))]^2 ⋅ mm(I_(N_d)).
$

Similarly, by applying the treatment used to obtain @ref:course_book[eq. 35], the resolution matrix, $mm(R)_vvh(m)$, may be computed using
$
mm(R)_vvh(m) ≈ mm(G)_(α p)^(-1) mm(G) = tilde(mm(G))_(α p)^(-1) mm(G)^TT mm(W) mm(G)
$ <eq:res-matrix>
where we recall that the limiting case $mm(R)_vvh(m) = mm(I)$ identifies the optimal case where the data is fully resolved by the inversion.

Finally, we may estimate the covariance matrix of our reconstruction at the Core-Mantle Boundary, $vvh(d)_"CMB" = mm(G)_"CMB" vvh(m)$ as
$
mm(Σ)_(vvh(d)_"CMB")
&≈ mm(G)_"CMB" mm(Σ)_vvh(m) mm(G)_"CMB"^TT\
&= mm(G)_"CMB" tilde(mm(G))_(α p)^(-1) mm(G)^TT mm(W) Cov(vvh(d)) mm(W) mm(G) tilde(mm(G))_(α p)^(-1) mm(G)_"CMB"^TT\
&= [hat(σ)_"MAD" (vvh(d) - mm(G) vvh(m))]^2 mm(G)_"CMB" tilde(mm(G))_(α p)^(-1) mm(G)^TT mm(W)^2 mm(G) tilde(mm(G))_(α p)^(-1) mm(G)_"CMB"^TT.
$

where the scalar data covariance and diagonal weight matrix are exploited in the last step to significantly reduce the computational burden.

// #pagebreak(weak: true)
// = Results
// We will investigate the performance of three different inversions starting out with the unregularised least-squares solution, which can be expected to perform poorly due to the illposedness of the problem. To promote stability in the inversion, we subsequently present and analyse the results of both regularisation in both the $L_2$ and $L_1$ norms.

// == Unregularised Least Squares Solution <sec:results-ls>

// #let fig-residual-ls = [
//   #box(width: 50%)[
//     #figure(
//       image("export/residual_ls.png"),
//       caption: [
//         Residual of naïve least squares solution.
//       ]
//     ) <fig:residual-ls>
//   ]
// ]

// #let fig-sol-ls-cmb = [
//   #box(width: 50%)[
//     #figure(
//       image("export/sol_ls_cmb.png"),
//       caption: [
//         Bad reconstruction of radial magnetic field at CMB due to unstable inversion.
//       ]
//     ) <fig:sol-ls-cmb>
//   ]
// ]

// #meander.reflow({
//   import meander: *

//   placed(top + right, stack(fig-residual-ls, v(1em), fig-sol-ls-cmb))
//   container()
//   content([
// If we disregard the analysis carried out in @sec:theory-math and attempt to construct a model by inversion using the naïve, unregularised least-squares solution given by @eq:m-hat-ls, we recover the reconstruction of the radial magnetic field at the CMB shown in @fig:sol-ls-cmb.

// Immediately, we note the rapid oscillations exhibited by the field towards the poles, which does not agree with our expectations. Further, the magnitude of the field oscillations is much larger than that found in the literature @ref:hammer, which prescribes an amplitude of approximately #qty("1", "mT"). This is a consequence of the aforementioned ill-posedness of the inverse problem. Further, inspection of the residuals $vv(ε) = vvh(d) - mm(G) vvh(m)$ as shown in @fig:residual-ls reveals that these are far from normally distributed and instead appear to follow a long-tailed distribution such as the Huber distribution, which has been fitted to the data alongside the Laplace distribution.

// It is for this reason that we have chosen the Huber weights given in @eq:huber-weights as the M-estimator for our IRLS algorithm, which has the effect of reducing the influence of observations further from the centre of the data distribution.
// ])})

// // #pagebreak(weak: true)
// == Non-Robust $L_2$-Regularised Solution
// We consider the specific case of a non-robust $L_2$-regularised solution, which by the treatment in @sec:theory-math and identifying
// $
// p = 2
// , wider
// mm(W) = mm(W)_"m2" = mm(I)
// $

// gives
// $
// vv(m)_(α 2)
// &= (mm(G)^TT mm(W) mm(G) + α^2 mm(H)^TT mm(W)_"m2" mm(H))^(-1) mm(G)^TT mm(W) vv(d)\
// &= (mm(G)^TT mm(G) + α^2 mm(H)^TT mm(H))^(-1) mm(G)^TT vv(d) = tilde(mm(G))_(α 2)^(-1) mm(G)^TT vv(d) = mm(G)_(α 2)^(-1) vv(d)\
// $

// where $mm(G)_(α 2)^(-1) = (mm(G)^TT mm(G) + α^2 mm(H)^TT mm(H))^(-1) mm(G)^TT$ may be used to find the model covariance and resolution matrices through @eq:cov-matrix and @eq:res-matrix respectively.
// Notably, since neither of the iteratively constructed weight matrices feature, this regularization can be carried out without resorting to @alg:irls-r.

// #pagebreak(weak: true)
// #let fig-L2-curve = [
//   #box(width: 50%)[
//     #figure(
//       image("export/L_curve_L2.png"),
//       caption: [L-curve for the $L_2$ least-squares solution]
//     ) <fig:l2-curve>
//   ]
// ]
// #let fig-l2-field-map = [
//   #box(width: 50%)[
//     #figure(
//       image("export/L2_field_map.png"),
//       caption: [Projection of $L_2$-regularised solution onto CMB]
//     ) <fig:l2-field-map>
//   ]
// ]
// #let fig-l2-residual = [
//   #box(width: 50%)[
//     #figure(
//       image("export/L2_residual.png"),
//       caption: [Distribution of residuals for the $L_2$-regularised model]
//     ) <fig:l2-residual>
//   ]
// ]
// #let fig-l2-different_power_spectra = [
//   #box(width: 50%)[
//     #figure(
//       image("export/L2_different_power_spectra.png"),
//       caption: [Power spectra for different $L_2$ regularization strengths]
//     ) <fig:l2-different-power-spectra>
//   ]
// ]
// #meander.reflow({
//     import meander: *
//   placed(top + right, stack(fig-L2-curve))
//   container()
//   content([
// To select an appropriate regularization parameter, the L-curve method is applied. @fig:l2-curve shows the L-curve as a scatter plot, where the corner identifies the optimal trade-off between data misfit and model norm. The corresponding regularization parameter is
// $
// alpha^2 = 1.60 times 10^(-8).
// $
// The non-robust $L_2$-regularized solution is obtained from the closed-form expression derived above, which arises from minimizing a quadratic objective function.

// This value provides a reasonable balance between fitting the observations and suppressing instability that would otherwise produce dominant and noisy high-frequency spatial oscillations.
// ])})

// #subpar.grid(
//   align: top,
//   figure(
//     image("export/L2_comparison.png"),
//     caption: [
//       Under-regularized, $alpha^2 = 5 times 10^(-11)$.
//     ]
//   ),
//   figure(
//     image("export/L2_comparison_2.png"),
//     caption: [
//       Optimal regularized, $alpha_*^2 = 1.60 times 10^(-8)$
//     ]
//   ),
//   figure(
//     image("export/L2_comparison_3.png"),
//     caption: [
//       Over-regularized, $alpha^2 = 5 times 10^(-7)$
//     ]
//   ),


//   columns: (1fr, 1fr, 1fr),
//   caption: [
//     Solutions to the non-robust $L_2$-regularized problem for various regularization parameters.
//   ],
//   label: <fig:l2-solutions>,
// )
// To illustrate the effect of regularization, three representative solutions are considered: an under-regularized solution, the optimal solution, and an over-regularized solution, shown in @fig:l2-solutions. For small values of $alpha^2$, the regularization is too weak and the inversion fits the data too closely. This leads to overly oscillatory behavior and amplification of high-degree spherical harmonic components, producing unrealistic small-scale features in the CMB field. For large values of $alpha^2$, the regularization becomes too strong and suppresses regions with high spatial frequencies model parameters excessively, resulting in an overly smooth field. The optimal solution, as shown in @fig:l2-field-map, lies between these two extremes and notably preserves the large reversed-flux region over the South Atlantic Ocean while reducing noise.

// #pagebreak(weak: true)
// #meander.reflow({
//     import meander: *
//   placed(top + right, stack(fig-l2-field-map, v(1em), fig-l2-different_power_spectra, v(1em), fig-l2-residual))
//   container()
//   content([
// @fig:l2-different-power-spectra shows the power spectra for different values of $alpha^2$. When $alpha^2$ is small, the spectrum retains relatively strong high-degree contributions, indicating that high-frequency structure and noise are insufficiently damped. When $alpha^2$ is large, the high-degree components are strongly suppressed and the spectrum is dominated by low-degree terms.

// Among the examined solutions, the case with $alpha_*^2 = 1.60 times 10^(-8)$ shown in green provides a good balance between retaining the structure of the reconstruction while suppressing high-frequency spatial fluctuations, as also reflected by its intermediate decay behavior in the power spectrum.

// The residual distribution for the optimal $L_2$-regularized solution is shown in @fig:l2-residual. It is strongly peaked around zero and approximately symmetric, but is clearly not normally distributed. However, the residuals exhibit heavier tails than a Gaussian distribution and are more closely matched by a Laplace distribution. This suggests that the errors are not purely Gaussian and may include modelling errors, outliers, or unresolved physical processes. This behavior highlights a limitation of the non-robust $L_2$ formulation and motivates the use of robust $L_1$ regularization.
// ])
// })


// #pagebreak(weak: true)
// == Robust $L_1$-Regularised Solution <sec:results-l1>
// #let alpha_sq_l1_star = num("7.54e-5")

// #let fig-l-curve-l1 = [
//   #box(width: 50%)[
//     #figure(
//       image("export/L_curve_L1.png"),
//       caption: [
//         L-curve for the robust $L_1$-regularised solutions.
//         The corner is determined using the orthogonal chord distance.
//       ]
//     ) <fig:l-curve-l1>
//   ]
// ]

// #meander.reflow({
//   import meander: *

//   placed(top + right, stack(fig-l-curve-l1))
//   container()
//   content([
// For the $L_1$-regularised solution with robust weights, we require the full machinery described in @sec:math-regularisation and must resort to iteratively approaching the correct weights $mm(W), mm(W)_"m1"$ using the IRLS algorithm given in @alg:irls-r with the norm order given by $p=1$.

// We again seek to find an optimal regularisation parameter $α^2$ by employing the L-curve heuristic, which must notably be carried out with the appropriately weighted norms. That is, we must use $norm(mm(G) vvh(m)_(r α 1) - vvh(d))^2_mm(W)$ and $norm(mm(H) vvh(m))_1$ to produce the L-curve on which we seek the corner by finding the point with the largest orthogonal distance to the chord between the extremal values of $α$. This method offers a more robust method of estimating the L-curve corner than those seeking the maximal curvature, which were observed to be inaccurate in some cases due to their reliance on estimating the second derivative.
// ])
// })

// The regularisation parameter with optimal trade-off between the regularisation and data fidelity norms is found to be $α_*^2 = #alpha_sq_l1_star$ and gives the solution seen in @fig:l1-solutions alongside under-regularised and over-regularised solutions, where @fig:l1-under-regularised and @fig:l1-solution are computed to tolerance $τ=10^(-9)$ and @fig:l1-over-regularised is computed for $τ=num("5E-7")$.

// We note that insufficient regularisation pushes the solution towards the oscillatory, unstable regime first seen in @sec:results-ls, whereas the over-regularised solution promotes sparsity to the extent that important features are removed. A larger plot of the optimal solution is available in @app:l1-plot. For the remainder of the analysis, we select the parameter $α_*^2 = #alpha_sq_l1_star$.

// #subpar.grid(
//   align: top,
//   figure(image("export/sol_irls_L1_cmb_alpha_5.00e-06.png"), caption: [
//     Under-regularised, $α^2 = num("5e-6")$
//   ]), <fig:l1-under-regularised>,
//   figure(image("export/sol_irls_L1_cmb.png"), caption: [
//     Optimal regularisation, $α^2_* = #alpha_sq_l1_star$
//   ]), <fig:l1-solution>,
//   figure(image("export/sol_irls_L1_cmb_alpha_1.00e+00.png"), caption: [
//     Over-regularised, $α^2 = 1$
//   ]), <fig:l1-over-regularised>,
//   columns: (1fr, 1fr, 1fr),
//   caption: [Solutions to the robust $L_1$-regularised problem for various regularisations parameters.],
//   label: <fig:l1-solutions>,
// )

// #let fig-spectrum-l1 = [
//   #box(width: 50%)[
//     #figure(
//       image("export/spectrum_L1.png"),
//       caption: [
//         Power spectrum of the robust $L_1$-regularised solution. Note in particular that the power decreases for higher degrees.
//       ]
//     ) <fig:spectrum-l1>
//   ]
// ]

// #let fig-covariance-l1 = [
//   #box(width: 50%)[
//     #figure(
//       image("export/covariance_L1.png"),
//       caption: [
//         The model covariance matrix and its associated standard deviations. Saw-tooth pattern arises from increments in degree of spherical harmonics.
//       ]
//     ) <fig:covariance-l1>
//   ]
// ]

// #meander.reflow({
//   import meander: *

//   placed(top + right, stack(fig-spectrum-l1))
//   container()
//   content([

// If we inspect the power spectrum of the solution shown in @fig:spectrum-l1, we find that the contribution to the total power spectrum decreases for higher-order components of the model — this is encouraging, as a diverging power spectrum would indicate that the truncated orders are significant and the model derived in @sec:theory-physics would need to be reconsidered in this case. The fact that most of the power is contained in the first degree can be understood through the largely dipolar shape of the obtained field, where the field lines mostly leave the core in the southern hemisphere and enter it in the northern hemisphere. A perfect dipole would, of course, be captured wholly by the first degree, as evident by @eq:semi-discrete-model.

// It is also worthwhile to consider the rate of decay in the components of @fig:spectrum-l1, which has largely leveled off at a high constant value of $≈ qty("4e9", "nT^2")$ at higher degrees. This is somewhat concerning as it implies that the truncated degrees are not entirely insignificant to the  prediction of the field at the CMB. We will contemplate this further below and return to this point in the discussion.
// ])})

// // #pagebreak(weak: true)
// === Covariance and Resolution Characteristics <sec:results-l1-stats>
// #meander.reflow({
//   import meander: *

//   placed(top + right, stack(fig-covariance-l1))
//   container()
//   content([
// When also considering higher degrees, we recover a region of flux reversed relative to the simple dipole model that may be considered as a null hypothesis. In order to give credence to its existence, we first compute an estimate of  model covariance matrix, $mm(Σ)_vvh(m)$ by the method described in @sec:math-cov-res-matrix. We note again that the estimate is produced under a number of assumptions about the data, in particular that the data errors are uncorrelated and have a common variance, which may not be the case if the dominant source of errors are local influences from the ionosphere and magnetosphere.

// The first 50 indices of the model covariance are shown in @fig:covariance-l1 where the diagonal elements can be seen to dominate until approximately the 14th degree. This is at least in part a result of having used the _plug-in_ estimate for the data covariance $Σ_vvh(d)$ as forewarned. In the lower figure the ordering of the Gauss coefficients becomes apparent by the jagged nature of the parameter standard deviations. The flattened parameter index $i$ follows the natural ordering of @eq:semi-discrete-model and each saw-tooth of the plot is attributed to an increment of the degree $0 < n ≤ N$ with the associated orders $0≤m≤n$ immediately following it. The mapping from the first Gauss coefficients of each degree, $g_n^0$, and the model parameter index $i$ can be found in @table:degree-index-map. In conclusion, higher indices generally belong to higher frequency components of the model.
// ])
// })

// #pad(0em)[
//   #let (ns, is_) = {
//     let ns = ()
//     let is_ = ()
//     let N = 20

//     let i = 0
//     for n in range(1, N+1) {
//       ns.push(n)
//       is_.push(i)

//       for m in range(0, n+1) {
//         i += 1  // g^m_n
//         if m == 0 { continue }  // Skip h^0_n
//         i += 1  // g^m_n
//       }
//     }

//     (ns, is_)
//   }
//   #set text(size: 10pt)
//   #let big = 1pt
//   #let small = 0.5pt
//   #let padding = 1pt
//   #figure(
//     table(
//       stroke: small,
//       inset: 3pt,
//       columns: ns.len() + 1,
//       table.vline(stroke: big),
//       table.hline(stroke: big, start: 0, end:1),
//       pad([$g^0_n$], padding),
//       table.vline(stroke: big),
//       table.hline(stroke: big, start: 0, end:1),
//       ..ns.map(x => [$#x$]),
//       pad([$i$], padding),
//       ..is_.map(x => [$#x$]),
//       table.hline(stroke: big, start: 0, end:1),
//     ),
//     caption: [Mapping from the first Gauss coefficient of each degree $n$ to the model parameter index $i$.]
//   ) <table:degree-index-map>
// ]


// // #pagebreak(weak: true)

// #let fig-resolution-l1 = [
//   #box(width: 50%)[
//     #figure(
//       image("export/resolution_L1.png"),
//       caption: [
//         Model resolution matrix and its diagonal. Jagged edges in diagonal arise from increments in the degree of the spherical harmonics.
//       ]
//     ) <fig:resolution-l1>
//   ]
// ]
// It may be initially surprising that the standard deviation decreases for higher-frequency Gauss coefficients, but this arises from the bias introduced by the regularisation. As the harmonic degree $n$ increases, the signal at the satellite becomes increasingly attenuated by the $(a\/r)^(n+2)$ factor of @eq:semi-discrete-model since $a\/r_"sat" ≈ 0.925$. This is what leads to instability in the unregularised problem as discussed in @sec:math-illposedness, but in the regularised problem these terms are instead dominated by the regularisation scheme, which penalises the size of the projection onto the CMB in the $L_1$ norm and promotes sparsity. This drives the high-frequency coefficients towards zero and explains the phenomenon observed in @fig:covariance-l1 while also giving a theoretical understanding of the over-regularised case shown in @fig:l1-over-regularised.

// #pagebreak(weak: true)
// #meander.reflow({
//   import meander: *

//   placed(top + right, stack(fig-resolution-l1))
//   container()
//   content([
// To give credence to this hypothesis, we compute the model resolution matrix, $mm(R)_vvh(m)$, which is depicted in @fig:resolution-l1. It largely follows the same shape as @fig:covariance-l1, but the interpretation is different — where the covariance matrix contains information about the variance of the model parameters, the resolution matrix offers insights into the _bias_ of the model. This can be understood by its construction from the pseudo-inverse and forward projection, which probes the null-space of the inversion. Diagonal entries close to unity correspond with perfectly resolved parameters, whereas values close to zero indicate that only very little information about those parts of the model are recovered. Off-diagonal elements are understood to capture the model's ability to discern how information of the measurement data is distributed onto the parameters, where higher values indicate a larger influence of one parameter on the other.

// We observe that higher frequency components have lower resolutions and are almost entirely lost during the inversion as a result of the strong attenuation discussed in @sec:math-illposedness. This also suggests that even if higher degrees were included in the model, we would not be able to derive significant information from them given the same dataset, $hat(D)$ — the higher degrees are already dominated by the bias introduced through the regulariser.

// In conclusion, we find that the model has poor resolution of smaller spatial features and that the reconstruction of those are dominated by the a priori choice of regulariser.

// Lastly, we may compute estimates of the variance in the magnetic field predictions at the CMB, which are the diagonal entries of $Σ_vvh(d)_"CMB"$ described in @sec:math-cov-res-matrix. From these we can produce a map over the estimated standard errors of our prediction, as shown in @fig:l1-solution-pred-std, where we note that the region with reversed flux has standard errors lower than #qty("0.1", "mT"). From this, we can produce the binary mask shown in @fig:l1-solution-credibility, which demarcates regions of higher and lower confidence, here given by the conservative threshold that the prediction must lie further than 5 standard deviations from zero to be considered credible.

// #subpar.grid(
//   align: top,
//   figure(image("export/sol_irls_L1_cmb.png"), caption: [
//     The $L_1$ regularised solution
//   ]), <fig:l1-solution-2>,
//   figure(image("export/pred_std_dev_L1.png"), caption: [
//     Estimated standard error of predicted field at CMB
//   ]), <fig:l1-solution-pred-std>,
//   figure(image("export/pred_credible_L1.png"), caption: [
//     Regions where prediction is larger than 5 standard deviations
//   ]), <fig:l1-solution-credibility>,
//   columns: (1fr, 1fr, 1fr),
//   caption: [Investigation into credibility of $L_1$-regularised solution.],
//   label: <fig:l1-credibility>,
// )

// We find that our prediction satisfies this threshold for most of the surface of the CMB, but struggles to confidently recover the high-frequency edges of some patches, which agrees with the previous analysis of the model covariance and resolution matrices. In particular, we note that the region of reversed flux beneath the southern Atlantic satisfies the confidence threshold, indicating that the observed phenomenon is not merely an artefact of the inversion.
// ])
// })

// = Discussion <sec:discussion>

// Comparing the predicted field at the CMB produced by the $L_2$ and robust $L_1$ model, best viewed in @app:l2-plot and @app:l1-plot respectively, we clearly observe the effect of the different regularisations schemes for the high-frequency components in particular, which we found to have high variance and low resolution in @sec:results-l1. Where the $L_2$-regularisation steers the solution towards a smoothly varying surface, the $L_1$-regulariser promotes sparsity instead, which manifests as regions with only minute radial magnetic fields.

// However, for those parts of the model that are well-represented in the data, we find remarkable agreement between the two models. In particular, we note that the region of reversed flux in the southern Atlantic is found in both models, which further supports a rejection of the dipole model taken as the null-hypothesis. The region is large, and is thus not likely to be a feature of the regularisation bias, given the relatively low spatial frequencies associated with it. By the analysis in @sec:results-l1, we can thus expect the region to be well-resolved by at least the model regularised in the $L_1$-norm.

// The analysis of the prediction error carried out in @sec:results-l1 provides very high confidence in the existence of a region of reversed flux given the $5σ$ confidence threshold map of @fig:l1-solution-credibility, though as noted this result relies on both a plug-in estimate obtained from the Median Absolute Deviation of the residuals in the final iteration of @alg:irls-r as well as the assumption that the data errors are both univariate and uncorrelated. Given our understanding that these errors are likely the result of charges and currents in the ionosphere and magnetosphere — which also violate the interior sources assumption made in @sec:theory-physics — it seems unreasonable that they would share a common variance and be entirely uncorrelated. For this reason, it would be unwise to rely too heavily on the obtained estimate of the prediction error, hence the conservative $5σ$ threshold value.

// Additionally, we may return to the concerning finding in @fig:spectrum-l1 that the power carried by the truncated model parameters may be significant. Thankfully, those concerns are largely alleviated by the subsequent analysis regarding the regularisation bias presented in @sec:results-l1-stats, which explains that the Gauss coefficients associated with harmonic degrees larger than approximately $n=14$ are determined largely by the a priori choice of regulariser. This was also discussed in @sec:math-illposedness and arises due to the attenuation of high-frequency information in the model. Consequently, the truncation of the model does not noticeably diminish the amount of credible information about the geomagnetic behaviour at the CMB made available by our method.

// Lastly, we propose that the model may be significantly improved given a larger dataset, or by relaxing the assumption of interior sources in the forward projection to produce a more complex model that can directly incorporate dynamics or data from the ionosphere and magnetosphere. The latter would reduce the model residuals considerably and could be carried out either as a separate contribution to the field or by generalising the interior sources assumption.

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

== Large Plot of Solution to $L_2$-regularised Problem <app:l2-plot>
#figure(
  image("export/L2_field_map_v2.png"),
  caption: [
    // Large version of the plot shown in @fig:l2-solutions (b)
  ]
)
