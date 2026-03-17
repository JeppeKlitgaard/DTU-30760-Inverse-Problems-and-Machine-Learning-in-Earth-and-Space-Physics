#import "@preview/codly:1.3.0": *
#import "@preview/cetz:0.4.2"
#import "@preview/unify:0.7.1": num, qty
#import "@preview/algorithmic:0.1.0"
#import "@preview/physica:0.9.8": curl, grad, tensor, pdv, dv, TT

#import "@preview/meander:0.4.1"
#import "@preview/algorithmic:1.0.7"
#import algorithmic: style-algorithm, algorithm-figure, algorithm

#import "preamble_funcs.typ": mref, mathformatter
#import "style.typ": style

#show: style-algorithm
#show: style
#set page(numbering: "i")

#let vv = mathformatter(underbars: 0, bold: true, upright: false)
#let mm = mathformatter(underbars: 0, bold: false, upright: true)
#let vvh = x => math.hat(vv(x))

#set math.mat(delim:"[")

#let Cov = math.op("Cov")
#let diag = math.op("diag")
#let median = math.op("median")

#let wider = h(3em)

#let HH = math.sans(math.upright($H$))


= TODO

- ??: Abstract
- JK: Style report
// - ??: Discussion on ill-posedness (Hadamard sense)
//   - $1/r$ term
// - JK: Merge math theory sections
// - JK: Huber and MAD
- ??: Some recap/discussion of the assumptions used and how reasonable they are
- ML: Write-up L2, (i)
- JK: Write-up L1/Robust, (ii)
- BOTH: Covariance and resolution (iii) computation
- ML: (iv)
- ML: Introduction of mission and experimental design
- EXTRA: Resolve frequencies from modes (how far are we from seeing 60Hz?)
- ??: Discussion around patch size related to discretisation
- ??: Effect of large $α^2$
- ??: Comment on inadequacy of provided L curve estimate
  - Mention that it is a heuristic anyway
- EXTRA: Cholesky solver
- Typeset loss function

= Questions for TA/Chris
- We do not have access to a good estimate of the dispersion of the "untainted" Gaussian in this case, so presumably we need to iteratively find this using the residuals of the model?
- Lecture notes mention MAD as Mean Absolute Deviation and give formula using mean, but literature exclusively suggests MAD, Median Absolute Deviation.
- How are we suppose to approximate the data covariance matrix $Cov(vv(d))$? Do we assume the scale parameter may be found using MAD and that the samples are uncorrelated to get a diagonal matrix, or do we somehow use the robust weights as the covariance matrix?
- Will the characteristic blot size not be biased by $N$?
  - A: Can be investigated, but notice shape of power spectrum (peak is not necessarily maxed out at highest frequencies)
- Check understanding: Discrete jumps in covariance stems from increments of $n$ (order of spherical harmonic)
- Check understanding: When to use weights, when not to
- Check understanding: Cholesky OK, since PD

= Abstract

#outline()

= Artificial Intelligence Declaration

The body of work presented below is solely authored by the listed authors of this document, Jeppe Klitgaard and Ming-Ming Li, and none of the presented material has been produced by means of _Generative AI_ in the manner requiring citation as described by @ref:dtu-genai.

Generative AI has, however, been used selectively in the generation of limited parts of the code. GitHub Copilot tab-completion has been used during the development of the attached code and in the generation of an initial version of the following, which is also clearly sign-posted by comments in the attached code:
- Initial version of the plotting code for the residual plot
- Alternative L-curve corner finding algorithm: Normalised Menger Curvature
- Alternative L-curve corner finding algorithm: Chord Distance

The generated code output has in all instances been checked and revised by the authors and is _not_ considered to be a _significant_ contribution to total intellectual work undertaken and presented. 

#pagebreak()
#counter(page).update(1)
#set page(numbering: "1/1")
= Introduction

This report methods and implementations from the field of inverse problems and regularisation theory as applied to the determination of the radial magnetic field of Earth's surface and subsurface using data from ESA's _Swarm_ satellite mission recorded between XXX - YYY.

= Experimental Design

TODO

= Theory — Physics and Forward Model

Before an inverse problem can be solved, we must present the corresponding forward model of the system, which may be formulated in the linear, discretised case as
$
  vvh(d) = vv(d) + vv(e) = mm(G) vv(m) + vv(e),
$ <eq:forward-model>

where $mm(G) : cal(M) → cal(D)$ is the _system matrix_ which maps model parameters $vv(m)$ from the discrete model space, $cal(M) = ℝ^(N_m)$, to a vector $vv(d) = mm(G) vv(m)$ in the observation space, $cal(D) = ℝ^(N_d)$. Here, $vv(d)$ is understood to be the _deterministic_ component of the forward projection, and is colloquially — and somewhat inaccurately — referred to as the _noise-free data_. Naturally, we do not observe the forward projection of some idealised physical model, but rather observations $vvh(d)$ which also contain contributions from both the realisations of influencing stochastic processes and any systematic biases, both of which are captured by $vv(e)$. As we will later discover, $vv(e)$ cannot reasonably be assumed to be normally distributed.

In order to determine the system matrix $mm(G)$, we must consider carefully the physical model of the system. To this end, we start from the Ampère-Maxwell law and the assumption that the region of measurement — the orbits upon which the radial component $B_r$ is measured using the _Swarm_ satellites — is devoid of free charges and currents. That is, we assert that the magnetic field is curl-free:
$
  ∇ × vv(B) = μ_0 (vv(J) + ε_0 pdv(vv(E), t)) = vv(0).
$

Thus, we may define a _scalar magnetic potential_, $V$, such that
$
  vv(B) = -∇V.
$ <eq:scalar_potential>

Deeming it unlikely that any magnetic monopoles will present themselves, we may conclude that the field is also divergence-free, $∇ ⋅ vv(B) = 0$, which yields the Laplace equation for the scalar potential:
$
  ∇ ⋅ vv(B) = -∇^2 V = 0.
$

If we make the simplifying assumption that a spherical boundary surface of radius $a$ exists which delineates the electromagnetic sources from empty space, we find that the spatial domain of our model is the sphere, $S^2$. While such a spherical boundary is reasonable for the case of the Earth, assuming that no sources — here understood to be magnetic moments and electrical currents — exist outside of it is more difficult to justify. In particular, this neglects the charge carriers in the ionosphere and the influence of the magnetosphere, a discussion which we will return to later. 

With this in mind, we may choose the spherical harmonics $Y_n^m (ϑ, λ)$ as a convenient, orthogonal basis for the model space and continue our analysis in spherical coordinates $(r, ϑ, λ)$, in which we recall Laplace's equation to be of the form:
$
  ∇^2 V = 1/r^2 pdv(,r)(r^2 pdv(V, r)) + 1/(r^2 sin ϑ) pdv(, ϑ)(sin ϑ pdv(V, ϑ)) + 1/(r^2 sin^2 ϑ) pdv(V, λ, 2) = 0.
$

Separation of variables and the assumption of interior sources yield the solution @ref:geomath[p. 507, eq. 13]
$
  V(r, ϑ, λ)
  &=a ∑_(n=1)^∞ (a/r)^(n+1) ∑_(m=0)^n (g_n^m cos m λ + h_n^m sin m λ ) P_n^m (cos ϑ)\
$ <eq:lapl_sol>

in which $P_n^m$ are the associated Legendre polynomials and we have used the spherical harmonics as given by:
$
Y_(n,"even")^m = cos (m λ) P_n^m (cos ϑ)
wider
Y_(n,"odd")^m = sin (m λ) P_n^m (cos ϑ).
$

The parameters $g_n^m, h_n^m$ of @eq:lapl_sol are referred to as the _Gauss coefficients_ and the inverse problem is thus "simply" to determine these given the field measurements.

We do not, however, have access to the potential $V$ directly, which is introduced more as a mathematical convenience than a physical quantity. Instead, we measure the radial component of its gradient as given by @eq:scalar_potential,
$
  B_r (r, ϑ, λ)
  = - pdv(V, r)
  = ∑_(n=1)^∞ (n+1) (a/r)^(n+2) ∑_(m=0)^n (g_n^m cos m λ + h_n^m sin m λ ) P_n^m (cos ϑ).
$ <eq:B_r>

By identifying this as forward projection and truncating the spherical harmonic expansion to degree $N$, we obtain the semi-discrete model in terms of coordinates:
$
  B_r (r, ϑ, λ) = ∑_(n=1)^N (n+1) (a/r)^(n+2) ∑_(m=0)^n (g_n^m cos m λ + h_n^m sin m λ ) P_n^m (cos ϑ),
$ <eq:semi-discrete-model>
while the finite-dimensional model parameters are consequently given by
$
vv(m) = (g_1^0, g_1^1, h_1^1, g_2^0, ..., g_N^N, h_N^N)^T  ∈ cal(M)
$
where notably the $h_n^0$ coefficients are omitted as their terms vanish due to the $sin(m λ)$ factor in @eq:lapl_sol. This truncation is sensible under the assumption that the variations of the magnetic fields are smooth and that the features of interest are of sufficiently low frequency to be captured by the remaining basis functions.

By constructing a system of equations using @eq:semi-discrete-model and the associated discrete model parameters $vv(m)$ over an idealised, noise-free set $D = {d_i = B_r (r_i, ϑ_i, λ_i)}_(i=1)^(N_d)$, we produce the discrete forward operator as the system matrix, $mm(G)$, which is built as follows:
$
  vv(d) = mm(G) vv(m)\
  vec(
    B_r (vv(x)_1),
    B_r (vv(x)_2),
    ⋮,
    B_r (vv(x)_(N_d)),
  )
  =
  mat(
    // g=cos, n=1, m=0
    Ξ_1^0 (vv(x)_1),
    // h=sin, n=1, m=0, this vanishes
    // g=cos, n=1, m=1
    Ξ_1^1 (vv(x)_1),
    // h=sin, n=1, m=1
    ξ_1^1 (vv(x)_1),
    …,
    // g=cos, n=N, m=N
    Ξ_N^N (vv(x)_1),
    // h=sin, n=N, m=N
    ξ_N^N (vv(x)_1)
    ;
    Ξ_1^0 (vv(x)_2),
    Ξ_1^1 (vv(x)_2),
    ξ_1^1 (vv(x)_2),
    …,
    Ξ_N^N (vv(x)_2),
    ξ_N^N (vv(x)_2)
    ;
    ⋮,⋮,⋮,⋮,⋮,⋮
    ;
    Ξ_1^0 (vv(x)_(N_d)),
    Ξ_1^1 (vv(x)_(N_d)),
    ξ_1^1 (vv(x)_(N_d)),
    …,
    Ξ_N^N (vv(x)_(N_d)),
    ξ_N^N (vv(x)_(N_d))
    ;

  )
  vec(
    g_1^0,
    g_1^1,
    h_1^1,
    ⋮,
    g_N^N,
    h_N^N,
  )
$ <eq:forward-model-big>

where we have introduced some shorthand notation
$
vv(x)_i &= (r_i, ϑ_i, λ_i)\
Ξ_n^m (r, ϑ, λ) &= (n+1) (a/r)^(n+2) cos (m λ) P_n^m (cos ϑ)\
ξ_n^m (r, ϑ, λ) &= (n+1) (a/r)^(n+2) sin (m λ) P_n^m (cos ϑ).
$

In conclusion, we can build the discrete forward operator $mm(G) ∈ ℝ^(N_d × N_m), N_m = N(N+2)$ simply by evaluating the basis functions $Ξ_n^m$ and $ξ_n^m$ at the measurement locations, ${vv(x)_i}_(i=1)^(N_d)$. We note that a final assumption made is that the temporal variations in the field are locally much slower than the period over which the data has been collected, or equivalently that $B_r (r, ϑ, λ)$ is constant with respect to time.

= Theory — Mathematical Treatment of the Inverse Problem <sec:theory-math>

Having now produced a discrete forward model that may be used with our experimental data, we turn our attention to the associated _inverse problem_, that is, the estimation of model parameters $vv(m)$ given the observations $hat(D) = {hat(d)_i = B_r (r_i, ϑ_i, λ_i) + e_i}_(i=1)^(N_d)$. 

== Least Squares Solution <sec:math-ls>
Naïvely, we may obtain a solution to the inverse problem by the computing the inverse of $mm(G)$ with which we would left-multiply to obtain an estimate
$
vvh(m) = mm(G)^(-1) vvh(d),
$
except of course such an inversion does not exist since $mm(G)$ is not generally square. Instead, we seek the model parameters which minimise the misfit, $ϕ$, between the observed data $vvh(d)$ and the forward projection in a norm, here taken to be the $L_2$-norm:
$
vvh(m)
= arg min_vv(m) ϕ =  
arg min_vv(m) 1/2 norm(mm(G) vv(m) - vvh(d))_2^2
$

By inspection, the misfit $ϕ$ is convex, differentiable with respect to $vv(m)$, and coercive provided $mm(G)$ is full-rank, and thus the (global) minimum may be found by solving
$
∇_vv(m) ϕ = 0 &= ∇_vv(m) [1/2 (mm(G) vv(m) - vvh(d))^TT (mm(G) vv(m) - vvh(d))]\
&= 1/2 ∇_vv(m) [-vv(m)^TT mm(G)^TT vvh(d) + vv(m)^TT mm(G)^TT mm(G) vv(m) - vvh(d)^TT mm(G) vv(m) + vvh(d)^TT vvh(d)]\
&= 1/2 ∇_vv(m) [-2 vvh(d)^TT mm(G) vv(m) + vv(m)^TT mm(G)^TT mm(G) vv(m)] = - vvh(d)^TT mm(G) + vv(m)^TT mm(G)^TT mm(G) ⇔ mm(G)^TT mm(G) vvh(m) = mm(G)^TT vvh(d)\
vvh(m) &= (mm(G)^TT mm(G))^(-1) mm(G)^TT vvh(d) = (mm(G)^TT mm(G))^(-1) mm(G)^TT (vv(d) + vv(e)) = (mm(G)^TT mm(G))^(-1) mm(G)^TT (mm(G) vv(m) + vv(e)) = vv(m) + (mm(G)^TT mm(G))^(-1) mm(G)^TT vv(e)\
&= vv(m) + mm(G)^† vv(e)
$ <eq:m-hat-ls>

Clearly, our estimate will be rather poor unless $norm(mm(G)^† vv(e)) ≪ norm(vv(m))$. We can more clearly see how errors may enter our estimate if we consider using the constructed model to compute the radial magnetic field at the Earth's Core-Mantle Boundary (CMB).

== Stability Analysis <sec:math-illposedness>
Since the model is fully determined by the expansion coefficients in the angular coordinates, this is done simply by constructing a new system matrix $mm(G)_"CMB"$ as before, but replacing the ${r_i}_(i=1)^(N_d)$ radial coordinates with the radius of the CMB, $r_"CMB" = qty("3480.0", "km")$. By inspection of @eq:forward-model-big, this will an estimate of the field at the boundary, $vvh(d)_"CMB"$:
$
vvh(d)_"CMB" = mm(G)_"CMB" vvh(m) = mm(G)_"CMB" vv(m) + mm(G)_"CMB" mm(G)^† vv(e) = vv(d)_"CMB" + vvh(e)
$

Where $vvh(e) = mm(G)_"CMB" mm(G)^† vv(e)$ is the error of our prediction $vvh(d)_"CMB"$ given the realisation $vv(e)$. If we consider a specific spherical degree $n$ and @eq:semi-discrete-model, we find that the radial dependence factors out as $(n+1) (a\/r)^(n+2)$ which by identity $mm(G)^† mm(G) = mm(I)$ implies that $mm(G)^†$ must scale by the reciprocal, $(n+1)^(-1) (a\/r)^(-n-2)$. If we model all the measurements as having been taken at a single, mean radius, $r_"sat"$, we thus find
$
vvh(e)_n ∝ (n+1)^(-1) (a/r_"sat")^(-n-2) (n+1) (a/r_"CMB")^(n+2) vv(e) = (r_"sat" / r_"CMB")^(n+2) vv(e)
$

In our data, the average radius is $r_"sat" = qty("6891.51", "km")$, giving a spectral noise amplification factor $η_n=1.9803^(n+2), η_20 ≈ num("860000")$. Clearly, even minute errors in the high-frequency components of the measurements may dominate our estimate at the core.

This sensitivity to small changes causes the inverse problem to be severely _ill-posed_ in the sense of Hadamard and requires us to regularise our inversion.

== Regularisation Theory <sec:math-regularisation>
A common regularisation scheme is to introduce a penalty on the norm of model parameters, or some linear transformation thereof. In a general case, such an $L_p$-regularised solution may be found using the variational formulation
$
arg min_vv(m) 1/2 norm(mm(G) vv(m) - vvh(d))_mm(W)^2 + α^2/2 norm(mm(H) vv(m))_p^p,
$ <eq:robust-variational>

where $p$ is the _norm order_, $α^2$ is the _regularisation parameter_, and $mm(H) ∈ ℝ^(N_k × N_m)$ is a linear regularisation operator. The norm $norm(vv(v))_mm(W)^2 ≔ vv(v)^TT mm(W) vv(v)$ indicates the _weighted_ $L_2$ norm where $mm(W)$ is a diagonal matrix of _robust weights_ that must be determined jointly in the absence of an _a priori_ estimate.

To ensure differentiability at zero, the $L_p$ norm may be approximated using the relaxation scheme due to Ekblom:
$
norm(mm(H) vv(m))_p^p ≈ sum_(i=1)^(N_k) [((mm(H) vv(m))_i)^2 + δ^2]^(p\/2),
$
which yields an approximate solution to @eq:robust-variational of the form @ref:course_book[eq. 51]
$
vvh(m)_(r α p) = (mm(G)^TT mm(W) mm(G) + α^2 mm(H)^TT mm(W)_(m p) mm(H))^(-1) mm(G)^TT mm(W) vvh(d) = tilde(mm(G))_(α p)^(-1) mm(G)^TT mm(W) vvh(d) = mm(G)_(α p)^(-1) vvh(d).
$

Here, both the robust data weights $mm(W)$ and the (diagonal) regularisation weights may be obtained through the Iteratively Reweighted Least Squares (IRLS) algorithm, as described by @alg:irls-r. Specifically, the diagonal of $mm(W)_(m p)$ is populated by the vector $vv(w)_m$, which is updated element-wise at each iteration by
$
vv(w)_m = p (vv(γ)^2 + δ^2)^(p\/2 - 1), quad "where" quad vv(γ) = mm(H) vv(m),
$
with $δ$ denoting the _Ekblom perturbation_ which must satisfy $δ^2 ≪ norm(vv(γ))_2^2 \/ N_k.$ 

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
          Maximal iterations $N_"iter"$, convergence threshold $τ$
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

Where the functions used for the M-estimator and deviation scale estimation in @alg:irls-r, $cal(S)(vv(ε);vv(θ)_cal(S))$ and $cal(V)(vv(ε);vv(θ)_cal(V))$ respectively, may be chosen as the Huber weights @ref:course_book[p.9] and normalised Median Absolute Deviation (MAD) @ref:robust-stats[p.36], as given by:
$
w_"Huber" (ε; c=1.345) = cases(
  c abs(ε) &wide abs(ε) > c,
  1 &wide abs(ε) ≤ c,
)
$ <eq:huber-weights>

$
  hat(σ)_"MAD" (vv(ε)) = median(sum_(i=1)^(N_d) abs(ε_i - median(vv(ε))) ) / (Φ^(-1)(0.75))
$
where $Φ(p)$ is the CDF of a standard Gaussian distribution such that $Φ^(-1)(0.75) ≈ 0.675$ @ref:robust-stats[p.36].

// TODO better intro phrasing
== Covariance and Resolution Matrices <sec:math-cov-res-matrix>
We wish to be able to compute the covariance and resolution matrices of our model, which will give insight into respectively the uncertainty and bias of our obtained model, $vvh(m)$.

We note that $tilde(mm(G))_(α p)^(-1) = (mm(G)^TT mm(W) mm(G) + α^2 mm(H)^TT mm(W)_(m p) mm(H))^(-1)$ is symmetric by the diagonality of $mm(W), mm(W)_(m p)$.
Thus, by identity $Cov(mm(A) vv(b)) = mm(A) Cov(vv(b)) mm(A)^TT$, we obtain an estimated model covariance matrix
$
mm(Σ)_vvh(m) = Cov(vvh(m)_(r α p))
&≈ mm(G)_(α p)^(-1) Cov(vvh(d)) (mm(G)_(α p)^(-1))^TT\
&= tilde(mm(G))_(α p)^(-1) mm(G)^TT mm(W) Cov(vvh(d)) (tilde(mm(G))_(α p)^(-1) mm(G)^TT mm(W))^TT\
&= tilde(mm(G))_(α p)^(-1) mm(G)^TT mm(W) Cov(vvh(d)) mm(W) mm(G) tilde(mm(G))_(α p)^(-1).\
$ <eq:cov-matrix>

Here $mm(W)$ denotes the _frozen_ robust weights obtained at the final iteration of @alg:irls-r, which provide a reasonable, linear approximation in this case.

The variance of the measurements are not known a priori and we resort to a _plug-in_ estimate in which we assume the data to be uncorrelated and homoscedastic with a variance approximated by the Median Absolute Deviation of the residuals:
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

#pagebreak(weak: true)
= Results
== Unregularised Least Squares Solution

#let fig-residual-ls = [
  #box(width: 50%)[
    #figure(
      image("export/residual_ls.png"),
      caption: [
        Residual of naïve least squares solution.
      ]
    ) <fig:residual-ls>
  ]
]

#let fig-sol-ls-cmb = [
  #box(width: 50%)[
    #figure(
      image("export/sol_ls_cmb.png"),
      caption: [
        Bad reconstruction of radial magnetic field at CMB due to unstable inversion.
      ]
    ) <fig:sol-ls-cmb>
  ]
]

#meander.reflow({
  import meander: *

  placed(top + right, stack(fig-residual-ls, v(1em), fig-sol-ls-cmb))
  container()
  content([
If we disregard the analysis carried out in @sec:theory-math and attempt to construct a model by inversion using the naïve, unregularised least-squares solution given by @eq:m-hat-ls, we recover the projection of the radial magnetic field at the CMB shown in @fig:sol-ls-cmb.

Note that magnitude of field is much larger than the expected #qty("1", "mT") and field exhibits rapid oscillations towards the poles.
We note in particular the rapid oscilliations exhibited by the field towards the poles, which does not agree with our expectations. Additionally, the magnitude of the field oscillations is much larger than that found in the literature @ref:hammer, which prescribes an amplitude of approximately #qty("1", "mT"). This is a consequence of the aforementioned ill-posedness of the inverse problem. Further, inspection of the residuals $vv(ε) = vvh(d) - mm(G) vvh(m)$ as shown in @fig:residual-ls reveals that these are far from normally distributed and instead appear to follow a long-tailed distribution such as the Huber distribution.

For this reason, we choose the Huber weights given in @eq:huber-weights as the M-estimator for our IRLS algorithm, which has the effect of reducing the influence of observations further from the centre of the data distribution.
  ])
})

== Non-Robust $L_2$-Regularised Solution

We consider the specific case of a non-robust $L_2$-regularised solution, which by the treatment in @sec:theory-math and identifying
$
p = 2
, wider
mm(W) = mm(W)_"m2" = mm(I)
$

gives
$
vv(m)_(α 2)
&= (mm(G)^TT mm(W) mm(G) + α^2 mm(H)^TT mm(W)_"m2" mm(H))^(-1) mm(G)^TT mm(W) vv(d)\
&= (mm(G)^TT mm(G) + α^2 mm(H)^TT mm(H))^(-1) mm(G)^TT vv(d) = tilde(mm(G))_(α 2)^(-1) mm(G)^TT vv(d) = mm(G)_(α 2)^(-1) vv(d)\
$

where $mm(G)_(α 2)^(-1) = (mm(G)^TT mm(G) + α^2 mm(H)^TT mm(H))^(-1) mm(G)^TT$ may be used to find the covariance and resolution matrices through @eq:cov-matrix and @eq:res-matrix respectively.

Notably, since neither of the iteratively constructed weight matrices feature, this regularisation can be carried out without resorting to the IRLS algorithm.

Suggested TODO MML:
- Move and merge L2 section here
- Add plot
// - Consider reporting as $α$ instead of $α^2$ to match writing above.



== Robust $L_1$-Regularised Solution

#let fig-l-curve-l1 = [
  #box(width: 50%)[
    #figure(
      image("export/L_curve_L1.png"),
      caption: [
        L-curve for the robust $L_1$-regularised solutions.
        The corner is determined using the orthogonal chord distance.
      ]
    ) <fig:l-curve-l1>
  ]
]

#meander.reflow({
  import meander: *
  
  // placed(top + right, stack(fig-residual-ls, v(1em), fig-sol-ls-cmb))
  placed(top + right, stack(fig-l-curve-l1))
  
  container()
  content([
For the $L_1$-regularised solution with robust weights, we require the full machinery described in @sec:math-regularisation and must resort to iteratively approaching the correct weights $mm(W), mm(W)_"m1"$ using the IRLS algorithm given in @alg:irls-r with the norm order given by $p=1$.

We again seek to find an optimal regularisation parameter $α^2$ by employing the L-curve heuristic, which must notably be carried out with the appropriately weighted norms. That is, we must use $norm(mm(G) vvh(m)_(r α 1))^2_mm(W)$
])
})

= Finding an $L_2$ least-squares solution
#let fig-l2-curve = [
  #box(width: 50%)[
    #figure(
      image("export/L2_curve.png"),
      caption: [L curve for the L2 least-squares solution]
    ) <fig:l2-curve>
  ]
]
#meander.reflow({
    import meander: *
  placed(top + right, stack(fig-l2-curve))
  container()
  content([
The L-curve obtained for the non-robust L2-regularized inversion shows a clear corner, and thus the $alpha^2$ here is 
$
alpha^2 = 4.642 times 10^(-9).
$

This value gives a reasonable balance between the data misfit and the model norm. As the $alpha^2$ gets smaller, the regularization becomes weaker and the inversion starts to fit short-wavelength structure more strongly. As the $alpha^2$ gets larger, the regularization becomes stronger and overly smooth, meanwhile, the small-field is suppressed.
  ])
})



#let fig-l2-field-map = [
  #box(width: 50%)[
    #figure(
      image("export/L2_field_map.png"),
      caption: [L2 field map]
    ) <fig:l2-field-map>
  ]
]

#meander.reflow({
    import meander: *
  placed(top + right, stack(fig-l2-field-map))
  container()
  content([
In the core mantle boundary (CMB) field map obtained for the $alpha^2$ shows coherent large-scale flux patches. The main structures in this figure are broad and spatially smooth. This figure also shows the negative flux dominating parts of the northern hemisphere, especially in the high latitudes and positive flux dominating parts of the southern hemisphere, especially in the southern Atlantic, Africa and Indian Ocean areas. In the southern Atlantic, there is a circle feature shown positive flux in this map. In the power spectrum of estimated geomagnetic field at core mantle boundary (CMB) for the corresponding $alpha^2$ it shows that it is dominated by the lowest degrees and decreases toward higher degree. It means that the large-scale field is expected to be constrained better for the SWARM data and in the meanwhile, the smaller scales are more sensitive to instability and noise. The spectral decay is relatively smooth and the higher degree power is clearly dampened after 13 degree.
  ])
})
#let fig-l2-power-spectrum = [
  #box(width: 50%)[
    #figure(
      image("export/L2_power_spectrum.png"),
      caption: [L2 power spectrum]
    ) <fig:l2-power-spectrum>
  ]
]
#meander.reflow({
    import meander: *
  placed(top + right, stack(fig-l2-power-spectrum))
  container()
  content([
To explore how changing the value of $alpha^2$ affects solution, we produced the magnetic field maps and the power spectrum for  $alpha^2 = 4.642 times 10^(-11)$,
$4.642 times 10^(-10)$, $4.642 times 10^(-8)$, and
$4.642 times 10^(-7)$. When the $alpha^2$ becomes smaller the CMB field map becomes much more oscillatory and more alternating smaller patches appear. And for smaller $alpha^2$, the large power remains at higher degrees. The opposite trend occurs for the larger $alpha^2$. We have observed that in the CMB field map, the overall feature becomes smoother and the higher-degree part of the spectrum is strongly reduced, which may due to the regularization dominates the inversion too strongly, and the resolution of the field map is decreased. Among all the field maps and the power spectrums, the $alpha^2 = 4.642 times 10^(-9)$ present the best result.
  ])
})

#let fig-l2-residual = [
  #box(width: 50%)[
    #figure(
      image("export/L2_residual.png"),
      caption: [L2 residual]
    ) <fig:l2-residual>
  ]
]
#meander.reflow({
    import meander: *
  placed(top + right, stack(fig-l2-residual))
  container()
  content([
The residual distribution for the selected L2 solution is shown in Figure~@fig:l2-residual. The histogram is strongly peaked around zero and remains approximately symmetric, which shows that the inversion captures the main large-scale signal in the data. At the same time, the distribution exhibits noticeably heavier tails than a Gaussian distribution, and a small number of large residuals are still present. It clearly fits a Laplacian curve.
These remaining tails likely reflect non-Gaussian noise, modelling errors, or contributions from sources not fully represented by the forward model.
  ])
})

= model covariance matrix and resolution matrix solution
#let fig-covariance-matrix = [
  #box(width: 50%)[
    #figure(
      image("export/covariance_matrix.png"),
      caption: [covariance matrix]
    ) <fig:covariance-matrix>
  ]
]
#let fig-resolution-matrix = [
  #box(width: 50%)[
    #figure(
      image("export/resolution_matrix.png"),
      caption: [resolution matrix]
    ) <fig:resolution_matrix>
  ]
]
#meander.reflow({
    import meander: *
  placed(top + right, stack(fig-covariance-matrix, v(1em), fig-resolution-matrix))
  container()
  content([


To assess the reliability of the obtained model, we compute the model resolution matrix

$ R = (G^T G + alpha^2 H^T H)^(-1) G^T G. $

The diagonal elements of this matrix provide a measure of how well each model parameter is resolved. In the figure showing the diagonal of $R$, a clear trend is observed: the first parameters, corresponding to low spherical harmonic degrees, have values close to 1, while the values decrease rapidly for higher indices.

From the computed values, the mean resolution is approximately 0.389, with a maximum of about 1.003 and a minimum close to −0.008. The values close to 1 indicate that the large-scale components of the field are well resolved, while the rapid decrease towards zero shows that small-scale features are poorly constrained by the data. The small negative values are attributed to numerical effects.

The model covariance matrix is given by

$ C_m = (G^T G + alpha^2 H^T H)^(-1) G^T G (G^T G + alpha^2 H^T H)^(-1), $

and the standard deviations of the model parameters are obtained from

$ sigma_m = sqrt(d i a g(C_m)). $

In the corresponding figure, the uncertainties are relatively small for the low-degree coefficients and increase for higher-degree components. This reflects the stronger influence of the regularisation on small-scale features.

Together, the resolution and covariance figures show that the inversion reliably recovers the large-scale geomagnetic field, while the smaller-scale structure is less well constrained.
  ])
})
= comparison between L2 and L1 solution
The L2 solution minimizes the quadratic misfit $phi_(L 2) = sum_i r_i^2,$ which corresponds to assuming Gaussian-distributed residuals. In contrast, the L1 solution minimizes $phi_(L 1) = sum_i |r_i|, $ which corresponds to a Laplacian error distribution.

The L1 solution is obtained using an iteratively reweighted least-squares (IRLS) approach, where the weights are updated as

$ w_i = 1 / (|r_i|), $

leading to the weighted solution

$ m = (G^T W G + alpha^2 H^T H)^(-1) G^T W d. $

In the figures comparing the L2 and L1 field maps, both solutions recover the same large-scale structures, including the dominant hemispheric flux patterns. However, the L1 solution exhibits slightly sharper and more localized features, while the L2 solution appears smoother due to the quadratic penalisation of residuals.

The residual histograms further highlight the difference between the two approaches. The L2 residuals show a non-Gaussian distribution with heavy tails, while the L1 formulation is better suited to accommodate this behaviour. This indicates that the L1 solution is more robust to outliers in the data.
= convincing evidence for the reversed flux features at the core surface
The radial magnetic field at the core-mantle boundary is given by

$ B_r(r, theta, phi) = sum_(n=1)^N (n+1) (a/r)^(n+2)
  sum_(m=0)^n (g_n^m cos(m phi) + h_n^m sin(m phi)) P_n^m(cos theta). $

The L2-regularised field map clearly shows regions where $B_r$ is both positive and negative, forming coherent patches across the CMB. In the figure, negative flux dominates the northern hemisphere, while positive flux dominates large parts of the southern hemisphere, particularly in the South Atlantic region.

These patterns deviate from a simple dipole field and provide evidence for reversed flux patches. The stability of these features is supported by the smooth decay observed in the power spectrum and by the fact that the main structures remain consistent when varying $alpha^2$ within a reasonable range.

The resolution analysis shows that low-degree components are well resolved, while higher-degree components are not. Therefore, the large-scale reversed flux features observed in the figure can be considered robust, while smaller-scale details should be interpreted with caution.

Overall, the combination of the field map, power spectrum, residual analysis and resolution matrix provides convincing evidence that the reversed flux features at the core surface are physically meaningful and reliably recovered by the inversion.

= Data

The dataset used in this report consists of measured radial component of the Earth's magnetic field, obtained from the Swarm satellite mission. The observations of the radial component of the Earth's magnetic field have been documented at the satellite altitude and are also distributed globally. Each data point corresponds to a measurement result of the magnetic field with the spherical coordinates $(r, theta, lambda)$, where r represents the distance between the Earth's center and the satellite, $theta$ represents the colatitude of the observation and $lambda$ represents the longitude of the observation. 


= how one could build a better model of the core surface magnetic field


= Discussion <sec:discussion>

= Conclusion <sec:conclusion>

#bibliography("references.bib")