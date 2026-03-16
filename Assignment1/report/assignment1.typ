#import "@preview/codly:1.3.0": *
#import "@preview/cetz:0.4.2"
#import "@preview/unify:0.7.1": num, qty
#import "@preview/algorithmic:0.1.0"
#import "@preview/physica:0.9.8": curl, grad, tensor, pdv, dv, TT
// #import algorithmic: algorithm
#import "@preview/algorithmic:1.0.7"
#import algorithmic: style-algorithm, algorithm-figure, algorithm
#show: style-algorithm
#import "preamble_funcs.typ": mref, mathformatter
#import "style.typ": style

#show: style

#let vv = mathformatter(underbars: 0, bold: true, upright: false)
#let mm = mathformatter(underbars: 0, bold: false, upright: true)
#set math.mat(delim:"[")

#let Cov = math.op("Cov")
#let diag = math.op("diag")

#let wider = h(3em)

#let HH = math.sans(math.upright($H$))


= TODO

- ??: Abstract
- JK: Style report
// - ??: Discussion on ill-posedness (Hadamard sense)
//   - $1/r$ term
- JK: Merge math theory sections
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

= Introduction

This report methods and implementations from the field of inverse problems and regularisation theory as applied to the determination of the radial magnetic field of Earth's surface and subsurface using data from ESA's _Swarm_ satellite mission recorded between XXX - YYY.

= Experimental Design

TODO

= Theory

= The Forward Model

== Summary
- ✔ Assume: No currents or charges near instrument on satellite
- ✔ Assume: No magnetic monopoles exist
- Assume: The field is static
- Assume: Earth is spherical
- Choose: The correct orthogonal basis for spherical coordinates is spherical harmonics
- Assume: The field is smooth, so we can truncate the spherical harmonic expansion at a finite degree and order
- Thus: The forward problems goes from the spherical harmonic coefficients to the field measurements
- Thus: Inverse problem is to find the spherical harmonic coefficients from the field measurements.

==

Before an inverse problem can be solved, we must present the corresponding forward model of the system, which may be formulated in the linear, discretised case as
$
  vv(hat(d)) = vv(d) + vv(e) = mm(G) vv(m) + vv(e),
$ <eq:forward-model>

where $mm(G) : cal(M) → cal(D)$ is the _system matrix_ which maps model parameters $vv(m)$ from the discrete model space, $cal(M) = ℝ^(N_m)$, to a vector $vv(d) = mm(G) vv(m)$ in the observation space, $cal(D) = ℝ^(N_d)$. Here, $vv(d)$ is understood to be the _deterministic_ component of the forward projection, and is colloquially — and somewhat inaccurately — referred to as the _noise-free data_. Naturally, we do not observe the forward projection of some idealised physical model, but rather observations $vv(hat(d))$ which also contain contributions from both the realisations of influencing stochastic processes and any systematic biases, both of which are captured by $vv(e)$. As we will later discover, $vv(e)$ cannot reasonably be assumed to be normally distributed.

// here captured by $vv(e)$ 

// measurement vectors $vv(hat(d))$ in the data space, $cal(D) = ℝ^(N_d)$ subject to measurement noise $vv(e) ∈ cal(D)$, which notably may _not_ be assumed to be normally distributed. Here $vv(d) = mm(G) vv(m)$ is understood to be the deterministic — or colloquially, the _exact_ — component of the forward projection. Consequently, $vv(e)$ captures both the realisation of any stochastic influences and systematic biases.

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

Separation of variables and the assumption of interior sources yield the solution @geomath[p. 507, eq. 13]
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
where notably the $h_n^0$ coefficients are omitted as their terms vanish due to the $sin(m λ)$ factor in @eq:lapl_sol.

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

In conclusion, we can build the discrete forward operator $mm(G) ∈ ℝ^(N_d × N_m), N_m = N(N+2)$ simply by evaluating the basis functions $Ξ_n^m$ and $ξ_n^m$ at the measurement locations, ${vv(x)_i}_(i=1)^(N_d)$.

= Mathematical Treatment of the Inverse Problem <sec:theory-math>

Having now produced a discrete forward model that may be used with our experimental data, we turn our attention to the associated _inverse problem_, that is, the estimation of model parameters $vv(m)$ given the observations $hat(D) = {hat(d)_i = B_r (r_i, ϑ_i, λ_i) + e_i}_(i=1)^(N_d)$. Naïvely, we may invert the problem by the computing the inverse of $mm(G)$ with which we would left-multiply to obtain an estimate
$
vv(hat(m)) = mm(G)^(-1) vv(hat(d)),
$
except of course such an inversion does not exist since $mm(G)$ is not generally square. Instead, we seek the model parameters which minimise the misfit, $ϕ$, between the observed data $hat(vv(d))$ and the forward projection in a norm, here taken to be the $L_2$-norm:
$
hat(vv(m))
= arg min_vv(m) ϕ =  
arg min_vv(m) 1/2 norm(mm(G) vv(m) - vv(hat(d)))_2^2
$

By inspection, the misfit $ϕ$ is convex, differentiable with respect to $vv(m)$, and coercive provided $mm(G)$ is full-rank, and thus the (global) minimum may be found by solving
$
∇_vv(m) ϕ = 0 &= ∇_vv(m) [1/2 (mm(G) vv(m) - vv(hat(d)))^TT (mm(G) vv(m) - vv(hat(d)))]\
&= 1/2 ∇_vv(m) [-vv(m)^TT mm(G)^TT hat(vv(d)) + vv(m)^TT mm(G)^TT mm(G) vv(m) - vv(hat(d))^TT mm(G) vv(m) + hat(vv(d))^TT hat(vv(d))]\
&= 1/2 ∇_vv(m) [-2 hat(vv(d))^TT mm(G) vv(m) + vv(m)^TT mm(G)^TT mm(G) vv(m)] = - hat(vv(d))^TT mm(G) + vv(m)^TT mm(G)^TT mm(G)\
hat(vv(m)) &= (mm(G)^TT mm(G))^(-1) mm(G)^TT hat(vv(d)) = (mm(G)^TT mm(G))^(-1) mm(G)^TT (vv(d) - vv(e)) = (mm(G)^TT mm(G))^(-1) mm(G)^TT (mm(G) vv(m) - vv(e)) = vv(m) - (mm(G)^TT mm(G))^(-1) mm(G)^TT vv(e)\
&= vv(m) - mm(G)^† vv(e)
$

Clearly, our estimate will be rather poor unless $norm(mm(G)^† vv(e)) ≪ norm(vv(m))$. We can more clearly see how errors may enter our estimate if we consider using the constructed model to compute the radial magnetic field at the Earth's Core-Mantle Boundary (CMB).

Since the model is fully determined by the expansion coefficients in the angular coordinates, this is done simply by constructing a new system matrix $mm(G)_"CMB"$ as before, but replacing the ${r_i}_(i=1)^(N_d)$ radial coordinates with the radius of the CMB, $r_"CMB" = qty("3480.0", "km")$. By inspection of @eq:forward-model-big, this will an estimate of the field at the boundary, $hat(vv(d))_"CMB"$:
$
hat(vv(d))_"CMB" = mm(G)_"CMB" hat(vv(m)) = mm(G)_"CMB" vv(m) + mm(G)_"CMB" mm(G)^† vv(e) = vv(d)_"CMB" + hat(vv(e))
$

Where $hat(vv(e)) = mm(G)_"CMB" mm(G)^† vv(e)$ is the error of our prediction $hat(vv(d))_"CMB"$ given the realisation $vv(e)$. If we consider a specific spherical degree $n$ and @eq:semi-discrete-model, we find that the radial dependence factors out as $(n+1) (a\/r)^(n+2)$ which by identity $mm(G)^† mm(G) = mm(I)$ implies that $mm(G)^†$ must scale by the reciprocal, $(n+1)^(-1) (a\/r)^(-n-2)$. If we model all the measurements as having been taken at a single, mean radius, $r_"sat"$, we thus find
$
hat(vv(e))_n ∝ (n+1)^(-1) (a/r_"sat")^(-n-2) (n+1) (a/r_"CMB")^(n+2) vv(e) = (r_"sat" / r_"CMB")^(n+2) vv(e)
$

In our data, the average radius is $r_"sat" = qty("6891.51", "km")$, giving a spectral noise amplification factor $η_n=1.9803^(n+2), η_20 ≈ num("860000")$. Clearly, even minute errors in the high-frequency components of the measurements may dominate our estimate at the core.

This sensitivity to small changes causes the inverse problem to be severely _ill-posed_ in the sense of Hadamard and requires us to regularise our inversion.

A common regularisation scheme is to introduce a penality on the norm of model parameters, or some linear transformation thereof. In a general case, such an $L_p$-regularised solution may be found using the variational formulation
$
arg min_vv(m) 1/2 norm(mm(G) vv(m) - vv(d))_2^2 + α^2/2 norm(mm(H) vv(m))_p^p
$ <eq:robust-variational>

Where $p$ is the _norm order_, $α$ is the _regularisation parameter_, and $mm(H)$ is a linear regularisation operator.
We state with reference to @course_book[eq. 51] that the robust solution to @eq:robust-variational may be obtained iteratively and is of the form
$
vv(m)_(r α p) = (mm(G)^TT mm(W) mm(G) + α^2 mm(H)^TT mm(W)_"mp" mm(H))^(-1) mm(G)^TT mm(W) vv(d) = tilde(mm(G))_(α p)^(-1) mm(G)^TT mm(W) vv(d) = mm(G)_(α p)^(-1) vv(d)
$

Where $mm(W)$ and $mm(W)_"mp"$ are diagonal matrices obtained through the Iteratively Reweighted Least Squares (IRLS) algorithm associated with the robust reweighting and $L_p$-norm approximation respectively, as described by @alg:irls-r.

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
        Assign[$vv(γ)_k$][$bb(1)_(N_k) quad ∀ k ∈ cal(K)$]
        
        For(
          $i = 1 " to " N_"iter"$,
          {
            Comment[Construct system matrix and right-hand side, $mm(A) vv(m)_"new" = vv(b)$]
            LineComment(
              Assign[$vv(w)_(m k)$][$(vv(γ)_k^2 + δ_k^2)^(p_k\/2 -1) quad ∀ k$], [Compute Ekblom norm weights]
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
              Assign[$Δ_2$][$norm(vv(m) - vv(m)_"new")_2$],
              [Convergence criterion]
            )

            // Update model
            Comment[Update model]
            Assign[$vv(m)$][$vv(m)_"new"$]

            // Convergence check
            LineComment(
                If($Δ_2 < τ$, {
                Break
              }),
              [Stop if converged]
            )

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
          }
        )
        Return[$vv(m), vv(w), vv(γ)_(k∈cal(K))$]
      }
    )
  }
) <alg:irls-r>


We note that $tilde(mm(G))_(α p)^(-1)$ is symmetric by the diagonality of $mm(W), mm(W)_"mp"$.
Thus, by identity $Cov(mm(A) vv(b)) = mm(A) Cov(vv(b)) mm(A)^TT$ and the symmetry of  $tilde(mm(G))_(α p)^(-1)$, we obtain covariance matrix
$
Cov(vv(m)_(r α p))
&= mm(G)_(α p)^(-1) Cov(vv(d)) (mm(G)_(α p)^(-1))^TT\
&= tilde(mm(G))_(α p)^(-1) mm(G)^TT mm(W) Cov(vv(d)) (tilde(mm(G))_(α p)^(-1) mm(G)^TT mm(W))^TT\
&= tilde(mm(G))_(α p)^(-1) mm(G)^TT mm(W) Cov(vv(d)) mm(W) mm(G) tilde(mm(G))_(α p)^(-1).\
$ <eq:cov-matrix>

By applying the treatment used to obtain @course_book[eq. 35], the resolution matrix, $mm(R)_m$, may be computed using
$
mm(R)_m = mm(G)_(α p)^(-1) mm(G) = tilde(mm(G))_(α p)^(-1) mm(G)^TT mm(W) mm(G)
$ <eq:res-matrix>
where we recall that the limiting case $mm(R)_m = mm(I)$ identifies the optimal case where the data is fully resolved by the inversion.

== Unregularised Least Squares Solution



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

== Robust $L_1$-Regularised Solution

For the $L_1$-regularised solution with robust weights, we require the full machinery described in @sec:theory-math and must resort to iteratively approaching the correct weights $mm(W), mm(W)_"m1"$ using the IRLS algorithm given in @alg:irls-r with the norm order given by $p=1$.


= Finding an $L_2$ least-squares solution
The L-curve obtained for the non-robust L2-regularized inversion shows a clear corner, and thus the $alpha^2$ here is 
#align(center)[
  $ alpha^2 = 4.642 times 10^(-9). $
]
This value gives a reasonable balance between the data misfit and the model norm. As the $alpha^2$ gets smaller, the regularization becomes weaker and the inversion starts to fit short-wavelength structure more strongly. As the $alpha^2$ gets larger, the regularization becomes stronger and overly smooth, meanwhile, the small-field is suppressed.

In the core mantle boundary (CMB) field map obtained for the $alpha^2$ shows coherent large-scale flux patches. The main structures in this figure are broad and spatially smooth. This figure also shows the negative flux dominating parts of the northern hemisphere, especially in the high latitudes and positive flux dominating parts of the southern hemisphere, especially in the southern Atlantic, Africa and Indian Ocean areas. In the southern Atlantic, there is a circle feature shown positive flux in this map. In the power spectrum of estimated geomagnetic field at core mantle boundary (CMB) for the corresponding $alpha^2$ it shows that it is dominated by the lowest degrees and decreases toward higher degree. It means that the large-scale field is expected to be constrained better for the SWARM data and in the meanwhile, the smaller scales are more sensitive to instability and noise. The spectral decay is relatively smooth and the higher degree power is clearly dampened after 13 degree.

To explore how changing the value of $alpha^2$ affects solution, we produced the magnetic field maps and the power spectrum for  $alpha^2 = 4.642 times 10^(-11)$,
$4.642 times 10^(-10)$, $4.642 times 10^(-8)$, and
$4.642 times 10^(-7)$. When the $alpha^2$ becomes smaller the CMB field map becomes much more oscillatory and more alternating smaller patches appear. And for smaller $alpha^2$, the large power remains at higher degrees. The opposite trend occurs for the larger $alpha^2$. We have observed that in the CMB field map, the overall feature becomes smoother and the higher-degree part of the spectrum is strongly reduced, which may due to the regularization dominates the inversion too strongly, and the resolution of the field map is decreased. Among all the field maps and the power spectrums, the $alpha^2 = 4.642 times 10^(-9)$ present the best result.
= model covariance matrix and resolution matrix solution

= comparison between L2 and L1 solution

= convincing evidence for the reversed flux features at the core surface
#bibliography("references.bib")