#import "@preview/codly:1.3.0": *
#import "@preview/cetz:0.4.2"
#import "@preview/unify:0.7.1"
#import "@preview/algorithmic:0.1.0"
#import "@preview/physica:0.9.8": curl, grad, tensor, pdv, dv, TT
// #import algorithmic: algorithm

#import "preamble_funcs.typ": mref, mathformatter
#import "style.typ": style

#show: style

#let vv = mathformatter(underbars: 0, bold: true, upright: false)
#let mm = mathformatter(underbars: 0, bold: false, upright: true)
#set math.mat(delim:"[")

#let wider = h(3em)

= TODO

- Change from semi-discrete formulation to functional analysis formulation

= Introduction

This report describes methods and implementations from the field of inverse problems and regularisation theory as applied to the determination of the radial magnetic field of Earth's surface and subsurface using data from ESA's _Swarm_ satellite mission recorded between XXX - YYY.

= Experimental Design



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

Before an inverse problem can be solved, we must present the corresponding forward model of the system, which may be formulated in the linear case as
$
  vv(d) = cal(G) m + ε
$

where $cal(G) : cal(M) → cal(D)$ is the forward _linear operator_ which maps model parameters $m$ from the continuous model space, $cal(M)$, to measurement vectors $vv(d)$ in the data space, $cal(D) = ℝ^(N_cal(D))$ subject to measurement noise $ε ∈ cal(D)$, which notably may not be assumed to be Gaussian.

// TODO: MOVE BELOW TRUNCATION
Having started out with a continuous formulation of the forward problem, we will discretise the operator $cal(G)$ in order to obtain a system that may be solved numerically.

In order to determine the forward operator $cal(G)$, we must first consider carefully the physical model of the system, to which end we start from the Ampere-Maxwell law and the assumption that the field measurements, which are performed in orbit using the _Swarm_ satellites, are taken in the absence of free charges and currents nearby, that is, we assert that it the magnetic field is curl-free in the region of interest:
$
  ∇ × vv(B) = μ_0 (vv(J) + ε_0 pdv(vv(E), t)) = vv(0)
$

Thus we may define a _scalar magnetic potential_, $V$, such that
$
  vv(B) = -∇V
$ <eq:scalar_potential>

and by deeming it unlikely that we run into any magnetic monopoles, we conclude that the field is also divergence-free, $∇ ⋅ vv(B) = 0$, we obtain the Laplace equation for the scalar potential,
$
  ∇ ⋅ vv(B) = -∇^2 V = 0.
$

If we make the simplifying assumption that a spherical boundary field of radius $a$ exists which delineates the electromagnetic sources from empty space, we find that the continuous model space $cal(M)$ belongs to the Hilbert space defined on the unit sphere, that is $cal(M) ∈ L^2(S^2)$. It should be noted that the assumption of a spherical boundary field is reasonable, while the assumption that no sources, here understood to be magnetic moments and electrical currents, exist outside the boundary is more difficult to justify given the magnetisation of the Earth's crust and charge carriers in the ionosphere, a discussion which we will return to later.

With this in mind, we may choose the spherical harmonics $Y_n^m (ϑ, λ)$ as a convenient, orthogonal basis for the model space and continue our analysis in spherical coordinates $(r, ϑ, λ)$, in which we recall the Laplacian as
$
  ∇^2 V = 1/r^2 pdv(,r)(r^2 pdv(V, r)) + 1/(r^2 sin ϑ) pdv(, ϑ)(sin ϑ pdv(V, ϑ)) + 1/(r^2 sin^2 ϑ) pdv(V, λ, 2) = 0
$

Where separation of variables and the assumption of interior sources yields the solution @geomath[p. 507, eq. 13]
$
  V(r, ϑ, λ)
  &=a ∑_(n=1)^∞ (a/r)^(n+1) ∑_(m=0)^n (g_n^m cos m λ + h_n^m sin m λ ) P_n^m (cos ϑ)\
$ <eq:lapl_sol>

in which $P_n^m$ are the associated Legendre polynomials and we have used the spherical harmonics with
$
Y_(n,"even")^m = cos (m λ) P_n^m (cos ϑ)
wider
Y_(n,"odd")^m = sin (m λ) P_n^m (cos ϑ)
$

The parameters $g_n^m, h^n_m$ of @eq:lapl_sol are referred to as the _Gauss coefficients_ are determine the model, that is
$
m = {g_n^m, h_n^m} ∈ cal(M)
wider n ∈ ℕ_+, quad m ∈ {0, ..., n}
$
and the inverse problem is simply to determine these given the field measurements.

We do not, however, have access to the potential $V$ directly, which is introduced more as a mathematical convenience than a physical quantity. Instead, we measure the radial component of its gradient as given by @eq:scalar_potential,
$
  B_r (r, ϑ, λ)
  = - pdv(V, r)
  = ∑_(n=1)^∞ (n+1) (a/r)^(n+2) ∑_(m=0)^n (g_n^m cos m λ + h_n^m sin m λ ) P_n^m (cos ϑ)\
$ <eq:B_r>

By identifying this as the data and discretising the model and associated forward operator $cal(G)$ by truncation of the spherical harmonic expansion to degree $N$, we may obtain the discrete forward model
$
  B_r (r, ϑ, λ) = ∑_(n=1)^N (n+1) (a/r)^(n+2) ∑_(m=0)^n (g_n^m cos m λ + h_n^m sin m λ ) P_n^m (cos ϑ)\
$
which when expressed as a system of equations allows us to introduce the discrete forward operator as a matrix, $mm(G)$, which is built as follows:
$
  vv(d) = mm(G) vv(m)\
  vec(
    B_r (vv(x)_1),
    B_r (vv(x)_2),
    ⋮,
    B_r (vv(x)_(N_cal(D))),
  )
  =
  mat(
    // g=cos, n=1, m=0
    Ξ_1^0 (vv(x)_1),
    // h=sin, n=1, m=0
    ξ_1^0 (vv(x)_1),
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
    ξ_1^0 (vv(x)_2),
    Ξ_1^1 (vv(x)_2),
    ξ_1^1 (vv(x)_2),
    …,
    Ξ_N^N (vv(x)_2),
    ξ_N^N (vv(x)_2)
    ;
    ⋮,⋮,⋮,⋮,,⋮,⋮
    ;
    Ξ_1^0 (vv(x)_(N_cal(D))),
    ξ_1^0 (vv(x)_(N_cal(D))),
    Ξ_1^1 (vv(x)_(N_cal(D))),
    ξ_1^1 (vv(x)_(N_cal(D))),
    …,
    Ξ_N^N (vv(x)_(N_cal(D))),
    ξ_N^N (vv(x)_(N_cal(D)))
    ;

  )
  vec(
    g_1^0,
    h_1^0,
    g_1^1,
    h_1^1,
    ⋮,
    g_N^N,
    h_N^N,
  )
$

Where we have introduced some shorthand notation
$
vv(x)_i &= (r_i, ϑ_i, λ_i)\
Ξ_n^m (r, ϑ, λ) &= (n+1) (a/r)^(n+2) cos (m λ) P_n^m (cos ϑ)\
ξ_n^m (r, ϑ, λ) &= (n+1) (a/r)^(n+2) sin (m λ) P_n^m (cos ϑ).
$

In conclusion, we can build the discrete forward operator $mm(G) ∈ ℝ^(N_cal(D) × N_cal(M)), N_cal(M) = N(N+3)$ simply by evaluating the convenience functions $Ξ_n^m$ and $ξ_n^m$ at the measurement locations, ${vv(x_i)}_(i=1)^(N_cal(D))$.



= Mathematical Treatment of the Inverse Problem

We may summarise the previous section as a forward model in which the volume of the Earth has been been discretised into a finite number of elements each of which is assigned a value of the magnetic scalar potential, $V$.

---
We consider the general linear forward model
$
  vv(d) = mm(G) vv(m)
$

Where $vv(d)$ is a vector of radial magnetic field measurements and $vv(m)$ is a vector of the coefficients of the spherical harmonic expansion of the scalar potential, with $mm(G)$ being the


#set heading(numbering: "i.1)")
= Finding an $L_2$ least-squares solution

#bibliography("references.bib")