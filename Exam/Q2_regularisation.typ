#import "common/common.typ": *

#show: setup.with(2, "Regularisation")
#show: show-page1

#let questions = [
= Questions
- What is *meant by regularization* in the context of inverse problems and *why is often
necessary* in problems in Earth and Space physics?
- What is *L1-norm regularization*?
  - _Sparsity promoting, abs value_
- How can this be *implemented* and when should it be used?
  - _IRLS and Ekblom relaxation_.
]
#questions

= Approach
- Outline Hadamard conditions for well-posedness:
  - Existence, Uniqueness, Stability
- Regularisation can help reduce ill-posedness by biasing with prior information
  - Example — Assignment 1: Non-unique and unstable due to attenuation of higher frequencies

= Setup
- Loss function: $ϕ(vv(m))$
- Regulariser: $R(vv(m))$
- Regularisation parameter: $α^2$

= Tikhonov Regularisation
- Loss consists of _data fidelity_ and _regularisation_ terms: $ϕ(vv(m)) = E(vv(d), vv(m)) + α^2 R(vv(m))$
- For regulariser $R(vv(m)) = ||mm(H) vv(m)||_2$:
  - #inline-block-math(<eq:gen-tik>)[$hat(vv(m))_α = (mm(G)^T mm(G) + α^2 mm(H)^T mm(H))^(-1) mm(G)^T vv(d)$]
  - _Generalised Tikhonov_, ordinary if $mm(H)=mm(I)$
  - This gives solution through Moore-Penrose pseudoinverse:
  - Unregularised: $hat(vv(m)) = mm(G)^† vv(d) quad "where" quad mm(G)^† = (mm(G)^T mm(G))^(-1) mm(G)^T$
  - Regularised: $hat(vv(m))_α = mm(G)^(-α) vv(d) = (mm(G)^T mm(G) + α^2 mm(H)^T mm(H))^(-1) mm(G)^T vv(d)$

= Regularisation Generally
- Regularisation imposes _prior information_.
  - Using Gaussian assumption, can incorporate _a priori_ model mean and variance directly into closed-form expression. See Sec. 3.2 in lecture notes.
  - Assumption of homoskedasticity and independence yields @eq:gen-tik

= SVD
- Alternative method of regularisation is Truncated SVD (TSVD)
- Tikhonov has _smooth cut-off_, TSVD has _hard cut-off_ of singular values
- $vv(m)_α = ∑_(i=1)^k f_i (U_(⋅,i)^T vv(d))/(s_i) mm(V)_(⋅,i) quad "where" quad f_i = s_i^2/(s_i^2 + α^2)$

= Covariance and Resolution
- Model covariance: $mm(C)_hat(vv(m))_α = mm(G)^(-α) mm(C)_vv(d) (mm(G)^(-α))^T$
  - Diagonal gives variance of model parameters
- Model resolution: $mm(R)_hat(vv(m))_α = mm(G)^(-α) mm(G)$
  - Diagonal gives resolution of model parameters (1 fully resolved, 0 not resolved at all)
  - Trace is effective number of parameters

= Regularisation Parameter
+ Discrepancy Principle: choose $α$ s.t. $||vv(d) - mm(G) hat(vv(m))_α||_2^2 = δ σ^2$
  - $σ^2$ is _a priori_ data variance, $δ=1+ε$ is safety factor
  - Note: variance must be known _a priori_
+ L-Curve: Plot of $log R_α (vv(hat(m))_α)$ vs $log E (vv(d), vv(hat(m))_α)$
  - Choose $α$ at corner of L-curve
  - Note: Corner can be fairly smooth, especially if ill-posed
  - Plot looks like `\_`
+ Generalised Cross-Validation:
  - Intuition: Pick model that predicts omitted data point well
  - Generalised: assumes white, uncorrelated noise and that $(mm(G) mm(G)^(-α))_(k, k)$ is well-approximated by its average, $Tr(mm(G) mm(G)^(-α))\/N$
  - $g(α^2) = N ∑_(k=1)^N [ ((mm(G) vv(hat(m))_α)_k - vv(d)_k) / (Tr(mm(I) - mm(G) mm(G)^(-α))) ]^2$
  - Plot looks like `----v/`

= Sparse Regularisation
- We can impose sparsity through $p-norm$ regularisation
- The lower $p ∈ [0, 2]$ is, the more regularisation promotes sparsity
  - Typical choice: $L_1$-norm
- Often $L_p$-norms are implemented through _Ekblom relaxation_:
  - $hat(vv(m))_(m p) = (mm(G)^T mm(W) mm(G) + α^2 mm(H)^T mm(W)_(m p) mm(H))^(-1) mm(G)^T mm(W) vv(d)$
  - where $δ^2$ is positive and smaller than $γ_i^2 = ((mm(H) vv(m))_i)^2$
  - $mm(W)$ is _data weight matrix_ from robust estimation, $mm(W) = mm(I)$ if not robust weights
  - $mm(W)_(m p)$ is _Ekblom weight matrix_, diagonal with $w_i = (γ_i^2 + δ^2)^(p\/2 - 1)$

= IRLS
- Similarly to robust estimation:
  - Initialise $γ_i = 0$
  - Estimate $hat(vv(m))_(m p) = (mm(G)^T mm(W) mm(G) + α^2 mm(H)^T mm(W)_(m p) mm(H))^(-1) mm(G)^T mm(W) vv(d)$
  - Update $γ_i^2 = ((mm(H) vv(m))_i)^2$
  - Repeat until convergence

= Helpful intuition
- Draw some contour of constant misfit in a plot like @fig:lp-circle, This is disambiguated by the regularisation by where the expanding norm circle first intersects the misfit contour.


#let fig-lp-circle = [
  #figure(
    image("lp_norms.png"),
    caption: [Unit circles of various $L_p$ norms.],
  ) <fig:lp-circle>
]

#let fig-lp-penalty = [
  #figure(
    image("norm_penalties.png"),
    caption: [Penalty functions for various $L_p$ norms.],
  ) <fig:lp-penalty>
]

#v(1fr)
#grid(
  columns: 2,
  fig-lp-circle,
  fig-lp-penalty,
)

#show: show-page2

#questions
#LINE

#unwritten[
Ill-posedness: does not fit Hadamard definition of well-posedness:
- Existence
- Uniqueness
- Stability
]

#unwritten[
  Regularisation: Reducing ill-posedness by biasing solution with prior information.
]

#unwritten[Loss: ][$
  ϕ(vv(m))
  &= E(vv(d), vv(m)) + α^2 R(vv(m))\
  &= ||vv(d) - mm(G) vv(m)||_2^2 + α^2 ||mm(H) vv(m)||_p
$]

#unwritten[Tikhonov if $p=2$, ordinary if $mm(H) = mm(I)$, otherwise generalised]
#unwritten[Gen. Tikh, closed-form sol.: ][$
  hat(vv(m))_α = (mm(G)^T mm(G) + α^2 mm(H)^T mm(H))^(-1) mm(G)^T vv(d)
$]

#unwritten[Will skip Gaussian closed-form, TSVD, Cov+Res matrices]

#unwritten[
Reg. param selection:
- Discrepancy principle: $||vv(d) - mm(G) hat(vv(m))_α||_2^2 = δ σ^2$
  - Note: variance must be known _a priori_
- L-curve: corner of plot of $log R_α (vv(hat(m))_α)$ vs $log E (vv(d), vv(hat(m))_α)$
- GCV: $g(α^2) = N ∑_(k=1)^N [ ((mm(G) vv(hat(m))_α)_k - vv(d)_k) / (Tr(mm(I) - mm(G) mm(G)^(-α))) ]^2$
]

#unwritten[
  #box(width: 100%)[
    #fill-in-drawing(16em, 10em, caption: [Shape of regularisation parameter selection curves])
  ]
]

#colbreak(weak: true)
#unwritten[
= Sparse Regularisation
Lower $L_p$ norms give more sparsity due to 'pokiness' of unit circle, see below.
]

#unwritten[
  #box(width: 100%)[
    #fill-in-drawing(16em, 16em, caption: [
      Unit circles of various $L_p$ norms ($p ∈0.5, 1, 2, ∞$).\
      Draw misfit contour of 2D model.
      ])
  ]
]

#unwritten[
  Implemented via Ekblom relaxation and IRLS
]
$
  hat(vv(m))_(m p) &= (mm(G)^T mm(W) mm(G) + α^2 mm(H)^T mm(W)_(m p) mm(H))^(-1) mm(G)^T mm(W) vv(d)\
  w_i &= (γ_i^2 + δ^2)^(p\/2 - 1), quad γ_i^2 = ((mm(H) vv(m))_i)^2, quad 0 < δ^2 < γ_i^2
$

#unwritten[
IRLS:
+ Initialise $γ_i = 0$
+ Estimate $hat(vv(m))_(m p) = (mm(G)^T mm(W) mm(G) + α^2 mm(H)^T mm(W)_(m p) mm(H))^(-1) mm(G)^T mm(W) vv(d)$
+ Update $γ_i^2 = ((mm(H) vv(m))_i)^2$
+ Repeat until convergence
]
