#import "Exam/common/exports.typ": *

Some of these are very nit-picky, so feel very justified in shrugging and ignoring.

= Generally
- Having italicised differential is non-standard (ISO-80000 recommends upright), $d x → dif x$
- This might cause problems since it clashes a bit with the $d$ used for data, convention in document seems to be:
  - serif, italicised for scalars
  - serif, upright for elements of vectors
  - bold, serif, upright for vectors
- Maybe consistently use serif+italicised for scalars and serif+bold for vectors? Then operators like differential can get upright serif.

= Lecture notes, Chapter 2
- The notation $ M a x !$ and $M i n !$ seems non-standard (but maybe not in geophysics?) and should at least be in upright serifs, the kerning is a little jarring in the current form.
  - I would suggest something like $hat(vv(m)) = arg max_vv(m) p(vv(d)|vv(m))$ for the expressions in Section 2.2, for example.
- $pdv(ϕ, m_k) = 0$ condition is a little hazardeous without explicitly stating that $ϕ$ is at least differentiable and convex (formally I guess also lower-bounded and coercive is needed for gradient methods, but convexity I think is more commonly violated)
- ML estimator:
  - Currently states: "The ML estimator is, among all possible estimators, the one with the smallest variance. In addition the ML estimator is consistent and unbiased."
  - ML estimator is only lowest variance in asymptotic limit
  - ML estimator is only unbiased in asymptotic limit


= Lecture notes, Chapter 3

The notation $G^(-α)$ which is very easy to interpret as a matrix exponential, but I reckon it is also intuitive notation once you get it. Maybe elaborate a bit more on the notation.

= Lecture notes, Chapter 4
- $Δ m^1$ is easy to confuse with exponent, maybe $Δ m^((1))$?
  - I understand why you don't want it to be subscript since this would be easy to confuse with the indexing of the model parameters, alternatively express as a product: $p_i = p(vv(m)_i) ∏_(j=1)^M Δ m_j$ or with $m_j → m^((j))$, could even do $Δ cal(m)$ to denote that these are discretisations of the model parameter space, $cal(M)$.
- If adding `nanopinv`, consider briefly covering Preconditioned Crank-Nicolson (pCN)?
- Lecture notes don't feature anything about covariance functions and Gaussian processes, though this is covered well in lecture slides.

= Lecture notes, Chapter 5
- Mangled math near Eq (82) ("where j = 1, ...M"), lacks comma and mangled parentheses
- Near Eq (85): "where w is a vector containing all the weight and bias parameters"
  - Calling this a vector is a bit dodgy, at least stipulate that it is a vector of matrices.
  - Since the matrices don't have the same shape necessarily, they don't belong to the same vector space.
  - I think it would be better to refer to it as an _array_ or _tuple_ of parameter matrices, or something like that.
- This section is a bit long, consider maybe omitting Bayesian Neural Networks?

= Lecture: L7
- Curse of dimensionality: somewhat long derivation for something fairly conceptually simple and intuitive, though admittedly the way it was presented was very nice.

= Exercise 4.2
- Change: `np.NINF` to `-np.inf`, `np.Inf` to `np.inf` to be Numpy 2 compatible

= Lecture: L8
- S.3: Slides use the notation of "∫ dz INTEGRAND", which seems to be conventional in some parts of physics — I don't personally like it, but this is very much just a nit-pick
- S.4: $ρ$ used here, $d$ used on previous slide. I guess $ρ$ is already used for the prior
- S.3: Maybe add $x, z$ coordinate system
- I think in general the second part of the lecture would have benefitted from a shorter treatment. Particularly, I think both the conditional probability, sequential sampling and Gibbs sampling are fairly intuitive and can safely be covered fairly concisely.

= Advantages of switching to `nanopinv` over `pinv`
- `pinv` is very poorly packaged and distributed at the moment. This is much improved in `nanopinv`, which following Python packaging standards.
- `pinv` does many things that are very dangerous and non-idiomatic in Python. As an example, it very often uses mutable default arguments (lists, dicts), which gives rise to unexpected behaviour in Python.
- Further, it relies heavily on global state, which is prone to unexpected behaviour.
- It appears that `pinv` does not work with recent Python versions?
- The distributed files contain many temporary files (`.py~`, `.#pygstat.py`, `__pycache__/`) that are not suitable for distribution and anecdotally also caused problems when unzipping the files on Windows.
- `pinv` has incorrect use of docstrings
- `pinv` has imports within function calls


= Lecture: L10 (ML intro)
- I really like [Goodhart's law](https://en.wikipedia.org/wiki/Goodhart%27s_law) when thinking of train/test dataset splits.
  - I think it is mostly applied in social sciences and economics, but I think it is sufficiently profound that it could and should be mentioned in ML.

= Assignment 2
- Equation (1): $"exp"$ instead of $"Exp"$?
- Very nit-picky: $"TWIP"$ instead of $T W I P$
  - Kerning very bad otherwise, also upright indicates operator or single object.
  - Arguably $"TWIP"(bold(underline(ϕ)))$
- Mega nit-picky: $d x$ should be $dif x$ if differential, but here we have discretised integral so maybe $Δ x$ is the most appropriate notation?
- More information about where experimental data comes from

= `nanopinv` TODO
- I would like to add a few convenience methods and clean the code up a little more and provide a good example notebook. I will hopefully get a chance to have a peek at this during the 3-week June semester. Feel free to poke me.

= Distribute environment for course
- The course could benefit from distributing an environment for it (the current list on Learn works, but is complicated and could be much simplified by using `pixi` (if mixing Conda and Pypi packages), or just `uv` if all the used software is `pypi` distributed.)
- `nanopinv` could be included in the environment specification
- Installing the environment would be as simple as installing pixi and running `pixi install`.