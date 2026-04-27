= Lecture notes, Chapter 2

- The notation $ M a x !$ seems non-standard and should at least be in upright serifs, the kerning is quite strange in the current form.
  - I would suggest something like $hat(m) = arg max_m p(d | m)$ for the expressions in Section 2.2, for example.
= Lecture notes, Chapter 3

I really don't like the notation $G^(-α)$ which is very easy to interpret as a matrix exponential. Is there a (very) good reason why this notation is used?

= Lecture: L7
- Curse of dimensionality: somewhat long derivation for something fairly conceptually simple and intuitive, though admittedly the way it was presented was very nice.

= Exercise 4.2
- Change: `np.NINF` to `-np.inf`, `np.Inf` to `np.inf` to be Numpy 2 compatible

= Lecture: L8
- S.3: Slides use the, in my opinion, dreadful notation of "∫ dz INTEGRAND"
- S.4: $ρ$ used here, $d$ used on previous slide. I guess $ρ$ is already used for the prior
- S.3: Maybe add $x, z$ coordinate system
- I think in general the second part of the lecture would have benefitted from a shorter treatment. Particularly, I think both the conditional probability, sequential sampling and Gibbs sampling are fairly intuitive and can safely be covered fairly concisely.

= Use of `pinv`
- The distribution and source code of `pinv` is almost criminal
- It does many things that are considered very very bad practice, if not outright dangerous in Python. For example, it has mutable default arguments to functions (lists), which gives unexpected behaviour and is universally avoided
- It does not follow any packaging standards
- It does not work for any modern Python versions
- I think this course would be much better off avoiding any use of this package, perhaps implementing the necessary functionality ourselves or using larger, more widely-used packages that may be available.
- The distributed files contain many temporary files (`.py~`, `.#pygstat.py`, `__pycache__/`) that are not suitable for distribution and anecdotally also caused problems when unzipping the files on Windows.
- Incorrect use of docstrings
- Imports within function calls
-

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


= Question Assignment 2
- Does reconstruction look somewhat okay?
- Is it okay to have nugget in the covariance model?
- Is the TWIP number correct?
  - Seems to be quite far off when only considering the water table?
  - Have tested both SKFMM and FSM implementation
