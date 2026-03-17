= TODO

- ??: Abstract
- JK: Style report
  - Smaller caption text
  - Smaller caption gap
// - ??: Discussion on ill-posedness (Hadamard sense)
//   - $1/r$ term
// - JK: Merge math theory sections
// - JK: Huber and MAD
- ??: Some recap/discussion of the assumptions used and how reasonable they are
- ML: Write-up L2, (i)
  - How alpha changes solution (maybe 3 side-by-side plots?)
  - Power spectrum
- JK: Write-up L1/Robust, (ii)
  - Covariance/resolution
- ML: (iv)
- ML: Introduction of mission and experimental design
- ??: Discussion around patch size related to discretisation
- ??: Comment on inadequacy of provided L curve estimate
  - Mention that it is a heuristic anyway
- JK: TODO elaborate a bit more on iterative algorithm

= Questions for TA/Chris
- We do not have access to a good estimate of the dispersion of the "untainted" Gaussian in this case, so presumably we need to iteratively find this using the residuals of the model?
- Lecture notes mention MAD as Mean Absolute Deviation and give formula using mean, but literature exclusively suggests MAD, Median Absolute Deviation.
- How are we suppose to approximate the data covariance matrix $Cov(vv(d))$? Do we assume the scale parameter may be found using MAD and that the samples are uncorrelated to get a diagonal matrix, or do we somehow use the robust weights as the covariance matrix?
- Will the characteristic blot size not be biased by $N$?
  - A: Can be investigated, but notice shape of power spectrum (peak is not necessarily maxed out at highest frequencies)
- Check understanding: Discrete jumps in covariance stems from increments of $n$ (order of spherical harmonic)
- Check understanding: When to use weights, when not to
- Check understanding: Cholesky OK, since PD
