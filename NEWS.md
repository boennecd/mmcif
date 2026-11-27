# mmcif 0.1.2

* Fix errors when `mmcif_pd_univariate()`, `mmcif_pd_bivariate()`, and
  `mmcif_pd_cond()` receive left-truncation times.
* Require R 4.1.0 or later for the syntax used in the examples.
* Fix the link to `alabama::auglag()` in the help for `mmcif_fit()`.

# mmcif 0.1.1
* Issues using ATLAS (and possibly other Lapack versions) have been solved. The
  problem was a non-unique Q matrix from a QR decomposition (so the results 
  were not wrong just different).
* Functions are added to compute marginal figures post estimation. See
  `?mmcif_pd_univariate` and `?mmcif_pd_cond`.

# mmcif 0.1.0
* The first release on CRAN.
