# mmcif 0.1.2
* bugs are fixed in `mmcif_pd_univariate()`, `mmcif_pd_univariate()`, and
  `mmcif_pd_cond()` so they always work with `left_trunc`. They could give an 
  error in the previous version.  

# mmcif 0.1.1
* Issues using ATLAS (and possibly other Lapack versions) have been solved. The
  problem was a non-unique Q matrix from a QR decomposition (so the results 
  were not wrong just different).
* Functions are added to compute marginal figures post estimation. See
  `?mmcif_pd_univariate` and `?mmcif_pd_cond`.

# mmcif 0.1.0
* The first release on CRAN.
