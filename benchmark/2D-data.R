load("2D-data.RData")

# setup the object to do the computation
library(mmcif)
comp_obj <- mmcif_data(
  ~ a + b, dat, cause = cause, time = time, cluster_id = cluster_id,
  max_time = delta)

# assign a function to compute the log composite likelihood
library(fastGHQuad)
ghq_data <- with(gaussHermiteData(5L), list(node = x, weight = w))

ll_func <- \(par, n_threads = 1L)
mmcif:::mmcif_logLik(
  comp_obj$comp_obj, par = c(coef_risk, coef_traject, Sigma),
  ghq_data = ghq_data, n_threads = n_threads)

# the log composite likelihood at the true parameters
ll_func(c(coef_risk, coef_traject, Sigma))

# check the time to compute the log composite likelihood
bench::mark(
  `one thread` = ll_func(n_threads = 1L, c(coef_risk, coef_traject, Sigma)),
  `two threads` = ll_func(n_threads = 2L, c(coef_risk, coef_traject, Sigma)),
  `three threads` = ll_func(n_threads = 3L, c(coef_risk, coef_traject, Sigma)),
  `four threads` = ll_func(n_threads = 4L, c(coef_risk, coef_traject, Sigma)))
