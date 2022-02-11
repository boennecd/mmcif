load("2D-data.RData")

# setup the object to do the computation
library(mmcif)
comp_obj <- mmcif_data(
  ~ a + b, dat, cause = cause, time = time, cluster_id = cluster_id,
  max_time = delta)

# we need to find the combination of the spline bases that yield a straight
# line. You can skip this
comb_slope <- local({
  boundary_knots <- comp_obj$spline$boundary_knots
  pts <- seq(boundary_knots[1], boundary_knots[2], length.out = 1000)
  lm.fit(cbind(1, comp_obj$spline$expansion(pts)), pts)$coef
})

# assign a function to compute the log composite likelihood
library(fastGHQuad)
ghq_data <- with(gaussHermiteData(5L), list(node = x, weight = w))

ll_func <- \(par, n_threads = 1L)
mmcif:::mmcif_logLik(
  comp_obj$comp_obj, par = par, ghq_data = ghq_data, n_threads = n_threads)

# the log composite likelihood at the true parameters
coef_traject_spline <-
  rbind(t(coef_traject[1, ] %o% comb_slope[-1]),
        coef_traject[2, ] + comb_slope[1] * coef_traject[1, ],
        coef_traject[-(1:2), ])
true_values <- c(coef_risk, coef_traject_spline, Sigma)
ll_func(true_values)

# check the time to compute the log composite likelihood
bench::mark(
  `one thread` = ll_func(n_threads = 1L, true_values),
  `two threads` = ll_func(n_threads = 2L, true_values),
  `three threads` = ll_func(n_threads = 3L, true_values),
  `four threads` = ll_func(n_threads = 4L, true_values))
