
# MMCIF: Mixed Multivariate Cumulative Incidence Functions

[![R-CMD-check](https://github.com/boennecd/mmcif/workflows/R-CMD-check/badge.svg)](https://github.com/boennecd/mmcif/actions)

This package provides an implementation of the model introduced by
Cederkvist et al. (2018) to model within-cluster dependence of both risk
and timing in competing risk.

## Installation

The package can be installed from Github by calling

``` r
library(remotes)
install_github("boennecd/mmcif")
```

The code benefits from being build with automatic vectorization so
having e.g.  `-O3` in the `CXX17FLAGS` flags in your Makevars file may
be useful.

## The Model

The conditional cumulative incidence functions for cause
![k](https://render.githubusercontent.com/render/math?math=k "k") of
individual
![j](https://render.githubusercontent.com/render/math?math=j "j") in
cluster
![i](https://render.githubusercontent.com/render/math?math=i "i") is

![\\begin{align\*} F\_{kij}(t\\mid \\vec u\_i, \\vec\\eta\_i) &= \\pi\_k(\\vec z\_{ij}, \\vec u\_i) \\Phi(-\\vec x\_{ij}(t)^\\top\\vec\\gamma\_k - \\eta\_{ik}) \\\\ \\pi\_k(\\vec z\_{ij}, \\vec u\_i) &= \\frac{\\exp(\\vec z\_{ij}^\\top\\vec\\beta\_k + u\_{ik})}{1 + \\sum\_{l = 1}^K\\exp(\\vec z\_{ij}^\\top\\vec\\beta\_l + u\_{il})} \\\\ \\begin{pmatrix} \\vec U\_i \\\\ \\vec\\eta\_i \\end{pmatrix} &\\sim N^{(2K)}(\\vec 0;\\Sigma).\\end{align\*}](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign%2A%7D%20F_%7Bkij%7D%28t%5Cmid%20%5Cvec%20u_i%2C%20%5Cvec%5Ceta_i%29%20%26%3D%20%5Cpi_k%28%5Cvec%20z_%7Bij%7D%2C%20%5Cvec%20u_i%29%20%5CPhi%28-%5Cvec%20x_%7Bij%7D%28t%29%5E%5Ctop%5Cvec%5Cgamma_k%20-%20%5Ceta_%7Bik%7D%29%20%5C%5C%20%5Cpi_k%28%5Cvec%20z_%7Bij%7D%2C%20%5Cvec%20u_i%29%20%26%3D%20%5Cfrac%7B%5Cexp%28%5Cvec%20z_%7Bij%7D%5E%5Ctop%5Cvec%5Cbeta_k%20%2B%20u_%7Bik%7D%29%7D%7B1%20%2B%20%5Csum_%7Bl%20%3D%201%7D%5EK%5Cexp%28%5Cvec%20z_%7Bij%7D%5E%5Ctop%5Cvec%5Cbeta_l%20%2B%20u_%7Bil%7D%29%7D%20%5C%5C%20%5Cbegin%7Bpmatrix%7D%20%5Cvec%20U_i%20%5C%5C%20%5Cvec%5Ceta_i%20%5Cend%7Bpmatrix%7D%20%26%5Csim%20N%5E%7B%282K%29%7D%28%5Cvec%200%3B%5CSigma%29.%5Cend%7Balign%2A%7D "\begin{align*} F_{kij}(t\mid \vec u_i, \vec\eta_i) &= \pi_k(\vec z_{ij}, \vec u_i) \Phi(-\vec x_{ij}(t)^\top\vec\gamma_k - \eta_{ik}) \\ \pi_k(\vec z_{ij}, \vec u_i) &= \frac{\exp(\vec z_{ij}^\top\vec\beta_k + u_{ik})}{1 + \sum_{l = 1}^K\exp(\vec z_{ij}^\top\vec\beta_l + u_{il})} \\ \begin{pmatrix} \vec U_i \\ \vec\eta_i \end{pmatrix} &\sim N^{(2K)}(\vec 0;\Sigma).\end{align*}")

where there are
![K](https://render.githubusercontent.com/render/math?math=K "K")
competing risks. The
![\\vec x\_{ij}(t)^\\top\\vec\\gamma\_k](https://render.githubusercontent.com/render/math?math=%5Cvec%20x_%7Bij%7D%28t%29%5E%5Ctop%5Cvec%5Cgamma_k "\vec x_{ij}(t)^\top\vec\gamma_k")’s
for the trajectory must be constrained to be monotonically decreasing in
![t](https://render.githubusercontent.com/render/math?math=t "t").

## Example

We start with a simple example where there are
![K = 2](https://render.githubusercontent.com/render/math?math=K%20%3D%202 "K = 2")
competing risks and

![\\begin{align\*} \\vec x\_{ij}(t) &= \\left(\\text{arcthan}\\left(\\frac{t - \\delta/2}{\\delta/2}\\right), 1, a\_{ij}, b\_{ij}\\right) \\\\ a\_{ij} &\\sim N(0, 1) \\\\
 b\_{ij} &\\sim \\text{Unif}(-1, 1)\\\\ \\vec z\_{ij} &= (1, a\_{ij}, b\_{ij}) \\end{align\*}](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign%2A%7D%20%5Cvec%20x_%7Bij%7D%28t%29%20%26%3D%20%5Cleft%28%5Ctext%7Barcthan%7D%5Cleft%28%5Cfrac%7Bt%20-%20%5Cdelta%2F2%7D%7B%5Cdelta%2F2%7D%5Cright%29%2C%201%2C%20a_%7Bij%7D%2C%20b_%7Bij%7D%5Cright%29%20%5C%5C%20a_%7Bij%7D%20%26%5Csim%20N%280%2C%201%29%20%5C%5C%0A%20b_%7Bij%7D%20%26%5Csim%20%5Ctext%7BUnif%7D%28-1%2C%201%29%5C%5C%20%5Cvec%20z_%7Bij%7D%20%26%3D%20%281%2C%20a_%7Bij%7D%2C%20b_%7Bij%7D%29%20%5Cend%7Balign%2A%7D "\begin{align*} \vec x_{ij}(t) &= \left(\text{arcthan}\left(\frac{t - \delta/2}{\delta/2}\right), 1, a_{ij}, b_{ij}\right) \\ a_{ij} &\sim N(0, 1) \\
 b_{ij} &\sim \text{Unif}(-1, 1)\\ \vec z_{ij} &= (1, a_{ij}, b_{ij}) \end{align*}")

We set the parameters below and plot the conditional cumulative
incidences function when the random effects are zero and the covariates
are zero,
![a\_{ij} = b\_{ij} = 0](https://render.githubusercontent.com/render/math?math=a_%7Bij%7D%20%3D%20b_%7Bij%7D%20%3D%200 "a_{ij} = b_{ij} = 0").

``` r
# assign model parameters
n_causes <- 2L
delta <- 2

# set the betas
coef_risk <- c(.67, 1, .1, -.4, .25, .3) |> 
  matrix(ncol = n_causes)

# set the gammas
coef_traject <- c(-.8, -.45, .8, .4, -1.2, .15, .25, -.2) |> 
  matrix(ncol = n_causes)

# plot the conditional cumulative incidences when random effects and covariates
# are all zero
local({
  probs <- exp(coef_risk[1, ]) / (1 + sum(exp(coef_risk[1, ])))
  par(mar = c(5, 5, 1, 1), mfcol = c(1, 2))
  
  for(i in 1:2){
    plot(\(x) probs[i] * pnorm(
      -coef_traject[1, i] * atanh((x - delta / 2) / (delta / 2)) - 
        coef_traject[2, i]),
         xlim = c(1e-8, delta), ylim = c(0, 1), bty = "l",  xlab = "Time", 
         ylab = sprintf("Cumulative incidence; cause %d", i),
       yaxs = "i", xaxs = "i")
    grid()
  }
})
```

<img src="man/figures/README-assign_model_parameters-1.png" width="100%" />

``` r
# set the covariance matrix
Sigma <- c(0.306, 0.008, -0.138, 0.197, 0.008, 0.759, 0.251, 
-0.25, -0.138, 0.251, 0.756, -0.319, 0.197, -0.25, -0.319, 0.903) |> 
  matrix(2L * n_causes)
```

Next, we assign a function to simulate clusters. The cluster sizes are
uniformly sampled from one to the maximum size. The censoring times are
drawn from a uniform distribution from zero to
![3\\delta](https://render.githubusercontent.com/render/math?math=3%5Cdelta "3\delta").

``` r
library(mvtnorm)

# simulates a data set with a given number of clusters and maximum number of 
# observations per cluster
sim_dat <- \(n_clusters, max_cluster_size){
  stopifnot(max_cluster_size > 0,
            n_clusters > 0)
  
  cluster_id <- 0L
  apply(rmvnorm(n_clusters, sigma = Sigma), 1, \(rng_effects){
    U <- head(rng_effects, n_causes)
    eta <- tail(rng_effects, n_causes)
    
    n_obs <- sample.int(max_cluster_size, 1L)
    cluster_id <<- cluster_id + 1L
    
    # draw the cause
    covs <- cbind(a = rnorm(n_obs), b = runif(n_obs, -1))
    Z <- cbind(1, covs)
  
    cond_logits_exp <- exp(Z %*% coef_risk + rep(U, each = n_obs)) |> 
      cbind(1)
    cond_probs <- cond_logits_exp / rowSums(cond_logits_exp)
    cause <- apply(cond_probs, 1, 
                   \(prob) sample.int(n_causes + 1L, 1L, prob = prob))
    
    # compute the observed time if needed
    obs_time <- mapply(\(cause, idx){
      if(cause > n_causes)
        return(delta)
      
      # can likely be done smarter but this is more general
      coefs <- coef_traject[, cause]
      offset <- sum(Z[idx, ] * coefs[-1]) + eta[cause]
      rng <- runif(1)
      eps <- .Machine$double.eps
      root <- uniroot(
        \(x) rng - pnorm(
          -coefs[1] * atanh((x - delta / 2) / (delta / 2)) - offset), 
        c(eps^2, delta * (1 - eps)), tol = 1e-12)$root
    }, cause, 1:n_obs)
    
    cens <- runif(n_obs, max = 3 * delta)
    has_finite_trajectory_prob <- cause <= n_causes
    is_censored <- which(!has_finite_trajectory_prob | cens < obs_time)
    
    if(length(is_censored) > 0){
      obs_time[is_censored] <- pmin(delta, cens[is_censored])
      cause[is_censored] <- n_causes + 1L
    }
    
    data.frame(covs, cause = cause, time = obs_time, cluster_id)
  }, simplify = FALSE) |> 
    do.call(what = rbind)
}
```

We then sample a data set.

``` r
# sample a data set
set.seed(8401828)
n_clusters <- 10000L
max_cluster_size <- 5L
dat <- sim_dat(n_clusters, max_cluster_size = max_cluster_size)

# show some stats
NROW(dat) # number of individuals
#> [1] 29968
table(dat$cause) # distribution of causes (3 is censored)
#> 
#>     1     2     3 
#> 13168  5238 11562

# distribution of observed times by cause
tapply(dat$time, dat$cause, quantile, 
       probs = seq(0, 1, length.out = 11), na.rm = TRUE)
#> $`1`
#>           0%          10%          20%          30%          40%          50% 
#> 5.257375e-07 5.174792e-03 2.564810e-02 7.863649e-02 1.967725e-01 4.413597e-01 
#>          60%          70%          80%          90%         100% 
#> 8.252684e-01 1.310347e+00 1.734937e+00 1.949810e+00 1.999998e+00 
#> 
#> $`2`
#>           0%          10%          20%          30%          40%          50% 
#> 0.0002289942 0.0577883402 0.1415916659 0.2631213266 0.4419615800 0.6554170608 
#>          60%          70%          80%          90%         100% 
#> 0.9220373015 1.2436285180 1.5466109693 1.8220248514 1.9997238579 
#> 
#> $`3`
#>           0%          10%          20%          30%          40%          50% 
#> 0.0002596225 0.2815441000 0.6283755401 0.9823111117 1.3845409666 1.8905316826 
#>          60%          70%          80%          90%         100% 
#> 2.0000000000 2.0000000000 2.0000000000 2.0000000000 2.0000000000
```

Then we setup the C++ object to do the computation.

``` r
library(mmcif)
comp_obj <- mmcif_data(
  ~ a + b, dat, cause = cause, time = time, cluster_id = cluster_id, 
  max_time = delta, spline_df = 4L)
```

The `mmcif_data` function does not work with

![h(t) = \\text{arcthan}\\left(\\frac{t - \\delta/2}{\\delta/2}\\right)](https://render.githubusercontent.com/render/math?math=h%28t%29%20%3D%20%5Ctext%7Barcthan%7D%5Cleft%28%5Cfrac%7Bt%20-%20%5Cdelta%2F2%7D%7B%5Cdelta%2F2%7D%5Cright%29 "h(t) = \text{arcthan}\left(\frac{t - \delta/2}{\delta/2}\right)")

but instead with
![\\vec g(h(t))](https://render.githubusercontent.com/render/math?math=%5Cvec%20g%28h%28t%29%29 "\vec g(h(t))")
where
![\\vec g](https://render.githubusercontent.com/render/math?math=%5Cvec%20g "\vec g")
returns a natural cubic spline basis functions. The knots are based on
quantiles of
![h(t)](https://render.githubusercontent.com/render/math?math=h%28t%29 "h(t)")
evaluated on the event times. The degrees of freedom of the spline is
controlled with the `spline_df` argument. There is a `constraints`
element on the object returned by the `mmcif_data` function which
contains matrices that ensures that the
![\\vec x\_{ij}(t)^\\top\\vec\\gamma\_k](https://render.githubusercontent.com/render/math?math=%5Cvec%20x_%7Bij%7D%28t%29%5E%5Ctop%5Cvec%5Cgamma_k "\vec x_{ij}(t)^\top\vec\gamma_k")s
are monotonically decreasing if
![C\\vec\\zeta &gt; \\vec 0](https://render.githubusercontent.com/render/math?math=C%5Cvec%5Czeta%20%3E%20%5Cvec%200 "C\vec\zeta > \vec 0")
where ![C](https://render.githubusercontent.com/render/math?math=C "C")
is one of matrices and
![\\vec\\zeta](https://render.githubusercontent.com/render/math?math=%5Cvec%5Czeta "\vec\zeta")
is the concatenated vector of model parameters.

The time to compute the log composite likelihood is illustrated below.

``` r
NCOL(comp_obj$pair_indices) # the number of pairs in the composite likelihood
#> [1] 39726
length(comp_obj$singletons) # the number of clusters with one observation
#> [1] 1977

# we need to find the combination of the spline bases that yield a straight 
# line. You can skip this
comb_slope <- local({
  boundary_knots <- comp_obj$spline$boundary_knots
  pts <- seq(boundary_knots[1], boundary_knots[2], length.out = 1000)
  lm.fit(cbind(1, comp_obj$spline$expansion(pts)), pts)$coef
})

# assign a function to compute the log composite likelihood
library(fastGHQuad)
#> Loading required package: Rcpp
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
#> [1] -69864.07

# check the time to compute the log composite likelihood
bench::mark(
  `one thread` = ll_func(n_threads = 1L, true_values),
  `two threads` = ll_func(n_threads = 2L, true_values),
  `three threads` = ll_func(n_threads = 3L, true_values),
  `four threads` = ll_func(n_threads = 4L, true_values), 
  min_time = 4)
#> # A tibble: 4 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 one thread       617ms    618ms      1.62    19.8KB        0
#> 2 two threads      314ms    316ms      3.17        0B        0
#> 3 three threads    213ms    214ms      4.67        0B        0
#> 4 four threads     163ms    169ms      5.92        0B        0

# next, we compute the gradient of the log composite likelihood at the true 
# parameters. First we assign a few functions to verify the result. You can 
# skip these
upper_to_full <- \(x){
  dim <- (sqrt(8 * length(x) + 1) - 1) / 2
  out <- matrix(0, dim, dim)
  out[upper.tri(out, TRUE)] <- x
  out[lower.tri(out)] <- t(out)[lower.tri(out)]
  out
}
d_upper_to_full <- \(x){
  dim <- (sqrt(8 * length(x) + 1) - 1) / 2
  out <- matrix(0, dim, dim)
  out[upper.tri(out, TRUE)] <- x
  out[upper.tri(out)] <- out[upper.tri(out)] / 2
  out[lower.tri(out)] <- t(out)[lower.tri(out)]
  out
}

# then we can compute the gradient with the function from the package and with 
# numerical differentiation
gr_func <- \(par, n_threads = 1L)
  mmcif:::mmcif_logLik_grad(
    comp_obj$comp_obj, par = par, ghq_data = ghq_data, n_threads = n_threads)

gr_package <- gr_func(true_values)

true_values_upper <- 
  c(coef_risk, coef_traject_spline, Sigma[upper.tri(Sigma, TRUE)])
gr_num <- numDeriv::grad(
  \(x) ll_func(c(head(x, -10), upper_to_full(tail(x, 10)))), 
  true_values_upper, method = "simple")

# they are very close but not exactly equal as expected (this is due to the 
# adaptive quadrature)
rbind(
  `Numerical gradient` = 
    c(head(gr_num, -10), d_upper_to_full(tail(gr_num, 10))), 
  `Gradient package` = gr_package)
#>                        [,1]     [,2]      [,3]      [,4]      [,5]     [,6]
#> Numerical gradient 147.6561 288.3880 -95.08746 -72.04594 -248.7902 88.15532
#> Gradient package   148.4757 288.9434 -94.88415 -72.00976 -248.3909 88.32859
#>                        [,7]      [,8]      [,9]     [,10]     [,11]     [,12]
#> Numerical gradient 202.2185 -35.72080 -243.5082 -52.87079 -54.54853 -97.04450
#> Gradient package   202.6837 -35.36141 -243.3177 -52.62701 -53.72854 -95.98335
#>                        [,13]     [,14]     [,15]     [,16]    [,17]     [,18]
#> Numerical gradient -290.3414 -227.9222 -28.27648 -69.21042 77.55546 -355.0063
#> Gradient package   -289.9349 -227.7802 -28.19683 -69.18012 77.58975 -354.7077
#>                        [,19]     [,20]     [,21]    [,22]     [,23]    [,24]
#> Numerical gradient -124.8397 -157.9143 -113.1931 77.88803 -21.16448 73.71993
#> Gradient package   -124.4192 -157.7813 -113.3295 78.14345 -21.04247 73.82650
#>                       [,25]     [,26]    [,27]     [,28]     [,29]    [,30]
#> Numerical gradient 77.88803 -56.21169 18.07621 -118.9189 -21.16448 18.07621
#> Gradient package   78.14345 -51.78413 17.38095 -118.7676 -21.04247 17.38095
#>                        [,31]    [,32]    [,33]     [,34]    [,35]    [,36]
#> Numerical gradient -286.3637 36.95665 73.71993 -118.9189 36.95665 120.8844
#> Gradient package   -286.0021 36.98897 73.82650 -118.7676 36.98897 120.9828

# check the time to compute the gradient of the log composite likelihood
bench::mark(
  `one thread` = gr_func(n_threads = 1L, true_values),
  `two threads` = gr_func(n_threads = 2L, true_values),
  `three threads` = gr_func(n_threads = 3L, true_values),
  `four threads` = gr_func(n_threads = 4L, true_values), 
  min_time = 4)
#> # A tibble: 4 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 one thread       1.54s    1.54s     0.650    7.09KB        0
#> 2 two threads   786.01ms 794.73ms     1.26       336B        0
#> 3 three threads 549.67ms 556.88ms     1.78       336B        0
#> 4 four threads  430.52ms 441.56ms     2.24       336B        0
```

Then we optimize the parameters (TODO: there will be a wrapper to work
with the log Cholesky decomposition and an estimation function).

``` r
# find the starting values
system.time(start <- mmcif_start_values(comp_obj, n_threads = 4L))
#>    user  system elapsed 
#>   0.387   0.000   0.105

# the maximum likelihood without the random effects. Note that this is not 
# comparable with the composite likelihood
attr(start, "logLik")
#> [1] -26193.57

# computes the log Cholesky decomposition
log_chol <- \(x){
  x <- chol(x)
  diag(x) <- diag(x) |> log()
  x[upper.tri(x, TRUE)]
}
log_chol_inv <- \(x){
  dim <- (sqrt(8 * length(x) + 1) - 1) / 2
  out <- matrix(0, dim, dim)
  out[upper.tri(out, TRUE)] <- x
  diag(out) <- diag(out) |> exp()
  crossprod(out)
}

# examples of using log_chol and log_chol_inv
log_chol(Sigma)
#>  [1] -0.59208509  0.01446203 -0.13801455 -0.24947003  0.29228783 -0.24851681
#>  [7]  0.35612750 -0.29291060 -0.18532138 -0.21077242
stopifnot(all.equal(Sigma, log_chol(Sigma) |> log_chol_inv()))

# computes the log composite likelihood with a log Cholesky decomposition for
# the covariance matrix
ll_func_chol <- \(par, n_threads = 1L, ghq = ghq_data){
  n_vcov <- (2L * n_causes * (2L * n_causes + 1L)) %/% 2L
  par <- c(head(par, -n_vcov), tail(par, n_vcov) |> log_chol_inv())
  
  mmcif:::mmcif_logLik(
    comp_obj$comp_obj, par = par, ghq_data = ghq, n_threads = n_threads)
}

# the gradient of the log composite likelihood with a log Cholesky 
# decomposition for the covariance matrix
ll_func_chol_grad <- \(par, n_threads = 1L, ghq = ghq_data){
  n_vcov <- (2L * n_causes * (2L * n_causes + 1L)) %/% 2L
  vcov <- tail(par, n_vcov) |> log_chol_inv()
  par <- c(head(par, -n_vcov), vcov)
  
  gr <- mmcif:::mmcif_logLik_grad(
    comp_obj$comp_obj, par = par, ghq_data = ghq, n_threads = n_threads)
  
  # back propagate the gradients w.r.t. the random effects
  d_vcov <- tail(gr, 4L * n_causes * n_causes) |> matrix(2L * n_causes)
  C <- chol(vcov)
  d_vcov <- 2 * C %*% d_vcov 
  diag(d_vcov) <- diag(d_vcov) * diag(C)
  
  c(head(gr, -4L * n_causes * n_causes), d_vcov[upper.tri(d_vcov, TRUE)])
}

# set true value
truth <- c(coef_risk, coef_traject_spline, log_chol(Sigma))

# we can verify that the gradient is correct
gr_package <- ll_func_chol_grad(truth)
gr_num <- numDeriv::grad(ll_func_chol, truth, method = "simple")

rbind(`Numerical gradient` = gr_num, `Gradient package` = gr_package)
#>                        [,1]     [,2]      [,3]      [,4]      [,5]     [,6]
#> Numerical gradient 147.6561 288.3880 -95.08746 -72.04594 -248.7902 88.15532
#> Gradient package   148.4757 288.9434 -94.88415 -72.00976 -248.3909 88.32859
#>                        [,7]      [,8]      [,9]     [,10]     [,11]     [,12]
#> Numerical gradient 202.2185 -35.72080 -243.5082 -52.87079 -54.54853 -97.04450
#> Gradient package   202.6837 -35.36141 -243.3177 -52.62701 -53.72854 -95.98335
#>                        [,13]     [,14]     [,15]     [,16]    [,17]     [,18]
#> Numerical gradient -290.3414 -227.9222 -28.27648 -69.21042 77.55546 -355.0063
#> Gradient package   -289.9349 -227.7802 -28.19683 -69.18012 77.58975 -354.7077
#>                        [,19]     [,20]     [,21]     [,22]     [,23]    [,24]
#> Numerical gradient -124.8397 -157.9143 -33.15133 -9.135805 -15.47944 146.1701
#> Gradient package   -124.4192 -157.7813 -33.21198 -8.309043  -9.12881 146.2659
#>                        [,25]     [,26]    [,27]     [,28]    [,29]    [,30]
#> Numerical gradient -157.5147 -359.2097 145.8281 -256.4344 12.82710 158.5795
#> Gradient package   -158.5782 -358.6617 145.9577 -256.1650 12.85821 158.7373

# optimize the log composite likelihood
constraints <- comp_obj$constraints$vcov_lower
system.time(
  fit <- constrOptim(
    start$lower, \(par) -ll_func_chol(par, n_threads = 4L, ghq_data), 
    grad = \(par) -ll_func_chol_grad(par, n_threads = 4L, ghq_data), 
    method = "BFGS", ui = constraints, ci = rep(1e-8, NROW(constraints)),
    control = list(maxit = 10000L)))
#>    user  system elapsed 
#> 178.499   0.004  44.967

# the log composite likelihood at different points
ll_func_chol(truth, 4L)
#> [1] -69864.07
ll_func_chol(start$lower, 4L)
#> [1] -75100.48
-fit$value
#> [1] -69809.5
```

We show the estimated and true the conditional cumulative incidence
functions (the dashed curves are the estimates) when the random effects
are zero and the covariates are zero,
![a\_{ij} = b\_{ij} = 0](https://render.githubusercontent.com/render/math?math=a_%7Bij%7D%20%3D%20b_%7Bij%7D%20%3D%200 "a_{ij} = b_{ij} = 0").

``` r
local({
  # get the estimates
  coef_risk_est <- fit$par[comp_obj$indices$coef_risk] |> 
    matrix(ncol = n_causes)
  coef_traject_time_est <- fit$par[comp_obj$indices$coef_trajectory_time] |> 
    matrix(ncol = n_causes)
  coef_traject_est <- fit$par[comp_obj$indices$coef_trajectory] |> 
    matrix(ncol = n_causes)
  coef_traject_intercept_est <- coef_traject_est[5, ]
  
  # compute the risk probabilities  
  probs <- exp(coef_risk[1, ]) / (1 + sum(exp(coef_risk[1, ])))
  probs_est <- exp(coef_risk_est[1, ]) / (1 + sum(exp(coef_risk_est[1, ])))
  
  # plot the estimated and true conditional cumulative incidence functions. The
  # estimates are the dashed lines
  par(mar = c(5, 5, 1, 1), mfcol = c(1, 2))
  pts <- seq(1e-8, delta, length.out = 1000)
  
  for(i in 1:2){
    true_vals <- probs[i] * pnorm(
      -coef_traject[1, i] * atanh((pts - delta / 2) / (delta / 2)) - 
        coef_traject[2, i])
    
    estimates <- probs_est[i] * pnorm(
      -comp_obj$time_expansion(pts) %*% coef_traject_time_est[, i] - 
        coef_traject_intercept_est[i]) |> drop()
    
    matplot(pts, cbind(true_vals, estimates), xlim = c(1e-8, delta), 
            ylim = c(0, 1), bty = "l",  xlab = "Time", lty = c(1, 2),
            ylab = sprintf("Cumulative incidence; cause %d", i),
            yaxs = "i", xaxs = "i", type = "l", col = "black")
    grid()
  }
})
```

<img src="man/figures/README-compare_estimated_incidence_funcs-1.png" width="100%" />

Further illustrations of the estimated model are given below.

``` r
# the number of call we made
fit$counts
#> function gradient 
#>      219       42
fit$outer.iterations
#> [1] 2

# compare the estimates with the true values
rbind(`Estimate AGHQ` = head(fit$par, length(coef_risk)),
      Truth = c(coef_risk))
#>                    [,1]     [,2]       [,3]       [,4]      [,5]      [,6]
#> Estimate AGHQ 0.6897047 1.013095 0.08454711 -0.4051271 0.2324695 0.3197769
#> Truth         0.6700000 1.000000 0.10000000 -0.4000000 0.2500000 0.3000000
rbind(`Estimate AGHQ` = fit$par[comp_obj$indices$coef_trajectory],
      Truth = truth[comp_obj$indices$coef_trajectory])
#>                    [,1]      [,2]      [,3]      [,4]     [,5]     [,6]
#> Estimate AGHQ -2.720973 -3.359790 -6.254128 -4.666525 2.714610 0.799606
#> Truth         -2.716513 -3.335158 -6.159469 -4.649314 2.687941 0.800000
#>                    [,7]      [,8]      [,9]     [,10]     [,11]    [,12]
#> Estimate AGHQ 0.3647859 -4.068753 -4.983885 -9.147150 -6.847682 4.807038
#> Truth         0.4000000 -4.074770 -5.002737 -9.239203 -6.973971 4.856911
#>                   [,13]      [,14]
#> Estimate AGHQ 0.2179501 -0.2445392
#> Truth         0.2500000 -0.2000000

n_vcov <- (2L * n_causes * (2L * n_causes + 1L)) %/% 2L
Sigma
#>        [,1]   [,2]   [,3]   [,4]
#> [1,]  0.306  0.008 -0.138  0.197
#> [2,]  0.008  0.759  0.251 -0.250
#> [3,] -0.138  0.251  0.756 -0.319
#> [4,]  0.197 -0.250 -0.319  0.903
log_chol_inv(tail(fit$par, n_vcov))
#>             [,1]        [,2]       [,3]       [,4]
#> [1,]  0.27035089  0.07523342 -0.1343073  0.1920847
#> [2,]  0.07523342  0.84796892  0.2603790 -0.2853482
#> [3,] -0.13430731  0.26037902  0.7210608 -0.2790476
#> [4,]  0.19208466 -0.28534816 -0.2790476  0.9088556
```

## TODOs

The package is still under development. Here are a few TODOs:

-   Implement a function to compute the sandwich estimator.
-   Implement a function to do the estimation.
-   Support delayed entry.

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Cederkvist18" class="csl-entry">

Cederkvist, Luise, Klaus K Holst, Klaus K Andersen, and Thomas H
Scheike. 2018. “<span class="nocase">Modeling the cumulative incidence
function of multivariate competing risks data allowing for
within-cluster dependence of risk and timing</span>.” *Biostatistics* 20
(2): 199–217. <https://doi.org/10.1093/biostatistics/kxx072>.

</div>

</div>
