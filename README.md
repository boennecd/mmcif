
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
    
    data.frame(covs, cause, time = obs_time, cluster_id)
  }, simplify = FALSE) |> 
    do.call(what = rbind)
}
```

We then sample a data set.

``` r
# sample a data set
set.seed(8401828)
n_clusters <- 1000L
max_cluster_size <- 5L
dat <- sim_dat(n_clusters, max_cluster_size = max_cluster_size)

# show some stats
NROW(dat) # number of individuals
#> [1] 2962
table(dat$cause) # distribution of causes (3 is censored)
#> 
#>    1    2    3 
#> 1249  542 1171

# distribution of observed times by cause
tapply(dat$time, dat$cause, quantile, 
       probs = seq(0, 1, length.out = 11), na.rm = TRUE)
#> $`1`
#>           0%          10%          20%          30%          40%          50% 
#> 1.614943e-05 4.917920e-03 2.736760e-02 9.050192e-02 2.219183e-01 4.790993e-01 
#>          60%          70%          80%          90%         100% 
#> 8.505799e-01 1.358068e+00 1.743518e+00 1.953161e+00 1.999980e+00 
#> 
#> $`2`
#>          0%         10%         20%         30%         40%         50% 
#> 0.001447541 0.050092145 0.157118679 0.276071737 0.431094056 0.669010350 
#>         60%         70%         80%         90%        100% 
#> 0.964643472 1.336520397 1.607221087 1.863062755 1.993280013 
#> 
#> $`3`
#>          0%         10%         20%         30%         40%         50% 
#> 0.002122959 0.246899164 0.577581364 1.007698883 1.526068409 2.000000000 
#>         60%         70%         80%         90%        100% 
#> 2.000000000 2.000000000 2.000000000 2.000000000 2.000000000
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
evaluated on the event times. The knots differ for each type of
competing risk. The degrees of freedom of the splines is controlled with
the `spline_df` argument. There is a `constraints` element on the object
returned by the `mmcif_data` function which contains matrices that
ensures that the
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
#> [1] 3911
length(comp_obj$singletons) # the number of clusters with one observation
#> [1] 202

# we need to find the combination of the spline bases that yield a straight 
# line. You can skip this
comb_slope <- sapply(comp_obj$spline, \(spline){
  boundary_knots <- spline$boundary_knots
  pts <- seq(boundary_knots[1], boundary_knots[2], length.out = 1000)
  lm.fit(cbind(1, spline$expansion(pts)), pts)$coef
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
  rbind(comb_slope[-1, ] * rep(coef_traject[1, ], each = NROW(comb_slope) - 1), 
        coef_traject[2, ] + comb_slope[1, ] * coef_traject[1, ],
        coef_traject[-(1:2), ])
true_values <- c(coef_risk, coef_traject_spline, Sigma)
ll_func(true_values)
#> [1] -7087.145

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
#> 1 one thread      58.8ms   59.8ms      16.6    6.77KB        0
#> 2 two threads     30.8ms   31.3ms      31.8        0B        0
#> 3 three threads   20.4ms   20.6ms      48.3        0B        0
#> 4 four threads    15.9ms   16.3ms      61.2        0B        0

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
#>                         [,1]      [,2]      [,3]     [,4]      [,5]     [,6]
#> Numerical gradient -98.11703 -35.07656 -34.63039 48.14763 -5.095038 65.07167
#> Gradient package   -98.03479 -35.01920 -34.61153 48.14879 -5.050224 65.09122
#>                        [,7]      [,8]      [,9]     [,10]     [,11]     [,12]
#> Numerical gradient 54.21873 -43.88601 -27.00679 -25.23115 -36.50268 -69.74052
#> Gradient package   54.25802 -43.85578 -26.99131 -25.21205 -36.42770 -69.62866
#>                        [,13]     [,14]    [,15]    [,16]     [,17]     [,18]
#> Numerical gradient -60.25759 -66.81265 42.02254 14.85036 -7.649851 -2.324416
#> Gradient package   -60.21927 -66.79572 42.03519 14.85769 -7.640918 -2.294211
#>                        [,19]     [,20]    [,21]     [,22]    [,23]     [,24]
#> Numerical gradient -5.068106 -43.14671 10.28011 -1.492192 4.406297 0.2704534
#> Gradient package   -5.023948 -43.13251 10.26578 -1.464293 4.420291 0.2806485
#>                        [,25]    [,26]      [,27]     [,28]    [,29]      [,30]
#> Numerical gradient -1.492192 10.41117 -0.2873668 -4.417070 4.406297 -0.2873668
#> Gradient package   -1.464293 10.86468 -0.3569913 -4.402148 4.420291 -0.3569913
#>                        [,31]    [,32]     [,33]     [,34]    [,35]    [,36]
#> Numerical gradient -36.46242 6.307481 0.2704534 -4.417070 6.307481 8.955149
#> Gradient package   -36.42063 6.310710 0.2806485 -4.402148 6.310710 8.967326

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
#> 1 one thread     171.2ms  172.2ms      5.80    7.09KB        0
#> 2 two threads     92.5ms   94.2ms     10.5       336B        0
#> 3 three threads     62ms   63.5ms     15.4       336B        0
#> 4 four threads    48.9ms   52.5ms     18.0       336B        0
```

Then we optimize the parameters (TODO: there will be a wrapper to work
with the log Cholesky decomposition and an estimation function).

``` r
# find the starting values
system.time(start <- mmcif_start_values(comp_obj, n_threads = 4L))
#>    user  system elapsed 
#>   0.059   0.012   0.023

# the maximum likelihood without the random effects. Note that this is not 
# comparable with the composite likelihood
attr(start, "logLik")
#> [1] -2649.699

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

# set true values
truth <- c(coef_risk, coef_traject_spline, log_chol(Sigma))

# we can verify that the gradient is correct
gr_package <- ll_func_chol_grad(truth)
gr_num <- numDeriv::grad(ll_func_chol, truth, method = "simple")

rbind(`Numerical gradient` = gr_num, `Gradient package` = gr_package)
#>                         [,1]      [,2]      [,3]     [,4]      [,5]     [,6]
#> Numerical gradient -98.11703 -35.07656 -34.63039 48.14763 -5.095038 65.07167
#> Gradient package   -98.03479 -35.01920 -34.61153 48.14879 -5.050224 65.09122
#>                        [,7]      [,8]      [,9]     [,10]     [,11]     [,12]
#> Numerical gradient 54.21873 -43.88601 -27.00679 -25.23115 -36.50268 -69.74052
#> Gradient package   54.25802 -43.85578 -26.99131 -25.21205 -36.42770 -69.62866
#>                        [,13]     [,14]    [,15]    [,16]     [,17]     [,18]
#> Numerical gradient -60.25759 -66.81265 42.02254 14.85036 -7.649851 -2.324416
#> Gradient package   -60.21927 -66.79572 42.03519 14.85769 -7.640918 -2.294211
#>                        [,19]     [,20]    [,21]     [,22]    [,23]    [,24]
#> Numerical gradient -5.068106 -43.14671 5.158584 -4.345941 17.90292 27.53855
#> Gradient package   -5.023948 -43.13251 5.149804 -4.263097 18.55267 27.54659
#>                       [,25]     [,26]    [,27]     [,28]    [,29]    [,30]
#> Numerical gradient -25.5072 -46.19810 3.407793 -9.258671 6.517770 11.74634
#> Gradient package   -25.6095 -46.13604 3.421523 -9.233461 6.520487 11.76572

# optimize the log composite likelihood
constraints <- comp_obj$constraints$vcov_upper
system.time(
  fit <- constrOptim(
    start$lower, \(par) -ll_func_chol(par, n_threads = 4L, ghq_data), 
    grad = \(par) -ll_func_chol_grad(par, n_threads = 4L, ghq_data), 
    method = "BFGS", ui = constraints, ci = rep(1e-8, NROW(constraints)),
    control = list(maxit = 10000L)))
#>    user  system elapsed 
#>  43.429   0.000  10.980

# the log composite likelihood at different points
ll_func_chol(truth, 4L)
#> [1] -7087.145
ll_func_chol(start$lower, 4L)
#> [1] -7572.462
-fit$value
#> [1] -7050.352
```

Then we compute the sandwich estimator. The Hessian is currently
computed with numerical differentiation which is why it takes a while.

``` r
system.time(
  sandwich_est <- mmcif_sandwich(
    comp_obj, fit$par, n_threads = 4L, ghq_data = ghq_data))
#>    user  system elapsed 
#>  95.326   0.004  24.360
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
      -comp_obj$time_expansion(pts, cause = i) %*% coef_traject_time_est[, i] - 
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
#>      330      104
fit$outer.iterations
#> [1] 2

# compute the standard errors from the sandwich estimator
SEs <- diag(sandwich_est) |> sqrt()

# compare the estimates with the true values
rbind(`Estimate AGHQ` = fit$par[comp_obj$indices$coef_risk],
      `Standard errors` = SEs[comp_obj$indices$coef_risk],
      Truth = c(coef_risk))
#>                       [,1]       [,2]       [,3]        [,4]       [,5]
#> Estimate AGHQ   0.58704204 0.94945204 0.08989478 -0.40847151 0.20877780
#> Standard errors 0.07241647 0.06901773 0.10193329  0.09896056 0.07073334
#> Truth           0.67000000 1.00000000 0.10000000 -0.40000000 0.25000000
#>                      [,6]
#> Estimate AGHQ   0.5050836
#> Standard errors 0.1233503
#> Truth           0.3000000
rbind(`Estimate AGHQ` = fit$par[comp_obj$indices$coef_trajectory],
      `Standard errors` = SEs[comp_obj$indices$coef_trajectory],
      Truth = truth[comp_obj$indices$coef_trajectory])
#>                       [,1]       [,2]       [,3]       [,4]      [,5]
#> Estimate AGHQ   -2.7611877 -3.6360790 -6.5359546 -4.9824233 2.7889044
#> Standard errors  0.1115275  0.1320104  0.2121534  0.1572529 0.1048162
#> Truth           -2.8546006 -3.5847914 -6.5119295 -4.9573949 2.8655410
#>                       [,6]       [,7]       [,8]       [,9]      [,10]
#> Estimate AGHQ   0.78976445 0.32069232 -2.8899416 -3.3092790 -6.2134010
#> Standard errors 0.05135401 0.06088468  0.2250919  0.2291263  0.4477583
#> Truth           0.80000000 0.40000000 -2.5969205 -3.3415994 -6.0232223
#>                      [,11]     [,12]      [,13]     [,14]
#> Estimate AGHQ   -4.8802007 3.3642055 0.24694410 -0.347764
#> Standard errors  0.3178818 0.2575071 0.06955232  0.105923
#> Truth           -4.6611059 3.1144598 0.25000000 -0.200000

n_vcov <- (2L * n_causes * (2L * n_causes + 1L)) %/% 2L
Sigma
#>        [,1]   [,2]   [,3]   [,4]
#> [1,]  0.306  0.008 -0.138  0.197
#> [2,]  0.008  0.759  0.251 -0.250
#> [3,] -0.138  0.251  0.756 -0.319
#> [4,]  0.197 -0.250 -0.319  0.903
log_chol_inv(tail(fit$par, n_vcov))
#>             [,1]        [,2]        [,3]       [,4]
#> [1,]  0.38546325  0.05187605 -0.05642973  0.1599336
#> [2,]  0.05187605  0.81940111  0.24241750 -0.4245997
#> [3,] -0.05642973  0.24241750  0.68717057 -0.2424773
#> [4,]  0.15993363 -0.42459971 -0.24247726  1.0622478

# on the log Cholesky scale
rbind(`Estimate AGHQ` = fit$par[comp_obj$indices$vcov_upper],
      `Standard errors` = SEs[comp_obj$indices$vcov_upper],
      Truth = truth[comp_obj$indices$vcov_upper])
#>                       [,1]       [,2]       [,3]        [,4]      [,5]
#> Estimate AGHQ   -0.4766547 0.08355557 -0.1038692 -0.09089009 0.2773770
#> Standard errors  0.2077299 0.23538315  0.1575934  0.14488670 0.1153071
#> Truth           -0.5920851 0.01446203 -0.1380145 -0.24947003 0.2922878
#>                       [,6]      [,7]       [,8]       [,9]      [,10]
#> Estimate AGHQ   -0.2537725 0.2576015 -0.4949544 -0.1053979 -0.1506872
#> Standard errors  0.1064003 0.2399130  0.2008492  0.1502317  0.1871904
#> Truth           -0.2485168 0.3561275 -0.2929106 -0.1853214 -0.2107724
```

### Delayed Entry

We extend the previous example to the setting where there may be delayed
entry (left truncation). Thus, we assign a new simulation function. The
delayed entry is sampled by sampling a random variable from the uniform
distribution on -1 to 1 and taking the entry time as being the maximum
of this variable and zero.

``` r
library(mvtnorm)

# simulates a data set with a given number of clusters and maximum number of 
# observations per cluster
sim_dat <- \(n_clusters, max_cluster_size){
  stopifnot(max_cluster_size > 0,
            n_clusters > 0)
  
  cluster_id <- 0L
  replicate(n_clusters, simplify = FALSE, {
    n_obs <- sample.int(max_cluster_size, 1L)
    cluster_id <<- cluster_id + 1L
    
    # draw the covariates and the left truncation time
    covs <- cbind(a = rnorm(n_obs), b = runif(n_obs, -1))
    Z <- cbind(1, covs)
    
    delayed_entry <- pmax(runif(n_obs, -1), 0)
    cens <- rep(-Inf, n_obs)
    while(all(cens <= delayed_entry))
      cens <- runif(n_obs, max = 3 * delta)
    
    successful_sample <- FALSE
    while(!successful_sample){
      rng_effects <- rmvnorm(1, sigma = Sigma) |> drop()
      U <- head(rng_effects, n_causes)
      eta <- tail(rng_effects, n_causes)
      
      # draw the cause
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
      
      keep <- which(pmin(obs_time, cens) > delayed_entry)
      successful_sample <- length(keep) > 0
      if(!successful_sample)
        next
      
      has_finite_trajectory_prob <- cause <= n_causes
      is_censored <- which(!has_finite_trajectory_prob | cens < obs_time)
      
      if(length(is_censored) > 0){
        obs_time[is_censored] <- pmin(delta, cens[is_censored])
        cause[is_censored] <- n_causes + 1L
      }
    }
    
    data.frame(covs, cause, time = obs_time, cluster_id, delayed_entry)[keep, ]
  }) |> 
    do.call(what = rbind)
}
```

We sample a data set using the new simulation function.

``` r
# sample a data set
set.seed(51312406)
n_clusters <- 1000L
max_cluster_size <- 5L
dat <- sim_dat(n_clusters, max_cluster_size = max_cluster_size)

# show some stats
NROW(dat) # number of individuals
#> [1] 2524
table(dat$cause) # distribution of causes (3 is censored)
#> 
#>    1    2    3 
#>  976  435 1113

# distribution of observed times by cause
tapply(dat$time, dat$cause, quantile, 
       probs = seq(0, 1, length.out = 11), na.rm = TRUE)
#> $`1`
#>           0%          10%          20%          30%          40%          50% 
#> 1.388927e-06 1.155209e-02 6.301834e-02 2.279148e-01 5.695660e-01 9.795840e-01 
#>          60%          70%          80%          90%         100% 
#> 1.311657e+00 1.649644e+00 1.886742e+00 1.980891e+00 1.999995e+00 
#> 
#> $`2`
#>          0%         10%         20%         30%         40%         50% 
#> 0.002018665 0.090196750 0.280351482 0.429229049 0.658359909 0.906824403 
#>         60%         70%         80%         90%        100% 
#> 1.180596945 1.409365896 1.674830085 1.877513206 1.996200103 
#> 
#> $`3`
#>          0%         10%         20%         30%         40%         50% 
#> 0.005215882 0.462200984 0.836500967 1.188546420 1.599797465 2.000000000 
#>         60%         70%         80%         90%        100% 
#> 2.000000000 2.000000000 2.000000000 2.000000000 2.000000000

# distribution of the left truncation time
quantile(dat$delayed_entry, probs = seq(0, 1, length.out = 11))
#>           0%          10%          20%          30%          40%          50% 
#> 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 
#>          60%          70%          80%          90%         100% 
#> 7.270984e-05 2.151228e-01 4.462644e-01 7.109908e-01 9.989804e-01
```

Next, we fit the model as before but this time we pass the delayed entry
time.

``` r
library(mmcif)
comp_obj <- mmcif_data(
  ~ a + b, dat, cause = cause, time = time, cluster_id = cluster_id, 
  max_time = delta, spline_df = 4L, left_trunc = delayed_entry)

# we need to find the combination of the spline bases that yield a straight 
# line. You can skip this
comb_slope <- sapply(comp_obj$spline, \(spline){
  boundary_knots <- spline$boundary_knots
  pts <- seq(boundary_knots[1], boundary_knots[2], length.out = 1000)
  lm.fit(cbind(1, spline$expansion(pts)), pts)$coef
})

coef_traject_spline <- 
  rbind(comb_slope[-1, ] * rep(coef_traject[1, ], each = NROW(comb_slope) - 1), 
        coef_traject[2, ] + comb_slope[1, ] * coef_traject[1, ],
        coef_traject[-(1:2), ])
        
# set true values
truth <- c(coef_risk, coef_traject_spline, log_chol(Sigma))

# find the starting values
system.time(start <- mmcif_start_values(comp_obj, n_threads = 4L))
#>    user  system elapsed 
#>   0.055   0.000   0.018

# we can verify that the gradient is correct again
gr_package <- ll_func_chol_grad(truth)
gr_num <- numDeriv::grad(ll_func_chol, truth, method = "simple")

rbind(`Numerical gradient` = gr_num, `Gradient package` = gr_package)
#>                         [,1]      [,2]     [,3]     [,4]     [,5]     [,6]
#> Numerical gradient -47.70747 -8.790845 6.978009 7.569686 7.152367 6.220029
#> Gradient package   -47.65039 -8.753452 6.990982 7.570827 7.179427 6.232856
#>                        [,7]     [,8]      [,9]    [,10]     [,11]    [,12]
#> Numerical gradient 5.933652 8.550291 -28.04840 18.37015 -47.03097 86.43754
#> Gradient package   5.956640 8.573283 -28.04046 18.38455 -46.98865 86.49652
#>                       [,13]     [,14]     [,15]    [,16]    [,17]     [,18]
#> Numerical gradient 2.075114 -45.32322 -17.03091 13.92585 17.29465 -20.56575
#> Gradient package   2.097569 -45.31066 -17.02144 13.92966 17.30091 -20.54930
#>                       [,19]     [,20]    [,21]     [,22]     [,23]     [,24]
#> Numerical gradient 20.15351 -1.486694 6.759593 -5.759377 -2.593290 -14.53035
#> Gradient package   20.17827 -1.478664 6.752974 -5.687308 -2.036417 -14.52778
#>                       [,25]     [,26]    [,27]    [,28]     [,29]    [,30]
#> Numerical gradient 20.44250 -9.738749 5.922140 -10.9871 -14.59199 4.311828
#> Gradient package   20.36639 -9.701316 5.931385 -10.9683 -14.59037 4.324320

# optimize the log composite likelihood
constraints <- comp_obj$constraints$vcov_upper
system.time(
  fit <- constrOptim(
    start$lower, \(par) -ll_func_chol(par, n_threads = 4L, ghq_data), 
    grad = \(par) -ll_func_chol_grad(par, n_threads = 4L, ghq_data), 
    method = "BFGS", ui = constraints, ci = rep(1e-8, NROW(constraints)),
    control = list(maxit = 10000L)))
#>    user  system elapsed 
#>  51.646   0.004  13.388

# the log composite likelihood at different points
ll_func_chol(truth, 4L)
#> [1] -4744.71
ll_func_chol(start$lower, 4L)
#> [1] -5077.429
-fit$value
#> [1] -4723.988
```

Then we compute the sandwich estimator. The Hessian is currently
computed with numerical differentiation which is why it takes a while.

``` r
system.time(
  sandwich_est <- mmcif_sandwich(
    comp_obj, fit$par, n_threads = 4L, ghq_data = ghq_data))
#>    user  system elapsed 
#> 234.477   0.008  62.172
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
      -comp_obj$time_expansion(pts, cause = i) %*% coef_traject_time_est[, i] - 
        coef_traject_intercept_est[i]) |> drop()
    
    matplot(pts, cbind(true_vals, estimates), xlim = c(1e-8, delta), 
            ylim = c(0, 1), bty = "l",  xlab = "Time", lty = c(1, 2),
            ylab = sprintf("Cumulative incidence; cause %d", i),
            yaxs = "i", xaxs = "i", type = "l", col = "black")
    grid()
  }
})
```

<img src="man/figures/README-delayed_compare_estimated_incidence_funcs-1.png" width="100%" />

Further illustrations of the estimated model are given below.

``` r
# the number of call we made
fit$counts
#> function gradient 
#>      157       52
fit$outer.iterations
#> [1] 2

# compute the standard errors from the sandwich estimator
SEs <- diag(sandwich_est) |> sqrt()

# compare the estimates with the true values
rbind(`Estimate AGHQ` = fit$par[comp_obj$indices$coef_risk],
      `Standard errors` = SEs[comp_obj$indices$coef_risk],
      Truth = c(coef_risk))
#>                       [,1]       [,2]      [,3]       [,4]       [,5]      [,6]
#> Estimate AGHQ   0.57771870 0.98274466 0.1390372 -0.4127569 0.23022871 0.3438114
#> Standard errors 0.07591665 0.08423562 0.1052981  0.1034130 0.07867396 0.1175204
#> Truth           0.67000000 1.00000000 0.1000000 -0.4000000 0.25000000 0.3000000
rbind(`Estimate AGHQ` = fit$par[comp_obj$indices$coef_trajectory],
      `Standard errors` = SEs[comp_obj$indices$coef_trajectory],
      Truth = truth[comp_obj$indices$coef_trajectory])
#>                      [,1]       [,2]       [,3]       [,4]     [,5]       [,6]
#> Estimate AGHQ   -2.982989 -3.6259793 -6.6762595 -4.7860109 2.596630 0.88394488
#> Standard errors  0.164094  0.1649437  0.3373203  0.2230601 0.150313 0.06574697
#> Truth           -3.051344 -3.6657486 -6.6720198 -4.8559792 2.677752 0.80000000
#>                       [,7]       [,8]       [,9]      [,10]      [,11]
#> Estimate AGHQ   0.40170351 -2.6792465 -3.1392297 -5.6456604 -4.1521398
#> Standard errors 0.07498459  0.2104229  0.1883607  0.4000962  0.2553798
#> Truth           0.40000000 -2.7771001 -3.3481062 -6.2334375 -4.6449539
#>                     [,12]      [,13]      [,14]
#> Estimate AGHQ   2.7000387 0.24516624 -0.1692535
#> Standard errors 0.2252096 0.06476622  0.1199285
#> Truth           3.0259114 0.25000000 -0.2000000

n_vcov <- (2L * n_causes * (2L * n_causes + 1L)) %/% 2L
Sigma
#>        [,1]   [,2]   [,3]   [,4]
#> [1,]  0.306  0.008 -0.138  0.197
#> [2,]  0.008  0.759  0.251 -0.250
#> [3,] -0.138  0.251  0.756 -0.319
#> [4,]  0.197 -0.250 -0.319  0.903
log_chol_inv(tail(fit$par, n_vcov))
#>             [,1]       [,2]        [,3]       [,4]
#> [1,]  0.33704670 -0.1472958 -0.07178614  0.1340451
#> [2,] -0.14729577  0.3656632  0.41934150 -0.1207784
#> [3,] -0.07178614  0.4193415  0.72101419 -0.4180223
#> [4,]  0.13404511 -0.1207784 -0.41802234  0.5935066

# on the log Cholesky scale
rbind(`Estimate AGHQ` = fit$par[comp_obj$indices$vcov_upper],
      `Standard errors` = SEs[comp_obj$indices$vcov_upper],
      Truth = truth[comp_obj$indices$vcov_upper])
#>                       [,1]        [,2]       [,3]       [,4]      [,5]
#> Estimate AGHQ   -0.5437669 -0.25371446 -0.5998374 -0.1236504 0.7068118
#> Standard errors  0.2758097  0.23393357  0.3605871  0.1949168 0.1901816
#> Truth           -0.5920851  0.01446203 -0.1380145 -0.2494700 0.2922878
#>                       [,6]      [,7]       [,8]       [,9]      [,10]
#> Estimate AGHQ   -0.7895953 0.2308904 -0.1133140 -0.6814131 -1.3820540
#> Standard errors  0.5918266 0.2204139  0.3098877  0.1525837  0.8795477
#> Truth           -0.2485168 0.3561275 -0.2929106 -0.1853214 -0.2107724
```

## TODOs

The package is still under development. Here are a few TODOs:

-   Implement a function to do the estimation.

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
