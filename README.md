
# MMCIF: Mixed Multivariate Cumulative Incidence Functions

[![R-CMD-check](https://github.com/boennecd/mmcif/workflows/R-CMD-check/badge.svg)](https://github.com/boennecd/mmcif/actions)

This package provides an implementation of the model introduced by
Cederkvist et al. (2018) to model within-cluster dependence of both risk
and timing in competing risk. For interested readers, a vignette on
computational details can be found by calling
`vignette("mmcif-comp-details", "mmcif")`.

## Installation

The package can be installed from Github by calling

``` r
library(remotes)
install_github("boennecd/mmcif", build_vignettes = TRUE)
options(digits = 4)
```

The code benefits from being build with automatic vectorization so
having e.g.  `-O3` in the `CXX17FLAGS` flags in your Makevars file may
be useful.

## The Model

The conditional cumulative incidence functions for cause
![k](https://render.githubusercontent.com/render/math?math=k "k") of
individual ![j](https://render.githubusercontent.com/render/math?math=j
"j") in cluster
![i](https://render.githubusercontent.com/render/math?math=i "i") is

  
![\\begin{align\*} F\_{kij}(t\\mid \\vec u\_i, \\vec\\eta\_i) &=
\\pi\_k(\\vec z\_{ij}, \\vec u\_i) \\Phi(-\\vec
x\_{ij}(t)^\\top\\vec\\gamma\_k - \\eta\_{ik}) \\\\ \\pi\_k(\\vec
z\_{ij}, \\vec u\_i) &= \\frac{\\exp(\\vec z\_{ij}^\\top\\vec\\beta\_k +
u\_{ik})}{1 + \\sum\_{l = 1}^K\\exp(\\vec z\_{ij}^\\top\\vec\\beta\_l +
u\_{il})} \\\\ \\begin{pmatrix} \\vec U\_i \\\\ \\vec\\eta\_i
\\end{pmatrix} &\\sim
N^{(2K)}(\\vec 0;\\Sigma).\\end{align\*}](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign%2A%7D%20F_%7Bkij%7D%28t%5Cmid%20%5Cvec%20u_i%2C%20%5Cvec%5Ceta_i%29%20%26%3D%20%5Cpi_k%28%5Cvec%20z_%7Bij%7D%2C%20%5Cvec%20u_i%29%20%5CPhi%28-%5Cvec%20x_%7Bij%7D%28t%29%5E%5Ctop%5Cvec%5Cgamma_k%20-%20%5Ceta_%7Bik%7D%29%20%5C%5C%20%5Cpi_k%28%5Cvec%20z_%7Bij%7D%2C%20%5Cvec%20u_i%29%20%26%3D%20%5Cfrac%7B%5Cexp%28%5Cvec%20z_%7Bij%7D%5E%5Ctop%5Cvec%5Cbeta_k%20%2B%20u_%7Bik%7D%29%7D%7B1%20%2B%20%5Csum_%7Bl%20%3D%201%7D%5EK%5Cexp%28%5Cvec%20z_%7Bij%7D%5E%5Ctop%5Cvec%5Cbeta_l%20%2B%20u_%7Bil%7D%29%7D%20%5C%5C%20%5Cbegin%7Bpmatrix%7D%20%5Cvec%20U_i%20%5C%5C%20%5Cvec%5Ceta_i%20%5Cend%7Bpmatrix%7D%20%26%5Csim%20N%5E%7B%282K%29%7D%28%5Cvec%200%3B%5CSigma%29.%5Cend%7Balign%2A%7D
"\\begin{align*} F_{kij}(t\\mid \\vec u_i, \\vec\\eta_i) &= \\pi_k(\\vec z_{ij}, \\vec u_i) \\Phi(-\\vec x_{ij}(t)^\\top\\vec\\gamma_k - \\eta_{ik}) \\\\ \\pi_k(\\vec z_{ij}, \\vec u_i) &= \\frac{\\exp(\\vec z_{ij}^\\top\\vec\\beta_k + u_{ik})}{1 + \\sum_{l = 1}^K\\exp(\\vec z_{ij}^\\top\\vec\\beta_l + u_{il})} \\\\ \\begin{pmatrix} \\vec U_i \\\\ \\vec\\eta_i \\end{pmatrix} &\\sim N^{(2K)}(\\vec 0;\\Sigma).\\end{align*}")  

where there are
![K](https://render.githubusercontent.com/render/math?math=K "K")
competing risks. The ![\\vec
x\_{ij}(t)^\\top\\vec\\gamma\_k](https://render.githubusercontent.com/render/math?math=%5Cvec%20x_%7Bij%7D%28t%29%5E%5Ctop%5Cvec%5Cgamma_k
"\\vec x_{ij}(t)^\\top\\vec\\gamma_k")’s for the trajectory must be
constrained to be monotonically decreasing in
![t](https://render.githubusercontent.com/render/math?math=t "t"). The
covariates for the trajectory in this package are defined as

  
![\\vec x\_{ij}(t) = (\\vec h(t)^\\top, \\vec
z\_{ij}^\\top)^\\top](https://render.githubusercontent.com/render/math?math=%5Cvec%20x_%7Bij%7D%28t%29%20%3D%20%28%5Cvec%20h%28t%29%5E%5Ctop%2C%20%5Cvec%20z_%7Bij%7D%5E%5Ctop%29%5E%5Ctop
"\\vec x_{ij}(t) = (\\vec h(t)^\\top, \\vec z_{ij}^\\top)^\\top")  

for a spline basis ![\\vec
h](https://render.githubusercontent.com/render/math?math=%5Cvec%20h
"\\vec h") and known covariates ![\\vec
z\_{ij}](https://render.githubusercontent.com/render/math?math=%5Cvec%20z_%7Bij%7D
"\\vec z_{ij}") which are also used in the risk part of the model.

## Example

We start with a simple example where there are ![K
= 2](https://render.githubusercontent.com/render/math?math=K%20%3D%202
"K = 2") competing risks and

  
![\\begin{align\*} \\vec x\_{ij}(t) &=
\\left(\\text{arcthan}\\left(\\frac{t -
\\delta/2}{\\delta/2}\\right), 1, a\_{ij}, b\_{ij}\\right) \\\\ a\_{ij}
&\\sim N(0, 1) \\\\
b\_{ij} &\\sim \\text{Unif}(-1, 1)\\\\ \\vec z\_{ij} &= (1, a\_{ij},
b\_{ij})
\\end{align\*}](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign%2A%7D%20%5Cvec%20x_%7Bij%7D%28t%29%20%26%3D%20%5Cleft%28%5Ctext%7Barcthan%7D%5Cleft%28%5Cfrac%7Bt%20-%20%5Cdelta%2F2%7D%7B%5Cdelta%2F2%7D%5Cright%29%2C%201%2C%20a_%7Bij%7D%2C%20b_%7Bij%7D%5Cright%29%20%5C%5C%20a_%7Bij%7D%20%26%5Csim%20N%280%2C%201%29%20%5C%5C%0A%20b_%7Bij%7D%20%26%5Csim%20%5Ctext%7BUnif%7D%28-1%2C%201%29%5C%5C%20%5Cvec%20z_%7Bij%7D%20%26%3D%20%281%2C%20a_%7Bij%7D%2C%20b_%7Bij%7D%29%20%5Cend%7Balign%2A%7D
"\\begin{align*} \\vec x_{ij}(t) &= \\left(\\text{arcthan}\\left(\\frac{t - \\delta/2}{\\delta/2}\\right), 1, a_{ij}, b_{ij}\\right) \\\\ a_{ij} &\\sim N(0, 1) \\\\
 b_{ij} &\\sim \\text{Unif}(-1, 1)\\\\ \\vec z_{ij} &= (1, a_{ij}, b_{ij}) \\end{align*}")  

We set the parameters below and plot the conditional cumulative
incidences function when the random effects are zero and the covariates
are zero, ![a\_{ij} = b\_{ij}
= 0](https://render.githubusercontent.com/render/math?math=a_%7Bij%7D%20%3D%20b_%7Bij%7D%20%3D%200
"a_{ij} = b_{ij} = 0").

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
         xlim = c(0, delta), ylim = c(0, 1), bty = "l",  xlab = "Time", 
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
![3\\delta](https://render.githubusercontent.com/render/math?math=3%5Cdelta
"3\\delta").

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

  
![h(t) = \\text{arcthan}\\left(\\frac{t -
\\delta/2}{\\delta/2}\\right)](https://render.githubusercontent.com/render/math?math=h%28t%29%20%3D%20%5Ctext%7Barcthan%7D%5Cleft%28%5Cfrac%7Bt%20-%20%5Cdelta%2F2%7D%7B%5Cdelta%2F2%7D%5Cright%29
"h(t) = \\text{arcthan}\\left(\\frac{t - \\delta/2}{\\delta/2}\\right)")  

but instead with ![\\vec
g(h(t))](https://render.githubusercontent.com/render/math?math=%5Cvec%20g%28h%28t%29%29
"\\vec g(h(t))") where ![\\vec
g](https://render.githubusercontent.com/render/math?math=%5Cvec%20g
"\\vec g") returns a natural cubic spline basis functions. The knots are
based on quantiles of
![h(t)](https://render.githubusercontent.com/render/math?math=h%28t%29
"h(t)") evaluated on the event times. The knots differ for each type of
competing risk. The degrees of freedom of the splines is controlled with
the `spline_df` argument. There is a `constraints` element on the object
returned by the `mmcif_data` function which contains matrices that
ensures that the ![\\vec
x\_{ij}(t)^\\top\\vec\\gamma\_k](https://render.githubusercontent.com/render/math?math=%5Cvec%20x_%7Bij%7D%28t%29%5E%5Ctop%5Cvec%5Cgamma_k
"\\vec x_{ij}(t)^\\top\\vec\\gamma_k")s are monotonically decreasing if
![C\\vec\\zeta \>
\\vec 0](https://render.githubusercontent.com/render/math?math=C%5Cvec%5Czeta%20%3E%20%5Cvec%200
"C\\vec\\zeta \> \\vec 0") where
![C](https://render.githubusercontent.com/render/math?math=C "C") is one
of matrices and
![\\vec\\zeta](https://render.githubusercontent.com/render/math?math=%5Cvec%5Czeta
"\\vec\\zeta") is the concatenated vector of model parameters.

The time to compute the log composite likelihood is illustrated below.

``` r
NCOL(comp_obj$pair_indices) # the number of pairs in the composite likelihood
#> [1] 3911
length(comp_obj$singletons) # the number of clusters with one observation
#> [1] 202

# we need to find the combination of the spline bases that yield a straight 
# line to construct the true values using the splines. You can skip this
comb_slope <- sapply(comp_obj$spline, \(spline){
  boundary_knots <- spline$boundary_knots
  pts <- seq(boundary_knots[1], boundary_knots[2], length.out = 1000)
  lm.fit(cbind(1, spline$expansion(pts)), pts)$coef
})

# assign a function to compute the log composite likelihood
ll_func <- \(par, n_threads = 1L)
  mmcif_logLik(
    comp_obj, par = par, n_threads = n_threads, is_log_chol = FALSE)

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
#> 1 one thread      54.2ms   54.6ms      18.3        0B        0
#> 2 two threads       29ms   29.2ms      34.2        0B        0
#> 3 three threads   19.3ms   19.5ms      51.1        0B        0
#> 4 four threads    15.4ms   16.1ms      61.6        0B        0

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
gr_func <- function(par, n_threads = 1L)
  mmcif_logLik_grad(comp_obj, par, n_threads = n_threads, is_log_chol = FALSE)
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
#> 1 one thread     156.8ms  157.4ms      6.22      336B        0
#> 2 two threads       86ms   86.5ms     11.4       336B        0
#> 3 three threads   56.9ms   57.8ms     17.3       336B        0
#> 4 four threads      46ms     48ms     20.9       336B        0
```

Then we optimize the parameters.

``` r
# find the starting values
system.time(start <- mmcif_start_values(comp_obj, n_threads = 4L))
#>    user  system elapsed 
#>   0.056   0.000   0.018

# the maximum likelihood without the random effects. Note that this is not 
# comparable with the composite likelihood
attr(start, "logLik")
#> [1] -2649.699

# examples of using log_chol and log_chol_inv
log_chol(Sigma)
#>  [1] -0.59208509  0.01446203 -0.13801455 -0.24947003  0.29228783 -0.24851681
#>  [7]  0.35612750 -0.29291060 -0.18532138 -0.21077242
stopifnot(all.equal(Sigma, log_chol(Sigma) |> log_chol_inv()))

# set true values
truth <- c(coef_risk, coef_traject_spline, log_chol(Sigma))

# we can verify that the gradient is correct
gr_package <- mmcif_logLik_grad(
  comp_obj, truth, n_threads = 4L, is_log_chol = TRUE)
gr_num <- numDeriv::grad(
  mmcif_logLik, truth, object = comp_obj, n_threads = 4L, is_log_chol = TRUE, 
  method = "simple")

rbind(`Numerical gradient` = gr_num, `Gradient package` = gr_package)
#>                         [,1]      [,2]      [,3]     [,4]      [,5]     [,6]
#> Numerical gradient -98.11703 -35.07656 -34.63039 48.14763 -5.095038 65.07167
#> Gradient package   -98.03479 -35.01920 -34.61153 48.14879 -5.050224 65.09122
#>                        [,7]      [,8]      [,9]     [,10]     [,11]     [,12]
#> Numerical gradient 54.21873 -43.88601 -27.00679 -25.23115 -36.50268 -69.74052
#> Gradient package   54.25802 -43.85578 -26.99131 -25.21205 -36.42770 -69.62866
#>                        [,13]     [,14]    [,15]    [,16]     [,17]     [,18]
#> Numerical gradient -60.25759 -66.81265 42.02254 14.85036 -7.649850 -2.324416
#> Gradient package   -60.21927 -66.79572 42.03519 14.85769 -7.640918 -2.294211
#>                        [,19]     [,20]    [,21]     [,22]    [,23]    [,24]
#> Numerical gradient -5.068106 -43.14671 5.158584 -4.345941 17.90292 27.53855
#> Gradient package   -5.023948 -43.13251 5.149804 -4.263097 18.55267 27.54659
#>                       [,25]     [,26]    [,27]     [,28]    [,29]    [,30]
#> Numerical gradient -25.5072 -46.19810 3.407793 -9.258671 6.517770 11.74634
#> Gradient package   -25.6095 -46.13604 3.421523 -9.233461 6.520487 11.76572

# optimize the log composite likelihood
system.time(fit <- mmcif_fit(start$upper, comp_obj, n_threads = 4L))
#>    user  system elapsed 
#>  27.074   0.056   6.858

# the log composite likelihood at different points
mmcif_logLik(comp_obj, truth, n_threads = 4L, is_log_chol = TRUE)
#> [1] -7087.145
mmcif_logLik(comp_obj, start$upper, n_threads = 4L, is_log_chol = TRUE)
#> [1] -7572.462
-fit$value
#> [1] -7050.351
```

Then we compute the sandwich estimator. The Hessian is currently
computed with numerical differentiation which is why it takes a while.

``` r
system.time(sandwich_est <- mmcif_sandwich(comp_obj, fit$par, n_threads = 4L))
#>    user  system elapsed 
#>  88.983   0.004  22.606

# setting order equal to zero yield no Richardson extrapolation and just
# standard symmetric difference quotient. This is less precise but faster 
system.time(sandwich_est_simple <- 
              mmcif_sandwich(comp_obj, fit$par, n_threads = 4L, order = 0L))
#>    user  system elapsed 
#>  12.415   0.000   3.158
```

We show the estimated and true the conditional cumulative incidence
functions (the dashed curves are the estimates) when the random effects
are zero and the covariates are zero, ![a\_{ij} = b\_{ij}
= 0](https://render.githubusercontent.com/render/math?math=a_%7Bij%7D%20%3D%20b_%7Bij%7D%20%3D%200
"a_{ij} = b_{ij} = 0").

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
  pts <- seq(1e-8, delta * (1 - 1e-8), length.out = 1000)
  
  for(i in 1:2){
    true_vals <- probs[i] * pnorm(
      -coef_traject[1, i] * atanh((pts - delta / 2) / (delta / 2)) - 
        coef_traject[2, i])
    
    estimates <- probs_est[i] * pnorm(
      -comp_obj$time_expansion(pts, cause = i) %*% coef_traject_time_est[, i] - 
        coef_traject_intercept_est[i]) |> drop()
    
    matplot(pts, cbind(true_vals, estimates), xlim = c(0, delta), 
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
#>      110      104
fit$outer.iterations
#> [1] 2

# compute the standard errors from the sandwich estimator
SEs <- diag(sandwich_est) |> sqrt()
SEs_simple <- diag(sandwich_est_simple) |> sqrt()

# compare the estimates with the true values
rbind(`Estimate AGHQ` = fit$par[comp_obj$indices$coef_risk],
      `Standard errors` = SEs[comp_obj$indices$coef_risk],
      `Standard errors simple` = SEs_simple[comp_obj$indices$coef_risk],
      Truth = truth[comp_obj$indices$coef_risk])
#>                        cause1:risk:(Intercept) cause1:risk:a cause1:risk:b
#> Estimate AGHQ                       0.58620564    0.94954204    0.08950182
#> Standard errors                     0.07248608    0.06900108    0.10192677
#> Standard errors simple              0.07248630    0.06900042    0.10192677
#> Truth                               0.67000000    1.00000000    0.10000000
#>                        cause2:risk:(Intercept) cause2:risk:a cause2:risk:b
#> Estimate AGHQ                      -0.40865018    0.20842971     0.5046457
#> Standard errors                     0.09901495    0.07071248     0.1232734
#> Standard errors simple              0.09901363    0.07071228     0.1232732
#> Truth                              -0.40000000    0.25000000     0.3000000
rbind(`Estimate AGHQ` = fit$par[comp_obj$indices$coef_trajectory],
      `Standard errors` = SEs[comp_obj$indices$coef_trajectory],
      `Standard errors simple` = SEs_simple[comp_obj$indices$coef_trajectory],
      Truth = truth[comp_obj$indices$coef_trajectory])
#>                        cause1:spline1 cause1:spline2 cause1:spline3
#> Estimate AGHQ              -2.7627254     -3.6378274     -6.5395047
#> Standard errors             0.1116663      0.1321685      0.2124841
#> Standard errors simple      0.1116999      0.1322075      0.2124407
#> Truth                      -2.8546006     -3.5847914     -6.5119295
#>                        cause1:spline4 cause1:traject:(Intercept)
#> Estimate AGHQ              -4.9846541                  2.7904996
#> Standard errors             0.1574550                  0.1049764
#> Standard errors simple      0.1574231                  0.1049866
#> Truth                      -4.9573949                  2.8655410
#>                        cause1:traject:a cause1:traject:b cause2:spline1
#> Estimate AGHQ                0.79012484       0.32085836     -2.8832267
#> Standard errors              0.05139623       0.06090990      0.2243698
#> Standard errors simple       0.05139838       0.06091115      0.2242685
#> Truth                        0.80000000       0.40000000     -2.5969205
#>                        cause2:spline2 cause2:spline3 cause2:spline4
#> Estimate AGHQ              -3.3018386     -6.1991157     -4.8688943
#> Standard errors             0.2282964      0.4462072      0.3165857
#> Standard errors simple      0.2282202      0.4457886      0.3163094
#> Truth                      -3.3415994     -6.0232223     -4.6611059
#>                        cause2:traject:(Intercept) cause2:traject:a
#> Estimate AGHQ                           3.3577293       0.24607315
#> Standard errors                         0.2567064       0.06942385
#> Standard errors simple                  0.2565522       0.06942196
#> Truth                                   3.1144598       0.25000000
#>                        cause2:traject:b
#> Estimate AGHQ                -0.3471028
#> Standard errors               0.1057178
#> Standard errors simple        0.1057128
#> Truth                        -0.2000000

n_vcov <- (2L * n_causes * (2L * n_causes + 1L)) %/% 2L
Sigma
#>        [,1]   [,2]   [,3]   [,4]
#> [1,]  0.306  0.008 -0.138  0.197
#> [2,]  0.008  0.759  0.251 -0.250
#> [3,] -0.138  0.251  0.756 -0.319
#> [4,]  0.197 -0.250 -0.319  0.903
log_chol_inv(tail(fit$par, n_vcov))
#>             [,1]        [,2]        [,3]       [,4]
#> [1,]  0.38507267  0.04766207 -0.05680156  0.1582099
#> [2,]  0.04766207  0.81270528  0.24236467 -0.4260827
#> [3,] -0.05680156  0.24236467  0.68842695 -0.2425395
#> [4,]  0.15820990 -0.42608266 -0.24253949  1.0533736

# on the log Cholesky scale
rbind(`Estimate AGHQ` = fit$par[comp_obj$indices$vcov_upper],
      `Standard errors` = SEs[comp_obj$indices$vcov_upper],
      `Standard errors simple` = SEs_simple[comp_obj$indices$vcov_upper],
      Truth = truth[comp_obj$indices$vcov_upper])
#>                        vcov:risk1:risk1 vcov:risk1:risk2 vcov:risk2:risk2
#> Estimate AGHQ                -0.4771616       0.07680713       -0.1073361
#> Standard errors               0.2078095       0.23790573        0.1594426
#> Standard errors simple        0.2077903       0.23790612        0.1594369
#> Truth                        -0.5920851       0.01446203       -0.1380145
#>                        vcov:risk1:traject1 vcov:risk2:traject1
#> Estimate AGHQ                  -0.09153537           0.2776538
#> Standard errors                 0.14526369           0.1161468
#> Standard errors simple          0.14526139           0.1161449
#> Truth                          -0.24947003           0.2922878
#>                        vcov:traject1:traject1 vcov:risk1:traject2
#> Estimate AGHQ                      -0.2529550           0.2549543
#> Standard errors                     0.1064977           0.2403814
#> Standard errors simple              0.1064902           0.2403721
#> Truth                              -0.2485168           0.3561275
#>                        vcov:risk2:traject2 vcov:traject1:traject2
#> Estimate AGHQ                   -0.4961625             -0.1048815
#> Standard errors                  0.2021080              0.1505190
#> Standard errors simple           0.2021143              0.1505188
#> Truth                           -0.2929106             -0.1853214
#>                        vcov:traject2:traject2
#> Estimate AGHQ                      -0.1565378
#> Standard errors                     0.1889483
#> Standard errors simple              0.1888801
#> Truth                              -0.2107724

# on the original covariance matrix scale
vcov_est <- log_chol_inv(tail(fit$par, n_vcov))
vcov_est[lower.tri(vcov_est)] <- NA_real_
vcov_SE <- matrix(NA_real_, NROW(vcov_est), NCOL(vcov_est))
vcov_SE[upper.tri(vcov_SE, TRUE)] <- 
  attr(sandwich_est, "res vcov") |> diag() |> sqrt() |> 
  tail(n_vcov)

vcov_show <- cbind(Estimates = vcov_est, NA, SEs = vcov_SE) 
colnames(vcov_show) <- 
  c(rep("Est.", NCOL(vcov_est)), "", rep("SE", NCOL(vcov_est)))
print(vcov_show, na.print = "")
#>           Est.       Est.        Est.       Est.         SE        SE        SE
#> [1,] 0.3850727 0.04766207 -0.05680156  0.1582099  0.1600435 0.3039355 0.1788339
#> [2,]           0.81270528  0.24236467 -0.4260827            0.2729629 0.2206890
#> [3,]                       0.68842695 -0.2425395                      0.1159555
#> [4,]                                   1.0533736                               
#>             SE
#> [1,] 0.3015934
#> [2,] 0.3628394
#> [3,] 0.2063421
#> [4,] 0.2802714

Sigma # the true values
#>        [,1]   [,2]   [,3]   [,4]
#> [1,]  0.306  0.008 -0.138  0.197
#> [2,]  0.008  0.759  0.251 -0.250
#> [3,] -0.138  0.251  0.756 -0.319
#> [4,]  0.197 -0.250 -0.319  0.903
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
# line to construct the true values using the splines. You can skip this
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
#>   0.048   0.000   0.016

# we can verify that the gradient is correct again
gr_package <- mmcif_logLik_grad(
  comp_obj, truth, n_threads = 4L, is_log_chol = TRUE)
gr_num <- numDeriv::grad(
  mmcif_logLik, truth, object = comp_obj, n_threads = 4L, is_log_chol = TRUE, 
  method = "simple")

rbind(`Numerical gradient` = gr_num, `Gradient package` = gr_package)
#>                         [,1]      [,2]     [,3]     [,4]     [,5]     [,6]
#> Numerical gradient -47.70747 -8.790845 6.978009 7.569686 7.152368 6.220029
#> Gradient package   -47.65039 -8.753452 6.990982 7.570827 7.179427 6.232856
#>                        [,7]     [,8]      [,9]    [,10]     [,11]    [,12]
#> Numerical gradient 5.933652 8.550291 -28.04840 18.37015 -47.03097 86.43754
#> Gradient package   5.956640 8.573283 -28.04046 18.38455 -46.98865 86.49652
#>                       [,13]     [,14]     [,15]    [,16]    [,17]     [,18]
#> Numerical gradient 2.075114 -45.32322 -17.03091 13.92585 17.29465 -20.56575
#> Gradient package   2.097569 -45.31066 -17.02144 13.92966 17.30091 -20.54930
#>                       [,19]     [,20]    [,21]     [,22]     [,23]     [,24]
#> Numerical gradient 20.15351 -1.486694 6.759594 -5.759376 -2.593290 -14.53035
#> Gradient package   20.17827 -1.478664 6.752974 -5.687308 -2.036417 -14.52778
#>                       [,25]     [,26]    [,27]    [,28]     [,29]    [,30]
#> Numerical gradient 20.44250 -9.738749 5.922140 -10.9871 -14.59199 4.311828
#> Gradient package   20.36639 -9.701316 5.931385 -10.9683 -14.59037 4.324320

# optimize the log composite likelihood
system.time(fit <- mmcif_fit(start$upper, comp_obj, n_threads = 4L))
#>    user  system elapsed 
#>  78.293   0.000  20.330

# the log composite likelihood at different points
mmcif_logLik(comp_obj, truth, n_threads = 4L, is_log_chol = TRUE)
#> [1] -4744.71
mmcif_logLik(comp_obj, start$upper, n_threads = 4L, is_log_chol = TRUE)
#> [1] -5077.429
-fit$value
#> [1] -4723.969
```

Then we compute the sandwich estimator. The Hessian is currently
computed with numerical differentiation which is why it takes a while.

``` r
system.time(sandwich_est <- mmcif_sandwich(comp_obj, fit$par, n_threads = 4L))
#>    user  system elapsed 
#> 228.422   0.009  60.317

# setting order equal to zero yield no Richardson extrapolation and just
# standard symmetric difference quotient. This is less precise but faster 
system.time(sandwich_est_simple <- 
              mmcif_sandwich(comp_obj, fit$par, n_threads = 4L, order = 0L))
#>    user  system elapsed 
#>  31.326   0.000   8.309
```

We show the estimated and true the conditional cumulative incidence
functions (the dashed curves are the estimates) when the random effects
are zero and the covariates are zero, ![a\_{ij} = b\_{ij}
= 0](https://render.githubusercontent.com/render/math?math=a_%7Bij%7D%20%3D%20b_%7Bij%7D%20%3D%200
"a_{ij} = b_{ij} = 0").

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
  pts <- seq(1e-8, delta * (1 - 1e-8), length.out = 1000)
  
  for(i in 1:2){
    true_vals <- probs[i] * pnorm(
      -coef_traject[1, i] * atanh((pts - delta / 2) / (delta / 2)) - 
        coef_traject[2, i])
    
    estimates <- probs_est[i] * pnorm(
      -comp_obj$time_expansion(pts, cause = i) %*% coef_traject_time_est[, i] - 
        coef_traject_intercept_est[i]) |> drop()
    
    matplot(pts, cbind(true_vals, estimates), xlim = c(0, delta), 
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
#>      120      117
fit$outer.iterations
#> [1] 2

# compute the standard errors from the sandwich estimator
SEs <- diag(sandwich_est) |> sqrt()
SEs_simple <- diag(sandwich_est_simple) |> sqrt()

# compare the estimates with the true values
rbind(`Estimate AGHQ` = fit$par[comp_obj$indices$coef_risk],
      `Standard errors` = SEs[comp_obj$indices$coef_risk],
      `Standard errors simple` = SEs_simple[comp_obj$indices$coef_risk],
      Truth = truth[comp_obj$indices$coef_risk])
#>                        cause1:risk:(Intercept) cause1:risk:a cause1:risk:b
#> Estimate AGHQ                       0.57739741    0.98254944     0.1390546
#> Standard errors                     0.07589777    0.08422798     0.1052925
#> Standard errors simple              0.07589783    0.08422731     0.1052928
#> Truth                               0.67000000    1.00000000     0.1000000
#>                        cause2:risk:(Intercept) cause2:risk:a cause2:risk:b
#> Estimate AGHQ                       -0.4140296    0.23000878     0.3439262
#> Standard errors                      0.1033378    0.07871263     0.1175257
#> Standard errors simple               0.1033356    0.07871211     0.1175253
#> Truth                               -0.4000000    0.25000000     0.3000000
rbind(`Estimate AGHQ` = fit$par[comp_obj$indices$coef_trajectory],
      `Standard errors` = SEs[comp_obj$indices$coef_trajectory],
      `Standard errors simple` = SEs_simple[comp_obj$indices$coef_trajectory],
      Truth = truth[comp_obj$indices$coef_trajectory])
#>                        cause1:spline1 cause1:spline2 cause1:spline3
#> Estimate AGHQ              -2.9834789     -3.6264445     -6.6774679
#> Standard errors             0.1642122      0.1650987      0.3376146
#> Standard errors simple      0.1642446      0.1651228      0.3375251
#> Truth                      -3.0513444     -3.6657486     -6.6720198
#>                        cause1:spline4 cause1:traject:(Intercept)
#> Estimate AGHQ              -4.7865951                  2.5968313
#> Standard errors             0.2232754                  0.1504230
#> Standard errors simple      0.2231936                  0.1504262
#> Truth                      -4.8559792                  2.6777525
#>                        cause1:traject:a cause1:traject:b cause2:spline1
#> Estimate AGHQ                0.88419262       0.40167971     -2.6815559
#> Standard errors              0.06577916       0.07498845      0.2115471
#> Standard errors simple       0.06577133       0.07498862      0.2116242
#> Truth                        0.80000000       0.40000000     -2.7771001
#>                        cause2:spline2 cause2:spline3 cause2:spline4
#> Estimate AGHQ              -3.1416503     -5.6506497     -4.1551289
#> Standard errors             0.1897188      0.4024503      0.2574492
#> Standard errors simple      0.1898439      0.4023770      0.2573637
#> Truth                      -3.3481062     -6.2334375     -4.6449539
#>                        cause2:traject:(Intercept) cause2:traject:a
#> Estimate AGHQ                           2.6975075       0.24617100
#> Standard errors                         0.2257982       0.06481589
#> Standard errors simple                  0.2257853       0.06481115
#> Truth                                   3.0259114       0.25000000
#>                        cause2:traject:b
#> Estimate AGHQ                -0.1692573
#> Standard errors               0.1199820
#> Standard errors simple        0.1199801
#> Truth                        -0.2000000

n_vcov <- (2L * n_causes * (2L * n_causes + 1L)) %/% 2L
Sigma
#>        [,1]   [,2]   [,3]   [,4]
#> [1,]  0.306  0.008 -0.138  0.197
#> [2,]  0.008  0.759  0.251 -0.250
#> [3,] -0.138  0.251  0.756 -0.319
#> [4,]  0.197 -0.250 -0.319  0.903
log_chol_inv(tail(fit$par, n_vcov))
#>             [,1]       [,2]        [,3]       [,4]
#> [1,]  0.33610958 -0.1474731 -0.07154827  0.1427470
#> [2,] -0.14747306  0.3688796  0.41853415 -0.1078322
#> [3,] -0.07154827  0.4185341  0.72133123 -0.4205271
#> [4,]  0.14274700 -0.1078322 -0.42052705  0.5950032

# on the log Cholesky scale
rbind(`Estimate AGHQ` = fit$par[comp_obj$indices$vcov_upper],
      `Standard errors` = SEs[comp_obj$indices$vcov_upper],
      `Standard errors simple` = SEs_simple[comp_obj$indices$vcov_upper],
      Truth = truth[comp_obj$indices$vcov_upper])
#>                        vcov:risk1:risk1 vcov:risk1:risk2 vcov:risk2:risk2
#> Estimate AGHQ                -0.5451590      -0.25437372       -0.5950783
#> Standard errors               0.2756519       0.22890496        0.3423229
#> Standard errors simple        0.2756258       0.22884814        0.3420288
#> Truth                        -0.5920851       0.01446203       -0.1380145
#>                        vcov:risk1:traject1 vcov:risk2:traject1
#> Estimate AGHQ                   -0.1234124           0.7019541
#> Standard errors                  0.1956053           0.1751353
#> Standard errors simple           0.1956090           0.1749164
#> Truth                           -0.2494700           0.2922878
#>                        vcov:traject1:traject1 vcov:risk1:traject2
#> Estimate AGHQ                      -0.7723847           0.2462218
#> Standard errors                     0.5158879           0.2157852
#> Standard errors simple              0.5145938           0.2156687
#> Truth                              -0.2485168           0.3561275
#>                        vcov:risk2:traject2 vcov:traject1:traject2
#> Estimate AGHQ                  -0.08195509             -0.7200783
#> Standard errors                 0.28633866              0.1296342
#> Standard errors simple          0.28596812              0.1297394
#> Truth                          -0.29291060             -0.1853214
#>                        vcov:traject2:traject2
#> Estimate AGHQ                      -2.3470742
#> Standard errors                     1.4385800
#> Standard errors simple              1.4235827
#> Truth                              -0.2107724

# on the original covariance matrix scale
vcov_est <- log_chol_inv(tail(fit$par, n_vcov))
vcov_est[lower.tri(vcov_est)] <- NA_real_
vcov_SE <- matrix(NA_real_, NROW(vcov_est), NCOL(vcov_est))
vcov_SE[upper.tri(vcov_SE, TRUE)] <- 
  attr(sandwich_est, "res vcov") |> diag() |> sqrt() |> 
  tail(n_vcov)

vcov_show <- cbind(Estimates = vcov_est, NA, SEs = vcov_SE) 
colnames(vcov_show) <- 
  c(rep("Est.", NCOL(vcov_est)), "", rep("SE", NCOL(vcov_est)))
print(vcov_show, na.print = "")
#>           Est.       Est.        Est.       Est.         SE        SE        SE
#> [1,] 0.3361096 -0.1474731 -0.07154827  0.1427470  0.1852985 0.2343531 0.2271176
#> [2,]            0.3688796  0.41853415 -0.1078322            0.1655676 0.2287801
#> [3,]                       0.72133123 -0.4205271                      0.1473798
#> [4,]                                   0.5950032                               
#>             SE
#> [1,] 0.2639076
#> [2,] 0.2790536
#> [3,] 0.2561399
#> [4,] 0.1855851

Sigma # the true values
#>        [,1]   [,2]   [,3]   [,4]
#> [1,]  0.306  0.008 -0.138  0.197
#> [2,]  0.008  0.759  0.251 -0.250
#> [3,] -0.138  0.251  0.756 -0.319
#> [4,]  0.197 -0.250 -0.319  0.903
```

### Delayed Entry with Different Strata

We may allow for different transformations for groups of individuals.
Specifically, we can replace the covariates for the trajectory

  
![\\vec x\_{ij}(t) = (\\vec h(t)^\\top, \\vec
z\_{ij}^\\top)^\\top](https://render.githubusercontent.com/render/math?math=%5Cvec%20x_%7Bij%7D%28t%29%20%3D%20%28%5Cvec%20h%28t%29%5E%5Ctop%2C%20%5Cvec%20z_%7Bij%7D%5E%5Ctop%29%5E%5Ctop
"\\vec x_{ij}(t) = (\\vec h(t)^\\top, \\vec z_{ij}^\\top)^\\top")  

with

  
![\\vec x\_{ij}(t) = (\\vec h\_{l\_{ij}}(t)^\\top, \\vec
z\_{ij}^\\top)^\\top](https://render.githubusercontent.com/render/math?math=%5Cvec%20x_%7Bij%7D%28t%29%20%3D%20%28%5Cvec%20h_%7Bl_%7Bij%7D%7D%28t%29%5E%5Ctop%2C%20%5Cvec%20z_%7Bij%7D%5E%5Ctop%29%5E%5Ctop
"\\vec x_{ij}(t) = (\\vec h_{l_{ij}}(t)^\\top, \\vec z_{ij}^\\top)^\\top")  

where there are
![L](https://render.githubusercontent.com/render/math?math=L "L") strata
each having their own spline basis
![h\_{l}(t)](https://render.githubusercontent.com/render/math?math=h_%7Bl%7D%28t%29
"h_{l}(t)") and
![l\_{ij}](https://render.githubusercontent.com/render/math?math=l_%7Bij%7D
"l_{ij}") is the strata that observation
![j](https://render.githubusercontent.com/render/math?math=j "j") in
cluster ![i](https://render.githubusercontent.com/render/math?math=i
"i") is in. This is supported in the package using the `strata` argument
of `mmcif_data`. We illustrate this by extending the previous example.
First, we assign new model parameters and plot the cumulative incidence
functions as before but for each strata.

``` r
# assign model parameters
n_causes <- 2L
delta <- 2

# set the betas
coef_risk <- c(.9, 1, .1, -.2, .5, 0, 0, 0, .5,
               -.4, .25, .3, 0, .5, .25, 1.5, -.25, 0) |> 
  matrix(ncol = n_causes)

# set the gammas
coef_traject <- c(-.8, -.45, -1, -.1, -.5, -.4, 
                  .8, .4, 0, .4, 0, .4, 
                  -1.2, .15, -.4, -.15, -.5, -.25, 
                  .25, -.2, 0, -.2, .25, 0) |> 
  matrix(ncol = n_causes)

# plot the conditional cumulative incidences when random effects and covariates
# are all zero
local({
  for(strata in 1:3 - 1L){
    probs <- exp(coef_risk[1 + strata * 3, ]) / 
      (1 + sum(exp(coef_risk[1 + strata * 3, ])))
    par(mar = c(5, 5, 1, 1), mfcol = c(1, 2))
    
    for(i in 1:2){
      plot(\(x) probs[i] * pnorm(
            -coef_traject[1 + strata * 2, i] * 
              atanh((x - delta / 2) / (delta / 2)) - 
              coef_traject[2 + strata * 2, i]),
           xlim = c(0, delta), ylim = c(0, 1), bty = "l",  
           xlab = sprintf("Time; strata %d", strata + 1L), 
           ylab = sprintf("Cumulative incidence; cause %d", i),
         yaxs = "i", xaxs = "i")
      grid()
    }
  }
})
```

<img src="man/figures/README-strata_assign_model_parameters-1.png" width="100%" /><img src="man/figures/README-strata_assign_model_parameters-2.png" width="100%" /><img src="man/figures/README-strata_assign_model_parameters-3.png" width="100%" />

``` r
# the probabilities of each strata
strata_prob <- c(.2, .5, .3)

# set the covariance matrix
Sigma <- c(0.306, 0.008, -0.138, 0.197, 0.008, 0.759, 0.251, 
-0.25, -0.138, 0.251, 0.756, -0.319, 0.197, -0.25, -0.319, 0.903) |> 
  matrix(2L * n_causes)
```

Then we define a simulation function.

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
    strata <- sample.int(length(strata_prob), 1L, prob = strata_prob)
    
    # keep only the relevant parameters
    coef_risk <- coef_risk[1:3 + (strata - 1L) * 3, ]
    coef_traject <- coef_traject[
        c(1:2 + (strata - 1L) * 2L, 1:2 + 6L + (strata - 1L) * 2L), ]
    
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
    
    data.frame(covs, cause, time = obs_time, cluster_id, delayed_entry, 
               strata)[keep, ]
  }) |> 
    do.call(what = rbind) |>
    transform(strata = factor(sprintf("s%d", strata)))
}
```

We sample a data set using the new simulation function.

``` r
# sample a data set
set.seed(14712915)
n_clusters <- 1000L
max_cluster_size <- 5L
dat <- sim_dat(n_clusters, max_cluster_size = max_cluster_size)

# show some stats
NROW(dat) # number of individuals
#> [1] 2518
table(dat$cause) # distribution of causes (3 is censored)
#> 
#>    1    2    3 
#>  639  791 1088

# distribution of observed times by cause
tapply(dat$time, dat$cause, quantile, 
       probs = seq(0, 1, length.out = 11), na.rm = TRUE)
#> $`1`
#>           0%          10%          20%          30%          40%          50% 
#> 2.490829e-06 2.253828e-02 1.256518e-01 3.135076e-01 6.053037e-01 9.441454e-01 
#>          60%          70%          80%          90%         100% 
#> 1.228836e+00 1.553497e+00 1.808667e+00 1.948847e+00 1.999976e+00 
#> 
#> $`2`
#>           0%          10%          20%          30%          40%          50% 
#> 2.160791e-10 1.023409e-03 1.814844e-02 1.198223e-01 3.706364e-01 8.522511e-01 
#>          60%          70%          80%          90%         100% 
#> 1.431528e+00 1.802011e+00 1.959062e+00 1.997557e+00 2.000000e+00 
#> 
#> $`3`
#>           0%          10%          20%          30%          40%          50% 
#> 0.0001884811 0.4570256790 0.8623890468 1.2217385036 1.6003899319 2.0000000000 
#>          60%          70%          80%          90%         100% 
#> 2.0000000000 2.0000000000 2.0000000000 2.0000000000 2.0000000000

# within strata
tapply(dat$time, interaction(dat$cause, dat$strata), quantile, 
       probs = seq(0, 1, length.out = 11), na.rm = TRUE)
#> $`1.s1`
#>           0%          10%          20%          30%          40%          50% 
#> 0.0001021571 0.0239064960 0.1246552149 0.3120237157 0.6869977010 0.8928431571 
#>          60%          70%          80%          90%         100% 
#> 1.2835876919 1.6481639909 1.9030873600 1.9707152763 1.9998886078 
#> 
#> $`2.s1`
#>         0%        10%        20%        30%        40%        50%        60% 
#> 0.00606347 0.09269759 0.28663148 0.46865425 0.69802525 0.91903340 1.13146065 
#>        70%        80%        90%       100% 
#> 1.39654072 1.70728367 1.86849141 1.97797010 
#> 
#> $`3.s1`
#>          0%         10%         20%         30%         40%         50% 
#> 0.004762652 0.567088551 0.903437056 1.294738798 1.601903845 2.000000000 
#>         60%         70%         80%         90%        100% 
#> 2.000000000 2.000000000 2.000000000 2.000000000 2.000000000 
#> 
#> $`1.s2`
#>         0%        10%        20%        30%        40%        50%        60% 
#> 0.00151253 0.06265371 0.19412702 0.36749402 0.59876779 0.97118195 1.20354442 
#>        70%        80%        90%       100% 
#> 1.42461147 1.66548337 1.86809811 1.99503894 
#> 
#> $`2.s2`
#>           0%          10%          20%          30%          40%          50% 
#> 2.160791e-10 2.177026e-04 5.389362e-03 4.519772e-02 2.076902e-01 6.664729e-01 
#>          60%          70%          80%          90%         100% 
#> 1.505552e+00 1.865910e+00 1.978992e+00 1.998866e+00 2.000000e+00 
#> 
#> $`3.s2`
#>          0%         10%         20%         30%         40%         50% 
#> 0.008547627 0.474137327 0.889172671 1.314716811 1.692051326 2.000000000 
#>         60%         70%         80%         90%        100% 
#> 2.000000000 2.000000000 2.000000000 2.000000000 2.000000000 
#> 
#> $`1.s3`
#>           0%          10%          20%          30%          40%          50% 
#> 2.490829e-06 1.177193e-03 1.072074e-02 5.402910e-02 4.040542e-01 9.582478e-01 
#>          60%          70%          80%          90%         100% 
#> 1.250237e+00 1.783403e+00 1.973985e+00 1.993983e+00 1.999976e+00 
#> 
#> $`2.s3`
#>           0%          10%          20%          30%          40%          50% 
#> 7.133340e-07 1.677743e-03 1.950116e-02 1.115359e-01 4.179140e-01 9.792991e-01 
#>          60%          70%          80%          90%         100% 
#> 1.574613e+00 1.854558e+00 1.970326e+00 1.997426e+00 2.000000e+00 
#> 
#> $`3.s3`
#>           0%          10%          20%          30%          40%          50% 
#> 0.0001884811 0.4129234051 0.6988777092 1.0541778269 1.3786663743 1.7758420114 
#>          60%          70%          80%          90%         100% 
#> 2.0000000000 2.0000000000 2.0000000000 2.0000000000 2.0000000000

# distribution of strata
table(dat$strata)
#> 
#>   s1   s2   s3 
#>  577 1233  708

# distribution of the left truncation time
quantile(dat$delayed_entry, probs = seq(0, 1, length.out = 11))
#>        0%       10%       20%       30%       40%       50%       60%       70% 
#> 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.1751577 
#>       80%       90%      100% 
#> 0.4210481 0.6771705 0.9998299
```

Next, we fit the model as before but this time we with strata specific
fixed effects and transformations.

``` r
library(mmcif)
comp_obj <- mmcif_data(
  ~ strata + (a + b) : strata - 1, dat, cause = cause, time = time, 
  cluster_id = cluster_id, max_time = delta, spline_df = 4L, 
  left_trunc = delayed_entry, strata = strata)
```

``` r
# we need to find the combination of the spline bases that yield a straight 
# line to construct the true values using the splines. You can skip this
comb_slope <- sapply(comp_obj$spline, \(spline){
  boundary_knots <- spline$boundary_knots
  pts <- seq(boundary_knots[1], boundary_knots[2], length.out = 1000)
  lm.fit(cbind(1, spline$expansion(pts)), pts)$coef
})

coef_traject_spline <- lapply(1:length(unique(dat$strata)), function(strata){
  slopes <- coef_traject[1 + (strata - 1) * 2, ]
  comb_slope[-1, ] * rep(slopes, each = NROW(comb_slope) - 1)
}) |> 
  do.call(what = rbind)

coef_traject_spline_fixef <- lapply(1:length(unique(dat$strata)), function(strata){
  slopes <- coef_traject[1 + (strata - 1) * 2, ]
  intercepts <- coef_traject[2 + (strata - 1) * 2, ]
  
  fixefs <- coef_traject[7:8 + (strata - 1) * 2, ]
  
  rbind(intercepts + comb_slope[1, ] * slopes,
        fixefs)
}) |> 
  do.call(what = rbind)

coef_traject_spline <- rbind(coef_traject_spline, coef_traject_spline_fixef)

# handle that model.matrix in mmcif_data gives a different permutation of 
# the parameters
permu <- c(seq(1, 7, by = 3), seq(2, 8, by = 3), seq(3, 9, by = 3))

# set true values
truth <- c(coef_risk[permu, ], 
           coef_traject_spline[c(1:12, permu + 12L), ], 
           log_chol(Sigma))

# find the starting values
system.time(start <- mmcif_start_values(comp_obj, n_threads = 4L))
#>    user  system elapsed 
#>   0.257   0.000   0.067

# we can verify that the gradient is correct again
gr_package <- mmcif_logLik_grad(
  comp_obj, truth, n_threads = 4L, is_log_chol = TRUE)
gr_num <- numDeriv::grad(
  mmcif_logLik, truth, object = comp_obj, n_threads = 4L, is_log_chol = TRUE, 
  method = "simple")

rbind(`Numerical gradient` = gr_num, `Gradient package` = gr_package)
#>                         [,1]      [,2]     [,3]      [,4]      [,5]     [,6]
#> Numerical gradient -27.70091 0.5557241 3.166838 -4.274770 -7.888019 20.91197
#> Gradient package   -27.68741 0.5709704 3.176976 -4.267357 -7.870269 20.91979
#>                         [,7]      [,8]      [,9]    [,10]    [,11]     [,12]
#> Numerical gradient -15.66595 -9.175612 -22.87761 12.02118 7.578959 -9.006909
#> Gradient package   -15.66319 -9.169393 -22.87455 12.02100 7.604564 -9.003399
#>                       [,13]    [,14]     [,15]    [,16]     [,17]    [,18]
#> Numerical gradient 10.10185 24.54336 -36.42806 23.52133 -11.56557 18.89238
#> Gradient package   10.10739 24.56492 -36.41997 23.52472 -11.55781 18.89517
#>                        [,19]    [,20]    [,21]     [,22]    [,23]     [,24]
#> Numerical gradient -1.542221 12.25438 10.04058 -16.47070 15.98909 -18.40081
#> Gradient package   -1.535355 12.26080 10.04263 -16.46659 16.00106 -18.39514
#>                        [,25]     [,26]    [,27]    [,28]    [,29]     [,30]
#> Numerical gradient -10.09486 -3.850488 15.51596 10.29838 9.962551 0.2324760
#> Gradient package   -10.09144 -3.849274 15.52042 10.30535 9.965124 0.2378367
#>                       [,31]     [,32]    [,33]    [,34]     [,35]    [,36]
#> Numerical gradient 25.19498 -16.81884 48.56917 10.18288 -9.259465 7.030678
#> Gradient package   25.20503 -16.79895 48.57437 10.19814 -9.238732 7.035232
#>                       [,37]    [,38]    [,39]    [,40]    [,41]    [,42]
#> Numerical gradient 2.899077 2.785980 7.046566 6.839909 6.109974 2.291247
#> Gradient package   2.904297 2.792787 7.049303 6.841548 6.110747 2.291344
#>                        [,43]     [,44]   [,45]    [,46]     [,47]    [,48]
#> Numerical gradient -2.230388 -28.17980 37.4577 4.681710 -28.40913 27.85014
#> Gradient package   -2.230170 -28.16845 37.4678 4.685886 -28.40242 27.86207
#>                       [,49]     [,50]      [,51]    [,52]     [,53]    [,54]
#> Numerical gradient 15.30995 -1.697091 -0.3346146 13.86711 -4.902067 42.64469
#> Gradient package   15.31827 -1.692730 -0.3307758 13.87094 -4.887511 42.66410
#>                       [,55]    [,56]     [,57]    [,58]    [,59]     [,60]
#> Numerical gradient 22.99682 47.49255 -29.16648 13.91418 26.06323 -8.109339
#> Gradient package   23.00098 47.51655 -29.14291 13.91586 26.07051 -8.101688
#>                         [,61]    [,62]     [,63]    [,64]     [,65]     [,66]
#> Numerical gradient -0.2252433 12.53314 -5.329883 25.05756 -29.12199 -10.13247
#> Gradient package   -0.2215722 12.56738 -5.194179 25.06128 -29.13691 -10.11251
#>                         [,67]    [,68]     [,69]    [,70]
#> Numerical gradient -0.9177198 19.72815 -5.217907 17.82464
#> Gradient package   -0.9119268 19.72419 -5.214498 17.84479

# optimize the log composite likelihood
system.time(fit <- mmcif_fit(start$upper, comp_obj, n_threads = 4L))
#>    user  system elapsed 
#> 167.352   0.004  43.175

# the log composite likelihood at different points
mmcif_logLik(comp_obj, truth, n_threads = 4L, is_log_chol = TRUE)
#> [1] -3188.463
mmcif_logLik(comp_obj, start$upper, n_threads = 4L, is_log_chol = TRUE)
#> [1] -3453.665
-fit$value
#> [1] -3113.285
```

Then we compute the sandwich estimator. The Hessian is currently
computed with numerical differentiation which is why it takes a while.

``` r
system.time(sandwich_est <- mmcif_sandwich(comp_obj, fit$par, n_threads = 4L))
#>    user  system elapsed 
#> 414.968   0.012 107.909

# setting order equal to zero yield no Richardson extrapolation and just
# standard symmetric difference quotient. This is less precise but faster 
system.time(sandwich_est_simple <- 
              mmcif_sandwich(comp_obj, fit$par, n_threads = 4L, order = 0L))
#>    user  system elapsed 
#>  70.603   0.004  18.388
```

We show the estimated and true the conditional cumulative incidence
functions (the dashed curves are the estimates) when the random effects
are zero and the covariates are zero, ![a\_{ij} = b\_{ij}
= 0](https://render.githubusercontent.com/render/math?math=a_%7Bij%7D%20%3D%20b_%7Bij%7D%20%3D%200
"a_{ij} = b_{ij} = 0").

``` r
local({
  # get the estimates
  coef_risk_est <- fit$par[comp_obj$indices$coef_risk] |> 
    matrix(ncol = n_causes)
  coef_traject_time_est <- fit$par[comp_obj$indices$coef_trajectory_time] |> 
    matrix(ncol = n_causes)
  coef_traject_est <- fit$par[comp_obj$indices$coef_trajectory] |> 
    matrix(ncol = n_causes)
  
  for(strata in 1:3 - 1L){
    # compute the risk probabilities  
    probs <- exp(coef_risk[1 + strata * 3, ]) / 
      (1 + sum(exp(coef_risk[1 + strata * 3, ])))
    probs_est <- exp(coef_risk_est[1 + strata, ]) / 
      (1 + sum(exp(coef_risk_est[1 + strata, ])))
  
    # plot the estimated and true conditional cumulative incidence functions. The
    # estimates are the dashed lines
    par(mar = c(5, 5, 1, 1), mfcol = c(1, 2))
    pts <- seq(1e-8, delta * (1 - 1e-8), length.out = 1000)
    
    for(i in 1:2){
      true_vals <- probs[i] * pnorm(
        -coef_traject[1 + strata * 2, i] * 
          atanh((pts - delta / 2) / (delta / 2)) - 
          coef_traject[2 + strata * 2, i])
      
      estimates <- probs_est[i] * pnorm(
        -comp_obj$time_expansion(pts, cause = i, which_strata = strata + 1L) %*% 
          coef_traject_time_est[, i] - 
          coef_traject_est[13 + strata, i]) |> drop()
      
      matplot(pts, cbind(true_vals, estimates), xlim = c(0, delta), 
              ylim = c(0, 1), bty = "l",  
              xlab = sprintf("Time; strata %d", strata + 1L), lty = c(1, 2),
              ylab = sprintf("Cumulative incidence; cause %d", i),
              yaxs = "i", xaxs = "i", type = "l", col = "black")
      grid()
    }
  }
})
```

<img src="man/figures/README-strata_compare_estimated_incidence_funcs-1.png" width="100%" /><img src="man/figures/README-strata_compare_estimated_incidence_funcs-2.png" width="100%" /><img src="man/figures/README-strata_compare_estimated_incidence_funcs-3.png" width="100%" />

Further illustrations of the estimated model are given below.

``` r
# the number of call we made
fit$counts
#> function gradient 
#>      259      256
fit$outer.iterations
#> [1] 2

# compute the standard errors from the sandwich estimator
SEs <- diag(sandwich_est) |> sqrt()
SEs_simple <- diag(sandwich_est_simple) |> sqrt()

# compare the estimates with the true values
rbind(`Estimate AGHQ` = fit$par[comp_obj$indices$coef_risk],
      `Standard errors` = SEs[comp_obj$indices$coef_risk],
      `Standard errors simple` = SEs_simple[comp_obj$indices$coef_risk],
      Truth = truth[comp_obj$indices$coef_risk])
#>                        cause1:risk:stratas1 cause1:risk:stratas2
#> Estimate AGHQ                     0.7922891           -0.1756971
#> Standard errors                   0.1630274            0.1172098
#> Standard errors simple            0.1630266            0.1172090
#> Truth                             0.9000000           -0.2000000
#>                        cause1:risk:stratas3 cause1:risk:stratas1:a
#> Estimate AGHQ                   0.001345903              1.0592843
#> Standard errors                 0.197079775              0.1606144
#> Standard errors simple          0.197078348              0.1606098
#> Truth                           0.000000000              1.0000000
#>                        cause1:risk:stratas2:a cause1:risk:stratas3:a
#> Estimate AGHQ                      0.53808422             0.03497927
#> Standard errors                    0.09332375             0.14147083
#> Standard errors simple             0.09332368             0.14147068
#> Truth                              0.50000000             0.00000000
#>                        cause1:risk:stratas1:b cause1:risk:stratas2:b
#> Estimate AGHQ                      0.05853282             -0.1336840
#> Standard errors                    0.24079177              0.1514553
#> Standard errors simple             0.24079179              0.1514553
#> Truth                              0.10000000              0.0000000
#>                        cause1:risk:stratas3:b cause2:risk:stratas1
#> Estimate AGHQ                      0.09248032           -0.2900559
#> Standard errors                    0.30169043            0.1931723
#> Standard errors simple             0.30169003            0.1931741
#> Truth                              0.50000000           -0.4000000
#>                        cause2:risk:stratas2 cause2:risk:stratas3
#> Estimate AGHQ                    0.06880132            1.4504241
#> Standard errors                  0.11723453            0.1619753
#> Standard errors simple           0.11723389            0.1619755
#> Truth                            0.00000000            1.5000000
#>                        cause2:risk:stratas1:a cause2:risk:stratas2:a
#> Estimate AGHQ                       0.3515460              0.5659669
#> Standard errors                     0.1696191              0.1020672
#> Standard errors simple              0.1696033              0.1020673
#> Truth                               0.2500000              0.5000000
#>                        cause2:risk:stratas3:a cause2:risk:stratas1:b
#> Estimate AGHQ                      -0.3921941              0.7420469
#> Standard errors                     0.1273544              0.2507171
#> Standard errors simple              0.1273546              0.2507167
#> Truth                              -0.2500000              0.3000000
#>                        cause2:risk:stratas2:b cause2:risk:stratas3:b
#> Estimate AGHQ                       0.0720306             0.06860908
#> Standard errors                     0.1624276             0.23688057
#> Standard errors simple              0.1624277             0.23688021
#> Truth                               0.2500000             0.00000000
rbind(`Estimate AGHQ` = fit$par[comp_obj$indices$coef_trajectory],
      `Standard errors` = SEs[comp_obj$indices$coef_trajectory],
      `Standard errors simple` = SEs_simple[comp_obj$indices$coef_trajectory],
      Truth = truth[comp_obj$indices$coef_trajectory])
#>                        cause1:strata1:spline1 cause1:strata1:spline2
#> Estimate AGHQ                      -3.2015649             -3.5616728
#> Standard errors                     0.2427267              0.2601628
#> Standard errors simple              0.2427929              0.2602482
#> Truth                              -2.8429994             -3.3002000
#>                        cause1:strata1:spline3 cause1:strata1:spline4
#> Estimate AGHQ                      -6.7806828             -5.1295718
#> Standard errors                     0.4561244              0.2911589
#> Standard errors simple              0.4561008              0.2911978
#> Truth                              -6.2296492             -4.5237055
#>                        cause1:strata2:spline1 cause1:strata2:spline2
#> Estimate AGHQ                      -3.7185904             -4.5735116
#> Standard errors                     0.3607189              0.3233739
#> Standard errors simple              0.3614938              0.3239846
#> Truth                              -3.5537493             -4.1252500
#>                        cause1:strata2:spline3 cause1:strata2:spline4
#> Estimate AGHQ                      -8.5917576             -6.2961255
#> Standard errors                     0.7584787              0.4286720
#> Standard errors simple              0.7598372              0.4288698
#> Truth                              -7.7870616             -5.6546319
#>                        cause1:strata3:spline1 cause1:strata3:spline2
#> Estimate AGHQ                      -1.7822819             -2.0857876
#> Standard errors                     0.2583571              0.2725111
#> Standard errors simple              0.2584518              0.2726432
#> Truth                              -1.7768747             -2.0626250
#>                        cause1:strata3:spline3 cause1:strata3:spline4
#> Estimate AGHQ                      -4.0312462             -3.0266751
#> Standard errors                     0.3677312              0.2754555
#> Standard errors simple              0.3677587              0.2755583
#> Truth                              -3.8935308             -2.8273159
#>                        cause1:traject:stratas1 cause1:traject:stratas2
#> Estimate AGHQ                        2.8069771               3.6795408
#> Standard errors                      0.2388469               0.3802249
#> Standard errors simple               0.2388900               0.3809116
#> Truth                                2.4723583               3.5529479
#>                        cause1:traject:stratas3 cause1:traject:stratas1:a
#> Estimate AGHQ                        1.7970961                 0.8953815
#> Standard errors                      0.2394907                 0.1051754
#> Standard errors simple               0.2395300                 0.1051785
#> Truth                                1.4264739                 0.8000000
#>                        cause1:traject:stratas2:a cause1:traject:stratas3:a
#> Estimate AGHQ                        -0.01421898                -0.1032773
#> Standard errors                       0.10684325                 0.1968660
#> Standard errors simple                0.10684365                 0.1968652
#> Truth                                 0.00000000                 0.0000000
#>                        cause1:traject:stratas1:b cause1:traject:stratas2:b
#> Estimate AGHQ                          0.4838568                 0.4798645
#> Standard errors                        0.1813106                 0.1593105
#> Standard errors simple                 0.1813118                 0.1593235
#> Truth                                  0.4000000                 0.4000000
#>                        cause1:traject:stratas3:b cause2:strata1:spline1
#> Estimate AGHQ                          0.5819734              -5.533182
#> Standard errors                        0.2326716               2.017375
#> Standard errors simple                 0.2326633               2.062676
#> Truth                                  0.4000000              -7.489220
#>                        cause2:strata1:spline2 cause2:strata1:spline3
#> Estimate AGHQ                       -6.886268             -14.255027
#> Standard errors                      1.408168               4.522870
#> Standard errors simple               1.435869               4.605603
#> Truth                               -8.626866             -16.184799
#>                        cause2:strata1:spline4 cause2:strata2:spline1
#> Estimate AGHQ                      -16.349271             -2.8068700
#> Standard errors                      2.906072              0.2070903
#> Standard errors simple               2.905909              0.2070736
#> Truth                              -11.544659             -2.4964067
#>                        cause2:strata2:spline2 cause2:strata2:spline3
#> Estimate AGHQ                      -2.9653221             -5.7972382
#> Standard errors                     0.2000256              0.3980463
#> Standard errors simple              0.2000474              0.3979077
#> Truth                              -2.8756220             -5.3949331
#>                        cause2:strata2:spline4 cause2:strata3:spline1
#> Estimate AGHQ                      -4.3406259             -3.3754001
#> Standard errors                     0.2491784              0.2139228
#> Standard errors simple              0.2491216              0.2141860
#> Truth                              -3.8482196             -3.1205084
#>                        cause2:strata3:spline2 cause2:strata3:spline3
#> Estimate AGHQ                      -3.8418560             -7.5333498
#> Standard errors                     0.1949043              0.4045862
#> Standard errors simple              0.1951083              0.4049575
#> Truth                              -3.5945275             -6.7436663
#>                        cause2:strata3:spline4 cause2:traject:stratas1
#> Estimate AGHQ                      -5.1395682                5.611321
#> Standard errors                     0.2326657                2.068628
#> Standard errors simple              0.2327319                2.114993
#> Truth                              -4.8102745                7.839085
#>                        cause2:traject:stratas2 cause2:traject:stratas3
#> Estimate AGHQ                        2.4675290               3.3365015
#> Standard errors                      0.2025788               0.2052060
#> Standard errors simple               0.2025419               0.2053849
#> Truth                                2.4130284               2.9537855
#>                        cause2:traject:stratas1:a cause2:traject:stratas2:a
#> Estimate AGHQ                          0.7064067                0.12637470
#> Standard errors                        0.1862586                0.08221796
#> Standard errors simple                 0.1862453                0.08221864
#> Truth                                  0.2500000                0.00000000
#>                        cause2:traject:stratas3:a cause2:traject:stratas1:b
#> Estimate AGHQ                         0.20851914                 0.1158351
#> Standard errors                       0.07836823                 0.2438536
#> Standard errors simple                0.07836695                 0.2438453
#> Truth                                 0.25000000                -0.2000000
#>                        cause2:traject:stratas2:b cause2:traject:stratas3:b
#> Estimate AGHQ                       -0.008881347                -0.1020800
#> Standard errors                      0.139340228                 0.1220400
#> Standard errors simple               0.139340383                 0.1220394
#> Truth                               -0.200000000                 0.0000000

n_vcov <- (2L * n_causes * (2L * n_causes + 1L)) %/% 2L
Sigma
#>        [,1]   [,2]   [,3]   [,4]
#> [1,]  0.306  0.008 -0.138  0.197
#> [2,]  0.008  0.759  0.251 -0.250
#> [3,] -0.138  0.251  0.756 -0.319
#> [4,]  0.197 -0.250 -0.319  0.903
log_chol_inv(tail(fit$par, n_vcov))
#>              [,1]       [,2]         [,3]        [,4]
#> [1,]  0.542162433 0.18203635 -0.007640552  0.29844029
#> [2,]  0.182036346 0.72822129  0.209117678  0.04279961
#> [3,] -0.007640552 0.20911768  0.943807137 -0.32343247
#> [4,]  0.298440294 0.04279961 -0.323432472  1.14442299

# on the log Cholesky scale
rbind(`Estimate AGHQ` = fit$par[comp_obj$indices$vcov_upper],
      `Standard errors` = SEs[comp_obj$indices$vcov_upper],
      `Standard errors simple` = SEs_simple[comp_obj$indices$vcov_upper],
      Truth = truth[comp_obj$indices$vcov_upper])
#>                        vcov:risk1:risk1 vcov:risk1:risk2 vcov:risk2:risk2
#> Estimate AGHQ                -0.3060948       0.24722558       -0.2024071
#> Standard errors               0.1954875       0.24681870        0.1543118
#> Standard errors simple        0.1954778       0.24681589        0.1543098
#> Truth                        -0.5920851       0.01446203       -0.1380145
#>                        vcov:risk1:traject1 vcov:risk2:traject1
#> Estimate AGHQ                  -0.01037672           0.2591734
#> Standard errors                 0.20120826           0.1759879
#> Standard errors simple          0.20119992           0.1759819
#> Truth                          -0.24947003           0.2922878
#>                        vcov:traject1:traject1 vcov:risk1:traject2
#> Estimate AGHQ                     -0.06589296           0.4053151
#> Standard errors                    0.10725230           0.2187823
#> Standard errors simple             0.10729922           0.2187843
#> Truth                             -0.24851681           0.3561275
#>                        vcov:risk2:traject2 vcov:traject1:traject2
#> Estimate AGHQ                  -0.07028316             -0.3215137
#> Standard errors                 0.23690625              0.1868295
#> Standard errors simple          0.23689608              0.1868211
#> Truth                          -0.29291060             -0.1853214
#>                        vcov:traject2:traject2
#> Estimate AGHQ                     -0.06857932
#> Standard errors                    0.15194135
#> Standard errors simple             0.15194689
#> Truth                             -0.21077242

# on the original covariance matrix scale
vcov_est <- log_chol_inv(tail(fit$par, n_vcov))
vcov_est[lower.tri(vcov_est)] <- NA_real_
vcov_SE <- matrix(NA_real_, NROW(vcov_est), NCOL(vcov_est))
vcov_SE[upper.tri(vcov_SE, TRUE)] <- 
  attr(sandwich_est, "res vcov") |> diag() |> sqrt() |> 
  tail(n_vcov)

vcov_show <- cbind(Estimates = vcov_est, NA, SEs = vcov_SE) 
colnames(vcov_show) <- 
  c(rep("Est.", NCOL(vcov_est)), "", rep("SE", NCOL(vcov_est)))
print(vcov_show, na.print = "")
#>           Est.      Est.         Est.        Est.        SE        SE        SE
#> [1,] 0.5421624 0.1820363 -0.007640552  0.29844029  0.211972 0.3974736 0.2962258
#> [2,]           0.7282213  0.209117678  0.04279961           0.2558114 0.3178490
#> [3,]                      0.943807137 -0.32343247                     0.1868292
#> [4,]                                   1.14442299                              
#>             SE
#> [1,] 0.3123872
#> [2,] 0.3869999
#> [3,] 0.3320827
#> [4,] 0.2143212

Sigma # the true values
#>        [,1]   [,2]   [,3]   [,4]
#> [1,]  0.306  0.008 -0.138  0.197
#> [2,]  0.008  0.759  0.251 -0.250
#> [3,] -0.138  0.251  0.756 -0.319
#> [4,]  0.197 -0.250 -0.319  0.903
```

## References

<div id="refs" class="references">

<div id="ref-Cederkvist18">

Cederkvist, Luise, Klaus K Holst, Klaus K Andersen, and Thomas H
Scheike. 2018. “Modeling the cumulative incidence function of
multivariate competing risks data allowing for within-cluster dependence
of risk and timing.” *Biostatistics* 20 (2): 199–217.
<https://doi.org/10.1093/biostatistics/kxx072>.

</div>

</div>
