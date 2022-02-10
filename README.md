
# MMCIF: Mixed Multivariate Cumulative Incidence Functions

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
for the trajectory are subject to monotonically decreasing.

## Example

We start with a simple example where there are
![K = 2](https://render.githubusercontent.com/render/math?math=K%20%3D%202 "K = 2")
competing risks and

![\\begin{align\*} \\vec x\_{ij}(t) &= \\left(\\text{arcthan}\\left(\\frac{t - \\delta/2}{\\delta/2}\\right), 1, a\_{ij}, b\_{ij}\\right) \\\\ a\_{ij} &\\sim N(0, 1) \\\\
 b\_{ij} &\\sim \\text{Unif}(-1, 1)\\\\ \\vec z\_{ij} &= (1, a\_{ij}, b\_{ij}) \\end{align\*}](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Balign%2A%7D%20%5Cvec%20x_%7Bij%7D%28t%29%20%26%3D%20%5Cleft%28%5Ctext%7Barcthan%7D%5Cleft%28%5Cfrac%7Bt%20-%20%5Cdelta%2F2%7D%7B%5Cdelta%2F2%7D%5Cright%29%2C%201%2C%20a_%7Bij%7D%2C%20b_%7Bij%7D%5Cright%29%20%5C%5C%20a_%7Bij%7D%20%26%5Csim%20N%280%2C%201%29%20%5C%5C%0A%20b_%7Bij%7D%20%26%5Csim%20%5Ctext%7BUnif%7D%28-1%2C%201%29%5C%5C%20%5Cvec%20z_%7Bij%7D%20%26%3D%20%281%2C%20a_%7Bij%7D%2C%20b_%7Bij%7D%29%20%5Cend%7Balign%2A%7D "\begin{align*} \vec x_{ij}(t) &= \left(\text{arcthan}\left(\frac{t - \delta/2}{\delta/2}\right), 1, a_{ij}, b_{ij}\right) \\ a_{ij} &\sim N(0, 1) \\
 b_{ij} &\sim \text{Unif}(-1, 1)\\ \vec z_{ij} &= (1, a_{ij}, b_{ij}) \end{align*}")

We set the parameters below and plot the conditional cumulative
incidence function when the random effects are zero and
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

# plot the conditional cumulative incidence when random effects and covariates
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
uniformly sampled from two to the maximum size. The censoring times are
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
    
    n_obs <- sample.int(max_cluster_size - 1L, 1L) + 1L
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
n_clusters <- 1000L
max_cluster_size <- 5L
dat <- sim_dat(n_clusters, max_cluster_size = max_cluster_size)

# show some stats
NROW(dat) # number of individuals
#> [1] 3573
table(dat$cause) # distribution of causes (3 is censored)
#> 
#>    1    2    3 
#> 1527  623 1423

# distribution of observed times by cause
tapply(dat$time, dat$cause, quantile, 
       probs = seq(0, 1, length.out = 11), na.rm = TRUE)
#> $`1`
#>           0%          10%          20%          30%          40%          50% 
#> 1.037532e-06 4.393086e-03 2.549201e-02 8.553160e-02 1.954232e-01 4.183183e-01 
#>          60%          70%          80%          90%         100% 
#> 7.995275e-01 1.341151e+00 1.727572e+00 1.944052e+00 1.999995e+00 
#> 
#> $`2`
#>           0%          10%          20%          30%          40%          50% 
#> 0.0009666954 0.0601977831 0.1514903014 0.2832340702 0.4514909913 0.6938606243 
#>          60%          70%          80%          90%         100% 
#> 0.9512206696 1.3040148597 1.5707337955 1.8322967104 1.9983337883 
#> 
#> $`3`
#>           0%          10%          20%          30%          40%          50% 
#> 0.0005559647 0.3034055659 0.6145405803 0.9248298794 1.3095068850 1.8354668445 
#>          60%          70%          80%          90%         100% 
#> 2.0000000000 2.0000000000 2.0000000000 2.0000000000 2.0000000000
```

Then we setup the C++ object to do the computation.

``` r
library(mmcif)
comp_obj <- mmcif_data(
  ~ a + b, dat, cause = cause, time = time, cluster_id = cluster_id, 
  max_time = delta)
```

The time to compute the log composite likelihood is illustrated below.

``` r
NCOL(comp_obj$pair_indices) # the number of pairs in the composite likelihood
#> [1] 5209

# assign a function to compute the log composite likelihood
library(fastGHQuad)
#> Loading required package: Rcpp
ghq_data <- with(gaussHermiteData(5L), list(node = x, weight = w))

ll_func <- \(par, n_threads = 1L)
  mmcif:::mmcif_logLik(
    comp_obj$comp_obj, par = c(coef_risk, coef_traject, Sigma), 
    ghq_data = ghq_data, n_threads = n_threads)

# the log composite likelihood at the true parameters
ll_func(c(coef_risk, coef_traject, Sigma))
#> [1] -8865.391

# check the time to compute the log composite likelihood
bench::mark(
  `one thread` = ll_func(n_threads = 1L, c(coef_risk, coef_traject, Sigma)),
  `two threads` = ll_func(n_threads = 2L, c(coef_risk, coef_traject, Sigma)),
  `three threads` = ll_func(n_threads = 3L, c(coef_risk, coef_traject, Sigma)),
  `four threads` = ll_func(n_threads = 4L, c(coef_risk, coef_traject, Sigma)))
#> # A tibble: 4 × 6
#>   expression         min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 one thread      91.3ms   91.5ms      10.9    23.5KB        0
#> 2 two threads     47.4ms   48.1ms      20.8      288B        0
#> 3 three threads   31.6ms   32.2ms      31.1      288B        0
#> 4 four threads    24.6ms   24.9ms      40.0      288B        0
```

Then we optimize the parameters (TODO: there will be a wrapper to work
with the log Cholesky decomposition and an estimation function; the
estimation time will be much reduced when the gradient is implemented).

``` r
# find the starting values
system.time(start <- mmcif_start_values(comp_obj, n_threads = 4L))
#>    user  system elapsed 
#>   0.049   0.000   0.014

# the maximum likelihood without the random effects. Note that this is not 
# comparable with the composite likelihood
attr(start, "logLik")
#> [1] -3135.4

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
# the covariance matrix. We handle the monotonicity constraint by log 
# transforming the slopes
ll_func_chol <- \(par, n_threads = 1L, ghq = ghq_data){
  n_vcov <- (2L * n_causes * (2L * n_causes + 1L)) %/% 2L
  par <- c(head(par, -n_vcov), tail(par, n_vcov) |> log_chol_inv())
  coef_trajectory_time <- comp_obj$indices$coef_trajectory_time
  par[coef_trajectory_time] <- -exp(par[coef_trajectory_time])
  
  mmcif:::mmcif_logLik(
    comp_obj$comp_obj, par = par, ghq_data = ghq, n_threads = n_threads)
}

# set starting values and true value
truth <- c(coef_risk, coef_traject, log_chol(Sigma))
truth[comp_obj$indices$coef_trajectory_time] <-
  log(-truth[comp_obj$indices$coef_trajectory_time])

start$lower[comp_obj$indices$coef_trajectory_time] <- 
  log(-start$lower[comp_obj$indices$coef_trajectory_time])

# optimize the log composite likelihood. First with a Laplace approximation
ghq_data_Laplace <- with(gaussHermiteData(1L), list(node = x, weight = w))
system.time(
  fit_Laplace <- optim(
    start$lower, \(par) -ll_func_chol(par, n_threads = 4L, ghq_data_Laplace), 
    control = list(maxit = 10000L)))
#>    user  system elapsed 
#> 291.600   0.007  72.917

# then we improve by using Adaptive Gauss-Hermite quadrature
system.time(
  fit <- optim(
    fit_Laplace$par, \(par) -ll_func_chol(par, n_threads = 4L, ghq_data), 
    control = list(maxit = 10000L)))
#>    user  system elapsed 
#> 490.367   0.008 122.614

# the log composite likelihood at different points
ll_func_chol(truth)
#> [1] -8865.391
ll_func_chol(start$lower, 4L)
#> [1] -9502.087
ll_func_chol(fit_Laplace$par, 4L)
#> [1] -8850.741
-fit$value
#> [1] -8843.321
```

``` r
# It took quite a few iterations. This may be much smaller with a gradient 
# approximation
fit_Laplace$counts
#> function gradient 
#>     3513       NA
fit$counts
#> function gradient 
#>     4253       NA

# compare the estimates with the true values
rbind(`Estimate Laplace` = head(fit_Laplace$par, length(coef_risk)),
      `Estimate AGHQ` = head(fit$par, length(coef_risk)),
      Truth = c(coef_risk))
#>                       [,1]      [,2]      [,3]       [,4]      [,5]      [,6]
#> Estimate Laplace 0.6820169 0.9213513 0.2091015 -0.4297888 0.3127709 0.3116543
#> Estimate AGHQ    0.6678179 0.9624318 0.2472146 -0.4646059 0.2948015 0.3738843
#> Truth            0.6700000 1.0000000 0.1000000 -0.4000000 0.2500000 0.3000000
rbind(`Estimate Laplace` = fit_Laplace$par[comp_obj$indices$coef_trajectory],
      `Estimate AGHQ` = fit$par[comp_obj$indices$coef_trajectory],
      Truth = truth[comp_obj$indices$coef_trajectory])
#>                        [,1]       [,2]      [,3]      [,4]      [,5]       [,6]
#> Estimate Laplace -0.2328143 -0.4796474 0.8355155 0.4058797 0.1630115 0.03025068
#> Estimate AGHQ    -0.2255714 -0.4848914 0.8394861 0.4171070 0.1481107 0.09699559
#> Truth            -0.2231436 -0.4500000 0.8000000 0.4000000 0.1823216 0.15000000
#>                       [,7]       [,8]
#> Estimate Laplace 0.1585224 -0.3294909
#> Estimate AGHQ    0.1766491 -0.3056605
#> Truth            0.2500000 -0.2000000

n_vcov <- (2L * n_causes * (2L * n_causes + 1L)) %/% 2L
Sigma
#>        [,1]   [,2]   [,3]   [,4]
#> [1,]  0.306  0.008 -0.138  0.197
#> [2,]  0.008  0.759  0.251 -0.250
#> [3,] -0.138  0.251  0.756 -0.319
#> [4,]  0.197 -0.250 -0.319  0.903
log_chol_inv(tail(fit$par, n_vcov))
#>             [,1]       [,2]        [,3]        [,4]
#> [1,]  0.40337303  0.2713524 -0.07630755  0.04580343
#> [2,]  0.27135242  1.3752989  0.22911260 -0.26401829
#> [3,] -0.07630755  0.2291126  0.68733205 -0.31152421
#> [4,]  0.04580343 -0.2640183 -0.31152421  0.72142485
log_chol_inv(tail(fit_Laplace$par, n_vcov))
#>             [,1]       [,2]        [,3]       [,4]
#> [1,]  0.24202930  0.3596424 -0.05965436 -0.1187252
#> [2,]  0.35964241  1.2683645  0.19424410 -0.2862462
#> [3,] -0.05965436  0.1942441  0.66533894 -0.3472695
#> [4,] -0.11872517 -0.2862462 -0.34726954  0.7829966
```

## TODOs

The package is still under development. Here are a few TODOs:

-   Implement the gradient. Most of the work is done in already in the
    package at <https://github.com/boennecd/ghq-cpp/tree/main/ghqCpp>.
-   Implement a setup functions that creates an expansion using natural
    cubic splines to make the model more flexible. The log composite
    likelihood can then be optimized with an optimizer that allows for
    linear inequality constraints.
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
