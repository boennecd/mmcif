# assign model parameters
n_causes <- 2L
delta <- 2

# set the betas
coef_risk <- c(.67, 1, .1, -.4, .25, .3) |>
  matrix(ncol = n_causes)

# set the gammas
coef_traject <- c(-.8, -.45, .8, .4, -1.2, .15, .25, -.2) |>
  matrix(ncol = n_causes)

# set the covariance matrix
Sigma <- c(0.306, 0.008, -0.138, 0.197, 0.008, 0.759, 0.251,
           -0.25, -0.138, 0.251, 0.756, -0.319, 0.197, -0.25, -0.319, 0.903) |>
  matrix(2L * n_causes)

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

set.seed(8401828)
n_clusters <- 1000L
max_cluster_size <- 5L
dat <- sim_dat(n_clusters, max_cluster_size = max_cluster_size)

save.image("2D-data.RData")
