#include <testthat.h>
#include "mmcif-logLik.h"

namespace {
constexpr double ghq_nodes[]{-3.43615911883774, -2.53273167423279, -1.75668364929988, -1.03661082978951, -0.342901327223705, 0.342901327223705, 1.03661082978951, 1.75668364929988, 2.53273167423279, 3.43615911883774},
               ghq_weights[]{7.6404328552326e-06, 0.00134364574678124, 0.0338743944554811, 0.240138611082314, 0.610862633735326, 0.610862633735326, 0.240138611082315, 0.033874394455481, 0.00134364574678124, 7.64043285523265e-06};

const ghqCpp::ghq_data ghq_dat_use{ghq_nodes, ghq_weights, 10};

} // namespace

context("mmcif_logLik works as expected with singleton data") {
  test_that("mmcif_logLik works") {
    /*
     set.seed(321)
     n_cov_traject <- 4L
     n_cov_risk <- 2L
     n_causes <- 3L

     g <- \(x) cbind(x, log(x))
     d_g <- \(x) cbind(1, 1/x)

     dput(covs_risk <- runif(n_cov_risk, -1) |> round(3))
     covs_traject <- runif(n_cov_traject, -1) |> round(3)
     vcov <-
     drop(rWishart(1, 2 * n_causes, diag(1/2/n_causes, 2 * n_causes))) |>
     drop() |> round(3)

     coefs_risk <- runif(n_causes * n_cov_risk, -1) |>
     round(3) |> matrix(n_cov_risk)
     coefs_traject <- runif(n_causes * n_cov_traject, -1) |>
     round(3) |> matrix(n_cov_traject)
     coefs_traject <- rbind(-.5, -.25, coefs_traject)
     c(coefs_risk, coefs_traject, vcov) |> dput()

# sample a plausible value
     u <- mvtnorm::rmvnorm(1, sigma = vcov) |> drop()
     etas <- covs_risk %*% coefs_risk |> drop()
     etas <- etas + u[1:n_causes]
     phats <- c(exp(etas), 1) / (1 + sum(exp(etas)))
     dput(cause <- sample(n_causes + 1, 1, prob = phats))
     stopifnot(cause < n_causes + 1) # we did not implement this case

     g_w_covs <- \(ti)
     cbind(g(ti), rep(covs_traject, each = length(ti)) |> matrix(length(ti)))
     rng <- runif(1)
     obs_time <- uniroot(
     \(ti){
     lp <- -g_w_covs(ti) %*% coefs_traject[, cause] - u[cause + n_causes]
     pnorm(lp) - rng
     },
     c(1e-16, 1000))$root

     d_g_w_covs <- \(ti)
     cbind(d_g(ti), matrix(0., length(ti), length(covs_traject)))

# the matrices we will need
     dput(covs_traject_w_time <- g_w_covs(obs_time) |> drop())
     dput(d_covs_traject_w_time <- d_g_w_covs(obs_time) |> drop())

     dput(gl <- fastGHQuad::gaussHermiteData(10L))

# compute the log densities we need. Start with the censored case
# trajectory probability is one
     library(ghqCpp)
     local({
     etas <- covs_risk %*% coefs_risk |> drop()
     res <- mixed_mult_logit_term(
     eta = as.matrix(etas), Sigma = vcov[1:n_causes, 1:n_causes],
     which_category = 0L, weights = gl$w, nodes = gl$x)
     dput(log(res))
     })

# then with the observed outcome
     local({
     idx_cause <- cause + n_causes

     res <- log(-d_covs_traject_w_time %*% coefs_traject[, cause]) |> drop()
     res <- res + dnorm(-covs_traject_w_time %*% coefs_traject[, cause] |> drop(),
     sd = sqrt(1 + vcov[idx_cause, idx_cause]), log = TRUE)

     M <- solve(vcov)
     M[idx_cause, idx_cause] <- M[idx_cause, idx_cause] + 1
     M <- solve(M)

     mean_cond <- M[, idx_cause] *
     sum(-covs_traject_w_time %*% coefs_traject[, cause])
     etas <- covs_risk %*% coefs_risk |> drop()
     etas <- etas + mean_cond[1:n_causes]

     res_integral <- mixed_mult_logit_term(
     eta = as.matrix(etas), Sigma = M[1:n_causes, 1:n_causes],
     which_category = cause, weights = gl$w, nodes = gl$x)

     dput(res + log(res_integral))
     })

# the censored case
     local({
     integrand <- 1
     etas <- covs_risk %*% coefs_risk |> drop()
     vcov_sub <- vcov[1:n_causes, 1:n_causes]

     for(cause in 1:n_causes){
     idx_cause <- cause + n_causes
     z <- -solve(vcov_sub, vcov[1:n_causes, idx_cause]) |> drop()
     s <- sqrt(1 + vcov[idx_cause, idx_cause] -
     sum(vcov[idx_cause, 1:n_causes] * (-z)))

     offset <- -covs_traject_w_time %*% coefs_traject[, cause] |> drop()

     integrand <- integrand - mixed_mult_logit_n_probit_term(
     eta = as.matrix(etas),
     which_category = cause, s = s, eta_probit = offset,
     Sigma = vcov_sub, z = z, weights = gl$w, nodes = gl$x)
     }
     dput(log(integrand))
     })
     */

    constexpr size_t n_causes{3}, n_cov_risk{2}, n_cov_traject{6};
    constexpr double covs_risk[]{0.912, 0.875},
                  covs_traject[]{4.1936419665673, 1.43356956082597, -0.524, -0.49, -0.219, -0.318},
                d_covs_traject[]{1, 0.238456217286129, 0, 0, 0, 0},
                           par[]{0.355, -0.737, -0.475, -0.365, 0.353, -0.874, -0.5, -0.25, -0.871, 0.322, -0.791, 0.474, -0.5, -0.25, -0.85, 0.782, 0.78, -0.123, -0.5, -0.25, -0.81, 0.091, -0.198, -0.974, 0.771, -0.126, -0.198, 0.068, 0.879, -0.145, -0.126, 1.158, 0.184, 1.053, -0.211, 0.595, -0.198, 0.184, 1.456, -0.43, -1.171, -0.332, 0.068, 1.053, -0.43, 2.234, 0.431, 0.462, 0.879, -0.211, -1.171, 0.431, 1.673, -0.001, -0.145, 0.595, -0.332, 0.462, -0.001, 0.832};

    ghqCpp::simple_mem_stack<double> mem;
    param_indexer indexer{n_cov_risk, n_cov_traject, n_causes};

    {
      // the censored case with one as the probability of the trajectory
      mmcif_data dat{covs_traject, d_covs_traject, covs_risk, false, n_causes};
      double const res{mmcif_logLik(par, indexer, dat, mem, ghq_dat_use)};
      constexpr double truth{-1.21423993502317};
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);
    }

    {
      // the observed case works
      constexpr unsigned cause{2L};
      mmcif_data dat{covs_traject, d_covs_traject, covs_risk, true, cause};
      double const res{mmcif_logLik(par, indexer, dat, mem, ghq_dat_use)};
      constexpr double truth{-4.26186714415087};
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);
    }

    // the censored case with a probability of the trajectory in (0, 1) works
    mmcif_data dat{covs_traject, d_covs_traject, covs_risk, true, n_causes};
    double const res{mmcif_logLik(par, indexer, dat, mem, ghq_dat_use)};
    constexpr double truth{-1.04608338904514};
    expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);
  }
}

namespace {
/*
 set.seed(321)
 n_cov_traject <- 4L
 n_cov_risk <- 2L
 n_causes <- 3L

 g <- \(x) cbind(x, log(x))
 d_g <- \(x) cbind(1, 1/x)

 vcov <-
 drop(rWishart(1, 2 * n_causes, diag(1/2/n_causes, 2 * n_causes))) |>
 drop() |> round(3)

 coefs_risk <- runif(n_causes * n_cov_risk, -1) |>
 round(3) |> matrix(n_cov_risk)
 coefs_traject <- runif(n_causes * n_cov_traject, -1) |>
 round(3) |> matrix(n_cov_traject)
 coefs_traject <- rbind(-.5, -.25, coefs_traject)
 c(coefs_risk, coefs_traject, vcov) |> dput()

# sample a plausible value
 u <- mvtnorm::rmvnorm(1, sigma = vcov) |> drop()

 g_w_covs <- \(ti, covs_traject)
 cbind(g(ti), rep(covs_traject, each = length(ti)) |> matrix(length(ti)))
 d_g_w_covs <- \(ti, covs_traject)
 cbind(d_g(ti), matrix(0., length(ti), length(covs_traject)))

 sample_obs <- \(){
 covs_risk <- runif(n_cov_risk, -1) |> round(3)
 covs_traject <- runif(n_cov_traject, -1) |> round(3)

 etas <- covs_risk %*% coefs_risk |> drop()
 etas <- etas + u[1:n_causes]
 phats <- c(exp(etas), 1) / (1 + sum(exp(etas)))
 cause <- sample(n_causes + 1, 1, prob = phats)
 stopifnot(cause < n_causes + 1) # we did not implement this case

 rng <- runif(1)
 obs_time <- uniroot(
 \(ti){
 lp <- -g_w_covs(ti, covs_traject) %*% coefs_traject[, cause] -
 u[cause + n_causes]
 pnorm(lp) - rng
 },
 c(1e-16, 1000))$root

 list(covs_risk = covs_risk,
 covs_traject_w_time = g_w_covs(obs_time, covs_traject) |> drop(),
 d_covs_traject_w_time = d_g_w_covs(obs_time, covs_traject) |> drop(),
 cause = cause)
 }

 dput(obs <- replicate(2, sample_obs(), simplify = FALSE))

 dput(gl <- fastGHQuad::gaussHermiteData(10L))
 */

constexpr size_t n_cov_risk{2},
                 n_cov_traject{6},
                 n_causes{3},
                 cause[]{0, 2};

param_indexer const indexer{n_cov_risk, n_cov_traject, n_causes};

constexpr double covs_risk1[]{-0.605, -0.876},
                 covs_risk2[]{0.584, 0.576},
              covs_traject1[]{0.85860256800553, -0.152449132274265, -0.734, 0.062, 0.957, 0.198},
              covs_traject2[]{4.78877056139744, 1.5662737107113, 0.304, 0.267, -0.669, 0.204},
            d_covs_traject1[]{1, 1.16468321580138, 0, 0, 0, 0},
            d_covs_traject2[]{1, 0.208821865065129, 0, 0, 0, 0},
                        par[]{-0.647, 0.076, -0.439, -0.98, -0.806, 0.324, -0.5, -0.25, 0.355, -0.737, -0.475, -0.365, -0.5, -0.25, 0.353, -0.874, -0.871, 0.322, -0.5, -0.25, -0.791, 0.474, -0.85, 0.782, 1.974, -0.235, 0.154, 0.195, 0.108, 0.059, -0.235, 0.401, 0.163, -0.161, 0.597, -0.155, 0.154, 0.163, 0.4, 0.026, 0.047, 0.179, 0.195, -0.161, 0.026, 0.506, -0.493, 0.018, 0.108, 0.597, 0.047, -0.493, 2.031, 0.052, 0.059, -0.155, 0.179, 0.018, 0.052, 0.578};
} // namespace

context("mmcif_logLik works as expected with bivariate data") {
  test_that("mmcif_logLik works when both individuals are observed") {
    /*
     library(ghqCpp)
     local({
     etas <- sapply(obs, \(x) x$covs_risk %*% coefs_risk)
     covs_traject <- sapply(
     obs,
     \(x, cause) -x$covs_traject_w_time %*% coefs_traject[, x$cause])
     d_covs_traject <- sapply(
     obs,
     \(x) -x$d_covs_traject_w_time %*% coefs_traject[, x$cause])

     V <- matrix(0, 2, 2 * n_causes)
     V[1, obs[[1]]$cause + n_causes] <- 1
     V[2, obs[[2]]$cause + n_causes] <- 1
     out <- sum(log(d_covs_traject)) +
     mvtnorm::dmvnorm(
     covs_traject, sigma = diag(2) + V %*% tcrossprod(vcov, V), log = TRUE)

     vcov_cond <- solve(crossprod(V) + solve(vcov))
     mean_cond <- vcov_cond %*% crossprod(V, covs_traject) |> drop()

     integral <- mixed_mult_logit_term(
     eta = etas + mean_cond[1:n_causes],
     Sigma = vcov_cond[1:n_causes, 1:n_causes],
     which_category = sapply(obs, `[[`, "cause"),
     weights = gl$w, nodes = gl$x)

     dput(out + log(integral))
     })
     */
    ghqCpp::simple_mem_stack<double> mem;

    mmcif_data
      obs1{covs_traject1, d_covs_traject1, covs_risk1, true, cause[0]},
      obs2{covs_traject2, d_covs_traject2, covs_risk2, true, cause[1]};

    double const res{mmcif_logLik(par, indexer, obs1, obs2, mem, ghq_dat_use)};
    constexpr double truth{-7.77301433778719};
    expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);
  }

  test_that("mmcif_logLik works when one individual is observed and one is censored") {
    /*
     library(ghqCpp)
     local({
     etas <- sapply(obs, \(x) x$covs_risk %*% coefs_risk)
     covs_traject <- sapply(
     obs,
     \(x, cause) -x$covs_traject_w_time %*% coefs_traject[, x$cause])
     d_covs_traject <- sapply(
     obs,
     \(x) -x$d_covs_traject_w_time %*% coefs_traject[, x$cause])

     v <- numeric(2 * n_causes)
     v[n_causes + obs[[1]]$cause] <- 1
     out <- log(d_covs_traject[1]) + dnorm(
     covs_traject[1], sd = sqrt(1 + v %*% vcov %*% v), log = TRUE)

     vcov_cond <- solve(tcrossprod(v) + solve(vcov))
     mean_cond <- drop(vcov_cond %*% v) * covs_traject[1]

     integral <- mixed_mult_logit_term(
     eta = etas + mean_cond[1:n_causes],
     Sigma = vcov_cond[1:n_causes, 1:n_causes],
     which_category = c(obs[[1]]$cause, 0L),
     weights = gl$w, nodes = gl$x)

     dput(out + log(integral))
     })
     */

    ghqCpp::simple_mem_stack<double> mem;

    {
      mmcif_data
        obs1{covs_traject1, d_covs_traject1, covs_risk1, true, cause[0]},
        obs2{covs_traject2, d_covs_traject2, covs_risk2, false, n_causes};

      double const res
        {mmcif_logLik(par, indexer, obs1, obs2, mem, ghq_dat_use)};
      constexpr double truth{-4.5240025456571};
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);
    }

    /*
     local({
     etas <- sapply(obs, \(x) x$covs_risk %*% coefs_risk)
     covs_traject <- sapply(
     obs,
     \(x, cause) -x$covs_traject_w_time %*% coefs_traject[, x$cause])
     d_covs_traject <- sapply(
     obs,
     \(x) -x$d_covs_traject_w_time %*% coefs_traject[, x$cause])

     v <- numeric(2 * n_causes)
     v[n_causes + obs[[1]]$cause] <- 1
     out <- log(d_covs_traject[1]) + dnorm(
     covs_traject[1], sd = sqrt(1 + v %*% vcov %*% v), log = TRUE)

     vcov_cond <- solve(tcrossprod(v) + solve(vcov))
     mean_cond <- drop(vcov_cond %*% v) * covs_traject[1]

     integral <- mixed_mult_logit_term(
     eta = as.matrix(etas[, 1]) + mean_cond[1:n_causes],
     Sigma = vcov_cond[1:n_causes, 1:n_causes],
     which_category = obs[[1]]$cause, weights = gl$w, nodes = gl$x)

     integrals <- sapply(1:n_causes, \(cause){
     idx_cause <- cause + n_causes
     rng_coefs <- solve(vcov_cond[1:n_causes, 1:n_causes],
     vcov_cond[1:n_causes, idx_cause])

     s <- sqrt(1 + vcov_cond[idx_cause, idx_cause] -
     vcov_cond[idx_cause, 1:n_causes] %*% rng_coefs)
     shift <- -obs[[2]]$covs_traject_w_time %*% coefs_traject[, cause] -
     mean_cond[idx_cause]

     rng_coefs <- -rng_coefs

     mixed_mult_logit_n_probit_term(
     eta = etas + mean_cond[1:n_causes],
     which_category = c(obs[[1]]$cause, cause), s = s, eta_probit = shift,
     Sigma = vcov_cond[1:n_causes, 1:n_causes], z = rng_coefs,
     weights = gl$w, nodes = gl$x)
     })
     integral <- integral - sum(integrals)

     dput(out + log(integral))
     })
     */

    {
      mmcif_data
        obs1{covs_traject1, d_covs_traject1, covs_risk1, true, cause[0]},
        obs2{covs_traject2, d_covs_traject2, covs_risk2, true, n_causes};

      constexpr double truth{-4.38551973954109};
      {
        double const res
          {mmcif_logLik(par, indexer, obs1, obs2, mem, ghq_dat_use)};
        expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);
      }

      double const res
        {mmcif_logLik(par, indexer, obs2, obs1, mem, ghq_dat_use)};
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);
    }
  }

  test_that("mmcif_logLik works when both individuals are censored") {
    /*
     local({
     etas <- sapply(obs, \(x) x$covs_risk %*% coefs_risk)
     mixed_mult_logit_term(
     eta = etas, Sigma = vcov[1:n_causes, 1:n_causes],
     which_category = integer(2), weights = gl$w, nodes = gl$x) |>
     log() |> dput()
     })
     */

    ghqCpp::simple_mem_stack<double> mem;

    {
      mmcif_data
        obs1{covs_traject1, d_covs_traject1, covs_risk1, false, n_causes},
        obs2{covs_traject2, d_covs_traject2, covs_risk2, false, n_causes};

      double const res
        {mmcif_logLik(par, indexer, obs1, obs2, mem, ghq_dat_use)};
      constexpr double truth{-3.01479043790748};
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);
    }

    /*
     local({
     etas <- sapply(obs, \(x) x$covs_risk %*% coefs_risk)

     integral <- mixed_mult_logit_term(
     eta = as.matrix(etas[, 2]), Sigma = vcov[1:n_causes, 1:n_causes],
     which_category = 0L, weights = gl$w, nodes = gl$x)

     integrals <- sapply(1:n_causes, \(cause){
     lp_traject <- -obs[[1]]$covs_traject_w_time %*% coefs_traject[, cause] |>
     drop()

     idx_cause <- cause + n_causes
     rng_coefs <-
     solve(vcov[1:n_causes, 1:n_causes], vcov[1:n_causes, idx_cause])

     s <- sqrt(1 + vcov[idx_cause, idx_cause] -
     vcov[idx_cause, 1:n_causes] %*% rng_coefs)

     rng_coefs <- -rng_coefs

     mixed_mult_logit_n_probit_term(
     eta = etas, which_category = c(cause, 0L), s = s,
     eta_probit = lp_traject, Sigma = vcov[1:n_causes, 1:n_causes],
     z = rng_coefs, weights = gl$w, nodes = gl$x)
     })

     log(integral - sum(integrals)) |> dput()
     })
     */

    {
      mmcif_data
        obs1{covs_traject1, d_covs_traject1, covs_risk1, true, n_causes},
        obs2{covs_traject2, d_covs_traject2, covs_risk2, false, n_causes};

      constexpr double truth{-2.18417107496657};
      {
        double const res
          {mmcif_logLik(par, indexer, obs1, obs2, mem, ghq_dat_use)};
        expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);
      }

      double const res
        {mmcif_logLik(par, indexer, obs2, obs1, mem, ghq_dat_use)};
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);
    }

    /*
     library(ghqCpp)
     local({
     vcov_sub <- vcov[1:n_causes, 1:n_causes]

     integrals_substract <- sapply(obs, \(obs_i){
     sapply(1:n_causes, \(cause){
     etas <- obs_i$covs_risk %*% coefs_risk |> drop()

     idx_cause <- cause + n_causes
     z <- -solve(vcov_sub, vcov[1:n_causes, idx_cause]) |> drop()
     s <- sqrt(1 + vcov[idx_cause, idx_cause] -
     sum(vcov[idx_cause, 1:n_causes] * (-z)))

     offset <- -obs_i$covs_traject_w_time %*% coefs_traject[, cause] |>
     drop()

     mixed_mult_logit_n_probit_term(
     eta = as.matrix(etas),
     which_category = cause, s = s, eta_probit = offset,
     Sigma = vcov_sub, z = z, weights = gl$w, nodes = gl$x)
     })
     })

     etas <- sapply(obs, \(x) x$covs_risk %*% coefs_risk)
     rng_coefs <- solve(vcov_sub, vcov[1:n_causes, -(1:n_causes)]) |> t()
     vcov_cond <- vcov[1:n_causes + n_causes, 1:n_causes + n_causes] -
     vcov[1:n_causes + n_causes, 1:n_causes] %*% t(rng_coefs)

     integrals_adds <- sapply(
     1:n_causes,
     \(cause1) sapply(
     1:n_causes, \(cause2){
     V <- matrix(0, 2, n_causes)
     V[1, cause1] <- 1
     V[2, cause2] <- 1

     shift <- -c(obs[[1]]$covs_traject_w_time %*% coefs_traject[, cause1],
     obs[[2]]$covs_traject_w_time %*% coefs_traject[, cause2])

     rng_coefs <- -V %*% rng_coefs
     vcov_cond_pbvn <- diag(2) + V %*% tcrossprod(vcov_cond, V)

     mixed_mult_logit_n_cond_pbvn(
     eta = etas, which_category = c(cause1, cause2),
     eta_pbvn = -shift, Psi = vcov_cond_pbvn, V = -rng_coefs,
     Sigma = vcov_sub, weights = gl$w, nodes = gl$x)
     }))

     log(1 - sum(integrals_substract) + sum(integrals_adds)) |> dput()
     })
     */

    {
      mmcif_data
        obs1{covs_traject1, d_covs_traject1, covs_risk1, true, n_causes},
        obs2{covs_traject2, d_covs_traject2, covs_risk2, true, n_causes};

      double const res
        {mmcif_logLik(par, indexer, obs1, obs2, mem, ghq_dat_use)};
      constexpr double truth{-2.00247380641032};
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);
    }
  }
}
