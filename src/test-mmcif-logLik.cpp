#include <testthat.h>
#include "mmcif-logLik.h"

namespace {
constexpr double ghq_nodes[]{-3.43615911883774, -2.53273167423279, -1.75668364929988, -1.03661082978951, -0.342901327223705, 0.342901327223705, 1.03661082978951, 1.75668364929988, 2.53273167423279, 3.43615911883774},
               ghq_weights[]{7.6404328552326e-06, 0.00134364574678124, 0.0338743944554811, 0.240138611082314, 0.610862633735326, 0.610862633735326, 0.240138611082315, 0.033874394455481, 0.00134364574678124, 7.64043285523265e-06};

const ghqCpp::ghq_data ghq_dat_use{ghq_nodes, ghq_weights, 10};

} // namespace

context("mmcif_log_Lik works as expected") {
  test_that("mmcif_log_Lik works with singleton data") {
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
     dput(covs_traject_w_tine <- g_w_covs(obs_time) |> drop())
     dput(d_covs_traject_w_time <- d_g_w_covs(obs_time) |> drop())

     dput(gl <- fastGHQuad::gaussHermiteData(10L))

# compute the log densities we need. Start with the censored where the
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
     res <- res + dnorm(-covs_traject_w_tine %*% coefs_traject[, cause] |> drop(),
     sd = sqrt(1 + vcov[idx_cause, idx_cause]), log = TRUE)

     M <- solve(vcov)
     M[idx_cause, idx_cause] <- M[idx_cause, idx_cause] + 1
     M <- solve(M)

     mean_cond <- M[, idx_cause] *
     sum(-covs_traject_w_tine %*% coefs_traject[, cause])
     etas <- covs_risk %*% coefs_risk |> drop()
     etas <- etas - mean_cond[1:n_causes]

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

     offset <- -covs_traject_w_tine %*% coefs_traject[, cause] |> drop()

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
      double const res{mmcif_log_Lik(par, indexer, dat, mem, ghq_dat_use)};
      constexpr double truth{-1.21423993502317};
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);
    }

    {
      // the observed case works
      constexpr unsigned cause{2L};
      mmcif_data dat{covs_traject, d_covs_traject, covs_risk, true, cause};
      double const res{mmcif_log_Lik(par, indexer, dat, mem, ghq_dat_use)};
      constexpr double truth{-3.76328226753652};
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);
    }

    // the censored case with a probability of the trajectory in (0, 1) works
    mmcif_data dat{covs_traject, d_covs_traject, covs_risk, true, n_causes};
    double const res{mmcif_log_Lik(par, indexer, dat, mem, ghq_dat_use)};
    constexpr double truth{-1.04608338904514};
    expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);
  }
}
