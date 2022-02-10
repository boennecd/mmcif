#include "mcif-logLik.h"
#include <testthat.h>

namespace {

/*
 set.seed(123)
 n_cov_traject <- 4L
 n_cov_risk <- 2L
 n_causes <- 3L

 g <- \(x) cbind(x, log(x))
 d_g <- \(x) cbind(1, 1/x)

 dput(covs_risk <- runif(n_cov_risk, -1) |> round(3))
 covs_traject <- runif(n_cov_traject, -1) |> round(3)

 coefs_risk <- runif(n_causes * n_cov_risk, -1) |>
 round(3) |> matrix(n_cov_risk)
 coefs_traject <- runif(n_causes * n_cov_traject, -1) |>
 round(3) |> matrix(n_cov_traject)
 coefs_traject <- rbind(-.5, -.25, coefs_traject)
 c(coefs_risk, coefs_traject) |> dput()

# sample a plausible value
 etas <- covs_risk %*% coefs_risk |> drop()
 phats <- c(exp(etas), 1) / (1 + sum(exp(etas)))
 dput(cause <- sample(n_causes + 1, 1, prob = phats))
 stopifnot(cause < n_causes + 1) # we did not implement this case

 g_w_covs <- \(ti)
 cbind(g(ti), rep(covs_traject, each = length(ti)) |> matrix(length(ti)))
 rng <- runif(1)
 obs_time <- uniroot(
 \(ti){
 lp <- -g_w_covs(ti) %*% coefs_traject[, cause]
 pnorm(lp) - rng
 },
 c(1e-16, 1000))$root

 d_g_w_covs <- \(ti)
 cbind(d_g(ti), matrix(0., length(ti), length(covs_traject)))

# the matrices we will need
 dput(covs_traject_w_time <- g_w_covs(obs_time) |> drop())
 dput(d_covs_traject_w_time <- d_g_w_covs(obs_time) |> drop())

 get_n_remove <- \(x, n){
 out <- x[1:n]
 eval(substitute(out <- out[-(1:n)], list(out = substitute(x), n = n)),
 parent.frame())
 out
 }
  */

constexpr size_t n_causes{3}, n_cov_risk{2}, n_cov_traject{6},
                 cause{1};
constexpr double covs_risk[]{-0.425, 0.577},
              covs_traject[]{0.00846107419528661, -4.77227913998632, -0.182, 0.766, 0.881, -0.909},
            d_covs_traject[]{1, 118.188302917503, 0, 0, 0, 0},
                       par[]{0.056, 0.785, 0.103, -0.087, 0.914, -0.093, -0.5, -0.25, 0.355, 0.145, -0.794, 0.8, -0.5, -0.25, -0.508, -0.916, -0.344, 0.909, -0.5, -0.25, 0.779, 0.386, 0.281, 0.989};
constexpr double gr_eps{1e-5};
} // namespace


context("mcif functions work") {
  test_that("works with an observed cause") {
    /*
     f1 <- \(x){
     coefs_risk <- get_n_remove(x, length(coefs_risk)) |>
     matrix(n_cov_risk)
     coefs_traject <- matrix(x, 2L + n_cov_traject)

     log(-sum(d_covs_traject_w_time * coefs_traject[, cause])) +
     dnorm(-sum(covs_traject_w_time * coefs_traject[, cause]), log = TRUE)
     }

     par <- c(coefs_risk, coefs_traject)
     f1(par) |> dput()
     numDeriv::grad(f1, par) |> dput()

     f2 <- \(x){
     coefs_risk <- get_n_remove(x, length(coefs_risk)) |>
     matrix(n_cov_risk)
     coefs_traject <- matrix(x, 2L + n_cov_traject)

     out <- log(-sum(d_covs_traject_w_time * coefs_traject[, cause])) +
     dnorm(-sum(covs_traject_w_time * coefs_traject[, cause]), log = TRUE)

     etas <- covs_risk %*% coefs_risk |> drop()
     denom <- exp(etas) |> sum() + 1

     out + etas[cause] - log(denom)
     }
     f2(par) |> dput()
     numDeriv::grad(f2, par) |> dput()
     */

    ghqCpp::simple_mem_stack<double> mem;
    param_indexer indexer{n_cov_risk, n_cov_traject, n_causes};
    mmcif_data const dat{covs_traject, d_covs_traject, covs_risk, true, cause};

    {
      constexpr double truth{2.33273860259647},
                      offset{-1},
                 true_grad[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.0286300077774044, -6.55678707593064, -0.10004644690055, 0.421074606123595, 0.484290767579494, -0.499682528669878, 0, 0, 0, 0, 0, 0};

      double res = mcif_logLik<false>(par, indexer, dat, mem);
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

      std::vector<double> grad(indexer.n_par_wo_vcov(), offset);
      res = mcif_logLik_grad<false>(par, grad.data(), indexer, dat, mem);
      expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

      for(size_t i = 0; i < grad.size(); ++i)
        expect_true
          (std::abs(grad[i] - offset - true_grad[i]) <
            (std::abs(true_grad[i]) + gr_eps) * gr_eps);
    }

    constexpr double truth{0.830481658506154},
                    offset{2},
               true_grad[]{0.15964478037357, -0.216741266624087, -0.330383467378871, 0.448544142550531, 0.0667994752319878, -0.0906901110958889, 0, 0, 0, 0, 0, 0, -0.0286300077774044, -6.55678707593064, -0.10004644690055, 0.421074606123595, 0.484290767579494, -0.499682528669878, 0, 0, 0, 0, 0, 0};

    double res = mcif_logLik<true>(par, indexer, dat, mem);
    expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

    std::vector<double> grad(indexer.n_par_wo_vcov(), offset);
    res = mcif_logLik_grad<true>(par, grad.data(), indexer, dat, mem);
    expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

    for(size_t i = 0; i < grad.size(); ++i)
      expect_true
        (std::abs(grad[i] - offset - true_grad[i]) <
          (std::abs(true_grad[i]) + gr_eps) * gr_eps);
  }

  test_that("works with a censored observation with probablity one for the trajectory") {
    /*
     f1 <- \(x){
     coefs_risk <- get_n_remove(x, length(coefs_risk)) |>
     matrix(n_cov_risk)
     coefs_traject <- matrix(x, 2L + n_cov_traject)

     etas <- covs_risk %*% coefs_risk |> drop()
     denom <- exp(etas) |> sum() + 1

     - log(denom)
     }

     par <- c(coefs_risk, coefs_traject)
     f1(par) |> dput()
     numDeriv::grad(f1, par) |> dput()
     */

    ghqCpp::simple_mem_stack<double> mem;
    param_indexer indexer{n_cov_risk, n_cov_traject, n_causes};
    mmcif_data const dat
      {covs_traject, d_covs_traject, covs_risk, false, n_causes};

    {
      constexpr double offset{3};

      double res = mcif_logLik<false>(par, indexer, dat, mem);
      expect_true(res == 0);

      std::vector<double> grad(indexer.n_par_wo_vcov(), offset);
      res = mcif_logLik_grad<false>(par, grad.data(), indexer, dat, mem);
      expect_true(res == 0);

      for(size_t i = 0; i < grad.size(); ++i)
        expect_true(grad[i] == offset);
    }

    constexpr double truth{-1.40828294409031},
                    offset{2},
               true_grad[]{0.15964478037357, -0.216741266624087, 0.0946165325639635, -0.128455857365894, 0.0667994752319878, -0.0906901110958889, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    double res = mcif_logLik<true>(par, indexer, dat, mem);
    expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

    std::vector<double> grad(indexer.n_par_wo_vcov(), offset);
    res = mcif_logLik_grad<true>(par, grad.data(), indexer, dat, mem);
    expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

    for(size_t i = 0; i < grad.size(); ++i)
      expect_true
        (std::abs(grad[i] - offset - true_grad[i]) <
          (std::abs(true_grad[i]) + gr_eps) * gr_eps);
  }

  test_that("works with a censored observation with probablity for the trajectory in (0, 1)") {
    /*
     f1 <- \(x){
     coefs_risk <- get_n_remove(x, length(coefs_risk)) |>
     matrix(n_cov_risk)
     coefs_traject <- matrix(x, 2L + n_cov_traject)

     etas <- covs_risk %*% coefs_risk |> drop()
     denom <- exp(etas) |> sum() + 1

     integrand_terms <- sapply(1:n_causes, \(cause){
     lp_traject <- -sum(covs_traject_w_time * coefs_traject[, cause])
     pnorm(lp_traject) * exp(etas[cause]) / denom
     })

     log(1 - sum(integrand_terms))
     }

     par <- c(coefs_risk, coefs_traject)
     f1(par) |> dput()
     numDeriv::grad(f1, par) |> dput()
     */

    ghqCpp::simple_mem_stack<double> mem;
    param_indexer indexer{n_cov_risk, n_cov_traject, n_causes};
    mmcif_data const dat
      {covs_traject, d_covs_traject, covs_risk, true, n_causes};

    {
      constexpr double offset{3};

      double res = mcif_logLik<false>(par, indexer, dat, mem);
      expect_true(res == 0);

      std::vector<double> grad(indexer.n_par_wo_vcov(), offset);
      res = mcif_logLik_grad<false>(par, grad.data(), indexer, dat, mem);
      expect_true(res == 0);

      for(size_t i = 0; i < grad.size(); ++i)
        expect_true(grad[i] == offset);
    }

    constexpr double truth{-0.532001042543614},
                    offset{7},
               true_grad[]{0.0443874774426018, -0.0602625281392007, 0.0477034585344386, -0.0647644600256424, -0.0190907533895791, 0.0259185051915334, 0.00211928952787073, -1.19533772377672, -0.0455864921686501, 0.191864027517853, 0.220668679188212, -0.227681985710412, 0.00109987164955784, -0.620357940705621, -0.0236585375711419, 0.0995738448324722, 0.114522920758956, -0.118162695760887, 0.000711200267143155, -0.401136552678053, -0.0152981103003366, 0.0643865520743632, 0.0740529404153219, -0.0764064958610772};

    double res = mcif_logLik<true>(par, indexer, dat, mem);
    expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

    std::vector<double> grad(indexer.n_par_wo_vcov(), offset);
    res = mcif_logLik_grad<true>(par, grad.data(), indexer, dat, mem);
    expect_true(std::abs(res - truth) < std::abs(truth) * 1e-8);

    for(size_t i = 0; i < grad.size(); ++i)
      expect_true
      (std::abs(grad[i] - offset - true_grad[i]) <
        (std::abs(true_grad[i]) + gr_eps) * gr_eps);
  }
}
