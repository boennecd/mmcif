#include "mmcif-misc.h"
#include "testthat.h"

context("mmcif functions work") {
  test_that("mmcif_logLik works with singleton data") {
    /*
     set.seed(111)
     Sig <- rWishart(1, 3, diag(3)) |> round(3) |> drop()
     x <- c(-1, 0, 1)
     dput(mvtnorm::dmvnorm(x, sigma = Sig, log = TRUE))
     */
    arma::mat Sig{2.498, -0.326, 0.361, -0.326, 0.629, -1.326, 0.361, -1.326, 3.632};
    Sig.reshape(3, 3);
    arma::vec x{-1, 0, 1};
    constexpr double truth{-3.51260818319038};
    ghqCpp::simple_mem_stack<double> mem;

    expect_true
      (std::abs(log_dmvn(x, Sig, mem) - truth) < std::abs(truth) * 1e-8);
  }
}
