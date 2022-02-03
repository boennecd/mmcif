#include <testthat.h>
#include "param-indexer.h"

context("param_indexer") {
  test_that("param_indexer works"){
    param_indexer const index{3, 2, 4};

    expect_true(index.n_cov_risk() == 3);
    expect_true(index.n_cov_traject() == 2);
    expect_true(index.n_causes() == 4);

    expect_true(index.risk() == index.risk(0));
    expect_true(index.risk(0) == 0);
    expect_true(index.risk(1) == 3);
    expect_true(index.risk(2) == 6);
    expect_true(index.risk(3) == 9);

    expect_true(index.surv() == index.surv(0));
    expect_true(index.surv(0) == 12);
    expect_true(index.surv(1) == 14);
    expect_true(index.surv(2) == 16);
    expect_true(index.surv(3) == 18);

    expect_true(index.vcov() == 20);

    expect_true(index.n_par<false>() == 36);
    expect_true(index.n_par<true>() == 30);

    std::vector<std::string> const names_full
      {"beta[1,1]","beta[2,1]","beta[3,1]",
       "beta[1,2]","beta[2,2]","beta[3,2]",
       "beta[1,3]","beta[2,3]","beta[3,3]",
       "beta[1,4]","beta[2,4]","beta[3,4]",
       "gamma[1,1]", "gamma[2,1]",
       "gamma[1,2]", "gamma[2,2]",
       "gamma[1,3]", "gamma[2,3]",
       "gamma[1,4]", "gamma[2,4]",
       "vcov[1,1]","vcov[2,1]","vcov[3,1]","vcov[4,1]",
       "vcov[1,2]","vcov[2,2]","vcov[3,2]","vcov[4,2]",
       "vcov[1,3]","vcov[2,3]","vcov[3,3]","vcov[4,3]",
       "vcov[1,4]","vcov[2,4]","vcov[3,4]","vcov[4,4]"};

    std::vector<std::string> const names_lower
      {"beta[1,1]","beta[2,1]","beta[3,1]",
       "beta[1,2]","beta[2,2]","beta[3,2]",
       "beta[1,3]","beta[2,3]","beta[3,3]",
       "beta[1,4]","beta[2,4]","beta[3,4]",
       "gamma[1,1]", "gamma[2,1]",
       "gamma[1,2]", "gamma[2,2]",
       "gamma[1,3]", "gamma[2,3]",
       "gamma[1,4]", "gamma[2,4]",
       "vcov[1]", "vcov[2]", "vcov[3]", "vcov[4]", "vcov[5]",
       "vcov[6]", "vcov[7]", "vcov[8]", "vcov[9]", "vcov[10]"};

    {
       auto res = index.param_names<false>();
       expect_true(res.size() == names_full.size());
       for(size_t i = 0; i < names_full.size(); ++i)
         expect_true(res[i] == names_full[i]);
    }

    auto res = index.param_names<true>();
    expect_true(res.size() == names_lower.size());
    for(size_t i = 0; i < names_lower.size(); ++i)
      expect_true(res[i] == names_lower[i]);
  }
}
