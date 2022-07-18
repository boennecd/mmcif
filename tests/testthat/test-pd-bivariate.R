test_that("mmcif_pd_bivariate works", {
  skip_if_not_installed("mets")
  library(mets)
  data(prt)

  # truncate the time
  max_time <- 90
  prt <- within(prt, {
    status[time >= max_time] <- 0
    time <- pmin(time, max_time)
  })

  # select the DZ twins and re-code the status
  prt_use <- subset(prt, zyg == "DZ") |>
    transform(status = ifelse(status == 0, 3L, status))

  # Gauss Hermite quadrature nodes and weights from fastGHQuad::gaussHermiteData
  ghq_data <- list(
    node = c(-3.43615911883774, -2.53273167423279, -1.75668364929988, -1.03661082978951,
             -0.342901327223705, 0.342901327223705, 1.03661082978951, 1.75668364929988,
             2.53273167423279, 3.43615911883774),
    weight = c(7.6404328552326e-06, 0.00134364574678124, 0.0338743944554811, 0.240138611082314,
               0.610862633735326,0.610862633735326, 0.240138611082315, 0.033874394455481,
               0.00134364574678124, 7.64043285523265e-06))

  # setup the object for the computation
  mmcif_obj <- mmcif_data(
    ~ country - 1, prt_use, status, time, id, max_time,
    2L, strata = country, ghq_data = ghq_data)

  # previous estimates
  par <- c(0.727279974859164, 0.640534073288067, 0.429437766165371, 0.434367104339573,
           -2.4737847536253, -1.49576564624673, -1.89966050143904, -1.58881346649412,
           -5.5431198001029, -3.5328359024178, -5.82305147022587, -3.4531896212114,
           -5.29132887832377, -3.36106297109548, -6.03690322125729, -3.49516746825624,
           2.55000711185704, 2.71995985605891, 2.61971498736444, 3.05976391058032,
           -5.97173564860957, -3.37912051983482, -5.14324860374941, -3.36396780694965,
           -6.02337246348561, -3.03754644968859, -5.51267338700737, -3.01148582224673,
           2.69665543753264, 2.59359057553995, 2.7938341786374, 2.70689750644755,
           -0.362056555418564, 0.24088005091276, 0.124070380635372, -0.246152029808377,
           -0.0445628476462479, -0.911485513197845, -0.27911988106887, -0.359648419277058,
           -0.242711959678559, -6.84897302527358) |>
    setNames(mmcif_obj$param_names$upper)

  # the test data we will use
  test_dat <- data.frame(
    country = factor(c("Norway", "Norway"), levels(prt_use$country)),
    status = c(1L, 2L), time = c(60, 75))

  # deriv_manual <- function(cause){
  #   names <- mmcif_obj$param_names$upper
  #   strata <- as.integer(test_dat$country)
  #   ti <- test_dat$time
  #   sig <- tail(par, 10) |> log_chol_inv()
  #
  #   is_spline <-
  #     lapply(cause, \(cause) grep(sprintf("cause%d.+spline\\d", cause), names))
  #   deriv_term <- mapply(\(ti, cause, is_spline, strata){
  #     log(
  #       -mmcif_obj$d_time_expansion(ti, cause = cause, which_strata = strata) %*%
  #         par[is_spline])
  #   }, ti, cause, is_spline, strata) |>
  #     sum()
  #
  #   dnorm_term_lp <- mapply(\(ti, cause, is_spline, strata){
  #     inter_term <- grep(sprintf("cause%d:traject.+countryNorway$", cause), names)
  #     dnorm_term_lp <-
  #       -mmcif_obj$time_expansion(ti, cause = cause, which_strata = strata) %*%
  #       par[is_spline] - par[inter_term]
  #     drop(dnorm_term_lp)
  #   }, ti, cause, is_spline, strata)
  #
  #   V <- matrix(0, 2L, 4L)
  #   V[1, 2 + cause[1]] <- 1
  #   V[2, 2 + cause[2]] <- 1
  #   M <- tcrossprod(V %*% sig, V)[1:2, 1:2] + diag(2)
  #
  #   dnorm_term <- mvtnorm::dmvnorm(dnorm_term_lp, sigma = M, log = TRUE)
  #
  #   library(ghqCpp)
  #   M <- solve(solve(sig) + crossprod(V))
  #   mlogit_offset <- par[grepl("risk:countryNorway", names)]
  #
  #   mu <- M %*% crossprod(V, dnorm_term_lp)
  #   lp <- mlogit_offset + mu[1:2]
  #   lp <- matrix(rep(lp, 2), length(lp))
  #
  #   logit_term <- ghqCpp::mixed_mult_logit_term(
  #     eta = lp,
  #     Sigma = M[1:2, 1:2], which_category = cause,
  #     weights = ghq_data$weight, nodes = ghq_data$node, use_adaptive = TRUE)
  #
  #   logit_term * exp(dnorm_term + deriv_term)
  # }
  # deriv_manual(c(1, 1)) |> dput()
  # deriv_manual(c(1, 2)) |> dput()
  # deriv_manual(c(2, 1)) |> dput()
  # deriv_manual(c(2, 2)) |> dput()

  der_der <- function(cause){
    mmcif_pd_bivariate(
      par = par, object = mmcif_obj, newdata = test_dat, cause = cause,
      strata = country, ghq_data = ghq_data, time = time, type =
        c("derivative", "derivative"))
  }

  expect_equal(der_der(c(1, 1)), 8.58238984981538e-05, tolerance = 1e-5)
  expect_equal(der_der(c(1, 2)), 1.65939236340968e-05, tolerance = 1e-5)
  expect_equal(der_der(c(2, 1)), 5.46437706632941e-06, tolerance = 1e-5)
  expect_equal(der_der(c(2, 2)), 5.69028230635104e-06, tolerance = 1e-5)

  # deriv_manual_max_time <- function(cause){
  #   names <- mmcif_obj$param_names$upper
  #   strata <- as.integer(test_dat$country)
  #   ti <- test_dat$time
  #   ti[1] <- max_time
  #   sig <- tail(par, 10) |> log_chol_inv()
  #
  #   is_spline <-
  #     lapply(cause, \(cause) grep(sprintf("cause%d.+spline\\d", cause), names))
  #   deriv_term <- mapply(\(ti, cause, is_spline, strata){
  #     log(
  #       -mmcif_obj$d_time_expansion(ti, cause = cause, which_strata = strata) %*%
  #         par[is_spline])
  #   }, ti[2], cause[2], is_spline[2], strata[2])
  #
  #   dnorm_term_lp <- mapply(\(ti, cause, is_spline, strata){
  #     inter_term <- grep(sprintf("cause%d:traject.+countryNorway$", cause), names)
  #     dnorm_term_lp <-
  #       -mmcif_obj$time_expansion(ti, cause = cause, which_strata = strata) %*%
  #       par[is_spline] - par[inter_term]
  #     drop(dnorm_term_lp)
  #   }, ti[2], cause[2], is_spline[2], strata[2])
  #
  #   V <- matrix(0, 1L, 4L)
  #   V[1, 2 + cause[2]] <- 1
  #   M <- tcrossprod(V %*% sig, V)[1, 1] + diag(1)
  #
  #   dnorm_term <- mvtnorm::dmvnorm(dnorm_term_lp, sigma = M, log = TRUE)
  #
  #   library(ghqCpp)
  #   M <- solve(solve(sig) + crossprod(V))
  #   mlogit_offset <- par[grepl("risk:countryNorway", names)]
  #
  #   mu <- M %*% crossprod(V, dnorm_term_lp)
  #   lp <- mlogit_offset + mu[1:2]
  #   lp <- matrix(rep(lp, 2), length(lp))
  #
  #   logit_term <- ghqCpp::mixed_mult_logit_term(
  #     eta = lp,
  #     Sigma = M[1:2, 1:2], which_category = cause,
  #     weights = ghq_data$weight, nodes = ghq_data$node, use_adaptive = TRUE)
  #
  #   logit_term * exp(dnorm_term + deriv_term)
  # }
  # deriv_manual_max_time(c(1, 1)) |> dput()
  # deriv_manual_max_time(c(1, 2)) |> dput()
  # deriv_manual_max_time(c(2, 1)) |> dput()
  # deriv_manual_max_time(c(2, 2)) |> dput()

  max_time_one <- function(cause, type){
    mmcif_pd_bivariate(
      par = par, object = mmcif_obj, newdata = test_dat, cause = cause,
      strata = country, ghq_data = ghq_data, time = c(max_time, time[2]),
      type = type)
  }

  expect_equal(max_time_one(c(1, 1), c("derivative", "derivative")), 0)
  expect_equal(max_time_one(c(1, 2), c("derivative", "derivative")), 0)
  expect_equal(max_time_one(c(2, 1), c("derivative", "derivative")), 0)
  expect_equal(max_time_one(c(2, 2), c("derivative", "derivative")), 0)

  expect_equal(
    max_time_one(c(1, 1), c("cumulative", "derivative")),
    0.00825796361590734, tolerance = 1e-5)
  expect_equal(
    max_time_one(c(1, 2), c("cumulative", "derivative")),
    0.00179906547679925, tolerance = 1e-5)
  expect_equal(
    max_time_one(c(2, 1), c("cumulative", "derivative")),
    0.000998769041537346, tolerance = 1e-5)
  expect_equal(
    max_time_one(c(2, 2), c("cumulative", "derivative")),
    0.000728556337892357, tolerance = 1e-5)

  # deriv_manual_max_time <- function(cause){
  #   names <- mmcif_obj$param_names$upper
  #   strata <- as.integer(test_dat$country)
  #   ti <- test_dat$time
  #   ti[1] <- max_time
  #   sig <- tail(par, 10) |> log_chol_inv()
  #
  #   library(ghqCpp)
  #   mlogit_offset <- par[grepl("risk:countryNorway", names)]
  #   lp <-
  #
  #   ghqCpp::mixed_mult_logit_term(
  #     eta = matrix(rep(mlogit_offset, 2), length(mlogit_offset)),
  #     Sigma = sig[1:2, 1:2], which_category = cause,
  #     weights = ghq_data$weight, nodes = ghq_data$node, use_adaptive = TRUE)
  # }
  # deriv_manual_max_time(c(1, 1)) |> dput()
  # deriv_manual_max_time(c(1, 2)) |> dput()
  # deriv_manual_max_time(c(2, 1)) |> dput()
  # deriv_manual_max_time(c(2, 2)) |> dput()

  max_time_two <- function(cause, type){
    mmcif_pd_bivariate(
      par = par, object = mmcif_obj, newdata = test_dat, cause = cause,
      strata = country, ghq_data = ghq_data, time = c(max_time, max_time),
      type = type)
  }

  expect_equal(max_time_two(c(1, 1), c("derivative", "derivative")), 0)
  expect_equal(max_time_two(c(1, 2), c("derivative", "derivative")), 0)
  expect_equal(max_time_two(c(2, 1), c("derivative", "derivative")), 0)
  expect_equal(max_time_two(c(2, 2), c("derivative", "derivative")), 0)

  expect_equal(max_time_two(c(1, 1), c("cumulative", "derivative")), 0)
  expect_equal(max_time_two(c(1, 2), c("cumulative", "derivative")), 0)
  expect_equal(max_time_two(c(2, 1), c("cumulative", "derivative")), 0)
  expect_equal(max_time_two(c(2, 2), c("cumulative", "derivative")), 0)

  expect_equal(max_time_two(c(1, 1), c("derivative", "cumulative")), 0)
  expect_equal(max_time_two(c(1, 2), c("derivative", "cumulative")), 0)
  expect_equal(max_time_two(c(2, 1), c("derivative", "cumulative")), 0)
  expect_equal(max_time_two(c(2, 2), c("derivative", "cumulative")), 0)

  expect_equal(
    max_time_two(c(1, 1), c("cumulative", "cumulative")),
    0.323763299379266, tolerance = 1e-5)
  expect_equal(
    max_time_two(c(1, 2), c("cumulative", "cumulative")),
    0.0397343239300881, tolerance = 1e-5)
  expect_equal(
    max_time_two(c(2, 1), c("cumulative", "cumulative")),
    0.0397343239300881, tolerance = 1e-5)
  expect_equal(
    max_time_two(c(2, 2), c("cumulative", "cumulative")),
    0.0153188705493866, tolerance = 1e-5)

  # adds up to one minus the survival probability
  cum_cum <- function(cause){
    mmcif_pd_bivariate(
      par = par, object = mmcif_obj, newdata = test_dat, cause = cause,
      strata = country, ghq_data = ghq_data, time = time, type =
        c("cumulative", "cumulative"))
  }

  cause_combs <- cbind(rep(1:3, each = 3), rep(1:3, 3))
  cause_combs <- cause_combs[rowSums(cause_combs) < 6, ]
  expect_equal(cum_cum(c(3L, 3L)), 1 - sum(apply(cause_combs, 1L, cum_cum)),
               tolerance = 1e-5)

  # we can recover the univariate
  test_dat_uni <- test_dat[1, ]
  uni <- function(cause){
    mmcif_pd_univariate(
      par = par, object = mmcif_obj, newdata = test_dat_uni, cause = cause,
      strata = country, ghq_data = ghq_data, time = time, type = "cumulative")
  }

  cause_combs <- cbind(1L, 1:3)
  expect_equal(uni(1L), sum(apply(cause_combs, 1L, cum_cum)), tolerance = 1e-5)

  cause_combs <- cbind(2L, 1:3)
  expect_equal(uni(2L), sum(apply(cause_combs, 1L, cum_cum)), tolerance = 1e-5)

  # integrating in one dimension gives the cumulative
  der_der <- Vectorize(function(ti1, ti2, cause){
    mmcif_pd_bivariate(
      par = par, object = mmcif_obj, newdata = test_dat, cause = cause,
      strata = country, ghq_data = ghq_data, time = c(ti1, ti2), type =
        c("derivative", "derivative"))
  }, c("ti1", "ti2"))
  cum_der <- Vectorize(function(ti2, cause){
    mmcif_pd_bivariate(
      par = par, object = mmcif_obj, newdata = test_dat, cause = cause,
      strata = country, ghq_data = ghq_data, time = c(time[1], ti2), type =
        c("cumulative", "derivative"))
  }, "ti2")

  cause_combs <- cbind(rep(1:2, 3), rep(1:3, each = 2)) |>
    split(rep(1:6, 2))

  for(combs in cause_combs){
    if(all(combs < 3)){
      tol <- 10 * sqrt(.Machine$double.eps)
      int <- integrate(
        der_der, lower = 0, upper = test_dat$time[1], ti2 = test_dat$time[2],
        rel.tol = tol, cause = combs)
      expect_equal(int$value, cum_der(test_dat$time[2], cause = combs),
                   tolerance = 100 * tol)
    }

    if(combs[2] < 3){
      int <- integrate(
        cum_der, lower = 0, upper = test_dat$time[2], rel.tol = tol,
        cause = combs)
      expect_equal(int$value, cum_cum(cause = combs),
                   tolerance = 100 * tol)
    }
  }

  # the same but when one is at the maximum time follow up
  old_time <- test_dat$time
  test_dat$time[1] <- max_time

  cause_combs <- cbind(rep(1:2, 3), rep(1:3, each = 2)) |>
    split(rep(1:6, 2))
  for(combs in cause_combs){
    # TODO: get rid of the last check when the additional branch is implemented
    #       in the C++ code
    if(combs[2] < 3 && combs[1] < 3){
      int <- integrate(
        cum_der, lower = 0, upper = test_dat$time[2], rel.tol = tol,
        cause = combs)
      expect_equal(int$value, cum_cum(cause = combs),
                   tolerance = 100 * tol)
    }
  }

  test_dat$time <- old_time

  # TODO: implement this case and add it one of the loops above
  tmp_data <- data.frame(
    country = factor(c("Norway", "Norway"), levels(prt_use$country)),
    status = c(3L, 2L), time = c(max_time, 75))

  expect_error(
    mmcif_pd_bivariate(
      par = par, object = mmcif_obj, newdata = tmp_data, cause = status,
      strata = country, ghq_data = ghq_data, time = time, type =
        c("cumulative", "cumulative")),
    "the case where one is censored and at the maximum follow-up and the other is a cummulative is not implemented")

  # TODO: test branches all branches (e.g. when at maximum time, see the code).
  #       The brances are:
  #         a. the first one is density and the second is an observed cif.
  #            Check (also at maximum time).
  #         b. both are observed cifs.
  #            Check (also at maximum time).
  #         c. the first one is censored and the second one is an observed cif.
  #            Check except at the maximum follow-up time.
  # TODO: test the code with truncation
  # TODO: test mmcif_pd_cond
})
