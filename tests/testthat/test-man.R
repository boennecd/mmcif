test_that("mmcif_x function calls in the manual pages gives the same", {
  skip_if_not_installed("mets")

  # prepare the data
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

  # randomly sub-sample
  set.seed(1)
  prt_use <- subset(
    prt_use, id %in% sample(unique(id), length(unique(id)) %/% 10L))

  ghq_data <- list(
    node = c(-3.43615911883774, -2.53273167423279, -1.75668364929988, -1.03661082978951, -0.342901327223705, 0.342901327223705, 1.03661082978951, 1.75668364929988, 2.53273167423279, 3.43615911883774),
    weight = c(7.6404328552326e-06, 0.00134364574678124, 0.0338743944554811, 0.240138611082314, 0.610862633735326, 0.610862633735326, 0.240138611082315, 0.033874394455481, 0.00134364574678124, 7.64043285523265e-06))

  mmcif_obj <- mmcif_data(
    ~ country - 1, prt_use, status, time, id, max_time,
    2L, strata = country, ghq_data = ghq_data)

  # get the staring values
  n_threads <- 2L
  start_vals <- mmcif_start_values(mmcif_obj, n_threads = n_threads)
  expect_snapshot_value(
    start_vals, cran = TRUE, style = "json2", tolerance = 1e-4)

  # mmcif_logLik(
  #   mmcif_obj, start_vals$upper, is_log_chol = TRUE, n_threads = n_threads) |>
  #   dput()
  expect_equal(
    mmcif_logLik(
      mmcif_obj, start_vals$upper, is_log_chol = TRUE, n_threads = n_threads),
    -2389.42390562039, tolerance = 1e-6)

  expect_snapshot_value(
    mmcif_logLik_grad(
      mmcif_obj, start_vals$upper, is_log_chol = TRUE, n_threads = n_threads),
    cran = TRUE, style = "json2", tolerance = 1e-4)

  # estimate the parameters
  skip_on_cran()
  ests <- mmcif_fit(start_vals$upper, mmcif_obj, n_threads = n_threads)
  expect_snapshot_value(ests[c("par", "value", "convergence")],
                        style = "serialize", tolerance = 1e-4)

  # get the sandwich estimator
  vcov_est <- mmcif_sandwich(
    mmcif_obj, ests$par, n_threads = n_threads, order = 2L)
  expect_snapshot_value(vcov_est, style = "serialize", tolerance = 1e-4)
})
