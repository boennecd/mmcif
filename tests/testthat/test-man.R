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
    node = c(-2.02018287045609, -0.958572464613819, 2.73349801485409e-17, 0.958572464613819, 2.02018287045609),
    weight = c(0.0199532420590459, 0.393619323152241, 0.945308720482941, 0.393619323152241, 0.0199532420590459))

  mmcif_obj <- mmcif_data(
    ~ country - 1, prt_use, status, time, id, max_time,
    2L, strata = country, ghq_data = ghq_data)

  # get the staring values
  n_threads <- 2L
  start_vals <- mmcif_start_values(mmcif_obj, n_threads = n_threads)
  expect_snapshot_value(start_vals, cran = TRUE, style = "json2")

  # mmcif_logLik(
  #   mmcif_obj, start_vals$upper, is_log_chol = TRUE, n_threads = n_threads) |>
  #   dput()
  expect_equal(
    mmcif_logLik(
      mmcif_obj, start_vals$upper, is_log_chol = TRUE, n_threads = n_threads),
    -2389.42390562039)

  expect_snapshot_value(
    mmcif_logLik_grad(
      mmcif_obj, start_vals$upper, is_log_chol = TRUE, n_threads = n_threads),
    cran = TRUE, style = "json2")

  # estimate the parameters
  skip_on_cran()
  ests <- mmcif_fit(start_vals$upper, mmcif_obj, n_threads = n_threads)
  expect_snapshot_value(ests[c("par", "value", "convergence")],
                        style = "serialize")

  # get the sandwich estimator
  vcov_est <- mmcif_sandwich(
    mmcif_obj, ests$par, n_threads = n_threads, order = 2L)
  expect_snapshot_value(vcov_est, style = "serialize")
})
