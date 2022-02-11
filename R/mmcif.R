#' Setup Object to Compute the Log Composite Likelihood
#'
#' @param formula \code{formula} for covariates in the risk and trajectories.
#' @param data \code{data.frame} with the covariate and outcome information.
#' @param cause an integer vector with the cause of each outcome. If there are
#' \code{n_causes} of outcome, then the vector should have values in
#' \code{1:(n_causes + 1)} with \code{n_causes + 1} indicating censoring.
#' @param time a numeric vector with the observed times.
#' @param cluster_id an integer vector with the cluster id of each individual.
#' @param max_time the maximum time after which there are no observed events. It
#' is denoted by \eqn{\delta} in the original article.
#' @param spline_df degrees of freedom to use for each spline in the
#' cumulative incidence functions.
#'
#' @importFrom stats model.frame model.matrix terms
#' @export
mmcif_data <- function(formula, data, cause, time, cluster_id, max_time,
                       spline_df = 3L){
  stopifnot(inherits(formula, "formula"),
            is.data.frame(data),
            length(spline_df) == 1, is.finite(spline_df), spline_df > 0)

  # setup the design matrix
  mf <- model.frame(formula, data)
  covs_risk <- model.matrix(terms(mf), mf)

  cause <- eval(substitute(cause), data, parent.frame())
  time_observed <- eval(substitute(time), data, parent.frame())
  cluster_id <- eval(substitute(cluster_id), data, parent.frame())

  n_causes <- length(unique(cause)) - 1L

  stopifnot(
    length(cause) > 0L,
    NROW(covs_risk) == length(cause),
    length(time_observed) == length(cause),
    length(cluster_id) == length(cluster_id),
    n_causes > 0L, all(1:(n_causes + 1L) %in% cause),
    is.numeric(max_time), length(max_time) == 1, is.finite(max_time),
    max_time > 0,
    max(time_observed[cause %in% 1:n_causes]) < max_time)

  # add the time transformations
  time_trans <- function(x) atanh((x - max_time / 2) / (max_time / 2))
  d_time_trans <- function(x) {
    outer <- (x - max_time / 2) / (max_time / 2)
    1/((1 - outer^2) * (max_time / 2))
  }

  time_observed_trans <- suppressWarnings(
    ifelse(time_observed < max_time, time_trans(time_observed), NA))

  is_observed <- cause %in% 1:n_causes
  spline <- monotone_ns(time_observed_trans[is_observed], df = spline_df)

  time_expansion <- function(x){
    x <- suppressWarnings(time_trans(x))
    spline$expansion(x, ders = 0L)
  }

  d_time_expansion <- function(x){
    z <- suppressWarnings(time_trans(x))
    d_x <- suppressWarnings(d_time_trans(x))
    spline$expansion(z, ders = 1L) * d_x
  }

  covs_trajectory <-  cbind(time_expansion(time_observed), covs_risk)
  d_covs_trajectory <- cbind(
    d_time_expansion(time_observed),
    matrix(0, NROW(covs_risk), NCOL(covs_risk)))

  has_finite_trajectory_prob <- time_observed < max_time

  # find all permutation of indices in each cluster
  stopifnot(
    # TODO: support singleton clusters (the C++ code is there)
    max(table(cluster_id)) > 1L)

  row_id <- seq_len(length(cluster_id))
  pair_indices <- tapply(row_id, cluster_id, function(ids){
    # TODO: do this smarter
    out <- expand.grid(first = ids, second = ids)
    out <- as.matrix(out[out$first > out$second, ])
    t(out)
  }, simplify = FALSE)

  pair_indices <- do.call(cbind, pair_indices)

  # create the C++ object
  # TODO: this function should be exposed to advanced users
  comp_obj <- mmcif_data_holder(
    covs_risk = t(covs_risk), covs_trajectory = t(covs_trajectory),
    d_covs_trajectory = t(d_covs_trajectory),
    has_finite_trajectory_prob = has_finite_trajectory_prob,
    cause = cause - 1L, n_causes = n_causes, pair_indices = pair_indices - 1L,
    singletons = integer())

  # create an object to do index the parameters
  n_coef_risk <- NCOL(covs_risk) * n_causes
  n_coef_trajectory <- (spline_df + NCOL(covs_risk)) * n_causes

  coef_trajectory_time <- matrix(
    seq_len(n_coef_trajectory) + n_coef_risk, spline_df + NCOL(covs_risk))
  coef_trajectory_time <- coef_trajectory_time[1:spline_df, ]

  idx_coef_trajectory <- seq_len(n_coef_trajectory) + n_coef_risk

  vcov_dim <- 2L * n_causes
  vcov_full <- max(idx_coef_trajectory) + seq_len(vcov_dim^2)
  vcov_lower <- max(idx_coef_trajectory) +
    seq_len((vcov_dim * (vcov_dim + 1L)) %/% 2L)

  indices <- list(
    coef_risk = seq_len(n_coef_risk),
    coef_trajectory = idx_coef_trajectory,
    coef_trajectory_time = c(coef_trajectory_time),
    vcov_full = vcov_full,
    vcov_lower = vcov_lower,
    n_causes = n_causes,
    n_par = max(vcov_full),
    n_par_lower_tri = max(vcov_lower),
    n_par_wo_vcov = max(idx_coef_trajectory))

  # create the constrain matrices
  n_constraints_per_spline <- NROW(spline$constraints)
  constraints <- matrix(0., n_constraints_per_spline * n_causes,
                        indices$n_par_wo_vcov)
  for(i in 1:n_causes)
    constraints[
      1:n_constraints_per_spline + (i - 1L) * n_constraints_per_spline,
      matrix(indices$coef_trajectory_time, spline_df)[, i]] <-
      spline$constraints

  constraints <- list(
    wo_vcov = constraints,
    vcov_full =
      cbind(constraints, matrix(0., NROW(constraints), vcov_dim^2)),
    vcov_lower =
      cbind(constraints, matrix(0., NROW(constraints),
                                (vcov_dim * (vcov_dim + 1L)) %/% 2L)))

  structure(
    list(comp_obj = comp_obj, pair_indices = pair_indices,
         covs_risk = covs_risk, covs_trajectory = covs_trajectory,
         time_observed = time_observed, cause = cause,
         time_trans = time_trans, d_time_trans = d_time_trans,
         time_expansion = time_expansion, d_time_expansion = d_time_expansion,
         max_time = max_time, indices = indices, spline = spline,
         d_covs_trajectory = d_covs_trajectory, constraints = constraints),
    class = "mmcif")
}

.check_n_threads <- function(n_threads){
  stopifnot(length(n_threads) == 1, is.finite(n_threads), n_threads > 0)
}

.log_chol <- \(x){
  x <- chol(x)
  diag(x) <- diag(x) |> log()
  x[upper.tri(x, TRUE)]
}

.log_chol_inv <- \(x){
  dim <- (sqrt(8 * length(x) + 1) - 1) / 2
  out <- matrix(0, dim, dim)
  out[upper.tri(out, TRUE)] <- x
  diag(out) <- diag(out) |> exp()
  crossprod(out)
}

#' Finds Staring Values
#'
#' @param object an object from \code{\link{mmcif_data}}.
#' @param n_threads the number of threads to use.
#' @param vcov_start starting value for the covariance matrix of the random
#' effects. \code{NULL} yields the identity matrix.
#'
#' @importFrom stats lm.fit constrOptim na.omit sd
#' @export
mmcif_start_values <- function(object, n_threads = 4L, vcov_start = NULL){
  stopifnot(inherits(object, "mmcif"))
  .check_n_threads(n_threads)

  # find intercepts and slopes in a simple model
  n_causes <- object$indices$n_causes
  cause <- object$cause
  time_observed <- object$time_observed

  is_observed <- cause %in% 1:n_causes

  slope_n_time <- tapply(
    time_observed[is_observed], cause[is_observed], function(time_observed){
    time_observed_trans <- object$time_trans(time_observed)

    time_sd <- sd(time_observed_trans)

    c(slope = -1 / time_sd, intercept = mean(time_observed_trans) / time_sd)
  })

  # find the linear combination which gives a line and a slope
  comb_slope <- local({
    boundary_knots <- object$spline$boundary_knots
    pts <- seq(boundary_knots[1], boundary_knots[2], length.out = 1000)
    lm.fit(cbind(1, object$spline$expansion(pts)), pts)$coef
  })

  comb_line <- local({
    design_mat <- na.omit(object$covs_trajectory)
    lm.fit(design_mat, rep(1, NROW(design_mat)))$coef
  })

  # start the optimization. Start only with the trajectories
  par <- numeric(object$indices$n_par_wo_vcov)
  coef_trajectory_time <- object$indices$coef_trajectory_time
  coef_trajectory <- object$indices$coef_trajectory

  par[coef_trajectory_time] <-
    comb_slope[-1] %o% sapply(slope_n_time, `[[`, "slope")
  par[coef_trajectory] <- par[coef_trajectory] +
    comb_line %o%
      (sapply(slope_n_time, `[[`, "intercept") +
         comb_slope[1] * sapply(slope_n_time, `[[`, "slope"))

  fn <- function(x, with_risk){
    if(!with_risk)
      par[coef_trajectory] <- x
    else
      par <- x
    -mcif_logLik(data_ptr = object$comp_obj, par = par, n_threads = n_threads,
                 with_risk = with_risk)
  }
  gr <- function(x, with_risk){
    if(!with_risk)
      par[coef_trajectory] <- x
    else
      par <- x
    grs <- -mcif_logLik_grad(data_ptr = object$comp_obj, par = par,
                             n_threads = n_threads, with_risk = with_risk)
    if(with_risk) grs else grs[coef_trajectory]
  }

  constraints <- object$constraints$wo_vcov
  opt_trajectory <-
    constrOptim(
      par[coef_trajectory], fn, gr, method = "BFGS",
      ui = constraints[, coef_trajectory], ci = rep(1e-8, NROW(constraints)),
      control = list(maxit = 1000L), with_risk = FALSE)
  stopifnot(opt_trajectory$convergence == 0)

  par[coef_trajectory] <- opt_trajectory$par
  opt <- constrOptim(
    par, fn, gr, method = "BFGS", ui = constraints,
    ci = rep(1e-8, NROW(constraints)), control = list(maxit = 1000L),
    with_risk = TRUE)
  stopifnot(opt$convergence == 0)

  if(is.null(vcov_start))
    vcov_start <- diag(2L * n_causes)
  stopifnot(is.matrix(vcov_start), all(dim(vcov_start) == 2L * n_causes),
            all(is.finite(vcov_start)))

  structure(list(full = c(opt$par, c(vcov_start)),
                 lower = c(opt$par, .log_chol(vcov_start))),
            logLik = -opt$value)
}
