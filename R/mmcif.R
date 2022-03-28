.ghq_data_default_if_null <- function(ghq_data){
  if(is.null(ghq_data))
    ghq_data <- list(
      node = c(-2.02018287045609, -0.958572464613819, 2.73349801485409e-17, 0.958572464613819, 2.02018287045609),
      weight = c(0.0199532420590459,  0.393619323152241, 0.945308720482941, 0.393619323152241, 0.0199532420590459))
  ghq_data
}

.check_is_log_chol <- function(is_log_chol)
  stopifnot(is.logical(is_log_chol), length(is_log_chol) == 1,
            is_log_chol %in% c(FALSE, TRUE))

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
#' @param left_trunc numeric vector with left-truncation times. \code{NULL}
#' implies that there is not any individuals with left-truncation.
#' @param ghq_data the default Gauss-Hermite quadrature nodes and weights to
#' use. It should be stored as a list with two elements called \code{"node"}
#' and \code{"weight"}. A default is provided if \code{NULL} is passed.
#'
#' @seealso
#' \code{\link{mmcif_fit}}, \code{\link{mmcif_start_values}} and
#' \code{\link{mmcif_sandwich}}.
#'
#' @importFrom stats model.frame model.matrix terms ave
#' @export
mmcif_data <- function(formula, data, cause, time, cluster_id, max_time,
                       spline_df = 3L, left_trunc = NULL, ghq_data = NULL){
  stopifnot(inherits(formula, "formula"),
            is.data.frame(data),
            length(spline_df) == 1, is.finite(spline_df), spline_df > 0)

  # setup the design matrix
  mf <- model.frame(formula, data)
  covs_risk <- model.matrix(terms(mf), mf)

  cause <- eval(substitute(cause), data, parent.frame())
  time_observed <- eval(substitute(time), data, parent.frame())
  cluster_id <- eval(substitute(cluster_id), data, parent.frame())

  left_trunc <- eval(substitute(left_trunc), data, parent.frame())
  if(is.null(left_trunc))
    left_trunc <- numeric(length(time_observed))

  n_causes <- length(unique(cause)) - 1L

  ghq_data <- .ghq_data_default_if_null(ghq_data)

  stopifnot(
    length(cause) > 0L,
    NROW(covs_risk) == length(cause),
    length(time_observed) == length(cause),
    all(is.finite(time_observed) & time_observed > 0),
    length(cluster_id) == length(cluster_id),
    n_causes > 0L, all(1:(n_causes + 1L) %in% cause),
    is.numeric(max_time), length(max_time) == 1, is.finite(max_time),
    max_time > 0,
    max(time_observed[cause %in% 1:n_causes]) < max_time,
    length(left_trunc) == length(cause),
    all(is.finite(left_trunc) & left_trunc >= 0),
    all(time_observed > left_trunc),
    is.list(ghq_data), all(c("node", "weight") %in% names(ghq_data)),
    length(ghq_data$node) == length(ghq_data$weight),
    length(ghq_data$node) > 0,
    all(sapply(ghq_data, function(x) all(is.finite(x)))))

  # add the time transformations
  time_trans <- function(x) atanh((x - max_time / 2) / (max_time / 2))
  d_time_trans <- function(x) {
    outer <- (x - max_time / 2) / (max_time / 2)
    1/((1 - outer^2) * (max_time / 2))
  }

  time_observed_trans <- suppressWarnings(
    ifelse(time_observed < max_time, time_trans(time_observed), NA))

  is_observed <- cause %in% 1:n_causes
  splines <- tapply(
    time_observed_trans[is_observed], cause[is_observed],
    function(time_observed_trans)
      monotone_ns(time_observed_trans, df = spline_df), simplify = FALSE)

  time_expansion <- function(x, cause){
    x <- suppressWarnings(time_trans(x))
    splines[[cause]]$expansion(x, ders = 0L)
  }

  d_time_expansion <- function(x, cause){
    z <- suppressWarnings(time_trans(x))
    d_x <- suppressWarnings(d_time_trans(x))
    splines[[cause]]$expansion(z, ders = 1L) * d_x
  }

  eval_trajectory_covs <- function(x, FUN, covs){
    varying_covs <- lapply(1:n_causes, FUN, x = x)
    do.call(cbind, lapply(varying_covs, cbind, covs))
  }

  covs_trajectory <- eval_trajectory_covs(
    time_observed, time_expansion, covs_risk)

  d_covs_trajectory <- eval_trajectory_covs(
    time_observed, d_time_expansion,
    matrix(0, NROW(covs_risk), NCOL(covs_risk)))

  has_finite_trajectory_prob <- time_observed < max_time

  # compute the covariates for the left truncation
  covs_trajectory_delayed <-
    matrix(NaN, NROW(covs_trajectory), NCOL(covs_trajectory))
  has_delayed_entry <- which(left_trunc > 0)

  covs_trajectory_delayed[has_delayed_entry, ] <-
    eval_trajectory_covs(
      left_trunc[has_delayed_entry], time_expansion,
      covs_risk[has_delayed_entry, ])

  # find all permutation of indices in each cluster
  cluster_length <- ave(cluster_id, cluster_id, FUN = length)

  row_id <- seq_len(length(cluster_id))
  pair_indices_res <- create_pair_indices(
    cluster_id[cluster_length > 1], row_id[cluster_length > 1])
  pair_indices <- pair_indices_res$pair_indices
  pair_cluster_id <- pair_indices_res$pair_cluster_id

  singletons <- row_id[cluster_length < 2]

  # create the C++ object
  # TODO: this function should be exposed to advanced users
  comp_obj <- mmcif_data_holder(
    covs_risk = t(covs_risk), covs_trajectory = t(covs_trajectory),
    d_covs_trajectory = t(d_covs_trajectory),
    has_finite_trajectory_prob = has_finite_trajectory_prob,
    cause = cause - 1L, n_causes = n_causes, pair_indices = pair_indices - 1L,
    singletons = singletons - 1L,
    covs_trajectory_delayed = t(covs_trajectory_delayed),
    pair_cluster_id = pair_cluster_id)

  # create an object to do index the parameters
  n_coef_risk <- NCOL(covs_risk) * n_causes
  n_coef_trajectory <- (spline_df + NCOL(covs_risk)) * n_causes

  coef_trajectory_time <- matrix(
    seq_len(n_coef_trajectory) + n_coef_risk, spline_df + NCOL(covs_risk))
  coef_trajectory_time <- coef_trajectory_time[1:spline_df, ]

  idx_coef_trajectory <- seq_len(n_coef_trajectory) + n_coef_risk

  vcov_dim <- 2L * n_causes
  vcov_full <- max(idx_coef_trajectory) + seq_len(vcov_dim^2)
  vcov_upper <- max(idx_coef_trajectory) +
    seq_len((vcov_dim * (vcov_dim + 1L)) %/% 2L)

  indices <- list(
    coef_risk = seq_len(n_coef_risk),
    coef_trajectory = idx_coef_trajectory,
    coef_trajectory_time = c(coef_trajectory_time),
    vcov_full = vcov_full,
    vcov_upper = vcov_upper,
    n_causes = n_causes,
    n_par = max(vcov_full),
    n_par_upper_tri = max(vcov_upper),
    n_par_wo_vcov = max(idx_coef_trajectory))

  # create the constrain matrices
  n_constraints_per_spline <- NROW(splines[[1]]$constraints)
  constraints <- matrix(0., n_constraints_per_spline * n_causes,
                        indices$n_par_wo_vcov)
  for(i in 1:n_causes)
    constraints[
      1:n_constraints_per_spline + (i - 1L) * n_constraints_per_spline,
      matrix(indices$coef_trajectory_time, spline_df)[, i]] <-
      splines[[i]]$constraints

  constraints <- list(
    wo_vcov = constraints,
    vcov_full =
      cbind(constraints, matrix(0., NROW(constraints), vcov_dim^2)),
    vcov_upper =
      cbind(constraints, matrix(0., NROW(constraints),
                                (vcov_dim * (vcov_dim + 1L)) %/% 2L)))

  # clean up
  rm(list = setdiff(
    ls(),
    c("comp_obj", "pair_indices", "singletons", "covs_risk",
      "covs_trajectory", "time_observed", "cause", "time_trans",
      "d_time_trans", "max_time", "indices", "splines", "d_covs_trajectory",
      "constraints", "covs_trajectory_delayed", "time_expansion",
      "d_time_expansion", "pair_cluster_id", "ghq_data")))

  structure(
    list(comp_obj = comp_obj, pair_indices = pair_indices,
         singletons = singletons, covs_risk = covs_risk,
         covs_trajectory = covs_trajectory, time_observed = time_observed,
         cause = cause, time_trans = time_trans, d_time_trans = d_time_trans,
         time_expansion = time_expansion, d_time_expansion = d_time_expansion,
         max_time = max_time, indices = indices, splines = splines,
         d_covs_trajectory = d_covs_trajectory, constraints = constraints,
         covs_trajectory_delayed = covs_trajectory_delayed,
         pair_cluster_id = pair_cluster_id, ghq_data = ghq_data),
    class = "mmcif")
}

.check_n_threads <- function(n_threads){
  stopifnot(length(n_threads) == 1, is.finite(n_threads), n_threads > 0)
}

#' Evaluates the Log Composite Likelihood and its Gradient
#'
#' @param object an object from \code{\link{mmcif_data}}.
#' @param n_threads the number of threads to use.
#' @param par numeric vector with parameters. This is either using a log
#' Cholesky decomposition for the covariance matrix or the covariance matrix.
#' @param ghq_data the Gauss-Hermite quadrature nodes and weights to
#' use. It should be stored as a list with two elements called \code{"node"}
#' and \code{"weight"}. A default is provided if \code{NULL} is passed.
#' @param is_log_chol logical for whether a log Cholesky decomposition is used
#' for the covariance matrix or the full covariance matrix.
#'
#'
#' @importFrom utils tail
#' @export
mmcif_logLik <- function(
  object, par, ghq_data = object$ghq_data, n_threads = 1L,
  is_log_chol = FALSE){
  stopifnot(inherits(object, "mmcif"))
  .check_n_threads(n_threads)
  .check_is_log_chol(is_log_chol)

  if(is_log_chol){
    n_causes <- object$indices$n_causes
    n_vcov <- (2L * n_causes * (2L * n_causes + 1L)) %/% 2L
    par <- c(head(par, -n_vcov), log_chol_inv(tail(par, n_vcov)))
  }

  mmcif_logLik_cpp(
    object$comp_obj, par = par, ghq_data = ghq_data, n_threads = n_threads)
}

#' @rdname mmcif_logLik
#' @export
mmcif_logLik_grad <- function(
  object, par, ghq_data = object$ghq_data, n_threads = 1L,
  is_log_chol = FALSE){
  stopifnot(inherits(object, "mmcif"))
  .check_n_threads(n_threads)
  .check_is_log_chol(is_log_chol)


  if(is_log_chol){
    par_org <- par

    n_causes <- object$indices$n_causes
    n_vcov <- (2L * n_causes * (2L * n_causes + 1L)) %/% 2L
    vcov <- log_chol_inv(tail(par, n_vcov))
    par <- c(head(par, -n_vcov), vcov)
  }

  gr <- mmcif_logLik_grad_cpp(
    object$comp_obj, par = par, ghq_data = ghq_data, n_threads = n_threads)
  if(!is_log_chol)
    return(gr)

  # back propagate the gradients w.r.t. the random effects
  d_vcov <- matrix(tail(gr, 4L * n_causes * n_causes), 2L * n_causes)
  C <- chol(vcov)
  d_vcov <- 2 * C %*% d_vcov
  diag(d_vcov) <- diag(d_vcov) * diag(C)

  c(head(gr, -4L * n_causes * n_causes), d_vcov[upper.tri(d_vcov, TRUE)])
}

#' Finds Staring Values
#'
#' @inheritParams mmcif_logLik
#' @param vcov_start starting value for the covariance matrix of the random
#' effects. \code{NULL} yields the identity matrix.
#'
#' @importFrom stats lm.fit constrOptim na.omit sd
#' @export
mmcif_start_values <- function(object, n_threads = 1L, vcov_start = NULL){
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
  comb_slope <- sapply(object$splines, function(spline){
    boundary_knots <- spline$boundary_knots
    pts <- seq(boundary_knots[1], boundary_knots[2], length.out = 1000)
    lm.fit(cbind(1, spline$expansion(pts)), pts)$coef
  })

  n_causes <- object$indices$n_causes
  n_cov_traject <- NCOL(object$covs_trajectory) %/% n_causes
  traject_factor <- gl(2, n_cov_traject)

  comb_line <-
    tapply(1:NCOL(object$covs_trajectory), traject_factor, function(indices){
      design_mat <- na.omit(object$covs_trajectory[, indices])
      lm.fit(design_mat, rep(1, NROW(design_mat)))$coef
    }, simplify = FALSE)
  comb_line <- do.call(cbind, comb_line)

  # start the optimization. Start only with the trajectories
  par <- numeric(object$indices$n_par_wo_vcov)
  coef_trajectory_time <- object$indices$coef_trajectory_time
  coef_trajectory <- object$indices$coef_trajectory

  par[coef_trajectory_time] <-
    comb_slope[-1, ] * rep(sapply(slope_n_time, `[[`, "slope"),
                           each = NROW(comb_slope) - 1)
  par[coef_trajectory] <- par[coef_trajectory] +
    c(comb_line *
        rep(sapply(slope_n_time, `[[`, "intercept") +
              comb_slope[1] * sapply(slope_n_time, `[[`, "slope"),
            each = NROW(comb_line)))

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
                 upper = c(opt$par, log_chol(vcov_start))),
            logLik = -opt$value)
}


#' Fits a Mixed Competing Risk Model
#'
#' @inheritParams mmcif_logLik
#' @param par numeric vector with parameters. This is using a log
#' Cholesky decomposition for the covariance matrix.
#' @param control,method,... arguments passed to \code{\link{constrOptim}}.
#'
#' @seealso
#' \code{\link{mmcif_data}}, \code{\link{mmcif_start_values}} and
#' \code{\link{mmcif_sandwich}}.
#'
#' @importFrom stats constrOptim
#' @export
mmcif_fit <- function(
  par, object, n_threads = 1L, control = list(maxit = 10000L), method = "BFGS",
  ...){
  stopifnot(inherits(object, "mmcif"))
  constraints <- object$constraints$vcov_upper

  constrOptim(
    par, function(par) -mmcif_logLik(
      object, par, n_threads = n_threads, is_log_chol = TRUE),
    grad = function(par) -mmcif_logLik_grad(
      object, par, n_threads = n_threads, is_log_chol = TRUE),
    method = method, ui = constraints,
    ci = rep(1e-8, NROW(constraints)),
    control = control, ...)
}

#' Computes the Sandwich Estimator
#'
#' Computes the sandwich estimator of the covariance matrix. The parameter that
#' is passed is using the log Cholesky decomposition. The Hessian is computed
#' using numerical differentiation with Richardson extrapolation to refine the
#' estimate.
#'
#' @inheritParams mmcif_logLik
#' @param par numeric vector with the parameters to compute the sandwich
#' estimator at.
#' @param eps determines the step size in the numerical differentiation using
#'  \code{max(eps^2, |par[i]| * eps)}.
#' @param scale scaling factor in the Richardson extrapolation.
#' @param tol relative convergence criteria in the extrapolation given
#' by \code{max(tol, |g[j]| * tol)} with \code{g} being the gradient.
#' @param order maximum number of iteration of the Richardson extrapolation.
#'
#' @seealso
#' \code{\link{mmcif_fit}} and \code{\link{mmcif_data}}.
#'
#' @export
mmcif_sandwich <- function(
  object, par, ghq_data = object$ghq_data, n_threads = 1L, eps = .01,
  scale = 2., tol = .00000001, order = 6L){
  stopifnot(inherits(object, "mmcif"))
  .check_n_threads(n_threads)
  ghq_data <- .ghq_data_default_if_null(ghq_data)

  cpp_res <- mmcif_sandwich_cpp(
    data_ptr = object$comp_obj, par = par, ghq_data = ghq_data,
    n_threads = n_threads, eps = eps, scale = scale, tol = tol, order = order)

  hessian <- cpp_res$hessian
  hessian[lower.tri(hessian)] <- t(hessian)[lower.tri(hessian)]
  meat <- cpp_res$meat
  res <- solve(hessian, t(solve(hessian, meat)))

  # construct the covariance matrix estimate on the original covariance matrix
  # scale
  jac <- diag(length(par))
  jac[object$indices$vcov_upper, object$indices$vcov_upper] <-
    jac_log_chol_inv(par[object$indices$vcov_upper])
  res_org <- tcrossprod(jac %*% res, jac)

  structure(res, meat = meat, hessian = hessian, `res vcov` = res_org)
}
