% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mmcif.R
\name{mmcif_start_values}
\alias{mmcif_start_values}
\title{Finds Staring Values}
\usage{
mmcif_start_values(object, n_threads = 1L, vcov_start = NULL)
}
\arguments{
\item{object}{an object from \code{\link{mmcif_data}}.}

\item{n_threads}{the number of threads to use.}

\item{vcov_start}{starting value for the covariance matrix of the random
effects. \code{NULL} yields the identity matrix.}
}
\value{
A list with
\itemize{
  \item an element called \code{"full"} with the starting value where the
  last components are the covariance matrix.
  \item an element called \code{"upper"} the staring values where the
  covariance matrix is stored as a log Cholesky decomposition. This is used
  e.g. for optimization with \code{\link{mmcif_fit}}.
}
}
\description{
Fast heuristic for finding starting values for the mixed cumulative incidence
functions model.
}
\examples{
if(require(mets)){
  # prepare the data
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
    prt_use, id \%in\% sample(unique(id), length(unique(id)) \%/\% 10L))

  n_threads <- 2L
  mmcif_obj <- mmcif_data(
    ~ country - 1, prt_use, status, time, id, max_time,
    2L, strata = country)

  # get the staring values
  start_vals <- mmcif_start_values(mmcif_obj, n_threads = n_threads)

  # the starting values
  print(start_vals)
}

}
