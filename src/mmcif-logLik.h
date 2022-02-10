#ifndef MMCIF_LOGLIK_H
#define MMCIF_LOGLIK_H

#include "ghq.h"
#include "simple-mem-stack.h"
#include "param-indexer.h"

/// information about a given individual
struct mmcif_data {
  /// covariates for the trajectory
  double const * cov_trajectory;
  /// derivatives of x w.r.t. time
  double const * d_cov_trajectory;
  /// covariates for the risks
  double const * cov_risk;
  /// is the probability of the trajectory in (0, 1) (finite) or is it always 1
  bool has_finite_trajectory_prob;
  /**
   * the flag for the cause (zero-based). Censoring is indicated by letting this
   * be the number of causes.
   */
  unsigned cause;
};

/**
 * computes the log composite likelihood term of a given pair. The parameters
 * should be stored with the full covariance matrix.
 */
double mmcif_logLik
  (double const * __restrict__ par, param_indexer const &indexer,
   mmcif_data const &obs1, mmcif_data const &obs2,
   ghqCpp::simple_mem_stack<double> &mem, ghqCpp::ghq_data const &dat);

/**
 * overload of singletons.
 */
double mmcif_logLik
  (double const * __restrict__ par, param_indexer const &indexer,
   mmcif_data const &obs, ghqCpp::simple_mem_stack<double> &mem,
   ghqCpp::ghq_data const &dat);

#endif
