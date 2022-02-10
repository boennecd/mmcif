#ifndef MCIF_LOGLIK_H
#define MCIF_LOGLIK_H

#include "mmcif-logLik.h"

/// computes the log likelihood of the model without the random effects.
template<bool with_risk>
double mcif_logLik
  (double const * __restrict__ par, param_indexer const &indexer,
   mmcif_data const &obs, ghqCpp::simple_mem_stack<double> &mem);

/// computes the gradient of mcif_logLik
template<bool with_risk>
double mcif_logLik_grad
  (double const * __restrict__ par, double * __restrict__ grad,
   param_indexer const &indexer, mmcif_data const &obs,
   ghqCpp::simple_mem_stack<double> &mem);

#endif
