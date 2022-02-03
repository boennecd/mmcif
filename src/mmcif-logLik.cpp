#include "mmcif-logLik.h"
#include "integrand-mixed-mult-logit-term.h"
#include "integrand-probit-term.h"
#include "integrand-cond-pbvn.h"

namespace {
bool is_censored(mmcif_data const &obs, param_indexer const &indexer){
  return obs.cause == indexer.n_causes();
}

/// computes the log likelihood when both outcomes are observed
double mmcif_log_Lik_both_obs
  (double const * __restrict__ par, param_indexer const &indexer,
   mmcif_data const &obs1, mmcif_data const &obs2,
   ghqCpp::simple_mem_stack<double> &mem, ghqCpp::ghq_data const &dat){
  return 0; // TODO: implement
}

/**
 * computes the log likelihood when the first outcome is observed but the
 * second is not
 */
double mmcif_log_Lik_one_obs
  (double const * __restrict__ par, param_indexer const &indexer,
   mmcif_data const &obs1, mmcif_data const &obs2,
   ghqCpp::simple_mem_stack<double> &mem, ghqCpp::ghq_data const &dat){
  return 0; // TODO: implement
}

/// computes the log likelihood when both outcomes are censored
double mmcif_log_Lik_both_cens
  (double const * __restrict__ par, param_indexer const &indexer,
   mmcif_data const &obs1, mmcif_data const &obs2,
   ghqCpp::simple_mem_stack<double> &mem, ghqCpp::ghq_data const &dat){
  return 0; // TODO: implement
}

} // namespace

double mmcif_log_Lik
  (double const * __restrict__ par, param_indexer const &indexer,
   mmcif_data const &obs, ghqCpp::simple_mem_stack<double> &mem,
   ghqCpp::ghq_data const &dat){
  bool const is_cens{is_censored(obs, indexer)};
  // TODO: implement
  if(is_cens)
    return 0;
  return 0;
}

double mmcif_log_Lik
  (double const * __restrict__ par, param_indexer const &indexer,
   mmcif_data const &obs1, mmcif_data const &obs2,
   ghqCpp::simple_mem_stack<double> &mem, ghqCpp::ghq_data const &dat){
  bool const is_cens1{is_censored(obs1, indexer)},
             is_cens2{is_censored(obs2, indexer)};

  if(is_cens1 && is_cens2)
    return mmcif_log_Lik_both_cens(par, indexer, obs1, obs2, mem, dat);
  else if(!is_cens1 && is_cens2)
    return mmcif_log_Lik_one_obs(par, indexer, obs1, obs2, mem, dat);
  else if(is_cens1 && !is_cens2)
    return mmcif_log_Lik_one_obs(par, indexer, obs2, obs1, mem, dat);
  return mmcif_log_Lik_both_obs(par, indexer, obs1, obs2, mem, dat);
}
