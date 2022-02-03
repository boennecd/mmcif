#include <RcppArmadillo.h>
#include "simple-mat.h"
#include <vector>
#include "param-indexer.h"
#include "mmcif-logLik.h"
#include "ghq.h"
#include "wmem.h"

using Rcpp::IntegerVector;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;

namespace {
simple_mat<double> NumericMatrix_to_simple_mat(NumericMatrix const from){
  simple_mat<double> out(from.nrow(), from.ncol());
  std::copy(from.begin(), from.end(), out.begin());
  return out;
}

/// holds the data to compute the composite likelihood
struct mmcif_data_holder {
  // all matrices are stored as [dimension per obs] x [# Observations]

  /// covariates for the trajectory
  simple_mat<double> covs_trajectory;
  /// derivatives of covs_trajectory w.r.t. time
  simple_mat<double> d_covs_trajectory;
  /// covariates for the risks
  simple_mat<double> covs_risk;
  /// is the probability of the trajectory in (0, 1) (finite) or is it always 1
  std::vector<char> has_finite_trajectory_prob;
  /**
   * the flag for the cause (zero-based). Censoring is indicated by letting this
   * be the number of causes.
   */
  std::vector<unsigned> cause;
  /// 2 x [# pairs] matrix with zero-based indices for the pairs
  simple_mat<size_t> pair_indices;
  /// indices for singletons
  std::vector<size_t> singletons;
  /// the indexer for the parameters
  param_indexer indexer;

  mmcif_data_holder
    (NumericMatrix const covs_trajectory_in,
     NumericMatrix const d_covs_trajectory_in, NumericMatrix const covs_risk_in,
     IntegerVector const has_finite_trajectory_prob_in,
     IntegerVector const cause_in, size_t const n_causes,
     Rcpp::IntegerMatrix pair_indices_in, IntegerVector const singletons_in):
    covs_trajectory{NumericMatrix_to_simple_mat(covs_trajectory_in)},
    d_covs_trajectory{NumericMatrix_to_simple_mat(d_covs_trajectory_in)},
    covs_risk{NumericMatrix_to_simple_mat(covs_risk_in)},
    has_finite_trajectory_prob
    {
      ([&]{
        std::vector<char> out;
        out.reserve(has_finite_trajectory_prob_in.size());
        for(int val : has_finite_trajectory_prob_in)
          out.emplace_back(val == 1);
        return out;
      })()
    },
    cause
    {
      ([&]{
        std::vector<unsigned> out;
        out.reserve(cause_in.size());
        for(int i : cause_in)
          if(i < 0)
            throw std::invalid_argument("cause has values less than zero");
          else
            out.emplace_back(i);
        return out;
      })()
    },
    pair_indices
    {
      ([&]{
        if(pair_indices_in.nrow() != 2)
          throw std::invalid_argument("pair_indices.nrow() != 2");
        simple_mat<size_t> out(2, pair_indices_in.ncol());

        auto to = out.begin();
        for(auto from = pair_indices_in.begin(); from != pair_indices_in.end();
            ++from, ++to){
          if(*from < 0)
            throw std::invalid_argument("pair_indices has indices less than zero");
          else
            *to = *from;
        }

        return out;
      })()
    },
    singletons
    {
      ([&]{
        std::vector<size_t> out;
        out.reserve(singletons_in.size());
        for(int idx : singletons_in)
          if(idx < 0)
            throw std::invalid_argument("singletons has indices less than zero");
          else
            out.emplace_back(idx);
        return out;
      })()
    },
    indexer{covs_risk.n_rows(), covs_trajectory.n_rows(), n_causes}
    {
      throw_if_invalidt();
    }

private:
  void throw_if_invalidt(){
    if(covs_risk.n_rows() != indexer.n_cov_risk())
      throw std::invalid_argument("covs_risk.n_rows() != indexer.n_cov_risk()");
    if(covs_trajectory.n_rows() != indexer.n_cov_traject())
      throw std::invalid_argument("covs_trajectory.n_rows() != indexer.n_cov_traject()");
    if(d_covs_trajectory.n_rows() != indexer.n_cov_traject())
      throw std::invalid_argument("d_covs_trajectory.n_rows() != indexer.n_cov_traject()");

    size_t const n_obs{covs_risk.n_cols()};
    if(covs_trajectory.n_cols() != n_obs)
      throw std::invalid_argument("covs_trajectory.n_cols() != n_obs");
    if(d_covs_trajectory.n_cols() != n_obs)
      throw std::invalid_argument("d_covs_trajectory.n_cols() != n_obs");
    if(has_finite_trajectory_prob.size() != n_obs)
      throw std::invalid_argument("has_finite_trajectory_prob.size() != n_obs");
    if(cause.size() != n_obs)
      throw std::invalid_argument("cause.size() != n_obs");

    for(size_t idx : pair_indices)
      if(idx >= n_obs)
        throw std::invalid_argument("pair_indices has an index that is out of bounds");
    for(size_t idx : singletons)
      if(idx >= n_obs)
        throw std::invalid_argument("singletons has an index that is out of bounds");

    std::vector<int> any_with_cause(indexer.n_causes(), 0);
    auto finite_prob = has_finite_trajectory_prob.begin();
    for(unsigned cause_i : cause){
      if(cause_i > indexer.n_causes())
        throw std::invalid_argument("there is cause with a value greater than the number of causes");
      else if(cause_i < indexer.n_causes())
        any_with_cause[cause_i] = 1;

      if(*finite_prob++ == 0 && cause_i < indexer.n_causes())
        throw std::runtime_error("there is an observed outcome with a probability of the trajectory which is one");
    }

    for(int any_with_cause_i : any_with_cause)
      if(any_with_cause_i < 1)
        throw std::invalid_argument("there is cause with no observations");
  }
};

ghqCpp::ghq_data ghq_data_from_list(Rcpp::List dat){
  NumericVector nodes = dat["node"],
                    weigths = dat["weight"];
  if(nodes.size() != weigths.size())
    throw std::runtime_error("nodes.size() != weigths.size()");

  return { &nodes[0], &weigths[0], static_cast<size_t>(nodes.size()) };
}

void throw_if_invalid_par
  (mmcif_data_holder const &data, NumericVector const par){
  if(static_cast<size_t>(par.size()) != data.indexer.n_par<false>())
    throw std::invalid_argument("invalid length of parameter vector");
}

mmcif_data mmcif_data_from_idx
  (mmcif_data_holder const &data, size_t const idx){
  return { data.covs_trajectory.col(idx), data.d_covs_trajectory.col(idx),
           data.covs_risk.col(idx), data.has_finite_trajectory_prob[idx] == 1,
           data.cause[idx] };
}

} // namespace


// [[Rcpp::export("mmcif_data_holder", rng = false)]]
SEXP mmcif_data_holder_to_R
  (NumericMatrix const covs_trajectoryn,
   NumericMatrix const d_covs_trajectoryn,
   NumericMatrix const covs_risk,
   IntegerVector const has_finite_trajectory_prob,
   IntegerVector const cause, size_t const n_causes,
   Rcpp::IntegerMatrix pair_indices, IntegerVector const singletons_in){
  return Rcpp::XPtr<mmcif_data_holder const>
    (new mmcif_data_holder
       (covs_trajectoryn, d_covs_trajectoryn, covs_risk,
        has_finite_trajectory_prob, cause, n_causes, pair_indices,
        singletons_in));
}

// [[Rcpp::export("mmcif_logLik", rng = false)]]
double mmcif_logLik_to_R
  (SEXP data_ptr, NumericVector const par, Rcpp::List ghq_data,
   unsigned const n_threads){

  Rcpp::XPtr<mmcif_data_holder const> data(data_ptr);
  throw_if_invalid_par(*data, par);
  wmem::setup_working_memory(n_threads);
  auto ghq_data_pass = ghq_data_from_list(ghq_data);

  double out{};
  size_t const n_pairs{data->pair_indices.n_cols()};
  double const * const par_ptr{&par[0]};

#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads) reduction(+:out)
#endif
  for(size_t i = 0; i < n_pairs; ++i){
    auto mmcif_dat1 = mmcif_data_from_idx(*data, data->pair_indices.col(i)[0]);
    auto mmcif_dat2 = mmcif_data_from_idx(*data, data->pair_indices.col(i)[1]);

    wmem::mem_stack().reset();
    out += mmcif_log_Lik(par_ptr, data->indexer, mmcif_dat1, mmcif_dat2,
                         wmem::mem_stack(), ghq_data_pass);
  }

  size_t const n_singletons{data->singletons.size()};
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads) reduction(+:out)
#endif
  for(size_t i = 0; i < n_singletons; ++i){
    auto mmcif_dat = mmcif_data_from_idx(*data, data->singletons[i]);

    wmem::mem_stack().reset();
    out += mmcif_log_Lik
      (par_ptr,  data->indexer, mmcif_dat, wmem::mem_stack(), ghq_data_pass);
  }

  return out;
}

