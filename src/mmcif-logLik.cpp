#include "mmcif-logLik.h"
#include "integrand-mixed-mult-logit-term.h"
#include "integrand-probit-term.h"
#include "integrand-cond-pbvn.h"
#include <numeric>
#include <Rmath.h> // Rf_dnorm4
#include "mmcif-misc.h"

using namespace ghqCpp;

namespace {
constexpr size_t ghq_target_size{100};
constexpr double adaptive_rel_eps{1e-6};

class mmcif_comp_helper {
  param_indexer const &indexer;
  double const * const par;
  simple_mem_stack<double> &mem;

public:
  mmcif_comp_helper(param_indexer const &indexer, double const *par,
                    simple_mem_stack<double> &mem):
    indexer{indexer}, par{par}, mem{mem} { }

  bool is_censored(mmcif_data const &obs){
    return obs.cause == indexer.n_causes();
  }

  double comp_lp_traject
  (mmcif_data const &obs, unsigned const cause){
    return -std::inner_product
    (obs.cov_trajectory, obs.cov_trajectory + indexer.n_cov_traject(),
     par + indexer.traject(cause), 0.);
  }

  double comp_lp_traject(mmcif_data const &obs){
    return comp_lp_traject(obs, obs.cause);
  }

  double comp_d_lp_traject
  (mmcif_data const &obs, unsigned const cause){
    return -std::inner_product
    (obs.d_cov_trajectory, obs.d_cov_trajectory + indexer.n_cov_traject(),
     par + indexer.traject(cause), 0.);
  }

  double comp_d_lp_traject(mmcif_data const &obs){
    return comp_d_lp_traject(obs, obs.cause);
  }

  void fill_vcov(arma::mat &vcov){
    auto vcov_dim = 2 * indexer.n_causes();
    vcov = mat_no_alloc(vcov_dim, vcov_dim, mem);
    std::copy(par + indexer.vcov(), par + indexer.n_par<false>(),
              vcov.begin());
  }

  double comp_trajector_cond_dens_obs_one
  (double const lp_traject, unsigned const cause){
    auto const n_causes = indexer.n_causes();
    double const *vcov{par + indexer.vcov()};
    auto const idx = cause + n_causes;
    double var{1 + vcov[idx + idx * 2 * n_causes]};
    return Rf_dnorm4(lp_traject, 0, std::sqrt(var), 1);
  }

  void fill_logit_offsets
  (double *eta, mmcif_data const &obs){
    for(size_t j = 0; j < indexer.n_causes(); ++j)
      eta[j] = std::inner_product
        (obs.cov_risk, obs.cov_risk + indexer.n_cov_risk(),
         par + indexer.risk(j), 0.);
  }

  void fill_cond_vcov_one_obs(arma::mat &res, unsigned const cause){
    auto const n_causes = indexer.n_causes();
    res = mat_no_alloc(2 * n_causes, 2 * n_causes, mem);
    std::copy(par + indexer.vcov(), par + indexer.n_par<false>(),
              res.begin());

    arma::mat Sigma_inv{mat_no_alloc(n_causes, n_causes, mem)};
    Sigma_inv = arma::inv_sympd(res);

    auto const idx = cause + n_causes;
    Sigma_inv(idx, idx) += 1;
    res = arma::inv_sympd(Sigma_inv);
  }
};

/// computes the log likelihood when both outcomes are observed
double mmcif_log_Lik_both_obs
  (double const * __restrict__ par, param_indexer const &indexer,
   mmcif_data const &obs1, mmcif_data const &obs2,
   simple_mem_stack<double> &mem, ghq_data const &dat){
  mmcif_comp_helper helper{indexer, par, mem};

  double lp_traject[]{helper.comp_lp_traject(obs1),
                      helper.comp_lp_traject(obs2)};
  double d_lp_traject[]{helper.comp_d_lp_traject(obs1),
                        helper.comp_d_lp_traject(obs2)};

  double out{std::log(d_lp_traject[0]) + std::log(d_lp_traject[1])};

  auto const n_causes = indexer.n_causes();
  auto const idx_cause_1 = obs1.cause + n_causes;
  auto const idx_cause_2 = obs2.cause + n_causes;

  arma::mat vcov_cond_dens{mat_no_alloc(2, 2, mem)};
  arma::mat V{mat_no_alloc(2, 2 * n_causes, mem)};
  V.zeros();
  V(0, idx_cause_1) = 1;
  V(0, idx_cause_2) = 1;

  arma::mat vcov;
  helper.fill_vcov(vcov);

  arma::vec const dens_point{lp_traject[0], lp_traject[1]};

  {
    vcov_cond_dens = V * vcov * V.t();
    vcov_cond_dens.diag() += 1;
    out += log_dmvn(dens_point, vcov_cond_dens, mem);
  }

  arma::mat logit_offsets{mat_no_alloc(n_causes, 2, mem)};
  helper.fill_logit_offsets(logit_offsets.colptr(0), obs1);
  helper.fill_logit_offsets(logit_offsets.colptr(1), obs2);

  arma::mat vcov_cond_inv{mat_no_alloc(2 * n_causes, 2 * n_causes, mem)};
  arma::mat vcov_cond    {mat_no_alloc(2 * n_causes, 2 * n_causes, mem)};
  vcov_cond_inv = arma::inv_sympd(vcov);
  vcov_cond_inv += V * V.t();
  vcov_cond = arma::inv_sympd(vcov_cond_inv);

  arma::mat vcov_cond_sub{mat_no_alloc(n_causes, n_causes, mem)};
  vcov_cond_sub = vcov_cond.submat(0, 0, n_causes - 1, n_causes - 1);

  arma::vec cond_mean = vcov_cond * V.t() * dens_point;
  cond_mean = cond_mean.subvec(0, n_causes - 1);

  arma::uvec const which_cat{obs1.cause + 1, obs2.cause + 1};

  auto mem_mark = mem.set_mark_raii();
  mixed_mult_logit_term<false> mult_term(logit_offsets, which_cat);
  rescale_shift_problem<false>
    prob_use(vcov_cond_sub, cond_mean, mult_term);

  adaptive_problem prob(prob_use, mem, adaptive_rel_eps);
  double res{};
  ghq(&res, dat, prob, mem, ghq_target_size);
  out += std::log(res);

  return out;
}

/**
 * computes the log likelihood when the first outcome is observed but the
 * second is not
 */
double mmcif_log_Lik_one_obs
  (double const * __restrict__ par, param_indexer const &indexer,
   mmcif_data const &obs1, mmcif_data const &obs2,
   simple_mem_stack<double> &mem, ghq_data const &dat){
  return 0; // TODO: implement
}

/// computes the log likelihood when both outcomes are censored
double mmcif_log_Lik_both_cens
  (double const * __restrict__ par, param_indexer const &indexer,
   mmcif_data const &obs1, mmcif_data const &obs2,
   simple_mem_stack<double> &mem, ghq_data const &dat){
  return 0; // TODO: implement
}



} // namespace

double mmcif_log_Lik
  (double const * __restrict__ par, param_indexer const &indexer,
   mmcif_data const &obs, simple_mem_stack<double> &mem,
   ghq_data const &dat){
  mmcif_comp_helper helper{indexer, par, mem};

  bool const is_cens{helper.is_censored(obs)};
  if(is_cens){
    auto const n_causes = indexer.n_causes();
    arma::mat logit_offsets{mat_no_alloc(n_causes, 1, mem)};
    helper.fill_logit_offsets(logit_offsets.begin(), obs);
    arma::uvec which_cat(1);

    arma::mat vcov;
    helper.fill_vcov(vcov);
    arma::mat vcov_sub{mat_no_alloc(n_causes, n_causes, mem)};
    vcov_sub = vcov.submat(0, 0, n_causes - 1, n_causes - 1);

    if(!obs.has_finite_trajectory_prob){
      which_cat[0] = 0;

      auto mem_mark = mem.set_mark_raii();
      mixed_mult_logit_term<false> mult_term(logit_offsets, which_cat);
      rescale_problem<false> prob_use(vcov_sub, mult_term);
      adaptive_problem prob(prob_use, mem, adaptive_rel_eps);
      double res{};
      ghq(&res, dat, prob, mem, ghq_target_size);
      return std::log(res);
    }

    auto mem_mark = mem.set_mark_raii();
    double integrand{1};
    for(size_t cause = 0; cause < indexer.n_causes(); ++cause){
      // TODO: allocations
      size_t const idx_traject{n_causes + cause};
      arma::vec vcov_col = vcov.col(idx_traject).subvec(0, n_causes - 1);
      arma::vec rng_coefs = arma::solve(vcov_sub, vcov_col,
                                        arma::solve_opts::likely_sympd);
      double const var_cond
        {1 + vcov(idx_traject, idx_traject) - arma::dot(vcov_col, rng_coefs)};
      rng_coefs *= -1;
      double const lp_traject{helper.comp_lp_traject(obs, cause)};

      mixed_probit_term<false> probit_term
        (std::sqrt(var_cond), lp_traject, rng_coefs);
      which_cat[0] = cause + 1;
      mixed_mult_logit_term<false> mult_term(logit_offsets, which_cat);

      combined_problem comp_prob({&probit_term, &mult_term});
      rescale_problem<false> prob_use(vcov_sub, comp_prob);
      adaptive_problem prob(prob_use, mem, adaptive_rel_eps);

      double res{};
      ghq(&res, dat, prob, mem, ghq_target_size);

      integrand -= res;
    }

    return std::log(integrand);
  }

  size_t const cause{obs.cause};
  double const lp_traject{helper.comp_lp_traject(obs, cause)};
  double const d_lp_traject{helper.comp_d_lp_traject(obs, cause)};

  double out{std::log(d_lp_traject)};
  out += helper.comp_trajector_cond_dens_obs_one(lp_traject, cause);

  auto const n_causes = indexer.n_causes();
  arma::mat logit_offsets{mat_no_alloc(n_causes, 1, mem)};
  helper.fill_logit_offsets(logit_offsets.begin(), obs);

  arma::mat vcov_cond;
  helper.fill_cond_vcov_one_obs(vcov_cond, cause);

  arma::vec cond_mean;
  cond_mean = vcov_cond.col(cause + n_causes) * lp_traject;
  cond_mean = cond_mean.subvec(0, n_causes - 1);

  arma::mat vcov_sub{mat_no_alloc(n_causes, n_causes, mem)};
  vcov_sub = vcov_cond.submat(0, 0, n_causes - 1, n_causes - 1);

  arma::uvec which_cat{static_cast<arma::uword>(cause + 1)};

  auto mem_mark = mem.set_mark_raii();
  mixed_mult_logit_term<false> mult_term(logit_offsets, which_cat);
  rescale_shift_problem<false> prob_use(vcov_sub, cond_mean, mult_term);
  adaptive_problem prob(prob_use, mem, adaptive_rel_eps);
  double res{};
  ghq(&res, dat, prob, mem, ghq_target_size);
  out += std::log(res);

  return out;
}

double mmcif_log_Lik
  (double const * __restrict__ par, param_indexer const &indexer,
   mmcif_data const &obs1, mmcif_data const &obs2,
   simple_mem_stack<double> &mem, ghq_data const &dat){
  mmcif_comp_helper helper{indexer, par, mem};
  bool const is_cens1{helper.is_censored(obs1)},
             is_cens2{helper.is_censored(obs2)};

  if(is_cens1 && is_cens2)
    return mmcif_log_Lik_both_cens(par, indexer, obs1, obs2, mem, dat);
  else if(!is_cens1 && is_cens2)
    return mmcif_log_Lik_one_obs(par, indexer, obs1, obs2, mem, dat);
  else if(is_cens1 && !is_cens2)
    return mmcif_log_Lik_one_obs(par, indexer, obs2, obs1, mem, dat);
  return mmcif_log_Lik_both_obs(par, indexer, obs1, obs2, mem, dat);
}
