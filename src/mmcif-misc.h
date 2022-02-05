#ifndef MMCIF_MISC_H
#define MMCIF_MISC_H

#include "simple_mem_stack.h"
#include <RcppArmadillo.h>
#include <numeric>

inline arma::mat mat_no_alloc
  (arma::uword const n_rows, arma::uword const n_cols,
   ghqCpp::simple_mem_stack<double> &mem){
  return { mem.get(n_rows * n_cols), n_rows, n_cols, false };
}

inline arma::vec vec_no_alloc
  (arma::uword const n_ele, ghqCpp::simple_mem_stack<double> &mem){
  return { mem.get(n_ele), n_ele, false };
}

inline double log_dmvn
  (arma::vec const &x, arma::mat const &Sigma,
   ghqCpp::simple_mem_stack<double> &mem){
  arma::uword const dim{Sigma.n_rows};
  arma::mat Sigma_chol{mat_no_alloc(dim, dim, mem)};
  Sigma_chol = arma::chol(Sigma, "lower");

  arma::vec x_scaled{vec_no_alloc(dim, mem)};
  x_scaled = arma::solve(arma::trimatl(Sigma_chol), x);

  double out{};
  for(arma::uword i = 0; i < dim; ++i)
    out -= x_scaled[i] * x_scaled[i];
  for(arma::uword i = 0; i < dim; ++i)
    out -= 2 * std::log(Sigma_chol(i, i));

  constexpr double log_2_pi{1.83787706640935};
  out -= static_cast<double>(dim) * log_2_pi;

  return out / 2;
}

#endif
