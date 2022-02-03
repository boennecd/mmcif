#ifndef MMCIF_MISC_H
#define MMCIF_MISC_H

#include "simple_mem_stack.h"
#include <RcppArmadillo.h>

inline arma::mat mat_no_alloc
  (arma::uword const n_rows, arma::uword const n_cols,
   ghqCpp::simple_mem_stack<double> &mem){
  return { mem.get(n_rows * n_cols), n_rows, n_cols, false };
}

inline arma::vec vec_no_alloc
  (arma::uword const n_ele, ghqCpp::simple_mem_stack<double> &mem){
  return { mem.get(n_ele), n_ele, false };
}

#endif
