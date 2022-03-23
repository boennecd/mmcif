// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// mmcif_data_holder_to_R
SEXP mmcif_data_holder_to_R(NumericMatrix const covs_trajectory, NumericMatrix const d_covs_trajectory, NumericMatrix const covs_risk, IntegerVector const has_finite_trajectory_prob, IntegerVector const cause, size_t const n_causes, Rcpp::IntegerMatrix pair_indices, IntegerVector const singletons, NumericMatrix const covs_trajectory_delayed);
RcppExport SEXP _mmcif_mmcif_data_holder_to_R(SEXP covs_trajectorySEXP, SEXP d_covs_trajectorySEXP, SEXP covs_riskSEXP, SEXP has_finite_trajectory_probSEXP, SEXP causeSEXP, SEXP n_causesSEXP, SEXP pair_indicesSEXP, SEXP singletonsSEXP, SEXP covs_trajectory_delayedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericMatrix const >::type covs_trajectory(covs_trajectorySEXP);
    Rcpp::traits::input_parameter< NumericMatrix const >::type d_covs_trajectory(d_covs_trajectorySEXP);
    Rcpp::traits::input_parameter< NumericMatrix const >::type covs_risk(covs_riskSEXP);
    Rcpp::traits::input_parameter< IntegerVector const >::type has_finite_trajectory_prob(has_finite_trajectory_probSEXP);
    Rcpp::traits::input_parameter< IntegerVector const >::type cause(causeSEXP);
    Rcpp::traits::input_parameter< size_t const >::type n_causes(n_causesSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type pair_indices(pair_indicesSEXP);
    Rcpp::traits::input_parameter< IntegerVector const >::type singletons(singletonsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix const >::type covs_trajectory_delayed(covs_trajectory_delayedSEXP);
    rcpp_result_gen = Rcpp::wrap(mmcif_data_holder_to_R(covs_trajectory, d_covs_trajectory, covs_risk, has_finite_trajectory_prob, cause, n_causes, pair_indices, singletons, covs_trajectory_delayed));
    return rcpp_result_gen;
END_RCPP
}
// mmcif_logLik_to_R
double mmcif_logLik_to_R(SEXP data_ptr, NumericVector const par, Rcpp::List ghq_data, unsigned n_threads);
RcppExport SEXP _mmcif_mmcif_logLik_to_R(SEXP data_ptrSEXP, SEXP parSEXP, SEXP ghq_dataSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type data_ptr(data_ptrSEXP);
    Rcpp::traits::input_parameter< NumericVector const >::type par(parSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type ghq_data(ghq_dataSEXP);
    Rcpp::traits::input_parameter< unsigned >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(mmcif_logLik_to_R(data_ptr, par, ghq_data, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// mmcif_logLik_grad_to_R
Rcpp::NumericVector mmcif_logLik_grad_to_R(SEXP data_ptr, NumericVector const par, Rcpp::List ghq_data, unsigned n_threads);
RcppExport SEXP _mmcif_mmcif_logLik_grad_to_R(SEXP data_ptrSEXP, SEXP parSEXP, SEXP ghq_dataSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type data_ptr(data_ptrSEXP);
    Rcpp::traits::input_parameter< NumericVector const >::type par(parSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type ghq_data(ghq_dataSEXP);
    Rcpp::traits::input_parameter< unsigned >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(mmcif_logLik_grad_to_R(data_ptr, par, ghq_data, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// mcif_logLik_to_R
double mcif_logLik_to_R(SEXP data_ptr, NumericVector const par, unsigned n_threads, bool const with_risk);
RcppExport SEXP _mmcif_mcif_logLik_to_R(SEXP data_ptrSEXP, SEXP parSEXP, SEXP n_threadsSEXP, SEXP with_riskSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type data_ptr(data_ptrSEXP);
    Rcpp::traits::input_parameter< NumericVector const >::type par(parSEXP);
    Rcpp::traits::input_parameter< unsigned >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< bool const >::type with_risk(with_riskSEXP);
    rcpp_result_gen = Rcpp::wrap(mcif_logLik_to_R(data_ptr, par, n_threads, with_risk));
    return rcpp_result_gen;
END_RCPP
}
// mcif_logLik_grad_to_R
Rcpp::NumericVector mcif_logLik_grad_to_R(SEXP data_ptr, NumericVector const par, unsigned n_threads, bool const with_risk);
RcppExport SEXP _mmcif_mcif_logLik_grad_to_R(SEXP data_ptrSEXP, SEXP parSEXP, SEXP n_threadsSEXP, SEXP with_riskSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type data_ptr(data_ptrSEXP);
    Rcpp::traits::input_parameter< NumericVector const >::type par(parSEXP);
    Rcpp::traits::input_parameter< unsigned >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< bool const >::type with_risk(with_riskSEXP);
    rcpp_result_gen = Rcpp::wrap(mcif_logLik_grad_to_R(data_ptr, par, n_threads, with_risk));
    return rcpp_result_gen;
END_RCPP
}
// ns_ptr
SEXP ns_ptr(const arma::vec& boundary_knots, const arma::vec& interior_knots);
RcppExport SEXP _mmcif_ns_ptr(SEXP boundary_knotsSEXP, SEXP interior_knotsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type boundary_knots(boundary_knotsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type interior_knots(interior_knotsSEXP);
    rcpp_result_gen = Rcpp::wrap(ns_ptr(boundary_knots, interior_knots));
    return rcpp_result_gen;
END_RCPP
}
// ns_eval
Rcpp::NumericMatrix ns_eval(SEXP ptr, Rcpp::NumericVector const points, int const ders);
RcppExport SEXP _mmcif_ns_eval(SEXP ptrSEXP, SEXP pointsSEXP, SEXP dersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< SEXP >::type ptr(ptrSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector const >::type points(pointsSEXP);
    Rcpp::traits::input_parameter< int const >::type ders(dersSEXP);
    rcpp_result_gen = Rcpp::wrap(ns_eval(ptr, points, ders));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP run_testthat_tests(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_mmcif_mmcif_data_holder_to_R", (DL_FUNC) &_mmcif_mmcif_data_holder_to_R, 9},
    {"_mmcif_mmcif_logLik_to_R", (DL_FUNC) &_mmcif_mmcif_logLik_to_R, 4},
    {"_mmcif_mmcif_logLik_grad_to_R", (DL_FUNC) &_mmcif_mmcif_logLik_grad_to_R, 4},
    {"_mmcif_mcif_logLik_to_R", (DL_FUNC) &_mmcif_mcif_logLik_to_R, 4},
    {"_mmcif_mcif_logLik_grad_to_R", (DL_FUNC) &_mmcif_mcif_logLik_grad_to_R, 4},
    {"_mmcif_ns_ptr", (DL_FUNC) &_mmcif_ns_ptr, 2},
    {"_mmcif_ns_eval", (DL_FUNC) &_mmcif_ns_eval, 3},
    {"run_testthat_tests", (DL_FUNC) &run_testthat_tests, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_mmcif(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
