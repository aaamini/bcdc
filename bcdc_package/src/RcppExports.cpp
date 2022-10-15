// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sample_int
int sample_int(int N);
RcppExport SEXP _bcdc_sample_int(SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_int(N));
    return rcpp_result_gen;
END_RCPP
}
// sample_int_vec
arma::uvec sample_int_vec(int N, int size);
RcppExport SEXP _bcdc_sample_int_vec(SEXP NSEXP, SEXP sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_int_vec(N, size));
    return rcpp_result_gen;
END_RCPP
}
// stick_break
arma::vec stick_break(arma::vec x);
RcppExport SEXP _bcdc_stick_break(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(stick_break(x));
    return rcpp_result_gen;
END_RCPP
}
// symmat_rbeta
arma::mat symmat_rbeta(arma::mat alpha, arma::mat beta);
RcppExport SEXP _bcdc_symmat_rbeta(SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(symmat_rbeta(alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// safe_exp
arma::vec safe_exp(arma::vec log_prob);
RcppExport SEXP _bcdc_safe_exp(SEXP log_probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type log_prob(log_probSEXP);
    rcpp_result_gen = Rcpp::wrap(safe_exp(log_prob));
    return rcpp_result_gen;
END_RCPP
}
// my_sampler
int my_sampler(arma::vec prob_vec);
RcppExport SEXP _bcdc_my_sampler(SEXP prob_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type prob_vec(prob_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(my_sampler(prob_vec));
    return rcpp_result_gen;
END_RCPP
}
// sample_index
int sample_index(arma::vec prob);
RcppExport SEXP _bcdc_sample_index(SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_index(prob));
    return rcpp_result_gen;
END_RCPP
}
// rgamma_vec
arma::vec rgamma_vec(arma::vec shape);
RcppExport SEXP _bcdc_rgamma_vec(SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(rgamma_vec(shape));
    return rcpp_result_gen;
END_RCPP
}
// rbeta_vec
arma::vec rbeta_vec(arma::vec alpha, arma::vec beta);
RcppExport SEXP _bcdc_rbeta_vec(SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(rbeta_vec(alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// rdirichlet
arma::vec rdirichlet(arma::vec theta);
RcppExport SEXP _bcdc_rdirichlet(SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(rdirichlet(theta));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_sbm_module();

static const R_CallMethodDef CallEntries[] = {
    {"_bcdc_sample_int", (DL_FUNC) &_bcdc_sample_int, 1},
    {"_bcdc_sample_int_vec", (DL_FUNC) &_bcdc_sample_int_vec, 2},
    {"_bcdc_stick_break", (DL_FUNC) &_bcdc_stick_break, 1},
    {"_bcdc_symmat_rbeta", (DL_FUNC) &_bcdc_symmat_rbeta, 2},
    {"_bcdc_safe_exp", (DL_FUNC) &_bcdc_safe_exp, 1},
    {"_bcdc_my_sampler", (DL_FUNC) &_bcdc_my_sampler, 1},
    {"_bcdc_sample_index", (DL_FUNC) &_bcdc_sample_index, 1},
    {"_bcdc_rgamma_vec", (DL_FUNC) &_bcdc_rgamma_vec, 1},
    {"_bcdc_rbeta_vec", (DL_FUNC) &_bcdc_rbeta_vec, 2},
    {"_bcdc_rdirichlet", (DL_FUNC) &_bcdc_rdirichlet, 1},
    {"_rcpp_module_boot_sbm_module", (DL_FUNC) &_rcpp_module_boot_sbm_module, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_bcdc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
