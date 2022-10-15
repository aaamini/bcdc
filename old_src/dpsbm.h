#ifndef __DPSBM__
#define __DPSBM__

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// void sbm_update_labels(
//     const arma::sp_mat& A,
//     const int s,
//     arma::uvec& z,
//     const int K,
//     arma::mat& m, 
//     arma::mat& mbar, 
//     const arma::vec& pri, 
//     const double alpha, const double beta);


arma::umat fit_dpsbm(const arma::sp_mat & A, 
                const double gam0 = 1,  
                const double alpha = 1, const double beta = 1,
                const int niter = 50, const int Zcap = 20, 
                const bool verb = true, const bool slice = false);


#endif /* __MODELS__ */
