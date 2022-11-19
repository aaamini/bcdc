#ifndef C1C7BDA3_83BF_4219_A248_05E3AF8C4D76
#define C1C7BDA3_83BF_4219_A248_05E3AF8C4D76


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// int find_tunc(arma::vec beta, double threshold);
// arma::uvec get_freq(const arma::uvec& z, const int K);
// arma::uvec get_up_freq(arma::uvec freq);
// arma::umat get_freq_minus_self(arma::uvec z, int K);


// arma::mat comp_blk_sums(arma::sp_mat At, arma::uvec z, int Kcap);
// arma::mat sp_compress_col(const arma::sp_mat& At, const arma::uvec& z, const int& Kcap);

List comp_blk_sums_and_sizes(const arma::sp_mat& At, const arma::uvec& z, const int Kcap, const bool div_diag = true);

// arma::vec sp_single_col_compress(const arma::sp_mat& A, const int& col_idx, const arma::uvec& z, const int& Kcap);

void print_progress(int itr, int itr_max);                    



#endif /* C1C7BDA3_83BF_4219_A248_05E3AF8C4D76 */
