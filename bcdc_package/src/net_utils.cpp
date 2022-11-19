// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

void print_progress(int itr, int itr_max) {
  int width = ceil(log10(itr_max));
  if (itr % 10 == 0) Rcout << "("
                            << std::setw(width) << itr 
                            << " / " 
                            << std::setw(width) << itr_max << ")\r";
}

//' @export
// [[Rcpp::export]]
arma::uvec get_freq(const arma::uvec& z, const int K) { 
  int n = z.n_elem;
  // int K = max(z)+1;
  arma::uvec freq(K,  arma::fill::zeros);
  // Rcout << K << freq; 
  for (int i = 0; i < n; i++) {
    freq(z(i))++;
  }
  return freq;
}


//' @export
// [[Rcpp::export]]
arma::mat comp_blk_sums(arma::sp_mat At, arma::uvec z, int Kcap) {
    // Compute block sums of a sparse matrix At w.r.t. labels in "z". The labels
    // in z are in the interval [0,1,...,Kcap]

    arma::sp_mat::const_iterator it     = At.begin();
    arma::sp_mat::const_iterator it_end = At.end();

    arma::mat lambda(Kcap, Kcap, arma::fill::zeros); 
    for(; it != it_end; ++it) {
        lambda(z(it.row()), z(it.col())) += (*it);
    }

    return lambda;
}


//' @export
// [[Rcpp::export]]
List comp_blk_sums_and_sizes(const arma::sp_mat& At, const arma::uvec& z, const int Kcap, const bool div_diag = true) {
    // z is a Kcap x 1 vector
    // At is a sparse n x n matrix

    int n = At.n_rows;
    // arma::mat lambda(Kcap, Kcap, arma::fill::zeros);

    arma::uvec zc = get_freq(z, Kcap); // zcounts

    arma::mat lambda = comp_blk_sums(At, z, Kcap);
    arma::umat NN = zc * zc.t() - arma::diagmat(zc);

    if (div_diag) {
        lambda.diag() /= 2; 
        NN.diag() /= 2;     // assumes that the diagonal of At is zero, 
                            // otherwise have to remove diag. first
    }
    
    return List::create( 
        Rcpp::Named("lambda") = lambda,
        Rcpp::Named("NN") = NN
    );
}
