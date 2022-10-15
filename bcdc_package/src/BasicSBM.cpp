// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "sampling.h"
#include "BasicSBM.h"

using namespace Rcpp;

BasicSBM::BasicSBM(const arma::sp_mat& A_)
{
    A = A_;
    n = A.n_rows;
}

BasicSBM::BasicSBM(const arma::sp_mat& A_, const int K) : K{K}
{
    // initialize and allocate variables
    A = A_;
    n = A.n_rows;
    M = arma::mat(K, K, arma::fill::zeros);
    N = arma::umat(K, K, arma::fill::zeros);
    set_z_to_random_labels();

    // Rcpp::print(wrap( blk_compressions[1]));
}

void BasicSBM::set_z_to_random_labels()
{
    z = sample_int_vec(K, n);
}

arma::vec BasicSBM::col_compress(const int& col_idx)
{
    arma::vec b(K, arma::fill::zeros);
    for (arma::sp_mat::const_col_iterator it = A.begin_col(col_idx); it != A.end_col(col_idx); ++it) {
        b(z(it.row())) += (*it);
    }
    return b;
}

arma::uvec BasicSBM::get_label_freq()
{
    arma::uvec freq(K, arma::fill::zeros);
    for (int i = 0; i < n; i++) {
        freq(z(i))++;
    }
    return freq;
}

arma::uvec BasicSBM::get_label_freq_except(const int& i) {
    arma::uvec freq = get_label_freq();
    freq(z(i))--;
    return freq;
}


void BasicSBM::update_blk_sums_and_sizes() {
    // arma::uvec zc = get_freq(z, K); // zcounts
    arma::uvec zc = get_label_freq(); // zcounts

    arma::sp_mat::const_iterator it     = A.begin();
    arma::sp_mat::const_iterator it_end = A.end();

    M = arma::mat(K, K, arma::fill::zeros);
    for(; it != it_end; ++it) {
        M(z(it.row()), z(it.col())) += (*it);
    }

    N = zc * zc.t() - arma::diagmat(zc);

    // assumes that the diagonal of At is zero,  otherwise have to remove diag. first
    M.diag() /= 2;
    N.diag() /= 2;
}


// void BasicSBM::update_z_element(const int i) {

//     arma::vec taui = col_compress(i); // sp_single_col_compress(A, i, z, K);
//     // arma::uvec mmi = get_freq(z, K);
//     // arma::uvec mmi = get_label_freq();
//     // mmi(z(i))--;
//     arma::uvec mmi = get_label_freq_except(i);
//     arma::vec log_prob = u * taui + v * mmi + log(pri);

//     z(i) = sample_index(safe_exp(log_prob));
// }


// arma::umat BasicSBM::run_gibbs(const int niter) {
//     // Run full Gibbs updates for "niter" iterations and record label history

//     arma::umat z_hist(n, niter+1);
//     z_hist.col(0) = z + 1;

//     // comp_count_tensors();
//     for (int iter = 0; iter < niter; iter++) {
//         update_eta();
//         update_pri();
//         for (int i = 0; i < n; i++) {
//             update_z_element(i);
//         }
//         // Rcpp::print(wrap(z.t()));
//         z_hist.col(iter+1) = z + 1;
//     } // iter

//     return z_hist;
//     // return Rcpp::List::create(
//     //     Rcpp::Named("z") = z_hist,
//     //     Rcpp::Named("xi") = xi_hist
//     // );
// }


// class BasicSBM {
//     public:
//         arma::sp_mat A;
//         int K;
//         int n;
//         arma::uvec z;
//         arma::mat eta;
//         arma::mat u;
//         arma::mat v;
//         arma::vec pri;

//         arma::mat M;
//         arma::umat N;

//         BetaParameters beta_params;
//         // double w0;
//         // double pi0;

//         // arma::vec z_log_prob_record; // for diagnostics

//         BasicSBM(const arma::sp_mat& A_,
//             // const arma::uvec z_init,
//             const int K,
//             const double alpha_eta,
//             const double beta_eta) : K{K} {

//             // initialize and allocate variables
//             A = A_;
//             n = A.n_rows;
//             beta_params.alpha = alpha_eta;
//             beta_params.beta = beta_eta;
//             set_z_to_random_labels();

//             eta = arma::mat(K, K, arma::fill::zeros);
//             u = arma::mat(K, K, arma::fill::zeros);
//             v = arma::mat(K, K, arma::fill::zeros);
//             pri = arma::vec(K, arma::fill::zeros);

//             M = arma::mat(K, K, arma::fill::zeros);
//             N = arma::umat(K, K, arma::fill::zeros);

//             // Rcpp::print(wrap( blk_compressions[1]));
//         }

//         void set_beta_params(const double alpha_eta, const double beta_eta) {
//             beta_params.alpha = alpha_eta;
//             beta_params.beta = beta_eta;
//         }

//         void set_z_to_random_labels() {
//              z = sample_int_vec(K, n);
//         }

//         void update_eta() {
//             // Update the eta-related tensors
//             // List out = comp_blk_sums_and_sizes(A, z, K);
//             // arma::mat lambda = out["lambda"];
//             // arma::umat NN = out["NN"];
//             // eta = symmat_rbeta(lambda + beta_params.alpha, NN - lambda + beta_params.beta);
//             update_blk_sums_and_sizes();

//             eta = symmat_rbeta(M + beta_params.alpha, N - M + beta_params.beta);

//             u = log(eta/(1-eta) + perturb);
//             v = log(1-eta + perturb);
//         }

//         void update_pri() {
//             arma::vec nn =  arma::conv_to<arma::vec>::from(get_freq(z, K));
//             pri = rdirichlet(nn + 1);
//         }

//         arma::vec col_compress(const int& col_idx) {
//             arma::vec b(K, arma::fill::zeros);

//             for (arma::sp_mat::const_col_iterator it = A.begin_col(col_idx); it != A.end_col(col_idx); ++it) {
//                 b(z(it.row())) += (*it);
//             }
//             return b;
//         }

//         void update_blk_sums_and_sizes() {
//             arma::uvec zc = get_freq(z, K); // zcounts

//             arma::sp_mat::const_iterator it     = A.begin();
//             arma::sp_mat::const_iterator it_end = A.end();

//             M = arma::mat(K, K, arma::fill::zeros);
//             for(; it != it_end; ++it) {
//                 M(z(it.row()), z(it.col())) += (*it);
//             }

//             N = zc * zc.t() - arma::diagmat(zc);

//             // assumes that the diagonal of At is zero,  otherwise have to remove diag. first
//             M.diag() /= 2;
//             N.diag() /= 2;
//         }


//         void update_z_element(const int i) {

//             arma::vec taui = col_compress(i); // sp_single_col_compress(A, i, z, K);
//             arma::uvec mmi = get_freq(z, K);
//             mmi(z(i))--;
//             arma::vec log_prob = u * taui + v * mmi + log(pri);

//             z(i) = sample_index(safe_exp(log_prob));
//         }


//         arma::umat run_gibbs(const int niter) {
//             // Run full Gibbs updates for "niter" iterations and record label history

//             arma::umat z_hist(n, niter+1);
//             z_hist.col(0) = z + 1;

//             // comp_count_tensors();
//             for (int iter = 0; iter < niter; iter++) {
//                 update_eta();
//                 update_pri();
//                 for (int i = 0; i < n; i++) {
//                     update_z_element(i);
//                 }
//                 // Rcpp::print(wrap(z.t()));
//                 z_hist.col(iter+1) = z + 1;
//             } // iter

//             return z_hist;
//             // return Rcpp::List::create(
//             //     Rcpp::Named("z") = z_hist,
//             //     Rcpp::Named("xi") = xi_hist
//             // );
//         }

//     private:
//         const double perturb = 1e-11;
// };


// void test_nsbm_cpp(List A, const int K, const int L) {

//     BasicSBM mynsbm(A, K, L);

//     mynsbm.run_gibbs(100);
// }

RCPP_MODULE(basic_sbm_module) {
      class_<BasicSBM>("BasicSBM")
      .constructor<arma::sp_mat, int>()
      .field("A", &BasicSBM::A)
      .field("K", &BasicSBM::K)
      .field("n", &BasicSBM::n)
      .field("z", &BasicSBM::z)
      .field("M", &BasicSBM::M)
      .field("N", &BasicSBM::N)
      .method("set_z_to_random_labels", &BasicSBM::set_z_to_random_labels)
      .method("update_blk_sums_and_sizes", &BasicSBM::update_blk_sums_and_sizes)
      .method("get_label_freq", &BasicSBM::get_label_freq)
      .method("get_label_freq_except", &BasicSBM::get_label_freq_except)
      .method("col_compress", &BasicSBM::col_compress)
//      .method("print", &BasicSBM::print)
    //   .method("update_z_element", &BasicSBM::update_z_element)
    //   .method("run_gibbs", &BasicSBM::run_gibbs)
    //   .method("update_pri", &BasicSBM::update_pri)
      ;
};

RCPP_EXPOSED_CLASS(BasicSBM);
