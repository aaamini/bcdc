// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <random>
#include <Rcpp/Benchmark/Timer.h>

#include "sampling.h" 
#include "utils.h" 

using namespace Rcpp;


// [[Rcpp::export]]
arma::umat multsbm_gibbs_sampler_fast(arma::sp_mat A, const int K, 
                    const double alpha = 1, const double beta = 1,
                    const int niter = 100) {

    int n = A.n_rows;

    arma::umat z_hist(n, niter, arma::fill::zeros);
    arma::uvec z = sample_int_vec(K, n);
    arma::mat B(K, K, arma::fill::zeros);   
    arma::vec pri(K); 

    for (int iter = 0; iter < niter; iter++) {
        List out = comp_blk_sums_and_sizes(A, z, K);
        arma::mat lambda = out["lambda"];
        arma::umat NN = out["NN"]; 
        B = symmat_rbeta(lambda + alpha, NN - lambda + beta);
        arma::mat uu = log(B/(1-B) + 1e-11);
        arma::mat vv = log(1-B + 1e-11);

        arma::vec nn =  arma::conv_to<arma::vec>::from(get_freq(z, K));
        pri = rdirichlet(nn + 1);

        for (int i = 0; i < n; i++) {
            arma::vec taui = sp_single_col_compress(A, i, z, K);
            arma::uvec mmi = get_freq(z, K);
            mmi(z(i))--;
            arma::vec Hi = uu * taui + vv * mmi +  log(pri);

            z(i) = sample_index(safe_exp(Hi));
        }
        z_hist.col(iter) = z;
    }
    return z_hist;
}


// [[Rcpp::export]]
arma::umat multsbm_collapsed_gibbs_sampler(arma::sp_mat& A, const int K, 
                    const double alpha = 1, const double beta = 1,
                    const int niter = 10) {

    int n = A.n_rows;
    arma::umat z_hist(n, niter, arma::fill::zeros);
    
    // Randomly initialize labels
    arma::uvec z = sample_int_vec(K, n);
    arma::uvec z_new(n, arma::fill::zeros);
    arma::vec pri(K, arma::fill::zeros);
    arma::mat Bet(K, K);

    for (int iter = 0; iter < niter; iter++) {
        for (int s = 0; s < n; s++) {
            Bet = comp_beta_matrix(A, z, K, alpha, beta);
            arma::vec nn =  arma::conv_to<arma::vec>::from(get_freq(z, K));
            pri = rdirichlet(nn + 1);

            arma::vec prob(K, arma::fill::ones);
            for (int rp = 0; rp < K; rp++) { // rp is the potential new value of z(s)
                z_new = z;
                z_new(s) = rp;
                arma::mat Bet_new =  comp_beta_matrix(A, z_new, K, alpha, beta);;
                prob(rp) *= arma::prod(arma::prod(Bet_new / (Bet + DBL_MIN))) * pri(rp);
                // prob(rp) *= static_cast<double>((nn(rp) + 1)) / nn(z(s)); 
            }

            z(s) = sample_index(prob);; // update z
        }
        z_hist.col(iter) = z;
    }

    return z_hist;    
    
}


// [[Rcpp::export]]
arma::umat multsbm_collapsed_gibbs_sampler_v2(
    arma::sp_mat& A, const int K, 
    const double alpha = 1, const double beta = 1,
    const int niter = 10
){

    int n = A.n_rows;
    arma::umat z_hist(n, niter, arma::fill::zeros);
    
    // Randomly initialize labels
    arma::uvec z = sample_int_vec(K, n);
    arma::uvec z_new(n, arma::fill::zeros);
    arma::vec pri(K, arma::fill::zeros);

    // initialize m and mbar
    List out = comp_blk_sums_and_sizes(A, z, K);
    arma::mat m = out["lambda"];
    arma::umat NN = out["NN"];
    arma::mat mbar = NN - m; 


    for (int iter = 0; iter < niter; iter++) {
        for (int s = 0; s < n; s++) {
            // Bet = comp_beta_matrix(A, z, K, alpha, beta);
            arma::vec nn =  arma::conv_to<arma::vec>::from(get_freq(z, K));
            pri = rdirichlet(nn + 1);

            // arma::vec prob(K, arma::fill::ones);
            arma::vec U = sp_single_col_compress(A, s, z, K);
            arma::uvec V = get_freq(z, K);
            V(z(s))--;

            int zs_old = z(s);
            // prob is K x 1 vector`
            arma::vec temp = comp_beta_ratio_prods_v1(m, mbar, U, V, zs_old, alpha, beta);
            
            
            // if (temp.has_nan()) {
            //     m.save("m.csv", arma::csv_ascii);
            //     mbar.save("mbar.csv", arma::csv_ascii);
            //     U.save("U.csv", arma::csv_ascii);
            //     V.save("V.csv", arma::csv_ascii);
            //     z.save("z.csv", arma::csv_ascii);
            //     Rcout << "s = " << s << "zs_old = " << zs_old << "\n";
            //     print(wrap(temp));
            // }
            
            arma::vec prob = temp % pri;

            z(s) = sample_index(prob);; // update z(s) -- this the zs_new we pick
            
            if (z(s) != zs_old) { // check if we need to update counts
                // update m and mbar
                arma::mat D = comp_blk_sums_diff_v2(U, z(s), zs_old);
                arma::mat DN = comp_blk_sums_diff_v2(arma::conv_to<arma::vec>::from(V), z(s), zs_old);
                m += D;
                mbar += DN - D;                
            }

            // if (iter == 5 && s > 100) {
            //     List out = comp_blk_sums_and_sizes(A, z, K);
            //     arma::mat mtemp = out["lambda"];
            //     arma::mat NN_temp = out["NN"];
            //     arma::mat mbar_temp = NN_temp - mtemp; 
            //     Rcout << "s = " << s << "\n";
            //     // print(wrap(mbar));
            //     // print(wrap(mtemp));
            //     // print(wrap(m-mtemp));
            //     print(wrap(mbar));
            //     print(wrap(mbar_temp));
            //     print(wrap(mbar-mbar_temp));
            //     Rcout << "\n";
            //     // if (s > 50) return z_hist;
            //  }
            
        }
        z_hist.col(iter) = z;
    }

    return z_hist;    
    
}

// [[Rcpp::export]]
arma::umat multsbm_collapsed_gibbs_sampler_v3(
    arma::sp_mat& A, const int K, 
    const double alpha = 1, const double beta = 1,
    const int niter = 10
){

    int n = A.n_rows;
    arma::umat z_hist(n, niter, arma::fill::zeros);
    
    // Randomly initialize labels
    arma::uvec z = sample_int_vec(K, n);
    arma::uvec z_new(n, arma::fill::zeros);
    arma::vec pri(K, arma::fill::zeros);

    // initialize m and mbar
    List out = comp_blk_sums_and_sizes(A, z, K);
    arma::mat m = out["lambda"];
    arma::umat NN = out["NN"];
    arma::mat mbar = NN - m; 


    for (int iter = 0; iter < niter; iter++) {
        for (int s = 0; s < n; s++) {
            // Bet = comp_beta_matrix(A, z, K, alpha, beta);
            arma::vec nn =  arma::conv_to<arma::vec>::from(get_freq(z, K));
            pri = rdirichlet(nn + 1);

            // arma::vec prob(K, arma::fill::ones);
            arma::vec U = sp_single_col_compress(A, s, z, K);
            arma::uvec V = get_freq(z, K);
            V(z(s))--;

            int zs_old = z(s);
            // prob is K x 1 vector`
            arma::vec temp = comp_log_beta_ratio_sums(m, mbar, U, V, zs_old, alpha, beta);
            
            
            // if (temp.has_nan()) {
            //     m.save("m.csv", arma::csv_ascii);
            //     mbar.save("mbar.csv", arma::csv_ascii);
            //     U.save("U.csv", arma::csv_ascii);
            //     V.save("V.csv", arma::csv_ascii);
            //     z.save("z.csv", arma::csv_ascii);
            //     Rcout << "s = " << s << "zs_old = " << zs_old << "\n";
            //     print(wrap(temp));
            // }
            
            arma::vec log_prob = temp + log(pri);

            z(s) = sample_index(safe_exp(log_prob)); // update z(s) -- this the zs_new we pick

            if (z(s) != zs_old) { // check if we need to update counts
                // update m and mbar
                arma::mat D = comp_blk_sums_diff_v2(U, z(s), zs_old);
                arma::mat DN = comp_blk_sums_diff_v2(arma::conv_to<arma::vec>::from(V), z(s), zs_old);
                m += D;
                mbar += DN - D;
            }
            // if (iter == 5 && s > 100) {
            //     List out = comp_blk_sums_and_sizes(A, z, K);
            //     arma::mat mtemp = out["lambda"];
            //     arma::mat NN_temp = out["NN"];
            //     arma::mat mbar_temp = NN_temp - mtemp; 
            //     Rcout << "s = " << s << "\n";
            //     // print(wrap(mbar));
            //     // print(wrap(mtemp));
            //     // print(wrap(m-mtemp));
            //     print(wrap(mbar));
            //     print(wrap(mbar_temp));
            //     print(wrap(mbar-mbar_temp));
            //     Rcout << "\n";
            //     // if (s > 50) return z_hist;
            //  }
            
        }
        z_hist.col(iter) = z;
    }

    return z_hist;    
    
}


// [[Rcpp::export]]
arma::umat multsbm_collapsed_gibbs_sampler_v4(
    arma::sp_mat& A, const int K, 
    const double alpha = 1, const double beta = 1,
    const int niter = 10
){

    int n = A.n_rows;
    arma::umat z_hist(n, niter, arma::fill::zeros);
    
    // Randomly initialize labels
    arma::uvec z = sample_int_vec(K, n);
    arma::uvec z_new(n, arma::fill::zeros);
    arma::vec pri(K, arma::fill::zeros);

    // initialize m and mbar
    List out = comp_blk_sums_and_sizes(A, z, K);
    arma::mat m = out["lambda"];
    arma::umat NN = out["NN"];
    arma::mat mbar = NN - m; 


    for (int iter = 0; iter < niter; iter++) {
        for (int s = 0; s < n; s++) {

            arma::vec nn =  arma::conv_to<arma::vec>::from(get_freq(z, K));
            pri = rdirichlet(nn + 1);

            sbm_update_labels(A, s, z, K, m, mbar, pri, alpha, beta);
            
        }
        z_hist.col(iter) = z;
    }

    return z_hist;    
    
}


// arma::umat fit_multsbm(arma::sp_mat A, const int K, 
//                     const double alpha = 1, const double beta = 1,
//                     const int niter = 100,
//                     bool collapsed = false) {
//     if (collapsed) {
//         return multsbm_collapsed_gibbs_sampler(A, K, alpha, beta, niter);
//     } else {
//         return multsbm_gibbs_sampler_fast(A, K, alpha, beta, niter);
//     }
// }
