#ifndef __BASICPSBM__
#define __BASICPSBM__
#include <RcppArmadillo.h>


class BasicPSBM {
    public:
        arma::sp_mat A;         // sparse adjacancy matrix
        arma::sp_mat mask;      // sparse mask matrix
        int K;                  // number of communities (or maximum number of communities)
        int n;
        arma::uvec z;           // label vector
        arma::mat M;            // block sums
        arma::mat N;           // block sizes


        BasicPSBM(const arma::sp_mat& A_, const arma::sp_mat& mask_);
        BasicPSBM(const arma::sp_mat& A_, const arma::sp_mat& mask_, const int K);


        arma::vec col_adj_compress(const int& col_idx);
        arma::vec col_mask_compress(const int& col_idx);

        arma::uvec get_label_freq();
        arma::uvec get_label_freq_except(const int&);

        void set_beta_params(const double alpha_eta, const double beta_eta);
        void set_z_to_random_labels();

        // void update_eta();
        // void update_pri();
        void update_blk_sums_and_sizes();
        // void update_z_element(const int i);

        // arma::umat run_gibbs(const int niter);

    private:
        const double perturb = 1e-11;
};


#endif /* __BASICSBM__ */


