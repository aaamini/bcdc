#ifndef __BASICSBM__
#define __BASICSBM__
#include <RcppArmadillo.h>



class BasicSBM {
    public:
        arma::sp_mat A;   // sparse adjacancy matrix
        int K;            // number of communities (or maximum number of communities)
        int n;
        arma::uvec z;  // label vector
        // arma::vec pri; 
        arma::mat M;    // block sums 
        arma::umat N;   // block sizes
        
        // double w0;
        // double pi0;

        // arma::vec z_log_prob_record; // for diagnostics

        BasicSBM(const arma::sp_mat& A_);
        BasicSBM(const arma::sp_mat& A_, const int K);


        arma::vec col_compress(const int& col_idx);
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
