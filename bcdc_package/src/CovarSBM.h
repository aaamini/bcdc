#ifndef B1EA7B78_E870_42A5_9166_3BC22BF3327A
#define B1EA7B78_E870_42A5_9166_3BC22BF3327A
#ifndef __COVARSBM__
#define __COVARSBM__

#include <RcppArmadillo.h>

#include "BasicSBM.h"
#include "BetaParameters.h"


class CovarSBM : public BasicSBM {
    public: 
        arma::mat eta; // connectivity matrix
        arma::mat u;
        arma::mat v;
        BetaParameters beta_params;
        int n_comm = 1;   // number of communities; we always matain K >= n_comm + 1
        int next_label = 1;
        double dp_concent_param;

        bool has_cont_features = false;
        bool has_disc_features = false;

        arma::mat Xc;
        arma::umat Xd;
        arma::mat xic;
        std::vector<arma::mat> xid;
        arma::uvec n_cats; 
        int feature_dim = 0;
        int feature_dim_cont = 0;
        int feature_dim_disc = 0;

        double s2 = 1;
        double tau2 = 1;
        // arma::mat dist_mat;

        // double w0;
        // double pi0;


        CovarSBM(const arma::sp_mat& A_, 
            const double alpha_eta, const double beta_eta,
            const double dp_concent_param);

        // void reset_K();
        void resize();
        // void resize_mats();

        void set_discrete_features(const arma::umat&);
        void set_continuous_features(const arma::mat&);
        void repopulate_from_prior();
        // arma::mat update_xi();
        void update_xi();
        // void update_dist_mat();
        void set_gauss_param(const double, const double);
        void sample_crp();
        void set_beta_params(const double alpha_eta, const double beta_eta);
        void update_eta();
        void update_z_element(const int i);
        arma::umat run_gibbs(const int niter);

    private:
        const double perturb = 1e-11;
};


#endif /* __COVARSBM__ */


#endif /* B1EA7B78_E870_42A5_9166_3BC22BF3327A */
