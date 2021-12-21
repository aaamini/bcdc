#ifndef __COVARSBM__
#define __COVARSBM__

#include <RcppArmadillo.h>

#include "BasicSBM.h"

struct BetaParameters {
    double alpha;
    double beta;

    BetaParameters() :alpha{1}, beta{1} {};
    BetaParameters(const double a, const double b) :alpha{a}, beta{b} {};

    void set(const double a, const double b)
    {
        alpha = a;
        beta = b;
    }
};

class CovarSBM : public BasicSBM {
    public: 
        arma::mat eta; // connectivity matrix
        arma::mat u;
        arma::mat v;
        BetaParameters beta_params;
        int n_comm = 1;   // number of communities; we always matain K >= n_comm + 1
        int next_label = 1;
        double dp_concent_param;

        bool has_node_features = false;
        arma::mat X;
        arma::mat xi;
        int feature_dim = 0;

        double s2 = 1;
        double tau2 = 1;
        arma::mat dist_mat;

        // double w0;
        // double pi0;


        CovarSBM(const arma::sp_mat& A_, 
            const double alpha_eta, const double beta_eta,
            const double dp_concent_param);

        // void reset_K();
        void resize();
        // void resize_mats();

        void repopulate_from_prior();
        arma::mat update_xi();
        void update_dist_mat();
        void set_node_features(const arma::mat&);
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
