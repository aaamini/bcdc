#ifndef __SBM__
#define __SBM__
#include <RcppArmadillo.h>

#include "BasicSBM.h"

struct BetaParameters {
    double alpha;
    double beta;

    BetaParameters() :alpha{1}, beta{1} {};
    BetaParameters(const double a, const double b) :alpha{a}, beta{b} {};
};

class SBM : public BasicSBM {
    public: 
        arma::vec pri; 
        arma::mat eta; // connectivity matrix
        arma::mat u;
        arma::mat v;
        BetaParameters beta_params;
        // double w0;
        // double pi0;


        SBM(const arma::sp_mat& A_, const int K, const double alpha_eta, const double beta_eta);

        void set_beta_params(const double alpha_eta, const double beta_eta);
        void update_pri();
        void update_eta();
        void update_z_element(const int i);
        arma::umat run_gibbs(const int niter);

    private:
        const double perturb = 1e-11;
};



#endif /* __SBM__ */
