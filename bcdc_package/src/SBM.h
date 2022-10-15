#ifndef C6DA5422_E96D_41E0_9027_9D4C4F0AF4E1
#define C6DA5422_E96D_41E0_9027_9D4C4F0AF4E1
#ifndef __SBM__
#define __SBM__
#include <RcppArmadillo.h>

#include "BasicSBM.h"
#include "BetaParameters.h"

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


#endif /* C6DA5422_E96D_41E0_9027_9D4C4F0AF4E1 */
