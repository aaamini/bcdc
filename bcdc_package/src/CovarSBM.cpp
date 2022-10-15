// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


#include "sampling.h"
#include "CovarSBM.h"


using namespace Rcpp;


CovarSBM::CovarSBM(const arma::sp_mat& A_,
            const double alpha_eta, const double beta_eta,
            const double dp_concent_param)
            : BasicSBM::BasicSBM(A_), dp_concent_param{dp_concent_param}
{
     // initialize and allocate variables
    beta_params.alpha = alpha_eta;
    beta_params.beta = beta_eta;

    sample_crp(); // This also sets K

    M = arma::mat(K, K, arma::fill::zeros);
    N = arma::umat(K, K, arma::fill::zeros);
    eta = arma::mat(K, K, arma::fill::zeros);
    u = arma::mat(K, K, arma::fill::zeros);
    v = arma::mat(K, K, arma::fill::zeros);
}

void CovarSBM::resize()
{
    // K has already changed elsewhere. We just resize everything to match K.
    M.resize(K, K);
    N.resize(K, K);
    eta.resize(K, K); // the new entries has to filled with eta prior
    u.resize(K, K);
    v.resize(K, K);

    if (has_cont_features)
    {
        xic.resize(K, feature_dim_cont);
    }
    if (has_disc_features)
    {
        // xid.resize(K, feature_dim_disc);
        for (int r = 0; r < feature_dim_disc; r++)
        {
            xid[r].resize(n_cats(r), K);
        }
    }
}

// void CovarSBM::reset_K()
// {
//     int K_old = K;
//     K = std::max( std::round(growth_factor*max_label), max_label+1.0);

//     M.resize(K, K);
//     N.resize(K, K);
//     eta.resize(K, K); // the new entries has to filled with eta prior
//     u.resize(K, K);
//     v.resize(K, K);

//     if (has_node_features) {
//         xi.resize(K, feature_dim);
//         for (int k = K_old; k < K; k++)
//         {
//             xi.row(k) = sqrt(tau2) * arma::vec(feature_dim, arma::fill::randn);
//         }
//         dist_mat.resize(K, n);
//     }

// }

void CovarSBM::set_continuous_features(const arma::mat& X_)
{
    Xc = X_;
    feature_dim_cont = Xc.n_cols;
    // feature_dim += feature_dim_cont;
    has_cont_features = true;
    xic = arma::mat(K, feature_dim_cont, arma::fill::zeros);
    // dist_mat = arma::mat(K, n, arma::fill::zeros);
}

void CovarSBM::set_discrete_features(const arma::umat& X_)
{
    Xd = X_;
    feature_dim_disc = Xd.n_cols;
    // feature_dim += feature_dim_disc;
    has_disc_features = true;
    //  xid = arma::mat(K, feature_dim_disc, arma::fill::zeros);

    //n_cats = arma::uvec(feature_dim_disc);
    n_cats = arma::max(Xd, 0).t() + 1;
    xid.reserve(feature_dim_disc);
    for (int r = 0; r < feature_dim_disc; r++)
    {
        xid.push_back(arma::mat(n_cats(r), K, arma::fill::zeros));
    }
}

void CovarSBM::set_gauss_param(const double s2_, const double tau2_)
{
    s2 = s2_;
    tau2 = tau2_;
}

void CovarSBM::set_beta_params(const double alpha_eta, const double beta_eta)
{
    // beta_params.alpha = alpha_eta;
    // beta_params.beta = beta_eta;
    beta_params.set(alpha_eta, beta_eta);
}


// This also sets K = n_comm + 1
void CovarSBM::sample_crp()
{
    z = arma::uvec(n, arma::fill::zeros);
    n_comm = 1;
    arma::uvec freq(n_comm, arma::fill::zeros);
    freq(z(0)) = 1;

    // arma::vec prob(K, arma::fill::zeros);
    // Rcpp::print(wrap(z.t()));
    // Rcpp::print(wrap(freq.t()));
    for (int i = 1; i < n; i++) {
        auto prob = arma::vec(n_comm+1, arma::fill::zeros);
        auto idx = arma::span(0, n_comm-1); // K-1 is the number of clusters
        prob(idx) = arma::conv_to<arma::vec>::from(
            freq(idx)
        );
        prob(n_comm) = dp_concent_param;
        z(i) = sample_index(prob);
        if (z(i) == n_comm) {
            n_comm++;
            freq.resize(n_comm);
        }
        freq(z(i))++;
        // Rcpp::print(wrap(z.t()));
        // Rcpp::print(wrap(freq.t()));
    }
    K = n_comm + 1; // one more than the total number of communities in "z"
    next_label = n_comm; // current labels are 0,...,n_comm-1
}


// void CovarSBM::update_dist_mat() {
//     // for (int i = 0; i < n; i++) {
//     //     for (int k = 0; k <= max_label; k++) { // compute one more row
//     //         dist_mat(k, i) = pow(arma::norm(X.row(i) - xi.row(k)), 2) / (2*s2);
//     //     }
//     // }
// }

void CovarSBM::update_xi()
{
    if (has_cont_features)
    {
        arma::mat centers(K, feature_dim_cont, arma::fill::zeros);
        auto sds = arma::vec(K, arma::fill::ones);
        auto freq = get_label_freq();

        // compute the mean of features over z-clusters
        for (int i = 0; i < n; i++) {
            centers.row(z(i)) += Xc.row(i);
        }

        double temp = 0;
        for (int k = 0; k < K; k++) {
            if (freq(k) > 0 || k == next_label) // We only need to update these indices
            {
                temp = tau2 / (freq(k)*tau2 + s2); // freq(k) will be zero for k == next_label
                centers.row(k) *= temp; // centers.row(k) will be zero for k == next_label
                sds(k) = sqrt(s2*temp);

                xic.row(k) = centers.row(k) + sds(k)*arma::rowvec(feature_dim_cont, arma::fill::randn);
            }
        }
    }
    if (has_disc_features)
    {
        for (int r  = 0; r < feature_dim_disc; r++)
        {
            arma::mat dirch_param(n_cats(r), K, arma::fill::ones); // Assumes Dirichelt prior parameter is 1
            for(int i = 0; i < n; i++)
            {
                dirch_param(Xd(i,r), z(i))++;
            }
            for (int k = 0; k < K; k++)
            {
                xid[r].col(k) = rdirichlet(dirch_param.col(k));
            }
        }
    }



    // double temp = 0;
    // for (int k = 0; k < max_label; k++) {
    //     temp = tau2 / (freq(k)*tau2 + s2);
    //     centers.row(k) *= temp;
    //     sds(k) = sqrt(s2*temp);
    // }
    // for (int k = max_label; k < K; k++) {
    //     sds(k) = sqrt(tau2);
    // }

    // // generate the normal sample
    // xi = arma::mat(K, feature_dim, arma::fill::randn);
    // for (int k = 0; k < K; k++) {
    //     xi.row(k) *= sds(k);
    // }
    // xi += centers;

    // return centers;
}

void CovarSBM::update_eta()
{
    // Update the eta-related tensors
    // List out = comp_blk_sums_and_sizes(A, z, K);
    // arma::mat lambda = out["lambda"];
    // arma::umat NN = out["NN"];
    // eta = symmat_rbeta(lambda + beta_params.alpha, NN - lambda + beta_params.beta);
    update_blk_sums_and_sizes();

    eta = symmat_rbeta(M + beta_params.alpha, N - M + beta_params.beta);

    u = log(eta/(1-eta) + perturb);
    v = log(1-eta + perturb);

    // if (arma::any(arma::any(u == arma::datum::inf))) {
    //     Rcout << "update_eta" << u << "\n" << eta
    //             << "\nM = \n" << M
    //             << "\nN = \n" << N
    //             << "\nN - M = \n" << N-M
    //             << "freq = " << get_label_freq()
    //             << "K = " << K << " , next_label = " << next_label;
    //     z.save("z.csv", arma::csv_ascii);
    //     A.save("A.csv", arma::coord_ascii);
    // }
}

void CovarSBM::repopulate_from_prior(){
    // Repopulate eta, xi, etc. over next_label row/column from the prior
    for (int l = 0; l < K; l++)
    {
        double temp = R::rbeta(beta_params.alpha, beta_params.beta);
        eta(next_label, l) = temp;
        eta(l, next_label) = temp;

        u(next_label, l) = log(temp/(1-temp) + perturb);
        v(next_label, l) = log(1-temp + perturb);

        u(l, next_label) = u(next_label, l);
        v(l, next_label) = v(next_label, l);
    }
    if (has_cont_features)
    {
        xic.row(next_label) = sqrt(tau2)*arma::rowvec(feature_dim_cont, arma::fill::randn);
    }
    if (has_disc_features)
    {
        for (int r  = 0; r < feature_dim_disc; r++)
        {
            arma::mat dirch_param(n_cats(r), K, arma::fill::ones); // Assumes Dirichelt prior parameter is 1
            xid[r].col(next_label) = rdirichlet(dirch_param.col(next_label));
        }
    }
    // if (arma::any(arma::any(u == arma::datum::inf))) {
    //     Rcout << "update_eta" << u << "\n" << eta;
    // }
}


void CovarSBM::update_z_element(const int i) {

    // int zi_old = z(i);
    arma::vec taui = col_compress(i); // sp_single_col_compress(A, i, z, K);
    // arma::uvec mmi = get_freq(z, K);
    // arma::uvec mmi = get_label_freq();
    // mmi(z(i))--;
    arma::vec mmi =  arma::conv_to<arma::vec>::from(
        get_label_freq_except(i)
    );

    // log_prob will be K x 1 >= (n_comm + 1) x 1
    arma::vec log_prob(K, arma::fill::zeros);

    // Rcout << "1\n";
    for (int k = 0; k < K; k++) {
        if (mmi(k) > 0 || k == next_label)
        {
            // Rcout << "2";
            // taui and mmi will be nonzero only for the correct positions,
            // hence the sum excludes irrelvant indices over the columns of u and v
            for (int l = 0; l < K; l++)
            {
                log_prob(k) += u(k,l) * taui(l) + v(k,l) * mmi(l);
            }
            // if (arma::any(log_prob == arma::datum::inf))
            // {
            //     Rcout << "->Inf\n";
            //     Rcout << log_prob.t() << "\n";
            //     Rcout << u << "\n";
            //     Rcout << v << "\n";
            // }
            // Rcout << "3";
            if (has_cont_features)
            {
                // Rcout << pow(arma::norm(X.row(i) - xi.row(k)), 2) / (2*s2) << "  ";
                log_prob(k) -= pow(arma::norm(Xc.row(i) - xic.row(k)), 2) / (2*s2);
            }
            // if (arma::any(log_prob == arma::datum::inf)) Rcout << "->Inf\n";
            // Rcout << log_prob.t() << "\n";
            // Rcout << "4\n";
            if (has_disc_features) {
                for (int r  = 0; r < feature_dim_disc; r++)
                {
                    log_prob(k) += log( xid[r](Xd(i,r), k) );
                }
            }

        }
        else
        {   // these indices are not part of the support; should get zero probability
            log_prob(k) = -arma::datum::inf;
        }
        if (mmi(k) > 0) // add log(psi_k)
        {
            log_prob(k) += log(mmi(k));
        }
    }

    // if (log_prob.has_nan()) Rcout << "->NaN\n";
    // Rcout << "5\n";
    // arma::vec log_prob = u * taui + v * mmi ; //+ log(pri);
    // log_prob += log(mmi); // will add zero at "next_label" position

    log_prob(next_label) += log(dp_concent_param);  // add log(psi_k) for k == next_label


    // Rcpp::print(wrap(log_prob));
    // Rcout << "\n";
    // Rcpp::print(wrap(safe_exp(log_prob)));

    // auto idx1 = arma::span(0, n_comm);
    // auto idx2 = arma::span(0, n_comm-1);

    // auto log_prob = arma::vec(n_comm+1, arma::fill::zeros);
    // for (int l = 0; l <= n_comm; l++) {
    //     if (l != next_label) {
    //         log_prob += u(idx1,l)*taui(l) + v(idx1,l)*mmi(l);
    //     }
    // }
    // arma::vec log_prob = u(idx1,idx2)*taui(idx2) + v(idx1,idx2)*mmi(idx2) ; //+ log(pri);
    // log_prob(idx2) += log( mmi(idx2) );
    // log_prob(max_label) += log(dp_concent_param);



    z(i) = sample_index(safe_exp(log_prob));

    // Rcout << "6\n";
    // Handle next_label creation
    if (z(i) == next_label) // We need to create a new empty cluster for futrue use
    {
        // Check for empty clusters other than "next_label"
        int ki = 0;
        for (ki = 0; ki < K; ki++)
        {
            if (mmi(ki) == 0 && ki != next_label)
            {
                break;
            }
        }

        if (ki < K)
        {
            // There is another empty cluster within our current range of
            // 0,..,K-1 >= n_comm. Resue it for next_label
            next_label = ki;
        }
        else
        {
            // No empty cluster anymore. We have to resize.
            K++;
            next_label = K-1;
            resize();
        }
        repopulate_from_prior(); // Repopulate eta, etc. over next_label row/column from the prior
    }
    // Rcout << "7\n";
}


arma::umat CovarSBM::run_gibbs(const int niter) {
    // Run full Gibbs updates for "niter" iterations and record label history

    arma::umat z_hist(n, niter+1);
    z_hist.col(0) = z + 1;

    for (int iter = 0; iter < niter; iter++)
    {
        update_eta();
        if (has_cont_features || has_disc_features)
        {
            update_xi();
        }
        for (int i = 0; i < n; i++)
        {
            update_z_element(i);
        }
        // Rcpp::print(wrap(z.t()));
        z_hist.col(iter+1) = z + 1;
    } // iter

    return z_hist;
    // return Rcpp::List::create(
    //     Rcpp::Named("z") = z_hist,
    //     Rcpp::Named("xi") = xi_hist
    // );
}

RCPP_MODULE(covar_sbm_module) {
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
    ;

    class_<CovarSBM>("CovarSBM")
        .derives<BasicSBM>("BasicSBM")
        .constructor<arma::sp_mat, int, double, double>()
        .field("eta", &CovarSBM::eta)
        .field("feature_dim", &CovarSBM::feature_dim)
        .field("Xc", &CovarSBM::Xc)
        .field("xic", &CovarSBM::xic)
        .field("Xd", &CovarSBM::Xd)
        .field("xid", &CovarSBM::xid)
        .field("u", &CovarSBM::u)
        .field("v", &CovarSBM::v)
        .field("n_cats", &CovarSBM::n_cats)
        // .field("max_label", &CovarSBM::max_label)
        .field("next_label", &CovarSBM::next_label)
        .field("has_cont_features", &CovarSBM::has_cont_features)
        .field("has_disc_features", &CovarSBM::has_disc_features)
        .field("s2", &CovarSBM::s2)
        .field("tau2", &CovarSBM::tau2)
        .method("update_eta", &CovarSBM::update_eta)
        .method("update_xi", &CovarSBM::update_xi)
        .method("update_z_element", &CovarSBM::update_z_element)
        .method("run_gibbs", &CovarSBM::run_gibbs)
        .method("sample_crp", &CovarSBM::sample_crp)
        .method("set_continuous_features", &CovarSBM::set_continuous_features)
        .method("set_discrete_features", &CovarSBM::set_discrete_features)
        .method("set_gauss_param", &CovarSBM::set_gauss_param)
        .method("set_beta_params", &CovarSBM::set_beta_params)
        .method("resize", &CovarSBM::resize)
        .method("repopulate_from_prior", &CovarSBM::repopulate_from_prior)
//      .method("print", &CovarSBM::print)
    ;
};

RCPP_EXPOSED_CLASS(CovarSBM);
