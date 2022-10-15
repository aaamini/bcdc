#include <RcppArmadillo.h>

#include "CovarSBM.h"
#include "SBM.h"

RCPP_MODULE(sbm_module) {
  Rcpp::class_<BasicSBM>("BasicSBM")
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

  Rcpp::class_<SBM>("SBM")
        .derives<BasicSBM>("BasicSBM")
        .constructor<arma::sp_mat, int, double, double>()
        .field("eta", &SBM::eta)
        .method("update_eta", &SBM::update_eta)
        .method("set_beta_params", &SBM::set_beta_params)
    //      .method("print", &SBM::print)
        .method("update_z_element", &SBM::update_z_element)
        .method("run_gibbs", &SBM::run_gibbs)
        .method("update_eta", &SBM::update_eta)
        .method("update_pri", &SBM::update_pri)
    ;

  Rcpp::class_<CovarSBM>("CovarSBM")
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

// RCPP_EXPOSED_CLASS(CovarSBM);
