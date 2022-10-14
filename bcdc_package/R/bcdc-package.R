## usethis namespace: start
#' @useDynLib bcdc, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
# NULL
Rcpp::loadModule("basic_sbm_module", TRUE)
Rcpp::loadModule("sbm_module", TRUE)
Rcpp::loadModule("covar_sbm_module", TRUE)

