
# Rcpp::sourceCpp("src/BasicSBM.cpp", verbose = T)
Rcpp::sourceCpp("src/SBM.cpp", verbose = T)
Rcpp::sourceCpp("src/CovarSBM.cpp", verbose = T)


library(ggplot2)
# theme_set(theme_minimal(base_size = 18))
theme_set(theme_minimal())
library(dplyr)
library(tictoc)

set.seed(185)
n = 50 # try 125
lambda = 10 # SBM average expected degree 
K = 2
alpha = 1
beta = 1
Zcap = 5  # increasing this slows down the collapsed one significantly

niter = 50
nreps = 50
nreps_per_net = 3 # Try also 5
include_dpsbm_flag = T

mu = 3
Mu = rbind(c(mu, 0), c(-mu,0))

comp_agg_nmi_path = function(z_list) {
   sapply(seq_along(z_list), function(it) hsbm::get_agg_nmi(z_list[[it]], list(z_tru))) 
}

convert_cpp_label_matrix_to_list = function(zmat) {
  lapply(1:ncol(zmat), function(it) list(zmat[,it]+1))  # the list() is there to treat these  multi-layer labels with only a single layer 
}

methods = list()

# Rcpp::sourceCpp("src/multsbm.cpp", verbose = T)
# methods[["SBM-fun"]] =  function(A) {
#   multsbm_gibbs_sampler_fast(A, K, alpha, beta = beta, niter = niter)
# }

methods[["SBM-class"]] =  function(A, X) {
   model = new(SBM, A, K, alpha, beta)
   model$run_gibbs(niter)[,-1]
}


methods[["SBM-CRP (1)"]] =  function(A, X) {
  model = new(CovarSBM, A, alpha, beta, 1)
  model$run_gibbs(niter)[,-1]
}

methods[["SBM-CRP (10)"]] =  function(A, X) {
  model = new(CovarSBM, A, alpha, beta, 10)
  model$run_gibbs(niter)[,-1]
}


methods[["CovarSBM (1)"]] =  function(A, X) {
  model = new(CovarSBM, A, alpha, beta, 1)
  # model$set_gauss_param(1, 1)
  model$set_node_features(X)
  model$run_gibbs(niter)[,-1]
}

methods[["CovarSBM (10)"]] =  function(A, X) {
  model = new(CovarSBM, A, alpha, beta, 10)
  # model$set_gauss_param(1, 1)
  model$set_node_features(X)
  model$run_gibbs(niter)[,-1]
}


# z_hist = methods[["SBM-CRP"]](A)
# apply(z_hist, 2, function(z) nett::compute_mutual_info(z, z_tru))
# cbind()


if (include_dpsbm_flag) {
  # DP-SBM
  Rcpp::sourceCpp("src/dpsbm.cpp", verbose = T)

  # methods[["DP-SBM-collapsed-v1"]] =  function(A) {
  #   fit_dpsbm_collapsed(A, Zcap = Zcap, niter = niter)
  # }

  methods[["DP-SBM-trunc-fun"]] =  function(A, X) {
    fit_dpsbm(A, Zcap = Zcap, niter = niter, gam0 = 1)
  }
}


mtd_names = names(methods)

res = NULL
for (rep in 1:nreps) {
  if (!rep %% 5) cat('.')
  z_tru = sample(1:K, n, replace = T)
  B = nett::gen_rand_conn(n, K, lambda = lambda) #, gamma = 0.5)
  A = nett::fast_sbm(z_tru, B)
  
  X = do.call(cbind, lapply(1:ncol(Mu), function(j) rnorm(n, mean = Mu[z_tru,j])))
  # diag(A) = sample(0:1, n, replace = T)  # This will destroy the code! -- comp_blk_sums_and_sizes assumes that the diagonal of A is zero
  
  for (j in 1:nreps_per_net) {
    res_curr = do.call(rbind, lapply(seq_along(methods), function(j) {
      # dt = system.time( z_list <- convert_cpp_label_matrix_to_list( methods[[j]](A) ) )["elapsed"]
      tic()
      z_list <- convert_cpp_label_matrix_to_list( methods[[j]](A, X) )
      dt = toc(quiet = T)
       data.frame(iter = 1:niter,
          dt = as.numeric(dt$toc - dt$tic), 
          rep = rep,
          rep_per_net = j,
          nmi = comp_agg_nmi_path(z_list),
          method = mtd_names[j])
    }))

    res = rbind(res, res_curr)
  }
}


p = res %>% 
  group_by(iter, method) %>% summarise(nmi = mean(nmi)) %>% 
  ggplot(aes(x = iter, y = nmi, color = method)) + 
  geom_line(size = 1.2) +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.8, 0.2),
    # legend.text = ggplot2::element_text(size=18),
  ) + 
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  ylab("Average NMI") + xlab("Iteration") +
  labs(title = sprintf("n = %d,  K = %d,  lam = %2.1f, mu = %2.1f", n, K, lambda, mu))

print(p)
ggsave(sprintf("covar_sbm_n=%d_K=%d_lam=%d_mu=%d.png", n, K, round(lambda), round(mu)),   width = 5, height = 4)


res %>% group_by(method) %>% summarise(avg_time = mean(dt))

# Tests 
# hsbm::seq_nmi_plot(z_list)
# comp_beta_matrix(A, z-1, K, alpha1, beta1)
# comp_Bet(A, z, K, alpha1, beta1)

# # Speed test
# microbenchmark::microbenchmark(regular =  fit_dpsbm(A, Zcap = Zcap, niter = niter),
#                                 collapsed =  fit_dpsbm_collapsed(A, Zcap = Zcap, niter = niter), times = 20)
