
Rcpp::sourceCpp("src/SBM.cpp", verbose = T)
Rcpp::sourceCpp("src/CovarSBM.cpp", verbose = T)
source("R/inference.R")

library(ggplot2)
theme_set(theme_minimal(base_size = 12))
# theme_set(theme_minimal())
library(dplyr)
library(tictoc)
library(parallel)
library(foreach)
library(patchwork)
# library(doParallel)
library(pbmcapply)

# set.seed(185)
n = 150 # try 125
# lambda = 10 # SBM average expected degree 
K = 2
alpha = 1
beta = 1
dp_concent = 1

niter = 100
nreps = 32
n_cores = 32

p = 0.1

methods = list()

methods[["SBM"]] =  function(A, X) {
   model = new(SBM, A, K, alpha, beta)
   get_map_labels(model$run_gibbs(niter))$z_map
  #  model$z
}

methods[["k-means"]] =  function(A, X) {
   kmeans(X, K, nstart = 20)$cluster
}

methods[["SC"]] =  function(A, X) {
   nett::spec_clust(A, K)
}


methods[["BCDC"]] =  function(A, X) {
  model = new(CovarSBM, A, alpha, beta, dp_concent)
  # model$set_gauss_param(1, 1)
  # model$set_node_features(X)
  model$set_discrete_features(X)
  get_map_labels(model$run_gibbs(niter))$z_map
  # model$z
}

mtd_names = names(methods)

runs = expand.grid(r = seq(0.1, .8, by = 0.1), rep = 1:nreps)

# cl <- parallel::makeForkCluster(n_cores)
# registerDoParallel(cl)

# res = do.call(rbind, lapply(1:nrow(runs), function(ri) {
# res = do.call(rbind, mclapply(1:nrow(runs), function(ri) {
res = do.call(rbind, pbmclapply(1:nrow(runs), mc.cores = n_cores, FUN = function(ri) {
# res = do.call(rbind, 
# foreach(ri = 1:nrow(runs)) %dopar% {
  rep = runs[ri,"rep"]
  r = runs[ri,"r"]

  eta = matrix(c(p, r*p, r*p, p), nrow=2)
  z_tru = sample(1:K, n, replace = T, prob = c(1,2))
  # eta = nett::gen_rand_conn(n, K, lambda = lambda) #, gamma = 0.5)
  A = nett::fast_sbm(z_tru, eta)

  theta = cbind(rdirichlet(rep(1,4)), rdirichlet(rep(1,4)))
  # theta = cbind(c(1,0,0,0),c(0,0,1,0)) # for test
  X = cbind(
    apply(theta[, z_tru], 2, function(th) sample(1:4, 1, prob = th)), 
    sample(1:4, n, T)
  ) 
  
  do.call(rbind, lapply(seq_along(methods), function(j) { 
    tic()
    zh <- as.vector(methods[[j]](A, X)) + 1
    dt = toc(quiet = T)
    data.frame(
        dt = as.numeric(dt$toc - dt$tic), 
        rep = rep,
        r = r,
        # nmi = nett::compute_mutual_info(zh, z_tru),
        nmi = nmi_wrapper(zh, z_tru), 
        method = mtd_names[j])
  }))
#})
# }, mc.cores = n_cores))
}))    

# arallel::stopCluster(cl)  

plt1 = res %>% 
  # group_by(method, mu) %>% summarise(nmi = mean(nmi)) %>% 
  ggplot(aes(x = as.factor(r), y = nmi, fill = method)) + 
  geom_boxplot() +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.9, 0.2),
    # legend.text = ggplot2::element_text(size=18),
  ) + 
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  ylab("NMI") + xlab("r") +
   labs(title = sprintf("n = %d,  p = %2.1f", n, p))

plt2 = res %>% 
  group_by(method, r) %>% summarise(nmi = mean(nmi)) %>% 
  ggplot(aes(x = r, y = nmi, color = method)) + 
  geom_line(size = 1.2) +
  ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = c(0.9, 0.2),
    # legend.text = ggplot2::element_text(size=18),
  ) + 
  ggplot2::guides(colour = ggplot2::guide_legend(keywidth = 2, keyheight = .75)) +
  ylab("NMI") + xlab("r") 
  # labs(title = sprintf("n = %d,  p = %2.1f,  r = %2.1f", n, p, r))

print(plt1 + plt2)
ggsave(sprintf("dis_n=%d_p=%2.1f.png", n, p), width = 10, height = 5)


res %>% 
  group_by(method) %>% 
  summarise(dt = mean(dt)) %>% 
  knitr::kable() %>% 
  print()
