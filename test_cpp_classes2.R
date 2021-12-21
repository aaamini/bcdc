# Test classes
# Rcpp::sourceCpp("src/BasicSBM.cpp", verbose = T)
# Rcpp::sourceCpp("src/SBM.cpp", verbose = T)
Rcpp::sourceCpp("src/CovarSBM.cpp", verbose = T)

set.seed(185)
n = 50 # try 125
lambda = 10 # SBM average expected degree 
K = 2
alpha = 1
beta = 1
dp_concent = 5
niter = 50 # number of iterations of the Gibbs sampler

z_tru = sample(1:K, n, replace = T)
B = nett::gen_rand_conn(n, K, lambda = lambda) #, gamma = 0.5)
A = nett::fast_sbm(z_tru, B)

mu = 2
Mu = rbind(c(mu, 0), c(-mu,0))
X = do.call(cbind, lapply(1:ncol(Mu), function(j) rnorm(n, mean = Mu[z_tru,j])))

# run the model with alpha = 1, beta = 1 and DP concent = 5
model = new(CovarSBM, A, alpha, beta, dp_concent)
zout = model$run_gibbs(niter) # 
nett::compute_mutual_info(z_tru, model$z+1)

# CovarSBM with node features
model = new(CovarSBM, A, alpha, beta, dp_concent)
model$set_node_features(X) 
zout = model$run_gibbs(niter) 
nett::compute_mutual_info(z_tru, model$z+1)

# model is an S4 class (it is an exported C++ class). You can view its content as follows:
str(model)

# ---- This is just FYI, not needed ---
# You can call various methods of the class
model = new(CovarSBM, A, alpha, beta, 1)
model$sample_crp()  # sets the labels to CRP
table(model$z)
model$update_eta()
model$eta
model$z
model$update_z_element(10)
model$get_label_freq()
table(model$z)

model$update_blk_sums_and_sizes()
model$M
model$N

model$repopulate_from_prior()


# model$col_compress(10)
