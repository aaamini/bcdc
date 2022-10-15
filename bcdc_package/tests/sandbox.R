set.seed(12)

# Uses some function from https://github.com/aaamini/nett package
# Sample from a SBM
n <- 150
K <- 2
z_tru <- sample(1:K, n, replace = TRUE, prob = c(1, 2)) # True labels
eta <- 0.1*matrix(c(1, .3, .3, 1), nrow = 2) # SBM connectivity matrix
A <- nett::fast_sbm(z_tru, eta)

# Create covariates
Mu <- rbind(c(1.5, 0), c(-1.5,0))
X <- do.call(cbind, lapply(1:ncol(Mu), function(j) rnorm(n, mean = Mu[z_tru,j])))

alpha <- beta <- 1
dp_concent <- 10
bcdc_model <- new(CovarSBM, A, alpha, beta, dp_concent) # Create the bcdc_model model as an S4 object

# Inspect the model object; you also write bcdc_model$ followed by [TAB] in Rstudio to
# get code completion suggestions for the list to properties and methods the
# object has.
str(bcdc_model)

# add continuous features to the model. For discrete features use
# bcdc_model$set_discrete_features()
bcdc_model$set_continuous_features(X)

# Runs the Gibbs sampler for "niter" iterations
niter = 50
zout = bcdc_model$run_gibbs(niter) # zout is an n x (10+1) matrix. The jth column is the label vector at the (j-1)th iteration

nmi_seq = apply(zout, 2, function(z) nett::compute_mutual_info(z,z_tru))
plot(1:(niter+1), nmi_seq,
     ylab = "NMI", xlab = "iteration", pch= 16)

# # Run the same chain for "niter" more iterations; the state is persistent
# zout = bcdc_model$run_gibbs(niter)
# plot(1:(niter+1), apply(zout, 2, function(z) nett::compute_mutual_info(z,z_tru)),
#      ylab = "NMI", xlab = "iteration", pch= 16)

bcdc_model$z # current model labels

# sbm_model <- new(SBM, A, K, alpha, beta)
# str(sbm_model)
# zout = sbm_model$run_gibbs(niter) # z
# nmi_seq = apply(zout, 2, function(z) nett::compute_mutual_info(z,z_tru))
# plot(1:(niter+1), nmi_seq,
#      ylab = "NMI", xlab = "iteration", pch= 16)

