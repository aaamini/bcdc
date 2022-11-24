library(R.matlab)

train_ratio = 0.1
out = readMat(sprintf("Data/facebook-ego-%2.2f.mat", train_ratio))
B = out$B
X = out$feature
idx_train = out$idx.train
idx_test = out$idx.test
dim(B)


SymNARM(B,X, 10, idx_train, idx_test, Burnin=15, Collections=15) 
