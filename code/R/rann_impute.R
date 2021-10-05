# Calculate k Approximate Nearest Neighbor (k-ANN) matrix for Imputation
rm(list=ls())
library(RANN)
library(dplyr)

# Methylation rates
load(file.path("Results","betas_ba.RData"))

betas_t <- t(betas)

k <- 500 # Number of nearest neighbors to use
eps <- 6 # Approximation error bound

# RANN
betas_ann <- nn2(betas_t, k = k, treetype = "bd", searchtype = "priority", eps = eps)

rownames(betas_ann$nn.dists) <- rownames(betas_t)
rownames(betas_ann$nn.idx) <- rownames(betas_t)

save(betas_ann, file=file.path("Results", paste0("ANN_eps", eps, "_k", k,".RData")))
