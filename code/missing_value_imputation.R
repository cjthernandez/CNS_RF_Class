#-----------------------------------------------------------------------------------
# this script performs missing value imputation of data from Nanopore methylation
# calls from IntraEpiGliom samples. Three imputation methods were implemented: 
# naive (replace NAs with 0.5)
# mean (replace NAs with the mean across all non-missing probes)
# kNN mean (replace NAs with mean of k non-missing approximate NN probes)
#
#------------------------------------------------------------------------------------ 
# Impute NA in each coverage group
rm(list=ls())
library(dplyr)

source(file.path("code","R","na_impute.R"))

impute <- "mean" # c(naive, mean, kneighbor)
dir.create(file.path("Results",paste(impute,"imputation",sep="_")))
kn <- 100

# Random Forest used as reference for imputation #### CHANGE NAME ACCORDING TO RF TRAINED
rf_mod <- "Brain_Metastasis" # "Capper" "32k" "32k500trees" "500trees" "Brain_Metastasis"
load(file.path("Results",paste0("rf.pred_",rf_mod,".RData")))
probes_rf <- rownames(rf.pred$importance)

# Which distance matrix?
if(impute=="kneighbor"){
eps <- 6
k <- 1000
load(file.path("Results",paste0("ANN_eps", eps, "_k", k,".RData"))) ### Distance matrix

dist_mat <- betas_ann
}

# Load list of methylation profiles
load(file=file.path("Results","betas_new_samples.RData")) # lst_betas

imp_dim = "cpg" # "sample"   "cpg"

run_time <- c()
cov_time <- c()
na_time <- c()

lst_betas_imp <- list()

names_lst <- names(lst_betas)

for(w in 1:length(lst_betas)){
  new_nanop <- lst_betas[[w]]
  
  # In case there are probes from the RF not included in Nanopore
  n_not_rf <- sum(is.na(match(probes_rf, names(new_nanop))))
  
  if(n_not_rf!=0){
  na_lst <- rep(NA, n_not_rf)
  names(na_lst) <- probes_rf[which(is.na(match(probes_rf, names(new_nanop))))]
  
  new_nanop <- c(new_nanop, na_lst)
  }
  
  message("Imputing sample ",w," of ",length(lst_betas), " ...", Sys.time())
    if(impute=="kneighbor"){
      new_samp <- na_impute_nanopore(new_nanop, dist_mat=dist_mat, probes_rf=probes_rf, impute=impute, imp_dim=imp_dim, k=kn)# [jj,]
    }else{
      new_samp <- na_impute_nanopore(new_nanop, probes_rf=probes_rf, impute=impute, imp_dim=imp_dim, k=kn)# [jj,]
    }
  
  lst_betas_imp[[w]] <- new_samp
    rm(new_samp)
}

names(lst_betas_imp) <- names_lst

save(lst_betas_imp, file=file.path("Results", paste(impute,"imputation",sep="_"),
                                   paste0("betas_new_samples_imp_",rf_mod,"_",impute,"_",kn,"_",imp_dim,".RData")))
