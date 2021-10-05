#-----------------------------------------------------------------------------------
# nested cross-validation 
#                                                                     
# Martin Sill
# m.sill@dkfz.de                                                                  
# 
# 2018-03-14 UTC
#
# Modified by: Camilo Hernandez-Toro
#
#------------------------------------------------------------------------------------                   
#options(max.print = 1000)
#options(stringsAsFactors = FALSE)
#options(scipen = 999)

rm(list=ls())

library(randomForest)
library(parallel)
library(minfi)
library(limma)
library(lubridate)

source(file.path("code","R","makefolds.R"))
source(file.path("code","R","train.R"))
source(file.path("code","R","calculateCVfold.R"))
source(file.path("code","R","batchadjust.R"))

ntrees <- 500 # In paper 10000
cores <- 4
seed <- 180314
p <- 10000 # Number of probes to be used
folds <- 3

#### Which CV/down-sampling test? ####
tests <- c("Capper", # 10k trees, no SD filter, importance selection
           "Capper32k", # 10k trees, SD filter 32k, importance selection
           "Capper500trees", # 500 trees, no SD filter, importance selection
           "Capper32k500trees", # 500 trees, SD filter, importance selection
           "probe_randsamp_FixedFold", # 500 trees, random selection
           "probe_randsamp_FixedFold_StdDev", # 500 trees, SD filter, random selection
           "brain_Metastasis_FixedFold_offset", # 500 trees, random selection
           "brain_Metastasis_full_offset") # 1000 trees, random selection

test <- tests[1]

if(grepl("Capper", test)){
  perc <- 1
  ntrees <- 10000
}else{
  perc <- c(1/10, 1/4, 1/3, 1/2, 2/3, 3/4, 1) # Varying percentage of probes to be selected in classifier
}
if(grepl("500trees",test)){
  ntrees <- 500
}

if(grepl("Metastasis_full", test)){
  perc <- 1
  ntrees <- 1000
}

message("loading filtered Mset ...",Sys.time())

if(grepl("Metastasis", test)){
  load(file.path("Results", "betas_ba_brainMetastasis.RData"))
  y <- c(anno$`methylation class:ch1`, y_met)
  # Remove "Uncertain" metastasis samples
  i_uncert <- which(grepl("Uncertain", y))
  y <- y[-i_uncert]
  betas <- betas[-i_uncert,]
  # Merge MELAN and Melanoma brain Metastasis
  i_melan <- which(grepl("MELAN", y))
  y[i_melan] <- "Melanoma brain metastasis"
  y <- as.factor(y)
  folds_name <- "nfolds_brainMetastasis.RData"
}else{
  load(file.path("Results","Mset_filtered.RData")) ##
  y <- as.factor(anno$`methylation class:ch1`)
  batch <- as.factor(anno$`material:ch1`)
  folds_name <- "nfolds.RData"
}

if(!file.exists(file.path("Results","CV",folds_name))){
  dir.create(file.path("Results","CV"),showWarnings = FALSE) 
  nfolds <- makenestedfolds(y,folds)
  save(nfolds,file=file.path("Results","CV",folds_name))
}
# Directory for testing probe number effects on classifier
dir.create(file.path("Results", "CV", test), showWarnings = FALSE)
load(file.path("Results","CV",folds_name))

message("performing nested CV ...", Sys.time())
message("check minimal class sizes for inner training loops")

# check minimal class sizes for inner training loops
minclasssize <- matrix(0,ncol=length(nfolds),nrow=length(nfolds))
for(i in 1:length(nfolds)){
  for(j in 1:length(nfolds))
    minclasssize[i,j]  <- min(table(y[nfolds[[i]][[2]][[j]]$train]))
}
colnames(minclasssize) <- paste0("innfold",1:folds)
rownames(minclasssize) <- paste0("fold",1:folds)
print(minclasssize)

##### Number of replicates. FixedFold tests.
if(grepl("FixedFold", test)){
  nrep <- 10
}else{
  nrep <- 1
}

# Outer fold
for(K in 1:folds){
  # Inner fold
  for(k in 0:folds){
    # For FixedFold replicate tests
    missclass_perc <- list()
    i <- 0
    
    if(k>0){  
      message("calculating fold ",K,".",k,"  ...",Sys.time())
      fold <- nfolds[[K]][[2]][[k]]
    }else{
      message("calculating outer fold ",K,"  ...",Sys.time())
      fold <- nfolds[[K]][[1]][[1]]
    }
    
    # Batch adjustment or betas (brain Metastasis)
    if(grepl("Metastasis", test)){
      badj <- list(betas.train=betas[fold$train,],betas.test=betas[fold$test,])
    }else{
    message("adjusting for batch effects ...",Sys.time())
    badj <- batchadjust(Mset_filtered, batch, fold)
    }
    
    # Percentages
    for(pp in perc){
      i <- i+1
      message("CV classifier - ", test, " ...", Sys.time())
      
      # Calculate CV fold with pp% of p probes
      rf.fold <- calcultateCVfold(badj, y, fold, floor(p*pp), cores, ntrees, test, nreps = nrep)
      
      if(grepl("FixedFold", test)){
        # f1_scores_perc[[i]] <- rf.fold[[1]]
        missclass_perc[[i]] <- rf.fold[[4]]
      }
      
      rf.scores <- rf.fold[[1]]
      rf.feat <- rf.fold[[2]]
      rf.imp <- rf.fold[[3]]
      
      fname <- paste("CVfold",K,k,"probes",floor(pp*100),"scores","RData",sep=".") # class probabilities
      save(rf.scores,file=file.path("Results","CV", test, fname))
      fname <- paste("CVfold",K,k,"probes",floor(pp*100),"feat","RData",sep=".") # features used in RF
      save(rf.feat,file=file.path("Results","CV", test, fname))
      fname <- paste("CVfold",K,k,"probes",floor(pp*100),"imp","RData",sep=".") # feature importance
      save(rf.imp,file=file.path("Results","CV", test, fname))
      rm(rf.scores,rf.feat,rf.imp,rf.fold)
    }
    
    # Saving files if FixedFold replicate test
    if(grepl("FixedFold", test)){
      names(missclass_perc) <- paste0("perc",round(perc,2))
      fname2 <- paste("CVfold",K,k,"probes","nrep",nrep,"missclasserr","RData",sep=".")
      save(missclass_perc,file=file.path("Results","CV", test, fname2))
    }
    # gc()
  }
}
message("finished ...",Sys.time())
