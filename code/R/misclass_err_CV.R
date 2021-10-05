# Classification error from raw and calibrated scores
rm(list=ls())

source(file.path("mnp_training-master","R","type2upper.R"))

### Which sampling test? ###
tests <- c("Capper", # 10k trees, no SD filter, importance selection
           "Capper32k", # 10k trees, SD filter 32k, importance selection
           "Capper500trees", # 500 trees, no SD filter, importance selection
           "Capper32k500trees",
           "probe_randsamp_FixedFold",
           "probe_randsamp_FixedFold_StdDev",
           "brain_Metastasis_FixedFold_offset",
           "brain_Metastasis_full_offset")

test <- tests[1]

upper <- T

load(file.path("Results","annotations.RData"))

# Get true labels
if(grepl("Metastasis",test)){
  load(file.path("Results","betas","betas_brainMetastasis.RData"))
  load(file.path("Results","CV","nfolds_brainMetastasis.RData"))
  y <- c(anno$`methylation class:ch1`, y_met)
  i_melan <- which(grepl("MELAN", y))
  y[i_melan] <- "Melanoma brain metastasis"
  i_uncert <- which(grepl("Uncertain", y))
  y <- as.factor(y[-i_uncert])  
}else{
  load(file.path("Results","CV","nfolds.RData"))
  y <- as.factor(anno$`methylation class:ch1`)
  # names(y) <- anno$geo_accession
}

# Get samples present in current fold
folds <- 3

if(grepl("FixedFold",test)){
  nrep <- 10
}else{
  nrep <- 1
}
if(grepl("Capper", test)){
  perc <- 100
}else{
  perc <- c(10, 25, 33, 50, 66, 75, 100)
}

if(grepl("full",test)){
  perc <- 100
}

lst_err <- list()

# Random sampling replicates
for(repp in 1:nrep){
  message("Calculating replicate ", repp, " of ", nrep," ...", Sys.time())
  # For misclass error from raw RF scores
  out_fold <- c()
  inn_fold <- c()
  p_per <- c()
  err <- c()
  # For misclass error from calibrated scores
  out_fold_cal <- c()
  inn_fold_cal <- c()
  p_per_cal <- c()
  err_cal <- c()
  
  # Outer fold
  for(K in 1:folds){
    # Inner fold
    for(k in 0:folds){
      message("Calculating fold ",K,".",k," ...",Sys.time())
      if(k>0){  
        fold <- nfolds[[K]][[2]][[k]]
      }else{
        fold <- nfolds[[K]][[1]][[1]]
      }
      # Percentages
      for(pp in perc){
        message("Percent ",pp, " ...", Sys.time())
        
        # Access calibrated scores (CVprobs)
        if(k==0){
          pname <- paste("probsCVfold",K,0,"probes",pp,"RData", sep=".")
          load(file.path("Results","CV", test, pname))
          
          if(grepl("FixedFold", test)){
            scores <- probs[[repp]]
          }else{
            scores <- probs
          }
          
          # Get misclass error from calibrated scores
          if(upper){
            e <- sum(type2upper(colnames(scores)[apply(scores,1,which.max)]) != type2upper(y[fold$test])) / length(fold$test)
          }else{
            e <- sum(colnames(scores)[apply(scores,1,which.max)] != y[fold$test]) / length(fold$test)
          }
          
          out_fold_cal <- c(out_fold_cal, K)
          inn_fold_cal <- c(inn_fold_cal, k)
          p_per_cal <- c(p_per_cal, pp)
          err_cal <- c(err_cal, e)
        }
        
        # Access raw scores (CVfold rf.scores)
        fname <- paste("CVfold",K,k,"probes",pp,"scores","RData",sep=".")
        load(file.path("Results","CV",test,fname))
        
        if(grepl("FixedFold", test)){
          scores <- rf.scores[[repp]]
        }else{
          scores <- rf.scores
        }
        
        # Get misclass error from raw scores
        if(upper){
          e <- sum(type2upper(colnames(scores)[apply(scores,1,which.max)]) != type2upper(y[fold$test])) / length(fold$test)
        }else{
          e <- sum(colnames(scores)[apply(scores,1,which.max)] != y[fold$test]) / length(fold$test)
        }
        
        out_fold <- c(out_fold, K)
        inn_fold <- c(inn_fold, k)
        p_per <- c(p_per, pp)
        err <- c(err, e)
      }
    }
  }
  
  # Raw scores
  err_scores <- data.frame(out_fold, inn_fold, p_per, err)
  # Calibrated scores
  err_probs <- data.frame(out_fold_cal, inn_fold_cal, p_per_cal, err_cal)
  colnames(err_probs) <- colnames(err_scores)
  err_scores$type <- rep("raw", nrow(err_scores))
  err_probs$type <- rep("cal", nrow(err_probs))
  # All scores
  err_scores_all <- rbind(err_scores, err_probs)
  
  if(!grepl("FixedFold", test)){ # Save misclass error dataframe
    if(upper){
      fout <- "missclass_err_all_upper.RData"
    }else{
    fout <- "missclass_err_all.RData"
    }
    save(err_scores_all, file=file.path("Results","CV",test,fout))
  }else{
    lst_err[[repp]] <- err_scores_all
  }
}

if(grepl("FixedFold", test)){ # Save misclass error list of dataframes for all replicates
  names(lst_err) <- paste0("rep", 1:nrep)
  err_scores_all <- lst_err
  # Save scores from all replicates
  if(upper){
    fout <- "missclass_err_all_upper.RData"
  }else{
    fout <- "missclass_err_all.RData"
  }
  save(err_scores_all,file=file.path("Results","CV",test,fout))
}

# }

