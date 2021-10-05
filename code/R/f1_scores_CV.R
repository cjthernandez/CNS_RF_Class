# F1 Scores from raw and calibrated scores
rm(list=ls())

source(file.path("data","R","type2upper.R"))
source(file.path("data","R","f1_score_fold.R"))

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
if(grepl("full", test)){
  perc <- 100
}

lst_f1 <- list()

# Random sampling replicates
for(repp in 1:nrep){
  message("Calculating replicate ", repp, " of ", nrep," ...", Sys.time())
  # For F1 scores from raw RF scores
  out_fold <- c()
  inn_fold <- c()
  p_per <- c()
  f1_w <- c()
  f1_m <- c()
  # For F1 scores from calibrated scores
  out_fold_cal <- c()
  inn_fold_cal <- c()
  p_per_cal <- c()
  f1_w_cal <- c()
  f1_m_cal <- c()
  
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
          pname <- paste("probsCVfold",K,0,"probes",pp,"RData", sep=".")###### Change file name if needed
          load(file.path("Results","CV", test, pname))###### Change directory if needed
          
          if(grepl("FixedFold", test)){
            scores <- probs[[repp]]
          }else{
            scores <- probs
          }
          
          # Get average and weighted average F1 scores from calibrated scores
          probs.fold <- f1_score_fold(scores, y, fold, upper)
          
          out_fold_cal <- c(out_fold_cal, K)
          inn_fold_cal <- c(inn_fold_cal, k)
          p_per_cal <- c(p_per_cal, pp)
          f1_m_cal <- c(f1_m_cal, probs.fold[1])
          f1_w_cal <- c(f1_w_cal, probs.fold[2])
        }
        
        # Access raw scores (CVfold rf.scores)
        fname <- paste("CVfold",K,k,"probes",pp,"scores","RData",sep=".")
        load(file.path("Results","CV",test,fname))
        
        if(grepl("FixedFold", test)){
          scores <- rf.scores[[repp]]
        }else{
          scores <- rf.scores
        }
        
        # Get average and weighted average F1 scores from raw scores
        scores.fold <- f1_score_fold(scores, y, fold, upper)
        
        out_fold <- c(out_fold, K)
        inn_fold <- c(inn_fold, k)
        p_per <- c(p_per, pp)
        f1_m <- c(f1_m, scores.fold[1])
        f1_w <- c(f1_w, scores.fold[2])
        #rm(rf.scores)
      }
    }
  }
  
  # Raw scores
  f1_scores <- data.frame(out_fold, inn_fold, p_per, f1_m, f1_w)
  # Calibrated scores
  f1_probs <- data.frame(out_fold_cal, inn_fold_cal, p_per_cal, f1_m_cal, f1_w_cal)
  colnames(f1_probs) <- colnames(f1_scores)
  f1_scores$type <- rep("raw", nrow(f1_scores))
  f1_probs$type <- rep("cal", nrow(f1_probs))
  # All scores
  f1_scores_all <- rbind(f1_scores, f1_probs)
  
  if(!grepl("FixedFold", test)){ # Save F1 scores dataframe
    if(upper){
      fout <- "f1_scores_all_upper.RData"
    }else{
      fout <- "f1_scores_all.RData"
    }
    save(f1_scores_all, file=file.path("Results","CV",test,fout))
  }else{
    lst_f1[[repp]] <- f1_scores_all
  }
}

if(grepl("FixedFold", test)){ # Save F1 scores list of dataframes for all replicates
  names(lst_f1) <- paste0("rep", 1:nrep)
  f1_scores_all <- lst_f1
  # Save scores from all replicates
  if(upper){
    fout <- "f1_scores_all_upper.RData"
  }else{
    fout <- "f1_scores_all.RData"
  }
  save(f1_scores_all,file=file.path("Results","CV",test,fout))
}

# }
