# Precision-Recall curves
rm(list=ls())

library(ROCR)
library(ggplot2)

source(file.path("code","R","type2upper.R"))
load(file.path("Results","annotations.RData"))

#### Which sampling test? ####
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

folds <-3

if(grepl("Metastasis", test)){
  load(file.path("Results","betas","betas_brainMetastasis.RData"))
  load(file.path("Results","CV","nfolds_brainMetastasis.RData"))
  # True labels
  y <- c(anno$`methylation class:ch1`, y_met)
  i_melan <- which(grepl("MELAN", y))
  y[i_melan] <- "Melanoma brain metastasis"
  i_uncert <- which(grepl("Uncertain", y))
  y <- y[-i_uncert]
}else{
  load(file.path("Results","CV","nfolds.RData"))
  # True labels
  y <- anno$`methylation class:ch1`
}

# Methylation families
if(upper){
  y <- type2upper(y)
}
# List of all classes or families
class_lst <- unique(y)

if(grepl("Capper", test)){
  per<- c(100)
}else{
  per <- c(10, 25, 33, 50, 66, 75, 100)
}

if(grepl("full", test)){
  per <- 100
}

if(grepl("FixedFold",test)){
  nrep <- 10
}else{
  nrep <- 1
}

c_auc <- c()
cal_auc <- c()
c_perc <- c()
c_rep <- c()
c_class <- c()
# c_fold <- c()

# Percentages
for(p in per){
  # Replicates
  for(r in 1:nrep){
    # Iterate over all tumor classes
    for(tc in 1:length(class_lst)){
      message("calculating PR-curve: Percent: ",p," class ",tc," of ",length(class_lst)," rep ",r," of ",nrep," ...", Sys.time())
      
      # Class probabilities
      cc <- list()
      cc_cal <- list()
      # True labels
      labb <- list()
      
      i <- 0 # folds counter
      
      # Iterate over all folds
      for(K in 1:folds){
        k <- 0 # Only PR curve from outer folds
        # Fold counter
        i <- i+1
        fold <- nfolds[[K]][[1]][[1]]
        
        pname <- paste("probsCVfold",K,k,"probes",p,"RData",sep=".")
        rname <- paste("CVfold",K,k,"probes",p,"scores","RData",sep=".")
        load(file.path("Results","CV",test,pname))#############
        load(file.path("Results","CV",test,rname))#############
        
        if(grepl("FixedFold", test)){
          scores <- rf.scores[[r]]
          probb <- probs[[r]]
        }else{
          scores <- rf.scores
          probb <- probs
        }
        
        # Get class probabilities for a single class over all samples (columns in rf.scores)
        if(upper){
          if(grepl("Metastasis", test)){
            c1 <- scores[,upper2type(class_lst[tc], met=T)]
            c1_cal <- probb[,upper2type(class_lst[tc], met=T)]
          }else{
            c1 <- scores[,upper2type(class_lst[tc])]
            c1_cal <- probb[,upper2type(class_lst[tc])]
          }
          # Get max score between tumor classes in broader tumor class ??can this be done?? OR SUM OF raw scores??
          if(is.null(dim(c1))){
            c2 <- c1
            c2_cal <- c1_cal
          }else{
            c2 <- apply(c1,1,max)
            c2_cal <- apply(c1_cal,1,max)
          }
        }
        else{
          c2 <- scores[,class_lst[tc]]
          c2_cal <- probb[,class_lst[tc]]
        }
        
        # Label the observations as positive and negative classes. REAL labels
        lab1 <- as.integer(class_lst[tc] == y[fold$test])
        
        cc[[i]] <- c2
        cc_cal[[i]] <- c2_cal
        labb[[i]] <- lab1
      }
      
      # Calculate AUC of PR curve using ROCR for raw and calibrated scores
      pred <- prediction(cc, labb)
      pred_cal <- prediction(cc_cal, labb)
      
      perf <- performance(pred, "prec","rec")
      perf_cal <- performance(pred_cal, "prec", "rec")
      
      au <- performance(pred, "aucpr")
      au_cal <- performance(pred_cal, "aucpr")
      auc <- unlist(au@y.values)
      auc_cal <- unlist(au_cal@y.values)
      
      # Store AUC from PR curve
      l_auc <- length(auc)
      c_perc <- c(c_perc, rep(p, l_auc))
      c_class <- c(c_class, rep(class_lst[tc],l_auc))
      c_rep <- c(c_rep, rep(r,l_auc))
      c_auc <- c(c_auc, auc)
      cal_auc <- c(cal_auc, auc_cal)
    }
  }
}

# AUC-PR from raw scores
df_raw_aucpr <- data.frame(factor(c_perc), factor(c_rep), factor(c_class), c_auc)
colnames(df_raw_aucpr) <- c("perc_f", "reps_f", "class_f", "auc")
df_raw_aucpr$type <- rep("raw", nrow(df_raw_aucpr))
# AUC-PR from calibrated scores
df_cal_aucpr <- data.frame(factor(c_perc), factor(c_rep), factor(c_class), cal_auc)
colnames(df_cal_aucpr) <- c("perc_f", "reps_f", "class_f", "auc")
df_cal_aucpr$type <- rep("cal", nrow(df_cal_aucpr))

df_aucpr <- rbind(df_raw_aucpr, df_cal_aucpr)

if(upper){
  aucname <- paste("PR_AUC_upper", p,"RData", sep=".")
}else{
  aucname <- paste("PR_AUC", p,"RData", sep=".")
}

save(df_aucpr, file=file.path("Results","CV",test,aucname))#########
