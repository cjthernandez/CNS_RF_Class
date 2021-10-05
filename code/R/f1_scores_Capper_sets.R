# F1-score training and validation
rm(list=ls())
library(tibble)
library(dplyr)
library(stringr)

source("mnp_training-master/R/type2upper.R")

mods <- c("Capper", "Capper32k", "Capper500trees", "Capper32k500trees") # Brain_Metastasis
m <- "Brain_Metastasis"
sets <- c("training", "validation")

upper <- T

f1s <- c()
mm <- c()
ss <- c()
tt <- c()

for(m in mods){
  for(s in sets){
    
    if(s == "validation"){
      load(file.path("Results","annotations_validation.RData"))
      # load(file.path("Results","betas","betas.ba_validation_2021-03-15.RData"))
      y <- anno$`methylation class:ch1`
      ii <- which(y =="PIN T, PB B")
      y[ii] <- "PIN T,  PB B"
      # i_melan <- which(grepl("MELAN", y))
      # y[i_melan] <- "Melanoma brain metastasis"
      # rm(betas_val)
    }else{
      load(file.path("Results","annotations.RData"))
      # load(file.path("Results", "betas","betas.ba.brainMetastasis_offset.RData"))
      y <- c(anno$`methylation class:ch1`)#, y_met)
      # Remove "Uncertain" metastasis samples
      # i_uncert <- which(grepl("Uncertain", y))
      # y <- y[-i_uncert]
      # rm(betas)
      # i_melan <- which(grepl("MELAN", y))
      # y[i_melan] <- "Melanoma brain metastasis"
    }
    
    if(upper){
      y <- type2upper(y)
      load(file.path("Results","predictions",paste0("Capper_",m,"_",s),"confusion_matrix_upper.RData"))
    }else{
      load(file.path("Results","predictions",paste0("Capper_",m,"_",s),"confusion_matrix.RData"))
    }
    
    n_y_fold <- table(y)
    
    f1_cal <- as.data.frame(conf_mat_cal[[4]])
    rownames(f1_cal) <- sapply(strsplit(rownames(f1_cal),":"), function(x){str_trim(x[2])})
    f1_cal <- rownames_to_column(f1_cal)
    f1_raw <- as.data.frame(conf_mat_raw[[4]])
    rownames(f1_raw) <- sapply(strsplit(rownames(f1_raw),":"), function(x){str_trim(x[2])})
    f1_raw <- rownames_to_column(f1_raw)
    
    # F1 score correction when Precision and/or Recall is 0
    # from: https://github.com/dice-group/gerbil/wiki/Precision,-Recall-and-F1-measure
    f1_cal <- f1_cal %>% select(rowname, Precision, Recall, F1) %>%
      mutate(n_F1 = ifelse(is.na(Precision), ifelse(is.na(Recall), 1, 0),
                           ifelse(is.na(Recall), 0, ifelse(is.nan(F1), 0, F1)))) # Avoid NA when dividing by 0
    f1_raw <- f1_raw %>% select(rowname, Precision, Recall, F1) %>%
      mutate(n_F1 = ifelse(is.na(Precision), ifelse(is.na(Recall), 1, 0),
                           ifelse(is.na(Recall), 0, ifelse(is.nan(F1), 0, F1)))) # Avoid NA when dividing by 0
    
    # Matching F1 scores and n by class
    ii <- match(f1_cal$rowname, names(n_y_fold))
    f1_cal$n <- as.integer(n_y_fold[ii])
    
    f1 <- f1_cal$n_F1
    nn <- f1_cal$n
    f1 <- f1[which(!is.na(nn))]
    nn <- nn[which(!is.na(nn))]
    
    # Calculate scores
    f1_w_cal <- weighted.mean(f1, nn)
    
    #
    ii <- match(f1_raw$rowname, names(n_y_fold))
    f1_raw$n <- as.integer(n_y_fold[ii])
    
    f1 <- f1_raw$n_F1
    nn <- f1_raw$n
    f1 <- f1[which(!is.na(nn))]
    nn <- nn[which(!is.na(nn))]
    
    # Calculate scores
    f1_w_raw <- weighted.mean(f1, nn)
    
    f1s <- c(f1s, f1_w_cal, f1_w_raw)
    tt <- c(tt, "cal", "raw")
    mm <- c(mm, m, m)
    ss <- c(ss, s, s)
    
  }
}

df_f1 <- data.frame(f1_score=f1s, model=mm, dataset=ss, type=tt)

if(upper){
  fout<-"f1-score_ALL_models_upper.RData"
}else{
  fout<-"f1-score_ALL_models.RData"
}

save(df_f1, file=file.path("Results","predictions","f1-score_ALL_models.RData"))

