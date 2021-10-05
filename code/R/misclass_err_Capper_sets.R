# Misclassification error Capper model variations Training and Validation sets
rm(list=ls())
library(dplyr)
library(ggplot2)

source(file.path("code","R","type2upper.R"))

upper <- F

models <- c("Capper",
            "Capper32k",
            "Capper500trees",
            "Capper32k500trees")

val <- c("training", "test")

c_err <- c()
c_mod <- c()
c_val <- c()
c_cal <- c()
df_mis <- list()

# Datasets
for(v in val){
if(v=="test"){
  load(file.path("Results","annotations_test.RData"))
  ii <- which(anno$`methylation class:ch1`=="PIN T, PB B")
  anno$`methylation class:ch1`[ii] <- "PIN T,  PB B"
  y <- anno$`methylation class:ch1`
  i_melan <- which(grepl("MELAN", y))
  y[i_melan] <- "Melanoma brain metastasis"
}else{
  load(file.path("Results","annotations.RData"))
  load(file.path("Results","betas","betas_brainMetastasis.RData"))
  y <- c(anno$`methylation class:ch1`, y_met)
  i_melan <- which(grepl("MELAN", y))
  y[i_melan] <- "Melanoma brain metastasis"
  i_uncert <- which(grepl("Uncertain", y))
  y <- (y[-i_uncert])
}

  if(upper){
    y <- type2upper(y)
  }

# Capper models
  for(m in models){
  # m <- "Capper"
  load(file.path("Results","predictions",paste0("Capper_",m,"_",v),"predictionNewSamples.RData"))
  pred.class <- colnames(rf.scores)[apply(rf.scores,1,which.max)]
  pred.class.cal <- colnames(probs)[apply(probs,1,which.max)]
  
  
  if(upper){
    pred.class <- type2upper(pred.class)
    pred.class.cal <- type2upper(pred.class.cal)
  }
  
  # Classes misclassified and their incorrect prediction
  i_mis_raw <- (pred.class != y)
  pred_raw <- pred.class[i_mis_raw]
  y_raw <- y[i_mis_raw]
  id_raw <- rownames(rf.scores)[i_mis_raw]
  raw_score <- rf.scores[cbind(which(i_mis_raw), apply(rf.scores,1,which.max)[i_mis_raw])]
  
  i_mis_cal <- (pred.class.cal != y)
  pred_cal <- pred.class.cal[i_mis_cal]
  y_cal <- y[i_mis_cal]
  id_cal <-  rownames(probs)[i_mis_cal]
  cal_score <- probs[cbind(which(i_mis_cal), apply(probs,1,which.max)[i_mis_cal])]
  
  y_m <- c(y_raw, y_cal)
  pred_m <- c(pred_raw, pred_cal)
  score_m <- c(raw_score, cal_score)
  type <- c(rep("raw", length(y_raw)), rep("cal", length(y_cal)))
  mods <- rep(paste(m,v,sep = "."), length(y_cal) + length(y_raw))
  samp_id <- c(id_raw, id_cal)
  
  df <- data.frame(mods, type, y_m, pred_m, score_m, samp_id)
  
  if(is.null(dim(df_mis))){
    df_mis <- df
  }else{
    df_mis <- rbind(df_mis, df)
  }
  
  # Misclassification errors
  e_raw <- sum(pred.class != y) / length(pred.class)
  e_cal <- sum(pred.class.cal != y) / length(pred.class.cal)
  c_err <- c(c_err, e_raw, e_cal)
  c_mod <- c(c_mod, rep(m,2))
  c_val <- c(c_val, rep(v,2))
  c_cal <- c(c_cal, "raw","cal")

}
# }

plt_mis <- df_mis %>% group_by(y_m, mods, type, pred_m) %>% summarize(count = n()) %>% ungroup()

df_err <- data.frame(c_mod, c_val, c_cal, c_err)
colnames(df_err) <- c("model", "dataset", "type", "misclass_err")

save(df_err, df_mis, file=file.path("Results","predictions","misclass_models.RData"))
