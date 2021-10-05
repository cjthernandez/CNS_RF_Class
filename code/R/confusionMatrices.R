# Confusion matrix for training and validation sets
rm(list=ls())
library(caret)

source(file.path("code","R","type2upper.R"))

val <- c("training", "test") # datasets
mods <- c("Capper","Capper32k", "Capper32k500trees", "Capper500trees") # RF models

upper <- F # Methylation class or Methylation family

v <- "test" # "test" "training"

# levels from training set
load(file.path("Results","annotations.RData"))
# load(file.path("Results", "betas_brainMetastasis.RData")) # Uncomment fro Brain Metastasis set
y <- c(anno$`methylation class:ch1`)#, y_met)
## Remove "Uncertain" metastasis samples
# i_uncert <- which(grepl("Uncertain", y))
# y <- y[-i_uncert]
# rm(betas)
# i_melan <- which(grepl("MELAN", y))
# y[i_melan] <- "Melanoma brain metastasis"
# y <- as.factor(y)

y_level <- y

## Uncomment for test set
# if(v=="test"){
#   load(file.path("Results","annotations_test.RData"))
#   # load(file.path("Results","betas","betas_ba_test.RData"))
#   y <- anno$`methylation class:ch1`
#   ii <- which(y =="PIN T, PB B")
#   y[ii] <- "PIN T,  PB B"
#   # i_melan <- which(grepl("MELAN", y))
#   # y[i_melan] <- "Melanoma brain metastasis"
#   # rm(betas_val)
# }
fout <- "confusion_matrix.RData"
if(upper){
  y <- type2upper(y)
  y_level <- type2upper(y_level)
}
y_level <- factor(y_level)
y <- factor(y, levels=levels(y_level))

for(m in mods){
  load(file.path("Results","predictions",paste0("Capper_",m,"_",v), "predictionNewSamples.RData")) # Suggested directory nomenclature
  if(!is.list(rf.scores)){
    pred.class <- colnames(rf.scores)[apply(rf.scores, 1, which.max)]
    pred.class.cal <- colnames(probs)[apply(probs, 1, which.max)]
  }else{
    pred.class <- names(rf.scores[[1]])[sapply(rf.scores, which.max, simplify = "array")]
    pred.class.cal <- names(probs[[1]])[sapply(probs, which.max, simplify="array")]
  }
  if(upper){
    pred.class <- type2upper(pred.class)
    pred.class.cal <- type2upper(pred.class.cal)
  }
  
  y_raw <- factor(pred.class, levels=levels(y_level)) # 
  y_cal <- factor(pred.class.cal, levels=levels(y_level)) # 
  
  # Confusion Matrices
  conf_mat_raw <- confusionMatrix(y_raw, y, mode="prec_recall")
  conf_mat_cal <- confusionMatrix(y_cal, y, mode="prec_recall")
  
  save(conf_mat_cal, conf_mat_raw, file=file.path("Results","predictions",paste0("Capper_",m,"_",v),fout))
}

#### Plotting confusion matrices ####
library(dplyr)
library(ggplot2)
library(tidyr)

mods <- c("Capper","Capper32k", "Capper32k500trees", "Capper500trees") # Brain_Metastasis

set <- c("training", "validation")

upper <- F

if(upper){
  fout <- "confusion_matrix_upper.RData"
}else{
  fout <- "confusion_matrix.RData"
}

m <- "Brain_metastasis"

for(m in mods){
  for(s in set){
    
    load(file.path("Results","predictions",paste0("Capper_",m,"_",s),fout))
    
    ## reshape data (tidy/tall form)
    dat2 <- conf_mat_cal$table %>% as_tibble() %>% group_by(Reference) %>% mutate(n_norm = n/sum(n)) %>% ungroup()
    
    plot_conf_mat <- ggplot(dat2, aes(Prediction, Reference)) +
      geom_tile(aes(fill = n_norm)) + 
      geom_text(data=dat2 %>% filter(n_norm>0 & n_norm<1), aes(Prediction, Reference, label = round(n, 1)), size=2) +
      theme(axis.text.x = element_text(angle = 45, hjust=1))
    
    #####
    pdf(file=file.path("Rplots", paste0("ConfMat_",m,"_",s,".pdf")), width=8, height=8)
    print(plot_conf_mat)
    dev.off()
    
  }
}
