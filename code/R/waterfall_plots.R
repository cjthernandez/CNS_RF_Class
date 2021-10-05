# Plots for general results of RF models
rm(list=ls())
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)

source(file.path("code","R","type2upper.R"))

sets <- c("training", "test")
models <- c("Capper","32k","500trees","32k500trees") # BrainMetastasis
v <- sets[2]
m <- models[1]

load(file.path("Results","predictions",paste0("Capper_",m,"_",v), "predictionNewSamples.RData"))

if(v=="test"){
  load("Results/annotations_test.RData")
  y <- anno$`methylation class:ch1`
  ii <- which(y =="PIN T, PB B")
  y[ii] <- "PIN T,  PB B"
  # i_melan <- which(grepl("MELAN", y))######BRAIN MET
  # y[i_melan] <- "Melanoma brain metastasis" ############## BRAIN MET
}else{
  load("Results/annotations.RData")
  # load(file.path("Results", "betas","betas.ba.brainMetastasis_offset.RData")) ######### BRAIN METT
  y <- c(anno$`methylation class:ch1`)#, y_met)
  # Remove "Uncertain" metastasis samples
  # i_uncert <- which(grepl("Uncertain", y))
  # y <- y[-i_uncert]
  # betas <- betas[-i_uncert,]
  # i_melan <- which(grepl("MELAN", y))
  # y[i_melan] <- "Melanoma brain metastasis"
}

# Matches at Methylation CLASS and FAMILY levels
if(!is.list(rf.scores)){
  pred.class.cal <- colnames(probs)[apply(probs,1,which.max)] # Calibrated scores
}else{
  pred.class.cal <- names(probs[[1]])[sapply(probs, which.max, simplify="array")]
}
# mismatches in methylation class and family
mm <- y != pred.class.cal
i_mm <- which(mm)
mmff <- type2upper(y) != type2upper(pred.class.cal)
i_mmf <- which(mmff)

if(set=="training"){ # For training set only methylation class
  pp_mm <- as.data.frame(probs[i_mm,])
  preds <- pred.class.cal[i_mm]
  trues <- y[i_mm]
}else{ # For test set only methylation family
  pp_mm <- as.data.frame(probs[i_mmff,])
  preds <- pred.class.cal[i_mmff]
  trues <- y[i_mmff]
}

i_c <- match(y, colnames(probs))
i_r <- seq_along(y)

ss <- probs[cbind(i_r, i_c)]

pp <- apply(probs,1,max)

# Data for plotting
df_plot <- data.frame(score_class=ss, score_max=pp, pred_class=pred.class.cal, true_class=y, mismatch=!mm, mismatch_mf=!mmff) %>%
  arrange(desc(score_max)) %>% mutate(true_fam=type2upper(true_class), pred_fam=type2upper(pred_class))
df_plot$id <- factor(1:nrow(df_plot), ordered = T)

sum(df_plot$score_max>=0.9)/nrow(df_plot)

# Plot all samples with their maximum probs score
plot_all_samp <- ggplot(df_plot, aes(x=id, y=score_max)) + geom_bar(stat="identity", aes(fill=mismatch))
# Plot only misclassified samples with methylation family
df_plot_mis <- subset(df_plot, !mismatch)
plot_misclass <- ggplot(df_plot_mis, aes(x=id, y=score_max)) + geom_bar(stat="identity", aes(fill=true_fam))

# Waterfall plots
pdf(file=file.path("Rplots",paste0("Waterfall_plot_",m,"_",s,".pdf")), width=8, height=8)
print(plot_all_samp)
print(plot_misclass)
dev.off()
