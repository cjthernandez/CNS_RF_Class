### Misclass error, F1 score, and PR-AUC across all Capper models from CV
rm(list=ls())
library(ggplot2)
library(dplyr)

source(file.path("code","R","type2upper.R"))

### Weighted F1 score
load("Results/CV/Capper32k/f1_scores_all.RData")
f1_32k <- f1_scores_all
f1_32k$upper <- rep(F, nrow(f1_32k))
f1_32k$model <- rep("Capper32k", nrow(f1_32k))
load("Results/CV/Capper32k/f1_scores_all_upper.RData")
f1_up_32k <- f1_scores_all
f1_up_32k$upper <- rep(T, nrow(f1_up_32k))
f1_up_32k$model <- rep("Capper32k", nrow(f1_up_32k))

load("Results/CV/Capper/f1_scores_all.RData")
f1_Capper <- f1_scores_all
f1_Capper$upper <- rep(F, nrow(f1_Capper))
f1_Capper$model <- rep("Capper", nrow(f1_Capper))
load("Results/CV/Capper/f1_scores_all_upper.RData")
f1_up_Capper <- f1_scores_all
f1_up_Capper$upper <- rep(T, nrow(f1_up_Capper))
f1_up_Capper$model <- rep("Capper", nrow(f1_up_Capper))

load("Results/CV/Capper500trees/f1_scores_all.RData")
f1_500 <- f1_scores_all
f1_500$upper <- rep(F, nrow(f1_500))
f1_500$model <- rep("Capper500trees", nrow(f1_500))
load("Results/CV/Capper500trees/f1_scores_all_upper.RData")
f1_up_500 <- f1_scores_all
f1_up_500$upper <- rep(T, nrow(f1_up_500))
f1_up_500$model <- rep("Capper500trees", nrow(f1_up_500))

load("Results/CV/Capper32k500trees/f1_scores_all.RData")
f1_32k500 <- f1_scores_all
f1_32k500$upper <- rep(F, nrow(f1_32k500))
f1_32k500$model <- rep("Capper32k500trees", nrow(f1_32k500))
load("Results/CV/Capper32k500trees/f1_scores_all_upper.RData")
f1_up_32k500 <- f1_scores_all
f1_up_32k500$upper <- rep(T, nrow(f1_up_32k500))
f1_up_32k500$model <- rep("Capper32k500trees", nrow(f1_up_32k500))

# load("Results/CV/brain_Metastasis_full_offset/f1_scores_all.RData")
# f1_BM <- f1_scores_all
# f1_BM$upper <- rep(F, nrow(f1_BM))
# f1_BM$model <- rep("Brain_Metastasis", nrow(f1_BM))
# load("Results/CV/brain_Metastasis_full_offset/f1_scores_all_upper.RData")
# f1_up_BM <- f1_scores_all
# f1_up_BM$upper <- rep(T, nrow(f1_up_BM))
# f1_up_BM$model <- rep("Brain_Metastasis", nrow(f1_up_BM))

f1_all <- rbind(f1_Capper, f1_32k, f1_500, f1_32k500, f1_BM, f1_up_Capper, f1_up_32k, f1_up_500, f1_up_32k500, f1_up_BM)

f1_plot <- f1_all %>% filter(inn_fold==0) %>% group_by(model, type, upper) %>% mutate(f1_score=mean(f1_w)) %>%
  select(-c(f1_w,f1_m,out_fold,inn_fold,p_per)) %>% distinct() %>% ungroup()
f1_plot$dataset <- rep("CV", nrow(f1_plot))

# Filter only OUTER folds
plot_f1_w <- ggplot(f1_plot, aes(x=mod, y=f1_w, color=type, group=out_fold)) +  geom_point(position=position_dodge(width=0.8)) +
  facet_grid(~upper) + geom_vline(xintercept=seq(1.5, length(unique(f1_plot$mod)) - 0.5, 1), color="grey")

pdf(file=file.path("Rplots","Weighted_F1_Capper_models.pdf"), width=10, height=5)
print(plot_f1_w)
dev.off()


### AUC-PR of OUTER folds
load("Results/CV/Capper32k/PR_AUC.100.RData")
auc_32k <- df_aucpr
auc_32k$mod <- rep("Capper32k", nrow(auc_32k))
auc_32k$upper <- rep(F, nrow(auc_32k))
load("Results/CV/Capper32k/PR_AUC_upper.100.RData")
auc_up_32k <- df_aucpr
auc_up_32k$mod <- rep("Capper32k", nrow(auc_up_32k))
auc_up_32k$upper <- rep(T, nrow(auc_up_32k))

load("Results/CV/Capper/PR_AUC.100.RData")
auc_Capper <- df_aucpr
auc_Capper$mod <- rep("Capper", nrow(auc_Capper))
auc_Capper$upper <- rep(F, nrow(auc_Capper))
load("Results/CV/Capper/PR_AUC_upper.100.RData")
auc_up_Capper <- df_aucpr
auc_up_Capper$mod <- rep("Capper", nrow(auc_up_Capper))
auc_up_Capper$upper <- rep(T, nrow(auc_up_Capper))

load("Results/CV/Capper500trees/PR_AUC.100.RData")
auc_500 <- df_aucpr
auc_500$mod <- rep("Capper500t", nrow(auc_500))
auc_500$upper <- rep(F, nrow(auc_500))
load("Results/CV/Capper500trees/PR_AUC_upper.100.RData")
auc_up_500 <- df_aucpr
auc_up_500$mod <- rep("Capper500t", nrow(auc_up_500))
auc_up_500$upper <- rep(T, nrow(auc_up_500))

load("Results/CV/Capper32k500trees/PR_AUC.100.RData")
auc_32k500 <- df_aucpr
auc_32k500$mod <- rep("Capper32k500t", nrow(auc_32k500))
auc_32k500$upper <- rep(F, nrow(auc_32k500))
load("Results/CV/Capper32k500trees/PR_AUC_upper.100.RData")
auc_up_32k500 <- df_aucpr
auc_up_32k500$mod <- rep("Capper32k500t", nrow(auc_up_32k500))
auc_up_32k500$upper <- rep(T, nrow(auc_up_32k500))

# load("Results/CV/brain_Metastasis_full_offset/PR_AUC.100.RData")
# auc_BM <- df_aucpr
# auc_BM$mod <- rep("BrainMet", nrow(auc_BM))
# auc_BM$upper <- rep(F, nrow(auc_BM))
# load("Results/CV/brain_Metastasis_full_offset/PR_AUC_upper.100.RData")
# auc_up_BM <- df_aucpr
# auc_up_BM$mod <- rep("BrainMet", nrow(auc_up_BM))
# auc_up_BM$upper <- rep(T, nrow(auc_up_BM))

auc_all <- rbind(auc_32k, auc_Capper, auc_500, auc_32k500, auc_up_Capper, auc_up_32k, auc_up_500, auc_up_32k500)
auc_all <- auc_all %>% group_by(class_f, upper) %>% mutate(min_auc=min(auc)) %>% ungroup()
auc_all$class_up <- ifelse(auc_all$upper, as.character(auc_all$class_f), type2upper(auc_all$class_f))

auc_all_up <- subset(auc_all, upper)

plot_upper <- ggplot(auc_all_up, aes(x=class_f, y=auc)) + geom_point(aes(color=type), position=position_dodge(width=0.8)) +
  labs(color="Type") + xlab("CNS-MF") + ylab("AUC-PR") +
  facet_grid(mod~.) + theme(axis.text.x=element_text(angle=60, hjust=1)) +
  geom_vline(xintercept=seq(1.5, length(unique(auc_all_up$class_f)) - 0.5, 1), color="grey")

auc_all$class_f <- factor(auc_all$class_f, levels = c("Control", upper2type("Control"),
                          "Other glioma", upper2type("Other glioma"),
                          "Nerve", upper2type("Nerve"),
                          "Pineal", upper2type("Pineal"),
                          "Mesenchymal", upper2type("Mesenchymal"),
                          "Melanocytic", upper2type("Melanocytic"),
                          "Plexus", upper2type("Plexus"),
                          "Glioma IDH", upper2type("Glioma IDH"),
                          "Hematopoietic", upper2type("Hematopoietic"),
                          "Ependymal", upper2type("Ependymal"),
                          "Sella", upper2type("Sella"),
                          "Glio-neuronal", upper2type("Glio-neuronal"),
                          "Glioblastoma", upper2type("Glioblastoma"),
                          "Embryonal", upper2type("Embryonal")))

auc_all$class_up <- factor(auc_all$class_up, levels = c("Control", upper2type("Control"),
                                                      "Other glioma", upper2type("Other glioma"),
                                                      "Nerve", upper2type("Nerve"),
                                                      "Pineal", upper2type("Pineal"),
                                                      "Mesenchymal", upper2type("Mesenchymal"),
                                                      "Melanocytic", upper2type("Melanocytic"),
                                                      "Plexus", upper2type("Plexus"),
                                                      "Glioma IDH", upper2type("Glioma IDH"),
                                                      "Hematopoietic", upper2type("Hematopoietic"),
                                                      "Ependymal", upper2type("Ependymal"),
                                                      "Sella", upper2type("Sella"),
                                                      "Glio-neuronal", upper2type("Glio-neuronal"),
                                                      "Glioblastoma", upper2type("Glioblastoma"),
                                                      "Embryonal", upper2type("Embryonal")))

auc_all_prob <- subset(auc_all, !upper & min_auc<0.9)

plot_problem <- ggplot(auc_all_prob, aes(x=class_f, y=auc, color=class_up)) +
  labs(color="CNS-MF", shape="Type") + xlab("CNS-MC") + ylab("AUC-PR") +
  geom_point(aes(shape=type), position=position_dodge(width=0.8)) + theme(axis.text.x=element_text(angle=60, hjust=1)) +
  facet_grid(mod~.) + geom_vline(xintercept=seq(1.5, length(unique(auc_all_prob$class_f)) - 0.5, 1), color="grey")

pdf(file=file.path("Rplots","AUC_PR_Capper_models.pdf"), width=10, height=5)
print(plot_upper)
print(plot_problem)
dev.off()



### Misclassification error
library(dplyr)

load("Results/CV/Capper/missclass_err_all.RData")
miss_cap <- err_scores_all
miss_cap$model <- rep("Capper", nrow(miss_cap))
miss_cap$upper <- rep(F, nrow(miss_cap))
load("Results/CV/Capper/missclass_err_all_upper.RData")
miss_up_cap <- err_scores_all
miss_up_cap$model <- rep("Capper", nrow(miss_up_cap))
miss_up_cap$upper <- rep(T, nrow(miss_up_cap))

load("Results/CV/Capper32k/missclass_err_all.RData")
miss_32k <- err_scores_all
miss_32k$model <- rep("Capper32k", nrow(miss_32k))
miss_32k$upper <- rep(F, nrow(miss_32k))
load("Results/CV/Capper32k/missclass_err_all_upper.RData")
miss_up_32k <- err_scores_all
miss_up_32k$model <- rep("Capper32k", nrow(miss_up_32k))
miss_up_32k$upper <- rep(T, nrow(miss_up_32k))

load("Results/CV/Capper500trees/missclass_err_all.RData")
miss_500 <- err_scores_all
miss_500$model <- rep("Capper500trees", nrow(miss_500))
miss_500$upper <- rep(F, nrow(miss_500))
load("Results/CV/Capper500trees/missclass_err_all_upper.RData")
miss_up_500 <- err_scores_all
miss_up_500$model <- rep("Capper500trees", nrow(miss_up_500))
miss_up_500$upper <- rep(T, nrow(miss_up_500))

load("Results/CV/Capper32k500trees/missclass_err_all.RData")
miss_32k500 <- err_scores_all
miss_32k500$model <- rep("Capper32k500trees", nrow(miss_32k500))
miss_32k500$upper <- rep(F, nrow(miss_32k500))
load("Results/CV/Capper32k500trees/missclass_err_all_upper.RData")
miss_up_32k500 <- err_scores_all
miss_up_32k500$model <- rep("Capper32k500trees", nrow(miss_up_32k500))
miss_up_32k500$upper <- rep(T, nrow(miss_up_32k500))

# load("Results/CV/brain_Metastasis_full_offset/missclass_err_all.RData")
# miss_BM <- err_scores_all
# miss_BM$model <- rep("BrainMet", nrow(miss_BM))
# miss_BM$upper <- rep(F, nrow(miss_BM))
# load("Results/CV/brain_Metastasis_full_offset/missclass_err_all_upper.RData")
# miss_up_BM <- err_scores_all
# miss_up_BM$model <- rep("BrainMet", nrow(miss_up_BM))
# miss_up_BM$upper <- rep(T, nrow(miss_up_BM))

miss_all <- rbind(miss_cap, miss_32k, miss_500, miss_32k500, miss_BM, miss_up_cap, miss_up_32k, miss_up_500, miss_up_32k500, miss_up_BM)
miss_plot <- miss_all %>% filter(inn_fold==0) %>% group_by(model, type, upper) %>% mutate(misclass_err=mean(err)) %>%
  select(-c(err,out_fold,inn_fold,p_per)) %>% distinct() %>% ungroup()
miss_plot$dataset <- rep("CV", nrow(miss_plot))

pdf(file=file.path("Rplots","Misclasserr_Capper_models.pdf"), width=10, height=5)
print(plot_miss_err)
dev.off()

