# Pheatmap of methylation rates with genomic location of the probes
rm(list=ls())

library(rtracklayer)
library(GenomicDistributions)
library(GenomicDistributionsData)
library(dplyr)
library(tidyr)

source(file.path("code","R","type2upper.R"))

# Annotation of EPIC array probes
probe_gen_loc <- import(file.path("Data","EPIC_annotation.bed"), format="bed")
# CpG island locations
cpg <- import(file.path("Data","CpGIslands","cpgIslandExt-hg38.bed"), format="bed")

# Find overlaps between probe locations and CpG locations
overlap <- findOverlaps(probe_gen_loc, cpg) # EPIC probes in CpGs

# Distribution of EPIC probes across genomic regions
gp = calcPartitionsRef(probe_gen_loc, "hg38")
cpg_gp <- mean(mcols(probe_gen_loc)$name %in% mcols(probe_gen_loc)@listData[[1]][overlap@from])
# plotPartitions(gp)

# These RF need to be generated.
rf.preds <- c("32k500trees",
              "32k",
              "500trees",
              "Capper")

partition <- c()
freq <- c()
cla <- c()

i <- 1
for(rf in rf.preds){
  ##### Distribution of probes selected in RF model across genomic regions
  load(file.path("Results",paste0("rf.pred_",rf,".RData")))
  # Finding probes shared between RF models
  if(i==1){
    shared <- rownames(rf.pred$importance)
  }else{
    shared <- intersect(shared, rownames(rf.pred$importance))
  }
  i <- i + 1
  
  rf_name <- rf
  
  if(rf_name != "Capper"){
    rf_name <- paste0("Capper",rf_name)
  }
  
  rf_imp <- rf.pred$importance
  n_probe <- nrow(rf_imp)
  df_imp <- data.frame(rf_imp)
  colnames(df_imp) <- colnames(rf_imp)
  
  # Partitions Full RF model
  prob_names_full_rf <- rownames(df_imp)
  prob_ind_full <-match(prob_names_full_rf, mcols(probe_gen_loc)$name)
  probe_full_rf <- probe_gen_loc[na.omit(prob_ind_full)]
  gp_full <- calcPartitionsRef(probe_full_rf, "hg38")
  
  cpg_full <- mean(prob_names_full_rf %in% mcols(probe_gen_loc)@listData[[1]][overlap@from]) # Probes in CpG

  # For dataframe
  if(length(cla)==0){
    partition <- c(levels(gp$partition), "CpG", levels(gp_full$partition), "CpG")
    freq <- c(gp$Freq/sum(gp$Freq), cpg_gp, gp_full$Freq/sum(gp_full$Freq), cpg_full)
    cla <- c(rep("EPIC Array", length(gp$Freq)+1), rep(rf_name, length(gp_full$Freq)+1))
  }else{
    partition <- c(partition, levels(gp_full$partition), "CpG")
    freq <- c(freq, gp_full$Freq/sum(gp_full$Freq), cpg_full)
    cla <- c(cla, rep(rf_name, length(gp_full$Freq)+1))
  }
  
}

prob_i <-match(shared, mcols(probe_gen_loc)$name)
probe_shared_rf <- probe_gen_loc[na.omit(prob_i)]
gp_shared <- calcPartitionsRef(probe_shared_rf, "hg38")

cpg_shared <- mean(shared %in% mcols(probe_gen_loc)@listData[[1]][overlap@from])

partition <- c(partition, levels(gp_shared$partition), "CpG")
freq <- c(freq, gp_shared$Freq/sum(gp_shared$Freq), cpg_shared)
cla <- c(cla, rep("Shared", length(gp_shared$Freq)+1))

df_partitions <- data.frame(partition, freq, cla)
df_partitions$partition <- as.factor(as.character(df_partitions$partition))

###Plotting###
library(ggplot2)

# General view of partitions
df_general_partitions <- df_partitions %>% filter(partition != "CpG")
df_cpg_partitions <- df_partitions %>% filter(partition == "CpG")

plot_general <- ggplot(df_general_partitions, aes(x=cla, y=freq, fill=partition)) +
  geom_col() +
  coord_flip() +
  ylab("Frequency")+
  xlab("")+
  theme(legend.position="bottom")+
  ggtitle("Probe distributions across partitions")

plot_cpg <- ggplot(df_cpg_partitions, aes(x=cla, y=freq, fill=partition)) +
  geom_col() +
  coord_flip() +
  ylab("Frequency")+
  xlab("")+
  theme(legend.position="bottom")+
  ggtitle("Probes in CpG")

pdf(file=file.path("Rplots","probeDistribution_GenomicPartitions_Capper_Models.pdf"), width = 7, height = 5)
print(plot_general)
print(plot_cpg)
dev.off()
