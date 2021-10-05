#-----------------------------------------------------------------------------------
# tSNE analysis
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

library(Rtsne)
library(RSpectra)
library("ggsci")
library(ggplot2)
library(ggrepel)
library(dplyr)

source(file.path("code","R","RSpectra_pca.R"))
source(file.path("code","R","type2upper.R"))

message("loading preprocessed data ...",Sys.time())
load(file.path("Results","betas_ba.RData")) #####################
betas0 <- betas
n_train <- nrow(betas0)
anno0 <- anno
load(file.path("Results","betas_ba_test.RData"))#########
ii <- which(anno$`methylation class:ch1`=="PIN T, PB B")
anno$`methylation class:ch1`[ii] <- "PIN T,  PB B"
betas1 <- betas
n_test <- nrow(betas1)

load(file.path("Results","betas_brainMetastasis.RData"))
n_met <- length(y_met)
betas_met <- betas

i_b <- na.omit(match(colnames(betas_met), colnames(betas0)))
i_n <- na.omit(match(colnames(betas0), colnames(betas_met)))

betas0 <- betas0[,i_b]
betas_met <- betas_met[,i_n]

betas <- rbind(betas0, betas1, betas_met)

# methylation classes
y <- as.factor(c(anno0$`methylation class:ch1`, anno$`methylation class:ch1`, y_met))
names(y) <- c(anno0$geo_accession, anno$geo_accession, rownames(betas_met))

# SD FILTERING
p <- 32000
betas <- betas[,order(-apply(betas,2,sd))[1:p]]

# calculate first 94 PCs
pca <- prcomp_svds(betas,k=94)

# extract the GSM IDs from the PCs - take only the first entry, splitted by "_"
id <- unlist(strsplit(rownames(pca$x),"_"))
gsm.id <- c(id[seq(1, 3*nrow(betas0), by=3)], 
            id[seq(3*nrow(betas0) + 1, 3*nrow(betas0) + 3*nrow(betas1), by=3)],
            rownames(pca$x)[(n_train+n_test + 1):nrow(pca$x)])

# calculate tSNE
res <- Rtsne(pca$x,pca=F,max_iter=2500,theta=0,verbose=T)

# Dataframe for plotting - make sure that the correct methylation classes are assigned to the GEO IDs
res_df <- data.frame(res$Y, geo_accession=gsm.id)
res_df <- data.frame(res_df, methyl_class=y[as.character(res_df$geo_accession)])
res_df$upper <- type2upper(res_df$methyl_class)
res_df$upper[res_df$upper=="Brain metastasis"] <- as.character(res_df$methyl_class[res_df$upper=="Brain metastasis"])
res_df <- res_df %>% mutate(set=ifelse(grepl("GSM24", geo_accession), "training", ifelse(grepl("GSM29", geo_accession),"test", "metastasis")))

### Plotting tSNE
tsne <- ggplot(res_df, aes(x=X1, y=X2, color=upper)) + ## either dist_df or res_df
  geom_point(aes(shape=set),size=1.5, alpha=0.5) +
  scale_shape_manual(values=c(4,16,17)) +
  guides(colour=guide_legend(override.aes=list(size=5, shape=16, alpha=0.5))) +
  xlab("") + ylab("") +
  ggtitle("t-SNE") +
  theme_light(base_size=10) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.position = "right",
        legend.text = element_text(colour="grey30", size = 6)) +
  scale_color_d3("category20", name="Broad tumor class")

###
dir.create("Rplots", showWarnings = FALSE)

pdf(file=file.path("Rplots","tSNE_training_test_metastasis_classes.pdf")) #, width=8, height=5)
print(tsne)
dev.off()
