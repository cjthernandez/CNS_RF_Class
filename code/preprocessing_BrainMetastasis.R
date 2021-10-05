#-----------------------------------------------------------------------------------
# this script process the data from the GSE44661 and GSE108576 repositories into 
# methylation rates and filters out the probes that are not shared with the
# GSE90496 repository
#
#------------------------------------------------------------------------------------   
# Preprocessing Brain Metatatic samples
rm(list=ls())

library(GEOquery)
library(limma)
library(dplyr)

# p_val <- 0.05

# Get annotations from GEO
geo1 <- "GSE44661"
gse <- getGEO("GSE44661", GSEMatrix=TRUE, getGPL=FALSE)
anno <- pData(gse$GSE44661_series_matrix.txt.gz) # include only Tissue and brain metastasis (Paper excludes MBM-13 and MBM-1 ¿?)
anno1 <- anno %>% filter(description.2=="Tissue", grepl("brain metastasis",`tissue:ch1`))

samp_names1 <- sub("$","_",sub("-",".", anno1$description)) # Sample names

# Load signal file
df <- read.table(file = file.path("Data","GSE44661_methylated_unmethylated_signals.txt.gz"), header = T, sep = "\t")

coln <- as.character(df$ID_REF)
betas10 <- c()
for(n in samp_names1){
    col_i <- grep(n, colnames(df))
    # For p-value filtering, column matching per sample is needed (different probes will pass for each sample)
  signal_samp <- df[,col_i] %>% #filter(.[[3]] < 0.05) %>%
    mutate(beta=.[[2]]/(.[[1]] + .[[2]] + 100))
  
  if(length(betas10)==0){
    betas10 <- signal_samp$beta
  }else{
    betas10 <- cbind(betas10, signal_samp$beta)
  }
}

rownames(betas10) <- coln
betas1 <- t(betas10)
rownames(betas1) <- paste0(samp_names1,geo1)

# Get annotations from GEO
geo2 <- "GSE108576"
gse2 <- getGEO("GSE108576", GSEMatrix=TRUE, getGPL=FALSE)
anno2 <- pData(gse2$GSE108576_series_matrix.txt.gz)

samp_names2 <- paste0("^",sub("$",".",sub("-",".", anno2$description.1))) # Sample names

# Load signal file
df2 <- read.table(file = file.path("Data","GSE108576_hoond_12-22-17_signal.txt.gz"), header = T, sep = "\t")

coln <- as.character(df2$ID)
betas20 <- c()
for(n in samp_names2){
    col_i <- grep(n, colnames(df2))
    # For p-value filtering, column matching per sample is needed (different probes will pass for each sample)
    signal_samp <- df2[,col_i] %>% #filter(.[[3]] < 0.05) %>%
      mutate(beta=.[[2]]/(.[[1]] + .[[2]] + 100))
    
    if(length(betas20)==0){
      betas20 <- signal_samp$beta
    }else{
      betas20 <- cbind(betas20, signal_samp$beta)
    }
}

rownames(betas20) <- coln
betas2 <- t(betas20)
rownames(betas2) <- paste0(substr(sub(".$","_", samp_names2),2,10), geo2)

betas1 <- betas1[,order(colnames(betas1))]
betas2 <- betas2[,order(colnames(betas2))]

betas <- rbind(betas1, betas2)

y_met <- c(sub("melanoma", "Melanoma", anno1$`tissue:ch1`), anno2$`disease state:ch1`) # Metastasis labels

# Directory
dir.create("Results", showWarnings = FALSE)

save(betas, y_met, file=file.path("Results","betas_brainMetastasis.RData"))



# Selecting probes shared with Capper data
rm(list=ls())
library(dplyr)

# New samples Nanopore
load(file.path("Results","betas_brainMetastasis.RData")) ######
metas <- betas
load(file.path("Results","betas_ba.RData")) ## Capper Methylation rates
betas0 <- betas

# Matching probes
betas2 <- betas0[,na.omit(match(colnames(metas),colnames(betas0)))]
betas2 <- betas2[,order(colnames(betas2))]
metas2 <- metas[,na.omit(match(colnames(betas0),colnames(metas)))]
metas2 <- metas2[,order(colnames(metas2))]

betas <- rbind(betas2, metas2)

save(betas,anno, y_met,file=file.path("Results","betas_ba_brainMetastasis_merged.RData"))

