#-----------------------------------------------------------------------------------
# preprocessing
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

library(minfi)
library(GEOquery)
library(limma)
library(openxlsx)
library(dplyr)

source(file.path("code","R","MNPprocessIDAT_functions.R")) 

dir.create("Results", showWarnings = FALSE)

# get sample annotation from GEO TRAIN (GSE90496) or TEST (GSE109379) sets
gse <- getGEO("GSE109379", GSEMatrix=TRUE, getGPL=FALSE)
anno <- pData(gse$GSE109379_series_matrix.txt.gz)

# read raw data downloaded from series GSE90496 and GSE109379 GEO Repositories
filepath <- file.path("Data","GSE109379_RAW",gsub("_Grn.*","",gsub(".*suppl/","",anno$supplementary_file)))
RGset <- read.metharray(filepath,verbose=TRUE) #GEOset

# Illumina normalization
message("running normalization ...",Sys.time())
# Preprocessing: Normalization and background correction
Mset <- MNPpreprocessIllumina(RGset)

# Probe filtering
message("probe filtering ...",Sys.time())
amb.filter <- read.table(file.path("code","filter","amb_3965probes.vh20151030.txt"),header=F)
epic.filter <- read.table(file.path("code","filter","epicV1B2_32260probes.vh20160325.txt"),header=F)
snp.filter <- read.table(file.path("code","filter","snp_7998probes.vh20151030.txt"),header=F)
xy.filter <- read.table(file.path("code","filter","xy_11551probes.vh20151030.txt"),header=F)
rs.filter <- grep("rs",rownames(Mset))
ch.filter <- grep("ch",rownames(Mset))

# filter CpG probes
remove <- unique(c(match(amb.filter[,1], rownames(Mset)),
                   match(epic.filter[,1], rownames(Mset)),
                   match(snp.filter[,1], rownames(Mset)),
                   match(xy.filter[,1], rownames(Mset)),
                   rs.filter,
                   ch.filter))

Mset_filtered <- Mset[-remove,]

### Saving Filtered Data ###
fname <- "Mset_filtered.RData" # # "Mset_filtered_test.RData" for TEST set
save(Mset_filtered, anno,file=file.path("Results", fname)) 

rm(Mset)
rm(GEOset)
rm(NWset)
rm(RGset)
gc()

#batch adjustment
message("performing batchadjustment ...",Sys.time())


methy <- getMeth(Mset_filtered)
unmethy <- getUnmeth(Mset_filtered)

rm(Mset_filtered)
gc()

# get FFPE/Frozen type
ffpe <- anno$`material:ch1`
batch <- ifelse(ffpe == "FFPE", 2, 1)

# remove batch effects by linear model
methy.ba <- 2^removeBatchEffect(log2(methy +1), batch)
unmethy.ba <- 2^removeBatchEffect(log2(unmethy +1), batch)

# recalculate betas, illumina like (BETAS = METHYLATION RATES)
#betas <- methy.all / (methy.all +unmethy.all +100)
betas <- methy.ba / (methy.ba +unmethy.ba +100)

betas <- as.data.frame(t(betas))

fname <- "betas_ba.RData" # "betas_ba_test.RData" for TEST set
aname <- "annotations.RData" # annotations_test.RData
save(betas,anno,file=file.path("Results",fname))
save(anno, file = file.path("Results",aname))
message("preprocessing finished ...",Sys.time())
