# Create betas files from methylation Nanopore data of new samples and Training set MATCHING probes.
rm(list=ls())
library(dplyr)

# Remove probes with NAs
count_na <- function(x){
  sum(is.na(x))
}
na_in_array <- function(x){
  return(count_na(x)>0)
}

missing_in_betas <- function(x, array.probes){
  array.probes[which(is.na(match(array.probes,names(x))))]
}

missing_values <- function(x, miss=T){
  if(miss){
    ii <- which(is.na(x))
  }else{
    ii <- which(!is.na(x))
  }
  names(x)[ii]
}

dir.create(file.path("Results","probesAnalysis"))
message("Select probes using Naopore data", Sys.time())

# New samples Nanopore
load(file.path("Data","LN229_Nanopolish.RData")) ###### Change file with new data
# Capper training data
load(file.path("Results","betas_ba.RData")) ## betas
betas0 <- betas

# namerows <- rownames(df)
# namecols <- colnames(df)
## When data comes as character matrix and missing values are "."
# if(is.character(df)){
#   # Points "." are missing values (NA)
#   # df[df=="."] <- NA
#   df <-apply(df, 2, as.numeric)
#   rownames(df) <- namerows
#   colnames(df) <- namecols
# }

# new_nanop <- data.frame(df)
new_nanop <- new_samp

# Matching probes
# Select only probes in betas shared with Nanopore
betas2 <- betas0[,na.omit(match(colnames(new_nanop),colnames(betas0)))]
betas2 <- betas2[,order(names(betas2))]

# Select only probes in Nanopore shared with betas
new_nanop2 <- new_nanop[,na.omit(match(colnames(betas0),colnames(new_nanop)))]
new_nanop2 <- new_nanop2[,order(names(new_nanop2))]

# Probes in betas not present in Nanopore
no_nanop <- colnames(betas0)[is.na(match(colnames(betas0), colnames(new_nanop)))]

save(no_nanop, file=file.path("Results","probesAnalysis","probes_missing_Nanopore.RData"))

# Probes in betas (array) with NAs in Nanopore
ii <- which(apply(new_nanop2, 2, na_in_array))

if(length(ii)>0){
  probes_na <- new_nanop2[,ii]
  save(probes_na, file=file.path("Results", "probesAnalysis", "probes_NA_Nanopore.RData"))
  # Remove probes with NAs
  new_nanop_narm <- new_nanop2[,-ii]
  betas_narm <- betas2[,-ii]
}else{
  new_nanop_narm <- new_nanop2
  betas_narm <- betas2
}

# Merged dataframe for visualization
betas <- bind_rows(betas_narm, new_nanop_narm)
save(betas,anno,file=file.path("Results","betas_ba_newSamples_merged.RData"))

# Dataframe to train RF
betas <- betas_narm
# Dataframe with new test samples
new_samp <- new_nanop_narm

save(betas, anno, file=file.path("Results","betas_ba_newSamples.RData"))
save(new_samp, file=file.path("Results","new_samples.RData"))


# Probes missing from Full RF model
# Capper RF classifiers
mod <- "Capper" ### Change accordingly
load(file.path("Results",paste0("rf.pred_",mod,".RData")))
rf_probes <- rownames(rf.pred$importance)

probes_miss_rf <- rf_probes[is.na(match(rf_probes, colnames(new_samp)))]
save(probes_miss_rf, file = file.path("Results","probesAnalysis",paste0("probes_missing_",mod,"RF.RData")))
