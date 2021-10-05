#-----------------------------------------------------------------------------------
# training Random Forest classifier
#                                                                     
# Martin Sill
# m.sill@dkfz.de                                                                  
#
# Modified by: Camilo Hernandez-Toro
#
#------------------------------------------------------------------------------------                   

rm(list=ls())

library(randomForest)
library(parallel)

source(file.path("code","R","train.R"))##

ntrees <- 10000
cores <- 4
seed <- 180314
p <- 10000 # Number of probes to be included in RF

message("loading preprocessed data ...",Sys.time())

#### For Brain Metastasis classifier
load(file.path("Results","betas_brainMetastasis.RData")) ################
betas0 <- betas

load(file.path("Results","betas_ba.RData"))

betas0 <- betas0[,na.omit(match(colnames(betas), colnames(betas0)))]

betas <- rbind(betas,betas0)

message("performing variable selection ...",Sys.time())

y <- c(anno$`methylation class:ch1`,y_met) #### rep("New",3)))

# Merge MELAN and Melanoma brain metastasis classes
i_melan <- which(grepl("MELAN", y)) 
y[i_melan] <- "Melanoma brain metastasis" 

# Remove Uncertain primary tumor metastasis
i_uncert <- which(grepl("Uncertain", y))
y <- as.factor(y[-i_uncert])
betas <- betas[-i_uncert,]
####

## Variance filtering to 32k probes
# betas <- betas[,order(-apply(betas,2,sd))[1:32000]]

set.seed(seed,kind ="L'Ecuyer-CMRG") 
message("seed: ",seed)
message("cores: ",cores)
message("ntrees: ",ntrees)  
message("n: ",nrow(betas))
message("p: ",ncol(betas))  

rf.varsel <- rfp(betas,
                 y,
                 mc=cores,
                 ntree=ntrees,
                 sampsize=rep(min(table(y)),length(table(y))),
                 importance=TRUE)

# get permutation variable importance (MeanDecreaseAccuracy column)
imp.meandecrease <- rf.varsel$importance[,dim(rf.varsel$importance)[2]-1]

# save variable selection forest
save(rf.varsel,file=file.path("Results","varsel.RData")) ###
rm(rf.varsel)

# reduce data matrix. Include only top 'p' probes with highest importance
or <- order(imp.meandecrease,decreasing=T)
betasy <- betas[,or[1:p]]

gc()

message("finished ...",Sys.time())

message("training classifier ...",Sys.time())

message("single core")
message("ntrees: ",ntrees)  
message("n: ",nrow(betasy))
message("p: ",ncol(betasy))

# Capper: final RF classifier: 10000 trees, mtry=100, 10000 probes feature selection
rf.pred <- randomForest(betasy,
                        y,
                        mc=cores,
                        ntree=ntrees,
                        strata=y,
                        mtry=100,#sqrt(ncol(betasy)),
                        sampsize=rep(min(table(y)),length(table(y))),
                        proximity=TRUE,
                        oob.prox=TRUE,
                        importance=TRUE,
                        keep.inbag=TRUE,
                        do.trace=FALSE,
                        seed=seed
)

message("finished ...",Sys.time())

#### CHANGE NAME ACCORDING TO RF TRAINED
mod <- "Capper" # "Capper" "32k" "32k500trees" "500trees" "Brain_Metastasis"
save(rf.pred,file=file.path("Results", paste0("rf.pred_", mod, ".RData")))

message("finished ...",Sys.time())
