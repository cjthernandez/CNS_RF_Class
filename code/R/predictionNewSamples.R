# Predicting classes for new samples
rm(list=ls())
library(randomForest)
library(glmnet)

source(file.path("code","R","train-pred_per_Samp.R"))

# Load classifier
mod <- "Capper" ### Change RF model accordingly
load(file.path("Results",paste0("rf.pred_",mod,".RData")))
probes_rf <- rownames(rf.pred$importance)
# Load calibrating model
load(file.path("Results","CV",mod,"calfit.100.RData"))

# Example with training data set (GSE90496)
set <- "training"
load(file.path("Results","betas_ba.RData")) ### Change file accordingly

predictions <- predictPerSample(new_samp=new_samp, rf_model=rf.pred, cal_model=cv.calfit)
rf.scores <- predictions$rf.scores
probs <- predictions$probs

dir.create(file.path("Results","predictions",paste0("Capper_",mod,"_",set)), recursive = T)
### Change file/directory name accordingly (RF used, samples predicted)
save(rf.scores, pred.class, probs, pred.class.cal, file=file.path("Results","predictions",
                                                                  paste0("Capper_",mod,"_",set),"predictionNewSamples.RData"))
