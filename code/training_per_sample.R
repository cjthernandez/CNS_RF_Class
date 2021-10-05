#-----------------------------------------------------------------------------------
# this script performs RF training and prediction with the probes of each Nanopore 
# methylation profile.
#
#------------------------------------------------------------------------------------ 
# Per sample training and predictions

rm(list=ls())
# load Cloud samples
load(file.path("Results","betas_new_samples.RData"))

# load Training dataset ((GSE90496))
load(file.path("Results","annotations.RData"))
load(file.path("Results","betas_ba.RData"))
y <- as.factor(anno$`methylation class:ch1`)

# load Calibration model "Capper"
load(file.path("Results","CV","Capper","calfit.100.RData"))

source(file.path("code","R","train-pred_per_Samp.R"))

lst_rf.scores <- c()
lst_probs <- c()
models <- list()
varsels <- list()
time_train <- c()
time_pred <- c()

# Sample names
nams <- sapply(strsplit(names(lst_betas), "_"), function(x) paste0(x[1], x[2])) #

for(i in 1:length(lst_betas)){
  name <- nams[i]
  new_s <- lst_betas[[i]]
  
  # Train RF using probes represented in new_s
  t0 <- proc.time()
  model <- trainRFperSamp(betas=betas, new_samp=new_s, y=y, new_samp_name = name, timed=timed, ntrees = 500, sdfilt=32000)
  t1 <- proc.time()
  
  time_train <- c(time_train, (t1-t0)[3])
  rf.pred <- model$rf.pred
  
  # RF Models and variable selection RFs
  models[[i]] <- rf.pred
  varsels[[i]] <- model$varsel
  
  # Predict using fitted RF
  t0 <- proc.time()
  pred <- predictPerSample(new_s, rf.pred, cv.calfit, new_samp_name = name, timed=timed)
  t1 <- proc.time()
  
  time_pred <- c(time_pred, (t1-t0)[3])
  
  rf.scores <- pred$rf.scores
  probs <- pred$probs
  
  # Raw scores and Calibrated scores
  if(is.null(length(lst_probs))){
    lst_rf.scores <- rf.scores
    lst_probs <- probs
  }else{
    lst_rf.scores <- rbind(lst_rf.scores, rf.scores)
    lst_probs <- rbind(lst_probs, probs)
  }
}

rownames(lst_rf.scores) <- nams
rownames(lst_probs) <- nams
names(models) <- nams
names(varsels) <- nams


save(models, varsels, time_train, file=file.path("Results",paste0("perSampTrain.RData")))

save(lst_rf.scores, lst_probs, time_pred, file=file.path("Results",paste0("perSampPred.RData")))
