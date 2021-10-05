# Per sample Training and Predicition functions

library(randomForest)
library(parallel)
library(glmnet)

source(file.path("code","R","train.R"))##
load(file.path("Results","annotations.RData"))

y_anno <- as.factor(anno$`methylation class:ch1`)

# Model fitting with probes matching probes in nanopore sample. 
trainRFperSamp <- function(betas, new_samp, y=y_anno, new_samp_name="new_samp", 
                           p=10000, ntrees=10000, cores=4, sdfilt=50000, seed=180314){#Default values similar to Capper et al.
  # Only use probes represented in new_samp (present in data and non-missing)
  i_na <- which(is.na(new_samp))
  if(length(i_na)!=0){
    new_no_na <- new_samp[-i_na]
  }else{
    new_no_na <- new_samp
  }
  
  betas <- betas[,na.omit(match(names(new_no_na), colnames(betas)))]
  
  # Variance filtering
  if(!isFALSE(sdfilt)){
    betas <- betas[,order(-apply(betas,2,sd))[1:sdfilt]]
  }
  
  message("performing variable selection ...",Sys.time())
  
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
  
  # Capper: final RF classifier
  rf.pred <- randomForest(betasy,
                          y,
                          mc=cores,
                          ntree=ntrees,
                          strata=y,
                          mtry=100,#sqrt(ncol(betas)),
                          sampsize=rep(min(table(y)),length(table(y))),
                          proximity=TRUE,
                          oob.prox=TRUE,
                          importance=TRUE,
                          keep.inbag=TRUE,
                          do.trace=FALSE,
                          seed=seed
  )
  
  message("finished ...",Sys.time())
  return(list(rf.pred=rf.pred, varsel=rf.varsel))
}

# Predictions per sample
predictPerSample <- function(new_samp, rf_model, cal_model=F, new_samp_name="new_samp", timed=F){
  probes_rf <- rownames(rf_model$importance)
  if(is.list(new_samp)){
    # Get raw RF class probabilities
    rf.scores <- lapply(new_samp, function(x){predict(rf_model, x[probes_rf], type="prob")})
    names(rf.scores) <- names(new_samp)
    # Calibrate class probabilities
    if(!isFALSE(cal_model)){
    probs <- lapply(rf.scores, function(x){predict(cal_model$glmnet.fit, newx=x, type="response", s=cal_model$lambda.1se)[,,1]})
    names(probs) <- names(rf.scores)
    }
  }else{
    # Get raw RF class probabilities
    rf.scores <- predict(rf_model, new_samp[probes_rf], type="prob")
    # Calibrate class probabilities
    if(!isFALSE(cal_model)){
    probs <- predict(cal_model$glmnet.fit, newx=rf.scores, type="response", s=cal_model$lambda.1se)[,,1]
    }
    if(is.matrix(new_samp)){
      rownames(rf.scores) <- rownames(new_samp)
      if(!isFALSE(cal_model)){
      rownames(probs) <- rownames(rf.scores)
      }
    }
  }
  
  return(list(rf.scores=rf.scores, probs=probs))
}
