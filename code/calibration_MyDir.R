#-----------------------------------------------------------------------------------
# this script fits for each CV fold a multinomial L2-penalized glmnet model to calibrate RF scores
# and one final calibration model using RF scores generated in the out loop of the CV 
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

library(rmarkdown)
library(glmnet)
library(doParallel)
library(HandTill2001)


#### Which sampling test? ####
tests <- c("Capper", # 10k trees, no SD filter, importance selection
           "Capper32k", # 10k trees, SD filter 32k, importance selection
           "Capper500trees", # 500 trees, no SD filter, importance selection
           "Capper32k500trees", # 500 trees, SD filter, importance selection
           "probe_randsamp_FixedFold", # 500 trees, random selection
           "probe_randsamp_FixedFold_StdDev", # 500 trees, SD filter, random selection
           "brain_Metastasis_FixedFold_offset", # 500 trees, random selection
           "brain_Metastasis_full_offset") # 1000 trees, random selection

test <- tests[5] # Modify to perform CV experiments and down sampling experiments

cores <- 4

if(grepl("FixedFold",test)){
  nrep <- 10
}else{
  nrep <-1
}

if(grepl("Capper", test)){
  perc <- 100  
}else{
  perc <- c(10, 25, 33, 50, 66, 75, 100)
}

if(grepl("full", test)){
  perc <- 100
}

registerDoParallel(cores)

message("loading data ...",Sys.time())
load(file.path("Results","annotations.RData"))
if(grepl("Metastasis", test)){
  load(file.path("Results","betas_brainMetastasis.RData"))
  load(file.path("Results","CV","nfolds_brainMetastasis.RData"))#
  y <- c(anno$`methylation class:ch1`, y_met)
  i_melan <- which(grepl("MELAN", y))
  y[i_melan] <- "Melanoma brain metastasis"
  i_uncert <- which(grepl("Uncertain", y))
  y <- as.factor(y[-i_uncert])
}else{
  load(file.path("Results","CV","nfolds.RData"))#
  y <- anno$`methylation class:ch1`
}

# Percentages
for(pr in perc){
  message("Percentage ", pr," ...",Sys.time())
  # Outer folds
  for(i in 1:length(nfolds)){
    lst_probs <- list()
    # Random sampling replicates
    for(rep in 1:nrep){
      scores <- list() 
      idx <- list()
      message("Replicate ",rep," of ", nrep," ...",Sys.time())
      # Inner folds
      for(j in 1:length(nfolds)){
        # Load rf.scores to fit calibration model
        fname <- paste("CVfold",i,j,"probes",pr,"scores","RData", sep=".")#
        load(file.path("Results","CV",test,fname))#
        if(!grepl("FixedFold",test)){
          scores[[j]] <- rf.scores
        }else{
          # Use scores from same replicate
          scores[[j]] <- rf.scores[[rep]]
        }
        idx[[j]] <- nfolds[[i]][[2]][[j]]$test
      }
      
      scores <- do.call(rbind,scores)
      idx <- unlist(idx)
      y <- y[idx]
      
      message("fitting calbriation model fold ", i," ...",Sys.time())
      # fit multinomial logistic ridge (alpha=) regression model (10-fold CV is default)
      suppressWarnings(cv.calfit <- cv.glmnet(y=y,x=scores,family="multinomial",type.measure="mse",
                                              alpha=0,nlambda=100,lambda.min.ratio=10^-6,parallel=TRUE))
      
      # Load rf.scores to be calibrated
      fname <- paste("CVfold",i,0,"probes",pr,"scores","RData", sep=".")
      load(file.path("Results","CV",test,fname))
      
      message("calibrating raw scores fold ",i," ...",Sys.time())
      
      if(!grepl("FixedFold",test)){
        # Calibrate scores
        probs <- predict(cv.calfit$glmnet.fit,newx=rf.scores,type="response",
                         s=cv.calfit$lambda.1se)[,,1] # use lambda estimated by 10fold CVlambda
        
        fname_probs <- paste("probsCVfold",i,0,"probes",pr,"RData", sep=".")
        save(probs, file=file.path("Results", "CV", test, fname_probs))
      }else{
        # Calibrate scores
        probs <- predict(cv.calfit$glmnet.fit,newx=rf.scores[[rep]],type="response",
                         s=cv.calfit$lambda.1se)[,,1] # use lambda estimated by 10fold CVlambda
        # Append to list of calibrations from each replicate
        lst_probs[[rep]] <- probs
      }
      
      err <- sum(colnames(probs)[apply(probs,1,which.max)] != y[nfolds[[i]][[1]][[1]]$test]) / length(nfolds[[i]][[1]][[1]]$test)
      
      message("misclassification error: ",err)
      
    }
    if(grepl("FixedFold",test)){
      probs <- lst_probs
      names(probs) <- paste0("rep", 1:nrep)
      
      fname_probs <- paste("probsCVfold",i,0,"probes",pr,"RData", sep=".")
      save(probs, file = file.path("Results", "CV", test, fname_probs))
    }
  }
}


#### Making the report and the calibration model
## Not defined for "FixedFold" tests
if(!grepl("FixedFold",test)){
  for(pr in perc){
    scores <- list()
    idx <- list()
    for(i in 1:length(nfolds)){
      fname <- paste("CVfold",i,0,"probes",pr,"scores","RData",sep=".") #
      load(file.path("Results","CV",test,fname)) #
      scores[[i]] <- rf.scores
      idx[[i]] <- nfolds[[i]][[1]][[1]]$test
    }
    scores <- do.call(rbind,scores)
    
    probl <- list()
    for(i in 1:length(nfolds)){
      fname <- paste("probsCVfold",i,0,"probes",pr,"RData", sep=".") #
      load(file.path("Results","CV",test,fname)) #
      probl[[i]] <- probs
    }
    probs <- do.call(rbind,probl)
    
    
    idx <- unlist(idx)
    y <- y[idx]
    
    ys <- colnames(scores)[apply(scores,1,which.max)]
    yp <- colnames(probs)[apply(probs,1,which.max)]
    
    errs <- sum(y!=ys)/length(y)
    errp <- sum(y!=yp)/length(y)
    
    message("overall misclassification error scores: ",errs, " percent: ",pr)
    message("overall misclassification error calibrated: ",errp, " percent: ",pr)
    
    message("fitting final calibration model ...",Sys.time())
    
    suppressWarnings(cv.calfit <- cv.glmnet(y=y,x=scores,family="multinomial",type.measure="mse",
                                            alpha=0,nlambda=100,lambda.min.ratio=10^-6,parallel=TRUE))
    
    
    cfname <- paste("calfit",pr,"RData", sep = ".")
    save(cv.calfit,file=file.path("Results","CV",test,cfname))#
    
    cvname <- paste0("CVresults_",pr,".RData")
    save(scores,probs,y,ys,yp,file=file.path("Results","CV",test,cvname)) #
    
    message("finished ...",Sys.time())
  }
}
