#-----------------------------------------------------------------------------------
# script to train the classifier in each CV fold, including batch adjustment and 
# feature selection. Includes CV experiments and down sampling experiments.
# 
# Martin Sill
# m.sill@dkfz.de                                                                  
# 
# 2018-03-14 UTC
#
# Modified by: Camilo Hernandez-Toro
#
#------------------------------------------------------------------------------------                   

calcultateCVfold <- function(badj, y, fold, p, cores, ntrees, test, nreps=10){### CHANGE TO badj when using different percentages of probes
  
  # # SD pre filtering to 32k probes, to speed up the example ################ Comment to use ALL probes
  if(test %in% c("Capper32k", "Capper32k500trees", "probe_randsamp_FixedFold_StdDev")){
    badj$betas.train <- badj$betas.train[,order(-apply(badj$betas.train,2,sd))[1:32000]]
  }
  
  # Select probes based on RF importance
  if(test %in% c("Capper","Capper32k","Capper32k500trees","Capper500trees","brain_Metastasis_full_offset")){
    message("performing variable selection ...",Sys.time())
    message("cores: ",cores)
    message("ntrees: ",ntrees)  
    message("n: ",nrow(badj$betas.train))
    message("p: ",ncol(badj$betas.train))  
    
    rf.varsel <- rfp(badj$betas.train,
                     y=y[fold$train],
                     mc=cores,
                     ntree=ntrees,
                     sampsize=rep(min(table(y[fold$train])),length(table(y[fold$train]))),
                     importance=TRUE)
    
    # get permutation variable importance
    imp.meandecrease <- rf.varsel$importance[,dim(rf.varsel$importance)[2]-1]
    
    # reduce data matrix
    or <- order(imp.meandecrease,decreasing=T)
    badj$betas.train <- badj$betas.train[,or[1:p]]
  }
  
  message("training classifier ...",Sys.time())
  message("cores: ",cores)
  message("ntrees: ",ntrees)  
  message("n: ",nrow(badj$betas.train))
  message("p: ", p)
  
  # Tests without replicates c("Capper","Capper32k","Capper500trees","Capper32k500trees", "brain_Metastasis_full_offset")){
  if(!(grepl("FixedFold", test))){
    # RF using the top p most important variables
    rf <- rfp(badj$betas.train,
              y[fold$train],
              sampsize=rep(min(table(y[fold$train])),length(table(y[fold$train]))),
              mc=cores,
              ntree=ntrees,
              importance=TRUE) 
    
    rf_imp <- rf$importance
    
    rf_feat <- badj$betas.train
    
    
    message("predicting test set ...",Sys.time())
    
    rf.scores <- predict(rf,badj$betas.test[,match(rownames(rf$importance), colnames(badj$betas.test))],type="prob")
    
    err <- sum(colnames(rf.scores)[apply(rf.scores,1,which.max)]!=y[fold$test])/length(fold$test)
    message("misclassification error: ",err)
    
    
    ## Returning the RF predictions, features, and feature importance on the fold
    return(list(rf.scores, rf_feat, rf_imp))
  }
  
  # Tests with replicates
  if(grepl("FixedFold", test)){
    
    miss.class <- list()
    scores.reps <- list()
    imp.reps <- list()
    feat.reps <- list()
    
    pb <- txtProgressBar(min = 0, max = nreps, initial = 0)
    
    # Replicates of probe random sampling
    for(r in 1:nreps){
      setTxtProgressBar(pb, r)
      
      # Random sampling of probes
      l <- ncol(badj$betas.train)
      rs_i <- sample(1:l, p, replace=F)
      badj$betas.train <- badj$betas.train[,rs_i]
      
      # RF using variables sampled randomly
      rf <- rfp(badj$betas.train,
                y[fold$train],
                sampsize=rep(min(table(y[fold$train])),length(table(y[fold$train]))),
                mc=cores,
                ntree=ntrees,
                importance=TRUE) 
      
      rf.scores <- predict(rf,badj$betas.test[,match(rownames(rf$importance),
                                                     colnames(badj$betas.test))],type="prob")
      
      imp.reps[[r]] <- rf$importance
      
      feat.reps[[r]] <- badj$betas.train
      
      scores.reps[[r]] <- rf.scores
      
      err <- sum(colnames(rf.scores)[apply(rf.scores,1,which.max)]!=y[fold$test])/length(fold$test)
      miss.class[[r]] <- err
    }
    
    names(miss.class) <- paste0("rep",1:nreps)
    names(scores.reps) <- paste0("rep",1:nreps)
    names(imp.reps) <- paste0("rep",1:nreps)
    names(feat.reps) <- paste0("rep",1:nreps)
    
    return(list(scores.reps, feat.reps, imp.reps, miss.class))
  }
}
