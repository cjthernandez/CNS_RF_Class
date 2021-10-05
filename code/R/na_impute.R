# Get number of NAs
count_na <- function(x){
  sum(is.na(x))
}
# Check if array contains NAs
na_in_array <- function(x){
  return(count_na(x)>0)
}

# Imputation of missing values in Nanopore data
na_impute_nanopore <- function(new_nanop, probes_rf, dist_mat=F, impute="mean", imp_dim="cpg", k=5){
  # Find probes with missing values
  imp_nanop <- new_nanop
  
  if(imp_dim=="sample"){
    # Probes with NA
    array_na <-which(apply(new_nanop, 2, na_in_array))
  }
  
  if(imp_dim=="cpg"){
    # When new_nanop is a single sample
    if(is.null(dim(new_nanop))){
      array_na <- 1
    }else{
      # Samples with NA
      array_na <- which(apply(new_nanop, 1, na_in_array))####### new_nanop, NEEDED??? TAKES TIME
    }
    # Get names for matching with distance matrix
    if(impute == "kneighbor"){
      if(is.null(dim(new_nanop))){
        probes_nanop <- names(new_nanop)
      }else{
        probes_nanop <- colnames(new_nanop)
      }
      d_mat <- dist_mat$nn.dists # Distances
      i_mat <- dist_mat$nn.idx # Index of NN
      probes_dist <- rownames(d_mat) # Probes in distance matrix
      probes_shared <- probes_dist[na.omit(match(probes_nanop, probes_dist))]
    }
  }
  
  for(i in array_na){
    # Impute probe with data from that probe in other samples
    if(imp_dim=="sample"){
      x <- unlist(new_nanop[,i])
    }
    # Impute probe with data from other probes in the same sample
    if(imp_dim=="cpg"){
      if(is.null(dim(new_nanop))){
        x <- new_nanop
      }else{
        x <- unlist(new_nanop[i,])
      }
    }
    # NA in sample
    i_na <- which(is.na(x))
    
    # Impute with 0.5 naively
    if(impute == "naive"){
      x[i_na] <- 0.5
    }
    
    # Impute with mean of all other samples/CpGs (non-NA values)
    if(impute == "mean"){
      x_mean <- mean(x, na.rm = T)
      x[i_na] <- x_mean
    }
    
    # Nearest Neighbor Imputation
    if(impute == "kneighbor"){
      x_imp <- c()
      if(k > (length(x) - count_na(x))){
        k <- length(x) - count_na(x)
      }
      probes_nano_na <- probes_nanop[i_na]
      na_dist <- match(probes_nano_na, probes_dist)
      p_dist <- probes_dist[na.omit(na_dist)] # Probes with missing values in current sample shared with distance matrix

      if(length(p_dist)==0)
      {
        message("sample ", i, " could not be imputed. All probes with NAs are missing from distance matrix")
        next 
      }
      
      # Impute only the probes with missing values that are included in RF
      probes_2_imp <- na.omit(match(probes_rf, probes_nano_na)) 
      
      for(j in probes_2_imp){
        probe  <- probes_nano_na[j] 

        if(!(probe %in% p_dist)){ # What happens if probe is not in dist_mat?
          xj <- mean(x, na.rm = T)
        }else{
          d <- i_mat[probe,]
          d_probes <- probes_dist[d] # kNN of probe
          rem <- na.omit(match(p_dist, d_probes)) # Remove probes with NA from kNNs of current probe
          d_probes <- d_probes[-rem]
          i_prob <- na.omit(match(d_probes, probes_nanop))[1:k]# NNs in distance matrix are ordered
          i_p <- probes_nanop[i_prob]
          # Mean of k nearest neighbors
          xj <- mean(x[i_p])
        }
        x_imp <- c(x_imp, xj)
      }
      x[i_na[probes_2_imp]] <- x_imp
    }
    
    if(imp_dim=="sample"){
      imp_nanop[,i] <- x
    }
    if(imp_dim=="cpg"){
      if(is.null(dim(imp_nanop))){
        imp_nanop <- x
      }else{
        imp_nanop[i,] <- x
      }
    }
  }
  
  if(is.null(dim(imp_nanop))){
    names(imp_nanop) <- names(new_nanop)
  }else{
    colnames(imp_nanop) <- colnames(new_nanop)
    rownames(imp_nanop) <- rownames(new_nanop)
  }
  
  return(imp_nanop)
}