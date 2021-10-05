# Calculate F1 score from CV results
library(caret)
library(tibble)
library(dplyr)
library(stringr)

# F1 score calculation
f1score <- function(p,r){
  # p -> Precisio
  # r -> Recall
  2 * p * r / (p + r)
}

# Calculate F1 scores (Weighted and Macro) from CV fold
f1_score_fold <- function(scores, y, fold, upper=F){
  # Simplifying rownames
  rownames(scores) <- sapply(strsplit(rownames(scores),"_"), function(x){x[1]})
  
  # True classes
  y_fold <- y[fold$test]
  
  # Predicted classes (levels(y) and colnames(scores) are in the same order)
  y_pred <- factor(levels(y)[apply(scores, 1, which.max)], levels = levels(y))
  
  if(upper){
    y_fold <- factor(type2upper(y_fold))
    y_pred <- factor(type2upper(y_pred), levels=levels(y_fold))
  }
  
  n_y_fold <- table(y_fold)
  
  # Generate confusion matrix and F1 scores
  con_mat <- confusionMatrix(y_pred, y_fold, mode="prec_recall")
  f1_by_class <- as.data.frame(con_mat[[4]])
  rownames(f1_by_class) <- sapply(strsplit(rownames(f1_by_class),":"), function(x){str_trim(x[2])})
  f1_by_class <- rownames_to_column(f1_by_class)
  
  # F1 score correction when Precision and/or Recall is 0
  # from: https://github.com/dice-group/gerbil/wiki/Precision,-Recall-and-F1-measure
  f1_by_class <- f1_by_class %>% select(rowname, Precision, Recall, F1) %>%
    mutate(n_F1 = ifelse(is.na(Precision), ifelse(is.na(Recall), 1, 0),
                       ifelse(is.na(Recall), 0, ifelse(is.nan(F1), 0, F1)))) # Avoid NA when dividing by 0
  
  # Matching F1 scores and n by class
  ii <- match(f1_by_class$rowname, names(n_y_fold))
  f1_by_class$n <- as.integer(n_y_fold[ii])
  
  # Calculate scores
  f1_weighted <- weighted.mean(f1_by_class$n_F1, f1_by_class$n)
  f1_macro <- mean(f1_by_class$n_F1)
  
  return(c(f1_macro, f1_weighted))
}