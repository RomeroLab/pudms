#' Create (train,test) grouped data sets from a grouped data set
#' 
#' @param grouped_dat a grouped data set
#' @param test_idx an integer from 1 to nfolds
#' @param nfolds number of folds
#' @param seed seed for the reproducibility
#' @param folds NULL or a list of (cv_labeled,cv_unlabeled). If NULL, cv_labeled and cv_unlabeled are generated randomly, then train/test sets are created. If not, train/test data sets are created from previous cv_labeled and cv_unlabeled list
#' @param return.datasets a logical value; if TRUE, test_grouped_dat / train_grouped_dat will be returned
#' @return a list containing (train_grouped_dat, test_grouped_dat,folds = (cv_labeled, cv_unlabeled))
#'@export
cv_grouped_dat = function(grouped_dat,
                          test_idx,
                          nfolds = 10,
                          seed = 1,
                          folds = NULL,
                          return.datasets = TRUE) {
  set.seed(seed)
  if (is.null(folds)) {
    # (cv1, ..., cv10)
    # (a multinomial vector adding up to the number of sequences in the unlabeled/labeled set)
    cv_unlabeled = t(with(grouped_dat, sapply(1:length(unlabeled), function(x) {
      rmultinom(n = 1,
                size = unlabeled[x],
                prob = rep(1 / nfolds, nfolds))
    })))
    colnames(cv_unlabeled) = paste("cv",1:nfolds,sep="")
    rownames(cv_unlabeled) = grouped_dat$seqId
    
    cv_labeled = t(with(grouped_dat, sapply(1:length(labeled), function(x) {
      rmultinom(n = 1,
                size = labeled[x],
                prob = rep(1 / nfolds, nfolds))
    })))
    dimnames(cv_labeled) = dimnames(cv_unlabeled)
  }else{
    cv_unlabeled = folds$cv_unlabeled
    cv_labeled = folds$cv_labeled
  }
  # folds contain all info to reproduce test/train grouped data sets
  folds = list(cv_labeled = cv_labeled, cv_unlabeled = cv_unlabeled)
  
  if(return.datasets){
    # create test_grouped_dat / train_grouped_dat
    test_grouped_dat = grouped_dat
    test_grouped_dat$unlabeled = cv_unlabeled[,test_idx]
    test_grouped_dat$labeled   = cv_labeled  [,test_idx]
    ridx = which(rowSums(test_grouped_dat[,c("labeled","unlabeled")])==0)
    if(length(ridx)>0){test_grouped_dat = test_grouped_dat[-ridx,]}
    
    train_grouped_dat = grouped_dat
    train_grouped_dat$unlabeled = rowSums(cv_unlabeled[,-test_idx,drop=F])
    train_grouped_dat$labeled   = rowSums(cv_labeled  [,-test_idx,drop=F])
    ridx = which(rowSums(train_grouped_dat[,c("labeled","unlabeled")])==0)
    if(length(ridx)>0){train_grouped_dat = train_grouped_dat[-ridx,]}
    
    r = list(
      train_grouped_dat = train_grouped_dat,
      test_grouped_dat = test_grouped_dat,
      folds = folds
    )
  }else{
    r = list(folds = folds)
  }
  return(r)
}
