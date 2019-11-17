#'@export
cv_grouped_dat = function(grouped_dat,
                          test_idx,
                          nfolds = 10,
                          seed = 1,
                          folds = NULL) {
  set.seed(seed)
  if (is.null(folds)) {
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
  
  folds = list(cv_labeled = cv_labeled, cv_unlabeled = cv_unlabeled)
  
  list(
    train_grouped_dat = train_grouped_dat,
    test_grouped_dat = test_grouped_dat,
    folds = folds
  )
}
