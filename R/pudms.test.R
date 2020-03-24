#' test a PU fit on a test data set
#' 
#'@param protein_dat a (full) grouped data set
#'@param test_idx an integer from 1 to nfolds, which will be used for a test set
#'@param nfolds the number of folds
#'@param pvalue a logical value, a default value is TRUE
#'@param seed a seed number for reproducibility
#'@param py1 a prevalence P(y=1)
#'@param maxit max number of iterations
#'@return a list containing 
#'@examples
#'library(puDMSdata)
#'data(grouped_DXS)
#'r = pudms.test(protein_dat = grouped_DXS, test_idx = 1, seed = 20122020,py1 = 0.1,pvalue = F )
#'@export
pudms.test = function(protein_dat, test_idx, nfolds=10, pvalue=T, seed = round(runif(1,min = 1,max = 1000)), py1, maxit = 1000){
  
  # train/test data split
  cv.test_idx = cv_grouped_dat(grouped_dat = protein_dat,test_idx = test_idx,nfolds = nfolds,seed = seed)
  
  # fitting using a train set
  fit.test_idx = pudms(protein_dat = cv.test_idx$train_grouped_dat,py1 = py1, verbose = T, pvalue = pvalue, maxit = maxit)
  
  # roc/auc on a test set
  roc.test_idx = corrected_roc_curve(coef = coef(fit.test_idx),test_grouped_dat = cv.test_idx$test_grouped_dat,verbose = T,py1 = py1)
  
  structure(list( pudms.fit.train =fit.test_idx, roc.test= roc.test_idx,folds= cv.test_idx$folds) ,class = "pudms.test")
}
