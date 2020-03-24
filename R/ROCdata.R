#' Create a data table (seqId, labeled, unlabeled, pr_test) to create a ROC curve
#' 
#' @param coef a named vector for the fitted coefficient vector
#' @param test_grouped_dat a data table (sequence, labeled, unlabeled, seqId) containing validation examples
#' @param order a considered order of effects
#' @param verbose a logical value
#' @return a data table containing (seqId, labeled, unlabeled, pr_test)
#'@export
ROCdata = function(coef,
                   test_grouped_dat,
                   order = 1,
                   basestate = NULL,
                   verbose=T) {
  
  unique_X <- function(X,seqId){
    # return unique rows of X where each row corresponds to each seqId
    Xaug = cbind(seqId,X)
    Xaug_u = Xaug[order(seqId), ]
    Xaug_u = Xaug_u[!duplicated(Xaug_u[,"seqId"]),]
    Xaug_u[,-1]
  }
  
  # create (X,z,wei) for a validation set
  cat("create Xtest for validation examples \n")
  a_Xdat_test = create_model_frame(
    grouped_dat = test_grouped_dat,
    order = order,
    aggregate = T, 
    basestate = basestate,
    verbose=verbose)
  
  # keep only unique sequences in Xtest
  Xtest_unique = with(a_Xdat_test, unique_X(X, seqId))
  
  if(dim(coef)[2]==1){coef = coef[,1]}
  if(is.null(names(coef))){stop("coef should be a named vector")}
  if(is.null(colnames(Xtest_unique))){stop("Xtest should have column names")}
  
  # keep coefficients in coef in Xtest_unique (no contribution otherwise)
  # keep columns in Xtest_unique which exist in coef
  a0 = coef[1]
  coef_test = coef[names(coef) %in% colnames(Xtest_unique)]
  Xtest = Xtest_unique[, colnames(Xtest_unique) %in% names(coef)]
  
  eta = Xtest %*% coef_test + a0
  pr_test = as.numeric(1 / (1 + exp(-eta)))
  
  with(test_grouped_dat, data.table(seqId = sort(unique(a_Xdat_test$seqId)), labeled, unlabeled, pr_test))
}