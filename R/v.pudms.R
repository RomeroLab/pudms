#' test a PU fit on a test data set
#' 
#'@import parallel
#'@import PUlasso
#'@param protein_dat input data. A data table containing (sequence, labeled, unlabeled, seqId)
#'@param py1 a numeric value, a numeric vector or NULL; the prevalence of positives in the unlabeled data. If length(py1) >1, optimal py1 will be chosen based on auc values on a test data set. If NULL (default), a sequence of py1 values (of length nhyperparam)--ranging from 0.001 to 0.5 interpolated in a log scale--will be considered.
#'@param nhyperparam an integer for the length of the py1 sequence if py1 == NULL
#'@param nfolds the number of subsamples. (nfolds -1)/nfolds splits will be used for training, and the rest will be used for testing.
#'@param test_idx a vector of indices of cross-validation models which will be fitted. Default is to fit the model for each of the cross-validation fold.
#'@param seed a seed number for reproducibility
#'@param order an integer; 1= main effects, 2= main effects + pairwise effects
#'@param refstate a character which will be used for the common reference state; the default is to use the most frequent amino acid as the reference state for each of the position. 
#'@param verbose a logical value. The default is TRUE
#'@param nobs_thresh the number of minimum required mutations per position
#'@param lambda l1 penalty
#'@param pvalue a logial value; if TRUE, p-values based on the asymptotic distribution are obtained
#'@param n_eff_prop proportion of an effective sample size
#'@param intercept a logical value; if TRUE, an estimated intercept is reported together with other coefficients
#'@param maxit maximum number of iterations
#'@param eps convergence threshold for the outer loop
#'@param inner_eps convergence threshold for the inner loop
#'@param initial_coef a vector representing an initial point where we start PUlasso algorithm from.
#'@param p.adjust.method method for multiple comparison
#'@param is.smooth.roc a logical value; if TRUE (default), a smoothed roc curve will be used in checking whether the roc curve is contained by the maximal roc curve at each py1
#'@param tol a numeric value; if the roc estimated curve <= y+tol, the estimated roc curve is determined to be contained by the maximal curve.
#'@param nCores the number of threads for computing.
#'@param full.fit a logical value; if TRUE, the model will be fitted using a full data set and at a chosen py1.
#'@param full.fit.pvalue a logical value; if TRUE, p-values for the full fit will be returned
#'@param outfile NULL or a string; if a string is provided, an output with the name of the string will be exported in a working directory. 
#'@importFrom stats runif
#'@return a list containing v.dmsfit (all fits using training/test splits), roc_curves (average roc curve at each py1), dmsfit (pudms.fit using a full data set at the selected py1),  folds (test/training split information), py1 (a sequence of py1 values used for searching), py1.opt (the selected py1 value based on the predictive performance of the models)
#'@export
v.pudms = function(protein_dat,
                   py1=NULL,
                   nhyperparam = 10,
                   nfolds = 5,
                   test_idx =1:nfolds,
                   seed = round(runif(1, min = 1, max = 1000)),
                   order = 1,
                   refstate =NULL,
                   verbose = T,
                   nobs_thresh=10,
                   lambda = 0,
                   pvalue = FALSE,
                   n_eff_prop = 1,
                   intercept = F,
                   maxit = 1000,
                   eps = 1e-3,
                   inner_eps = 0.01,
                   initial_coef = NULL,
                   p.adjust.method = "BH",
                   is.smooth.roc = TRUE,
                   tol = 1e-5,
                   nCores =1,
                   full.fit = FALSE,
                   full.fit.pvalue = FALSE,
                   outfile = NULL) {
  
  # if nCores>1, set up a parallel environment
  if(nCores>1){
    if(verbose) cat("creating a parallel environment...\n")
    isParallel = TRUE
    nCores = min(nCores,detectCores())
    cl <- makeCluster(nCores)
    
    # export libPaths to the cluster
    invisible(clusterCall(cl, function(x) .libPaths(x), .libPaths()))
    # export libraries to the cluster
    invisible(clusterEvalQ(cl, c(library(pudms),library(pbapply),library(data.table),library(PUlasso))))
    
  }else{
    cl = NULL
    isParallel = FALSE
  }
  if(verbose){pboptions(type="txt")}
  
  # train/test data split
  cvfolds = cv_grouped_dat(grouped_dat = protein_dat,test_idx = 1,nfolds = nfolds,seed = seed,return.datasets = F)
  # create training/test datasets for test_idx
  cv.datasets = lapply(X = 1:length(test_idx),FUN = function(i) {cv_grouped_dat(grouped_dat = protein_dat,test_idx = i,folds = cvfolds$folds)})
  
  gc()
  
  # ----------------------------------------------------------------------------------------------------
  # for each py1, we fit the model and obtain the modified roc curve
  if(length(py1)>1 | is.null(py1)){
    searchPy1 = TRUE
    # if py1==NULL, use a default py1 seq
    if(is.null(py1)){py1 = exp(seq(log(1e-3),log(0.5),length.out = nhyperparam))}
    # check whether py1 is a valid input
    if(!all(py1>0&py1<1)) {stop("all py1 values need to be in (0,1)")}
    
  }else{
    searchPy1 = FALSE
  }
  if(length(lambda)!=1) {stop ("currently, this function is implemented for a single lambda value")}
  
  
  # i: a running index for 1,,,test_idx
  # j: a running index for hyperparam 1,..., nhyperparam
  
  v.dmsfit = list() # list of length length(text_idx); each list contains length(py1) lists
  
  v.1round = function(x,i){
    dmsfit = pudms(protein_dat = cv.datasets[[i]]$train_grouped_dat,
                   py1 = py1[x],
                   order = order, 
                   refstate = refstate,
                   verbose= FALSE,
                   nobs_thresh=nobs_thresh,
                   lambda=lambda,
                   pvalue = pvalue,
                   n_eff_prop = n_eff_prop,
                   intercept = intercept,
                   maxit=maxit,
                   eps = eps,
                   inner_eps = inner_eps,
                   initial_coef = initial_coef,
                   p.adjust.method = p.adjust.method)
    
    roc = adjusted_roc_curve(coef = coef(dmsfit$fit),
                             test_grouped_dat = cv.datasets[[i]]$test_grouped_dat,
                             verbose = F,
                             py1 = py1[x],
                             plot = F)
    
    list(train_fit=dmsfit,test_roc = roc,py1 = py1[x]) 
  }
  
  for(i in 1:length(test_idx)){
    if(verbose) cat("fitting with a fold",i,"...\n")
    v.dmsfit[[i]] = pblapply(X = 1:length(py1), i=i, cl = cl, FUN =v.1round)
  }
  
  # assign names to v.dmsfit 
  names(v.dmsfit) = paste("fold",1:length(test_idx),sep = "")
  for (x in 1:length(test_idx)){names(v.dmsfit[[x]]) = paste("py1_",1:length(py1), sep = "")}
  
  # if length(test_idx)>1, for each hyperparameter py1, we obtain an averaged ROC curve
  
  if(length(test_idx)>1){
    if(verbose) {cat("obtaining an average ROC curve for each py1...\n")}
    if(isParallel) clusterExport(cl,varlist = c("v.dmsfit"),envir = environment())  
    roc_curves = pblapply(1:length(py1),
                          FUN = function(j){
                            roc_curves = lapply(1:length(test_idx), function(i) v.dmsfit[[i]][[j]]$test_roc$roc_curve)
                            m_roc_curve = roc.average(roc_curves)
                            m_roc_curve
                          },cl = cl)
  }else{
    roc_curves = lapply(1:length(py1), FUN = function(j) v.dmsfit[[1]][[j]]$test_roc$roc_curve)
  }
  
  if(searchPy1) {
    py1.opt = optimal_py1(roc_curves = roc_curves,py1 = py1,verbose = verbose,is.smooth.roc = is.smooth.roc,tol = tol)
  }else{
    py1.opt = NULL
  }
  
  if(full.fit & !is.null(py1.opt)){
    if (verbose) cat("fitting with a full dataset...\n")
    dmsfit = pudms(protein_dat = protein_dat,
                   py1 = py1.opt,
                   order = order, 
                   refstate = refstate,
                   verbose= FALSE,
                   nobs_thresh=nobs_thresh,
                   lambda=lambda,
                   pvalue = full.fit.pvalue,
                   n_eff_prop = n_eff_prop,
                   intercept = intercept,
                   maxit=maxit,
                   eps = eps,
                   inner_eps = inner_eps,
                   initial_coef = initial_coef,
                   p.adjust.method = p.adjust.method,
                   outfile = outfile)
  } else{
    dmsfit = NULL
  }
  if(isParallel) on.exit(stopCluster(cl))
  
  structure(list(v.dmsfit = v.dmsfit, 
                 roc_curves = roc_curves,
                 dmsfit = dmsfit, 
                 folds = cvfolds$folds,
                 py1 = py1,
                 py1.opt = py1.opt,
                 call = match.call()),class = "vpudms.fit")
  
}

