#' DMS experiments analysis using PUlasso method
#'
#'@param protein_dat input data. A data table containing (sequence, labeled, unlabeled, seqId)
#'@param py1 a numeric value representing the prevalence of positives in the unlabeled data
#'@param order an integer; 1= main effects, 2= main effects + pairwise effects
#'@param refstate a character which will be used for the common reference state; the default is to use the most frequent amino acid as the reference state for each of the position. 
#'@param verbose a logical value. The default is TRUE
#'@param nobs_thresh the number of minimum required mutations per position
#'@param lambda l1 penalty
#'@param nlambda if lambda= NULL, a sequence of nlambda is created for fitting
#'@param pvalue a logial value; if TRUE, p-values based on the asymptotic distribution are obtained
#'@param n_eff_prop proportion of an effective sample size
#'@param intercept a logical value; if TRUE, an estimated intercept is reported together with other coefficients
#'@param maxit maximum number of iterations
#'@param eps convergence threshold for the outer loop
#'@param inner_eps convergence threshold for the inner loop
#'@param p.adjust.method method for multiple comparison
#'@param initial_coef a vector representing an initial point where we start PUlasso algorithm from.
#'@param outfile NULL or a string; if a string is provided, an output with the name of the string will be exported in a working directory. 
#'@return a list containing a fit (from grpPUlasso), result_table (data.frame), and nobs
#'@import PUlasso
#'@importFrom stats p.adjust
#'@importFrom utils write.csv
#'@export
pudms <- function (protein_dat,
                   py1 = NULL,
                   order = 1, 
                   refstate = NULL,
                   verbose=T,
                   nobs_thresh=10,
                   lambda=0,
                   nlambda = 10,
                   pvalue = T,
                   n_eff_prop = NULL,
                   intercept = F,
                   maxit=1000,
                   eps = 1e-3,
                   inner_eps = 0.01,
                   initial_coef = NULL,
                   p.adjust.method = "BH",
                   outfile = NULL){

  if(verbose) cat(" 1. create model frames from an aggregated dataset:\n")
  Xprotein = create_model_frame(grouped_dat = protein_dat,
                                order = order, 
                                aggregate = T,
                                refstate = refstate,
                                verbose = verbose)
  
  # remove sequences which contain a feature whose number of mutations < nobs_thresh
  # at the end, the function checks whether filtered X matrix is of full rank. 
  # If X is not of full rank, function throws a warning
  filtered_Xprotein = filter_mut_less_than_thresh(Xprotein,
                                                  thresh = nobs_thresh,
                                                  checkResFullRank = T,
                                                  verbose = verbose)
  
  if(verbose) cat("\n\n 2. fit a model\n")
  fit = grpPUlasso(X = filtered_Xprotein$X,
                   z = filtered_Xprotein$z,
                   weights = filtered_Xprotein$wei,
                   py1 = py1,
                   verbose = verbose,
                   lambda = lambda, 
                   nlambda = nlambda,
                   maxit = maxit,
                   eps = eps,
                   inner_eps = inner_eps,
                   initial_coef = initial_coef
                   )
  
  
  if(pvalue){
    if(verbose) cat("\n\n 3. compute p-values\n")
    if(lambda>0){stop("currently p-value computation is supported only for lambda == 0")}
    
    if(is.null(n_eff_prop)){n_eff_prop = 1}
    if(n_eff_prop>1 || n_eff_prop <= 0){stop ("n_eff_prop should be between 0 and 1")}
    
    pvalues <-
      pval_pu(
        X = filtered_Xprotein$X,
        z = filtered_Xprotein$z,
        theta = as.numeric(fit$coef),
        py1 = py1,
        weights = filtered_Xprotein$wei,
        effective_n_prop = n_eff_prop
      )
    
  }else{
    pvalues = NULL
  }
  
  # export results to the excel file
  nobs = c(0,colSums(Diagonal(x=filtered_Xprotein$wei)%*%filtered_Xprotein$X))
  
  # summary 
  if(pvalue){
    dat = data.frame(coef= as.numeric(fit$coef),
                     se = pvalues$se, 
                     zvalue = pvalues$zvalue, 
                     p = pvalues$pvalue,
                     p.adj = p.adjust(pvalues$pvalue, method = p.adjust.method),
                     nobs = nobs,
                     eff_nobs = n_eff_prop*nobs)
  }else{
    n_eff_prop = 1
    if(dim(fit$coef)[2]==1){
      dat = data.frame(coef= as.numeric(fit$coef),
                     nobs = nobs,
                     eff_nobs = n_eff_prop*nobs)
    }else{
      dat = data.frame(coef= fit$coef,
                     nobs = nobs,
                     eff_nobs = n_eff_prop*nobs)
    }
    
  }
  
  if(!intercept){dat = dat[-1,]} # remove an intercept
  res = structure(list(fit = fit, result_table = dat, call = match.call()),class="pudms.fit")
  
  if(!is.null(outfile)) {
    cat("saving results as",outfile,"...\n")
    write.csv(res$result_table, file= outfile)}
  return(res)
}

