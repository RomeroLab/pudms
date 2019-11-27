#'@import PUlasso
#'@export
pudms <- function (protein_dat,
                   py1,
                   order = 1, 
                   aggregate = T, 
                   basestate = NULL,
                   verbose=T,
                   nobs_thresh=10,
                   lambda=0,
                   pvalue = T,
                   maxit=1000,
                   p.adjust.method = "BH"){

  cat(" 1. create a model matrix X from an aggregated dataset:\n")
  Xprotein = create_model_frame(grouped_dat = protein_dat,
                                order = order, 
                                aggregate = aggregate,
                                basestate = basestate,
                                verbose = verbose)
  
  # remove sequences which contain a feature whose number of mutations < nobs_thresh
  # at the end, the function checks whether filtered X matrix is of full rank. 
  # If X is not of full rank, function throws a warning
  filtered_Xprotein = filter_mut_less_than_thresh(Xprotein,
                                                  thresh = nobs_thresh,
                                                  checkResFullRank = T)
  
  cat("\n\n 2. fit a model\n")
  fit = grpPUlasso(X = filtered_Xprotein$X,
                   z = filtered_Xprotein$z,
                   weights = filtered_Xprotein$wei,
                   py1 = py1,
                   verbose = T,
                   lambda = lambda, 
                   maxit = maxit)
  
  
  if(pvalue){
    cat("\n\n 3. compute p-values\n")
    if(lambda>0){stop("currently pvalue computation is supported for lambda = 0")}
    pvalues <-
      pval_pu(
        X = filtered_Xprotein$X,
        z = filtered_Xprotein$z,
        theta = as.numeric(fit$coef),
        py1 = py1,
        weights = filtered_Xprotein$wei
      )
  }else{
    pvalues = NULL
  }
  
  # export results to the excel file
  nobs = c(0,colSums(Diagonal(x=filtered_Xprotein$wei)%*%filtered_Xprotein$X))
  
  # summary 
  dat = data.frame(coef= as.numeric(fit$coef),
                   se = pvalues$se, 
                   zvalue = pvalues$zvalue, 
                   p = pvalues$pvalue,
                   p.adj = p.adjust(pvalues$pvalue, method = p.adjust.method),
                   nobs = nobs)
  dat = dat[-1,] # remove an intercept
  res = structure(list(fitted = fit, result_table = dat, nobs = nobs),class="pudms.fit")
  
  return(res)
}

