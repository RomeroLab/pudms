#' Calculate a Wald statistic
#'
#'@import magrittr
#'@param idx indices for a group which we want to extract a p-value for
#'@param coef an estimated parameter
#'@param vcov a variance-covariance matrix
#'@param verbose a logical value. 
#'@export
wald_pvalue<-function(coef,vcov, position, excludeStates = "*", order = 1, position2 = NULL, verbose = TRUE){
  idx= extractIndices(coef = coef,position = position,excludeStates = excludeStates,order = order,position2 = position2,verbose= FALSE)
  if(verbose) cat("p-value for ",paste(rownames(coef[idx,,drop=F]),collapse = ","),"\n")
  coef_sub = matrix(coef[idx],ncol = 1)
  vcov_inv_sub = chol2inv(chol(vcov[idx,idx,drop=F]))
  chi2value = as.numeric(t(coef_sub)%*%vcov_inv_sub%*%coef_sub)
  p.grp = pchisq(q = chi2value,df = length(idx),lower.tail = F)
  r=c(chi2value,p.grp)
  names(r) = c("chisq2","p")
  r
}

# wald_pvalue(coef = fit$coef,vcov = vcov,position = 12,excludeStates = "*",order = 1,position2 = NULL,verbose = T)
# 
# wald_pvalue(coef = fit$coef,vcov = vcov,position = c(1,2),excludeStates = "*",order = 2,position2 = c(3,4),verbose = T)
