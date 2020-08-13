#' Obtain an asymptotic variance matrix of theta, standard errors, zvalues, and pvalues
#'
#'@param X an n by p matrix
#'@param z an n by 1 vector
#'@param theta a p+1 by 1 vector
#'@param py1 a numeric value
#'@param weights an n by 1 vector
#'@param effective_n_prop a numeric value
#'@importFrom stats pnorm 
#'@return a list (I = expected Fisher matrix, invI = inverse of I, se = sqrt(diag(invI)), zvalue, pvalue)
#'@export
pval_pu <- function(X,z,theta,py1, weights = rep(1,nrow(X)),effective_n_prop = 1){
  
  wei = weights
  N = sum(wei)
  nl = sum(wei[z==1])
  nu = sum(wei[z==0])
  
  h<-function(eta,nl,nu,py1){
    log(nl/(py1*nu))+eta-log(1+exp(eta))
  }
  
  eta = X%*%theta[-1]+theta[1]
  h_eta = h(eta,nl,nu,py1)
  Xint= cbind(rep(1,nrow(X)),X)
  
  # In(theta) = sum_i (wei_i*effective_n_prop) A''(hi) hi'^2 xi xi'
  w1 = exp(h_eta)/((1+exp(h_eta))^2) #A''(hi)
  w2 = (1/(1+exp(eta)))^2 #hi'^2
  
  W = Matrix::Diagonal(x = as.numeric(wei*effective_n_prop*w1*w2), n=nrow(Xint))
  i = t(Xint)%*%W%*%Xint
  ii = chol2inv(chol(i)) # symmetric, positive definite
  se = sqrt(diag(ii))
  zvalue = theta/se
  pval = pnorm(abs(zvalue),lower.tail = F)*2
  return(invisible(list(invI = ii, se = se, zvalue = zvalue, pvalue = pval)))
}