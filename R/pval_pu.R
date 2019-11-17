pval_pu <- function(X,z,theta,py1, weights = rep(1,nrow(X))){
  
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
  
  # In(theta) = sum_i w_i A''(hi) hi'^2 xi xi'
  w1 = exp(h_eta)/((1+exp(h_eta))^2) #A''(hi)
  w2 = (1/(1+exp(eta)))^2 #hi'^2
  
  W = Matrix::Diagonal(x = as.numeric(wei*w1*w2), n=nrow(Xint))
  i = t(Xint)%*%W%*%Xint
  ii = chol2inv(chol(i))
  se = sqrt(diag(ii))
  zvalue = theta/se
  pval = pnorm(abs(zvalue),lower.tail = F)*2
  return(invisible(list(I = i, invI = ii, se = se, zvalue = zvalue, pvalue = pval)))
}