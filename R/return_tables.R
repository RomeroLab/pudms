#'@import doParallel
#'@importFrom stats pchisq
return_tables <- function(Xprotein, fit, pvalues, nobs, n_eff_prop, p.adjust.method, nCores, verbose){
  
  if(nCores>1){
    if(verbose) cat("creating a parallel environment...\n")
    isParallel = TRUE
    nCores = min(nCores,detectCores())
    clustertype = ifelse(.Platform$OS.type=="windows", 'PSOCK', 'FORK') 
    cl <- makeCluster(nCores,type = clustertype)
    # export libPaths to the cluster
    invisible(clusterCall(cl, function(x) .libPaths(x), .libPaths()))
    # export libraries to the cluster
    invisible(clusterEvalQ(cl, c(library(pudms),library(pbapply),library(data.table),library(PUlasso))))
    
  }else{
    cl = NULL
    isParallel = FALSE
    registerDoSEQ()
  }
  
  # individual effects
  
  dat.i = data.frame(coef= as.numeric(fit$coef),
                     se = pvalues$se, 
                     zvalue = pvalues$zvalue, 
                     p = pvalues$pvalue,
                     p.adj = p.adjust(pvalues$pvalue, method = p.adjust.method),
                     nobs = nobs,
                     eff_nobs = n_eff_prop*nobs)
  
  block  = c(-1, Xprotein$blockidx) # to include intercept
  
  i=1;
  dat.g = foreach(i=1:length(unique(block)),
                  .combine = "rbind",
                  .packages = c("Matrix"))%dopar%{
    
    ridx = which(block == unique(block)[i]) # indices coef to each block
    coef_sub = matrix(fit$coef[ridx],ncol = 1)
    invI_sub = chol2inv(chol(pvalues$invI[ridx,ridx,drop=F]))
    chi2value = as.numeric(t(coef_sub)%*%invI_sub%*%coef_sub)
    
    p.grp = pchisq(q = chi2value,df = length(ridx),lower.tail = F)
    nobs.grp = sum(nobs[ridx])
    c(chi2value,p.grp,nobs.grp)
    
  }
  dat.g= data.frame(dat.g)
  colnames(dat.g) = c("chi2value","p.grp","nobs.grp")
  rownames(dat.g) = unique(block)
  
  p.grp.adj = p.adjust(p = dat.g[,"p.grp"],method = p.adjust.method)
  
  dat_grp = data.frame(dat.g[,-3], p.grp.adj, dat.g$nobs.grp)
  
  # obtain group size
  block = as.factor(block)
  block = factor(block, levels = rownames(dat_grp))
  ntimes = table(block)
  
  dat_grp = data.frame(lapply(dat_grp, rep, ntimes))
  sep = rep("*",nrow(dat.i))
  dat = data.frame(dat.i,sep,group=block,dat_grp)
  rownames(dat) = rownames(fit$coef)
  
  if(isParallel) on.exit(stopCluster(cl))
  return(dat)
}