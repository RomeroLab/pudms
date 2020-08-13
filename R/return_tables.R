#'@import doParallel
#'@import dplyr
#'@importFrom stats pchisq
return_tables <- function(Xprotein, fit, pvalues, nobs, n_eff_prop, p.adjust.method, nCores, verbose, exclude_gap,order,intercept){
  
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
  dat.i = dat.i %>% tibble::rownames_to_column("mut") %>% dplyr::mutate(position =c(-1, Xprotein$blockidx))
  
  # if exclude_gap == TRUE
  if(exclude_gap){
    extract_pos = function(x){x[1]%>% gsub(pattern = "[^0-9]",replacement = "") %>% unlist}
    extract_aa = function(x){x[2]}
    colnames<- colnames(Xprotein$X)
    S = "*"
    if(order==1){
      blockidx = colnames %>% strsplit(split="\\.") %>% sapply(extract_pos) 
      aam = colnames %>% strsplit(split="\\.") %>% sapply(extract_aa) 
      gapInd = (aam==S)
      blockidx[gapInd] = -2
      
    }else{
      pos1 = colnames %>%  strsplit(split=":") %>% lapply(function(x){x[1]}) %>% unlist
      pos2 = colnames %>%  strsplit(split=":") %>% lapply(function(x){x[2]}) %>% unlist
      midx = which(is.na(pos2))
      
      aam = pos1[midx] %>% strsplit(split="\\.") %>% sapply(extract_aa) 
      aa1 = pos1[-midx] %>% strsplit(split="\\.") %>% sapply(extract_aa) 
      aa2 = pos2[-midx] %>% strsplit(split="\\.") %>% sapply(extract_aa) 
      
      posm = pos1[midx] %>% strsplit(split="\\.") %>% sapply(extract_pos) 
      pos1 = pos1[-midx] %>% strsplit(split="\\.") %>% sapply(extract_pos) 
      pos2 = pos2[-midx] %>% strsplit(split="\\.") %>% sapply(extract_pos) 
      
      blockidx = c(posm, paste(pos1,pos2,sep = ":"))
      if(!all(c(-1,blockidx)==block)){stop("error in blockidx")}
      gapInd = c(aam==S, (aa1==S | aa2==S)) # if mutation contains a gap, TRUE
      blockidx[gapInd] = -2
      
    }
    block = c(-1,blockidx)
    
  }
  # check
  # fit$coef[block==-1,,drop=F]
  dat.i = dat.i %>% dplyr::mutate(group = block)
  
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
  
  dat_grp = data.frame(dat.g[,-3], p.grp.adj, nobs.grp=dat.g$nobs.grp)
  dat_grp = dat_grp %>% tibble::rownames_to_column("group")
  
  dat = merge(dat.i,dat_grp,by="group")
  dat[dat$group==-2,c("chi2value","p.grp","p.grp.adj","nobs.grp")] = NA
  
  dat= dat[order(match(dat$mut, dat.i$mut)),]
  if(!intercept){dat = dat[dat$group!=-1,]} # remove an intercept
  dat$group = NULL
  rownames(dat) = 1:nrow(dat)
  # # obtain group size
  # block = as.factor(block)
  # block = factor(block, levels = rownames(dat_grp))
  # ntimes = table(block)
  # 
  # 
  # dat_grp = data.frame(lapply(dat_grp, rep, ntimes))
  # sep = rep("*",nrow(dat.i))
  # dat = data.frame(dat.i,sep,group=block,dat_grp)
  # rownames(dat) = rownames(fit$coef)
  
  if(isParallel) on.exit(stopCluster(cl))
  return(dat)
}