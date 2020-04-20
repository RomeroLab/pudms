#' Average ROC curve
#'@import foreach 
#'@importFrom stats approx sd
#'@param roc_curves roc_curves which need to be averaged
#'@export
roc.average = function(roc_curves){
  
  x = lapply(roc_curves, function(x){x$fpr})
  x_pu = lapply(roc_curves, function(x){x$fpr_pu})
  
  xseq = function(x,n){
    # generate a sequence at which we interpolated ROC curves
    
    xmin = max(sapply(1:length(x), function(i) min(x[[i]]) )) # common minimum 
    xmax = min(sapply(1:length(x), function(i) max(x[[i]]) )) # common maximum 
    
    # makea range from 1e-4 to 1
    xmin = max(c(xmin,1e-4))
    xmax = min(c(xmax,1))
    
    seq(from=xmin,to = xmax,length.out = n)
  }
  x    = xseq(x,   round(mean(sapply(1:length(x),    function(i) length(x[[i]]))     )))
  x_pu = xseq(x_pu,round(mean(sapply(1:length(x_pu), function(i) length(x_pu[[i]]))  )))
  
  dt_m    = sapply(1:length(roc_curves),function(test_idx) approx(x = roc_curves[[test_idx]]$fpr,y = roc_curves[[test_idx]]$tpr, xout =  x ,ties = "mean")$y )
  dt_m_pu = sapply(1:length(roc_curves), function(test_idx) approx(x = roc_curves[[test_idx]]$fpr_pu,y = roc_curves[[test_idx]]$tpr_pu, xout =  x_pu,ties = "mean")$y)
  
  roc_curve_m =list(fpr_pu = x_pu,
                    tpr_pu = apply(dt_m_pu,1,mean,na.rm=T),
                    tpr_pu_sd = apply(dt_m_pu,1,sd,na.rm=T),
                    fpr = x,
                    tpr = apply(dt_m,1,mean,na.rm=T),
                    tpr_sd = apply(dt_m,1,sd,na.rm=T),
                    auc = mean(sapply(roc_curves, function(x){unique(x$auc)})),
                    auc_pu = mean(sapply(roc_curves, function(x){unique(x$auc_pu)})))
  return(roc_curve_m)
}


