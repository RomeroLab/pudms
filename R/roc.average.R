#' Average ROC curve
#'@import foreach 
#'@importFrom stats approx
#'@param roc_curves roc_curves which need to be averaged
#'@export
roc.average = function(roc_curves){
  
  x = lapply(roc_curves, function(x){x$fpr})
  x_pu = lapply(roc_curves, function(x){x$fpr_pu})
  x = sort(unlist(x))
  x_pu = sort(unlist(x_pu))
  
  dt_m    = sapply(1:length(roc_curves),function(test_idx) approx(x = roc_curves[[test_idx]]$fpr,y = roc_curves[[test_idx]]$tpr, xout =  x ,ties = "mean")$y )
  dt_m_pu = sapply(1:length(roc_curves), function(test_idx) approx(x = roc_curves[[test_idx]]$fpr_pu,y = roc_curves[[test_idx]]$tpr_pu, xout =  x_pu,ties = "mean")$y)
  
  roc_curve_m =list(fpr_pu = x_pu,
                    tpr_pu = apply(dt_m_pu,1,mean,na.rm=T),
                    fpr = x,
                    tpr = apply(dt_m,1,mean,na.rm=T),
                    auc = mean(sapply(roc_curves, function(x){unique(x$auc)})),
                    auc_pu = mean(sapply(roc_curves, function(x){unique(x$auc_pu)})))
  return(roc_curve_m)
}


