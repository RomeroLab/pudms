#' Search for optimal py1 values using roc_curves
#'
#'@import pbapply
#'@importFrom stats smooth.spline predict
#'@param roc_curves a list containing roc_curves of length(py1)
#'@param py1 a sequence of py1
#'@param tol NULL or a numerical number to specify a margin in determining whether the estimated ROC curve is smaller than the maximal curve. If NULL, 1sd from K-fold ROC curves at each value is used.
#'@param verbose a logical value
#'@export
optimal_py1 = function(roc_curves, py1, tol = NULL, verbose=T){
  
  if(verbose) {cat("choosing an optimal py1 value...\n")}else{pboptions(type="none")}
  isROCcontained = pbsapply(X =  1:length(py1), FUN = function(j) {
      
      x = roc_curves[[j]]
      if(!is.null(x$tpr_pu_sd)){
        dat = data.frame(fpr_pu = x$fpr_pu, tpr_pu = x$tpr_pu, tpr_pu_sd = x$tpr_pu_sd)
      }else{
        dat = data.frame(fpr_pu = x$fpr_pu, tpr_pu = x$tpr_pu)
      }
      
      roc.pumax = function(x,py1) sapply(1:length(x), function(i) ifelse(x[i] <= py1, x[i]/ py1, 1))
      tpr_pu.max = roc.pumax(dat$fpr_pu, py1[j])
      
      if(is.null(tol)){
        # if tol == NULL,
        if( !(is.null(x$tpr_pu_sd) | any(is.na(dat$tpr_pu_sd))) ){
          tol = dat$tpr_pu_sd
        }else{
          tol = 1e-4
        }
      }
      
      r= all(tpr_pu.max +tol >= dat$tpr_pu)
      r
      #ggplot(dat,aes(fpr_pu,tpr_pu))+geom_line()+geom_line(aes(fpr_pu,tpr_pu.max))+geom_line(aes(fpr_pu,tpr_pu.s),color="red")
    })
  
  aucs = sapply(1:length(py1) ,function(j) roc_curves[[j]]$auc_pu )
  if(all(isROCcontained==FALSE)){
    warning("all curves are not contained in the max ROC curve defined by py1 value. Redo the search with a smaller py1 value")
    py1.opt = NULL
  }else{
      viable.idx = which(isROCcontained)
      py1.opt = py1[which(aucs[viable.idx] == max(aucs[viable.idx]))]
  }
  return(py1.opt)
}