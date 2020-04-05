#' Search for optimal py1 values using roc_curves
#'
#'@import pbapply
#'@importFrom stats smooth.spline predict
#'@param roc_curves a list containing roc_curves of length(py1)
#'@param py1 a sequence of py1
#'@param tol a numerical number which allows small margin in determining whether the estimated ROC curve is smaller than the maximal curve.
#'@param is.smooth.roc a logical value; if TRUE (default), a smoothed roc curve will be used in checking whether the roc curve is contained by the maximal roc curve at each py1
#'@param verbose a logical value
#'@export
optimal_py1 = function(roc_curves, py1, tol = 1e-6, is.smooth.roc = T,verbose=T){
  
  if(verbose) {cat("choosing an optimal py1 value...\n")}else{pboptions(type="none")}
  isROCcontained = pbsapply(X =  1:length(py1), FUN = function(j) {
      
      x = roc_curves[[j]]
      dat = data.frame(fpr_pu = x$fpr_pu, tpr_pu = x$tpr_pu)
      
      if(is.smooth.roc) {
        # cat("this is executed\n")
        smooth.roc  = with(dat,smooth.spline(x = fpr_pu, y= tpr_pu))
        dat.sroc = predict(smooth.roc,x = dat$fpr_pu)
        idx = which(0 < dat$fpr_pu & dat$fpr_pu <=1 & 0 < dat.sroc$y & dat.sroc$y <=1) # keep observations within a valid range
        dat = dat[idx,]
        tpr_pu.s = dat.sroc$y[idx]
      }else{
        idx = which(0 < dat$fpr_pu & dat$fpr_pu <=1) 
        dat = dat[idx,]
      }
      
      roc.pumax = function(x,py1) sapply(1:length(x), function(i) ifelse(x[i] <= py1, x[i]/ py1, 1))
      tpr_pu.max = roc.pumax(dat$fpr_pu, py1[j])
      
      if(is.smooth.roc){
        r = all(tpr_pu.max +tol > tpr_pu.s)
      }else{
        r= all(tpr_pu.max +tol > dat$tpr_pu)
      }
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