#' Delete columns such that colSums <= thresh, and rows which contain a mutation into the deleted columns
#' 
#'@param Xprotein list (X,z,wei, seqId, blockId)
#'@param thresh an integer threshold
#'@param checkResFullRank a logical value 
#'@return a list (X,z,wei, seqId, blockId)
#'@export
filter_mut_less_than_thresh<-function(Xprotein, thresh=10, checkResFullRank=T){
  wX=Diagonal(x = Xprotein$wei)%*%Xprotein$X
  cidx_mut_less_thresh = colSums(wX)<=thresh
  ridx_mut_less_thresh = rowSums(wX[,cidx_mut_less_thresh])>0
  res = list(X=Xprotein$X[!ridx_mut_less_thresh,!cidx_mut_less_thresh],
             z=Xprotein$z[!ridx_mut_less_thresh],
             wei=Xprotein$wei[!ridx_mut_less_thresh],
             seqId=Xprotein$seqId[!ridx_mut_less_thresh],
             blockidx=Xprotein$blockidx[!cidx_mut_less_thresh])
  nremoved_seq = Xprotein$seqId[ridx_mut_less_thresh] %>% unique %>% length
  cat(nremoved_seq,"unique sequences which contain a feature whose nmuts <",thresh," are removed\n")
  if(checkResFullRank){
    cat("check whether a filtered X is a full rank matrix\n")
    q=qr(as.matrix(t(res$X)%*%res$X))
    if(q$rank!=ncol(res$X)){warning("filtered X is not a full rank matrix")
      }else{cat("filtered X is a full rank matrix\n")}
  }
  
  return(res)
}