#' Delete columns such that colSums <= nobs_thresh, and rows which contain a mutation into the deleted columns
#' 
#'@param Xprotein list (X,z,wei, seqId, blockId, refstate)
#'@param order an integer; 1= main effects, 2= main effects + pairwise effects
#'@param nobs_thresh an integer threshold
#'@param checkResFullRank a logical value 
#'@param verbose a logical value. The default is TRUE
#'@import Matrix
#'@import doParallel
#'@return a list (X,z,wei, seqId, blockId, refstate)
#'@export
filter_mut_less_than_thresh<-function(Xprotein, order= 1, nobs_thresh=10, checkResFullRank=T,verbose = T){
  
  if(verbose){cat("\n* remove sequences which contain a feature whose nmuts <= ",nobs_thresh,"\n", sep="")}
  
  if(order==1){
    midx = 1:ncol(Xprotein$X)
  }else{
    pos1 = Xprotein$blockidx  %>%  strsplit(split=":") %>% lapply(function(x){x[1]}) %>% unlist
    pos2 = Xprotein$blockidx  %>%  strsplit(split=":") %>% lapply(function(x){x[2]}) %>% unlist
    midx = which(is.na(pos2))
  }
  
  # we delete columns such that colSums (main-effects) <= nobs_thresh and rows which contain a mutation into the deleted columns
  #  filter.r1 = function(Xprotein,midx,nobs_thresh){
  
  wX=Matrix::Diagonal(x = Xprotein$wei)%*%Xprotein$X
  
  # return column indices to be deleted, and row indices to keep
  filter.r1 = function(wX, midx){
    
    ndelCols = 0; nchangeCols = 1; ridx = 1:nrow(wX)
    
    while(nchangeCols>0){
      cidx = which(Matrix::colSums(wX[ridx,midx]) <= nobs_thresh) # column idices to be deleted
      if(length(cidx) == length(midx)) stop("all columns have nmuts <=nobs_nobs_thresh")
      
      ndelCols.new = length(cidx); # number of columns to be deleted
      nchangeCols = ndelCols.new - ndelCols # change from the previous ndelCols
      ndelCols = ndelCols.new 
      
      # print(nchangeCols)
      ridx = which(Matrix::rowSums(wX[,cidx,drop=F])==0) #ridx to keep
    }
    
    return(list(ridx=ridx,cidx=cidx))
    
  }
  
  indices = filter.r1(wX,midx)
  ridx = indices$ridx; 
  cidx = Matrix::colSums(wX[indices$ridx,]) > nobs_thresh
  
  res = list(X=Xprotein$X[ridx,cidx],
             z=Xprotein$z[ridx],
             wei=Xprotein$wei[ridx],
             seqId=Xprotein$seqId[ridx],
             blockidx=Xprotein$blockidx[cidx],
             refstate=Xprotein$refstate)
  
  nremoved_seq = (Xprotein$seqId %>% unique %>% length) - (res$seqId %>% unique %>% length)
  if(nremoved_seq >0 & verbose) {cat(nremoved_seq," unique sequences which contain a feature whose nmuts <= ",nobs_thresh," are removed\n",sep="")}
  
  # check duplicated columns
  
  if(verbose) cat("\n* check whether there exist any duplicated columns\n")
  check_dupcols = check_duplicated_columns(res$X)
  dup_cols = check_dupcols$duplicated_columns
  
  if(length(check_dupcols$duplicated_indices)>0){
    # remove duplicated columns
    res$X = res$X[,-check_dupcols$duplicated_indices,drop=F]
    res$blockidx = res$blockidx[-check_dupcols$duplicated_indices]
    
    if(verbose){
      cat(length(check_dupcols$duplicated_indices),"duplicated columns are found and removed in the filtered X \n")
      for(i in 1:length(dup_cols)){
        cat(names(dup_cols[i]),"equals to", dup_cols[[i]],"\n")
      }
    }
    
  }else{
    if(verbose) cat("no duplicated columns in the filtered X \n")
    dup_cols = NULL
  }
  
  if(checkResFullRank){
    if(nrow(res$X) > ncol(res$X)){ # tall matrix
      if (verbose) cat("\n* check whether a filtered X is a full rank matrix\n")
      Xrank = rankMatrix(x = Matrix::crossprod(res$X),method = 'qr')
      if(Xrank!=ncol(res$X)){
        fullrank = FALSE
        warning("- filtered X is not a full rank matrix")
      }else{
        fullrank = TRUE
        if(verbose) cat("filtered X is a full rank matrix\n")}
    }else{
      if(verbose) cat("ncol(X) > nrow(X) \n")
      fullrank = FALSE
    }
    
  }else{
    fullrank = FALSE
    
  }
  
  res$fullrank = fullrank
  res$duplicated_columns = dup_cols
  return(res)
}    

