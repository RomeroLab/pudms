#' Create a model frame list from a grouped data set
#' 
#' @param grouped_dat a grouped data set
#' @param order an integer; 1= main effects, 2= main effects + pairwise effects
#' @param aggregate a logical value
#' @param refstate a character which will be used for the common reference state; the default is to use the most frequent amino acid as the reference state for each of the position. 
#' @param nCores the number of threads for computing.
#' @param verbose a logical value
#' @return a list (X,z,wei,seqId,blockidx, refstate)
#' @import Matrix
#' @import magrittr
#' @import pbapply
#' @importFrom stats relevel
#' @export
create_model_frame <- function(grouped_dat,
                               order = 1,
                               aggregate = T,
                               refstate = NULL,
                               nCores = 1,
                               verbose = T) {
  
  
  # create a seqmat
  list.seqmat<- create_seqmat(
    grouped_dat=grouped_dat,
    aggregate = aggregate,
    refstate = refstate,
    verbose = verbose
  )
  
  seqmat =list.seqmat$seqmat
  refstate = list.seqmat$refstate
  
  if(verbose){cat("convert to the sparse one-hot-encoding model matrix\n")}
  if(order>2){stop("order>2 is not supported yet")}
  if(nCores>1){
    X = sparse.model.matrix.parallel(order = order,data = seqmat,nCores = nCores,verbose=verbose)
  }else{
    # convert to the model matrix
    if(order==1){
      X = sparse.model.matrix(~.,data = seqmat)
    }else if(order==2){
      X = sparse.model.matrix(~.^2,data = seqmat)
    }
  }
  
  X = X[,-1]
  colnames = colnames(X)
  extract_pos = function(x){x[1]%>% gsub(pattern = "[^0-9]",replacement = "") %>% unlist}
  # unique block
  if(order == 1){
    blockidx = colnames %>% strsplit(split="\\.") %>% sapply(extract_pos) 
  }else if(order==2){
    pos1 = colnames %>%  strsplit(split=":") %>% lapply(function(x){x[1]}) %>% unlist
    pos2 = colnames %>%  strsplit(split=":") %>% lapply(function(x){x[2]}) %>% unlist
    midx = which(is.na(pos2))
    posm = pos1[midx] %>% strsplit(split="\\.") %>% sapply(extract_pos) 
    pos1 = pos1[-midx] %>% strsplit(split="\\.") %>% sapply(extract_pos) 
    pos2 = pos2[-midx] %>% strsplit(split="\\.") %>% sapply(extract_pos) 
    blockidx = c(posm, paste(pos1,pos2,sep = ":"))
  }
  
  K = nrow(grouped_dat) # number of unique sequence
  
  if(aggregate){
    
    X = rbind(X,X)
    z = c(rep(1,K), rep(0,K))
    wei = with(grouped_dat , c(labeled,unlabeled))
    seqId = with(grouped_dat, c(seqId, seqId))
    
    # remove columns where wei == 0
    ridx = which(wei!=0)
    X = X[ridx,]
    z = z[ridx]
    seqId = seqId[ridx]
    wei = wei[ridx]
    
    # summary
    res = list(X=X,z=z,wei=wei,seqId =seqId, blockidx = blockidx, refstate=refstate)
    
  }else{
    z = with(grouped_dat, c(rep(1,sum(labeled)), rep(0,sum(unlabeled)) ))
    seqId = c(with(grouped_dat, rep(seqId,labeled)), 
              with(grouped_dat, rep(seqId,unlabeled)))
    res = list(X=X,z=z,seqId =seqId, blockidx = blockidx, refstate = refstate)
  }
  res
}


