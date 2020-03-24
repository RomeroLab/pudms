#' Create a model frame list from a grouped data set
#' 
#' @param grouped_dat a grouped data set
#' @param order an integer of 1 or 2; main effects (order=1) or interaction effects (order=2) in X
#' @param aggregate a logical value
#' @param basestate a letter
#' @param verbose a logical value
#' @return a list (X,z,wei,seqId,blockidx)
#'@import Matrix
#'@import magrittr
#'@import pbapply
#'@import foreach
#'@export
create_model_frame <- function(grouped_dat, order = 1, aggregate = T, basestate = NULL,verbose=T){
  
  # obtain sequences from grouped_dat
  if(aggregate){
    sequence = with(grouped_dat, sequence)
  }else{
    sequence = c(with(grouped_dat, rep(sequence,labeled)), 
                 with(grouped_dat, rep(sequence,unlabeled)))
  }
  
  seqlen = nchar(grouped_dat$sequence[1]) # number of positions
  K = nrow(grouped_dat) # number of unique sequence
  
  if(verbose){cat("create a sequence matrix\n")}
  # create seq mat (each row is a length seqlen vector)
  seqmat = data.frame(pblapply(1:seqlen,function(i){substr(sequence,i,i)}))
  colnames(seqmat) = paste("P",0:(seqlen-1),".",sep="")
  
  # check number of unique factors in each position
  if(verbose){cat("check number of unique factors in each position\n")}
  cidx = pblapply(seqmat, function(x){length(levels(x))}) %>% unlist>1
  seqmat = seqmat[,cidx]
  seqlen = ncol(seqmat)
  if(verbose){cat("obtain ``base`` amino-acid states\n")}
  if(is.null(basestate)){
    # Obtain "base" amino-acid states
    # The most frequent state
    basestates <- pblapply( seqmat, function(x){ tbl = table(x); names(which(tbl == max(tbl))) }) %>% unlist
  }else{
    all_states = Reduce(intersect, lapply(seqmat,levels)) # levels which exist in all columns
    if(!basestate %in% all_states){stop("wrong basestate")}
    basestates = rep(basestate, seqlen)
  }
  
  # relevel each column of the seqmat so that ref = basestate
  seqmat_colname = colnames(seqmat)
  seqmat_colname = paste(basestates, gsub(seqmat_colname,pattern = "P",replacement = ""), sep="")
  seqmat = data.frame (lapply(1:seqlen, function(i){relevel(seqmat[,i],ref = basestates[i])} ))
  colnames(seqmat) = seqmat_colname
  
  if(verbose){cat("convert to the sparse one-hot-encoding model matrix\n")}
  # convert to the model matrix
  if(order==1){
    X = sparse.model.matrix(~.,data = seqmat)
  }else if(order==2){
    X = sparse.model.matrix(~.^2,data = seqmat)
  }else{
    stop("order>2 is not supported yet")
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
    res = list(X=X,z=z,wei=wei,seqId =seqId, blockidx = blockidx)
    
  }else{
    z = with(grouped_dat, c(rep(1,sum(labeled)), rep(0,sum(unlabeled)) ))
    seqId = c(with(grouped_dat, rep(seqId,labeled)), 
              with(grouped_dat, rep(seqId,unlabeled)))
    res = list(X=X,z=z,seqId =seqId, blockidx = blockidx)
  }
  res
}


