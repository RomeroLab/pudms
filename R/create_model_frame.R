#' Create a model frame list from a grouped data set
#' 
#' @param grouped_dat a grouped data set
#' @param order an integer; 1= main effects, 2= main effects + pairwise effects
#' @param aggregate a logical value
#' @param refstate a character which will be used for the common reference state; the default is to use the most frequent amino acid as the reference state for each of the position. 
#' @param verbose a logical value
#' @return a list (X,z,wei,seqId,blockidx, refstate)
#' @import Matrix
#' @import magrittr
#' @import pbapply
#' @importFrom stats relevel
#' @export
create_model_frame <- function(grouped_dat, order = 1, aggregate = T, refstate = NULL,verbose=T){
  
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
  if(!verbose){pboptions(type ="none")} # suppress a progress bar if verbose == F
  seqmat = data.frame(pblapply(1:seqlen,function(i){substr(sequence,i,i)}))
  colnames(seqmat) = paste("P",0:(seqlen-1),".",sep="")
  
  # check number of unique factors in each position
  pbo = pboptions()
  if(verbose){cat("check number of unique factors in each position\n")}
  if(!verbose){pboptions(type ="none")}
  cidx = pblapply(seqmat, function(x){length(levels(x))}) %>% unlist>1
  seqmat = seqmat[,cidx]
  seqlen = ncol(seqmat)
  
  if(is.null(refstate)){
    if(verbose){cat("obtain ``reference`` amino-acid states\n")}
    # Obtain "reference" amino-acid states
    # The most frequent state
    refstates <- pblapply( seqmat, function(x){ tbl = table(x); names(which(tbl == max(tbl)))[1] }) %>% unlist
  }else{
    if(length(refstate)==1){
      all_states = Reduce(intersect, lapply(seqmat,levels)) # levels which exist in all columns
      if(!refstate %in% all_states){stop("refstate amino acid does not exist in all positions")}
      refstates = rep(refstate, seqlen)
    }else{
      if(length(refstate)!= seqlen) stop("length(refstate) should be equal to 1 or length of protein")
      all_states = lapply(seqmat,levels)
      check_ext = sapply(1:length(refstate), function(i) length(intersect(refstate[i],all_states[[i]])))
      if(!all(check_ext==1)) stop("refstate amino acid does not exist in all positions")
      refstates = refstate
    }
    
  }
  # give names to refstates
  names(refstates) = paste("P",1:seqlen,sep="")
  
  # relevel each column of the seqmat so that ref = refstate
  seqmat_colname = colnames(seqmat)
  seqmat_colname = paste(refstates, gsub(seqmat_colname,pattern = "P",replacement = ""), sep="")
  seqmat = data.frame (lapply(1:seqlen, function(i){relevel(seqmat[,i],ref = refstates[i])} ))
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
    res = list(X=X,z=z,wei=wei,seqId =seqId, blockidx = blockidx, refstate=refstates)
    
  }else{
    z = with(grouped_dat, c(rep(1,sum(labeled)), rep(0,sum(unlabeled)) ))
    seqId = c(with(grouped_dat, rep(seqId,labeled)), 
              with(grouped_dat, rep(seqId,unlabeled)))
    res = list(X=X,z=z,seqId =seqId, blockidx = blockidx, refstate = refstates)
  }
  res
}


