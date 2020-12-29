#' Create a model frame list from a grouped data set
#' 
#' @param grouped_dat a grouped data set
#' @param aggregate a logical value
#' @param refstate a character which will be used for the common reference state; the default is to use the most frequent amino acid as the reference state for each of the position. 
#' @param verbose a logical value
#' @return a list (seqmat, refstate)
#' @import Matrix
#' @import magrittr
#' @import pbapply
#' @importFrom stats relevel
#' @export

create_seqmat<-function(grouped_dat, aggregate = T, refstate = NULL,verbose=T){
  # obtain sequences from grouped_dat
  if(aggregate){
    sequence = with(grouped_dat, sequence)
  }else{
    sequence = c(with(grouped_dat, rep(sequence,labeled)), 
                 with(grouped_dat, rep(sequence,unlabeled)))
  }
  
  seqlen = nchar(grouped_dat$sequence[1]) # number of positions
  
  
  if(verbose){cat("* create a sequence matrix\n")}
  # create seq mat (each row is a length seqlen vector)
  if(!verbose){pboptions(type ="none")} # suppress a progress bar if verbose == F
  seqmat = data.frame(pblapply(1:seqlen,function(i){substr(sequence,i,i)}),stringsAsFactors = TRUE)
  colnames(seqmat) = paste("P",0:(seqlen-1),".",sep="")
  
  if(is.null(refstate)){
    if(verbose){cat("* obtain ``reference`` amino-acid states\n")}
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
  names(refstates) = colnames(seqmat)
  
  # change "P" in the column names into the reference states
  seqmat_colname = colnames(seqmat)
  seqmat_colname = paste(refstates, gsub(seqmat_colname,pattern = "P",replacement = ""), sep="")
  
  # relevel each column of the seqmat so that ref = refstate
  seqmat = data.frame (lapply(1:seqlen, function(i){relevel(seqmat[,i],ref = refstates[i])} ))
  colnames(seqmat) = seqmat_colname
  
  # check number of unique factors in each position
  pbo = pboptions()
  if(verbose){cat("* check number of unique factors in each position\n")}
  if(!verbose){pboptions(type ="none")}
  cidx = pblapply(seqmat, function(x){length(levels(x))}) %>% unlist>1
  
  # keep positions with varying states
  seqmat = seqmat[,cidx]
  
  return(list(seqmat=seqmat,refstate = refstates))
}