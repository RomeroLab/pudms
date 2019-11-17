#'@import data.table
#'@import progress
#'
mut_to_seq = function(txtdat,WT){
  dat =strsplit(txtdat,split = ",")
  WT = foreach(i=1:nchar(WT),.combine = "c")%do% {substr(WT,i,i)}
  
  pb = progress::progress_bar$new(total = length(dat))
  seqdat = foreach(x=1:length(dat),.combine = "rbind")%do%{
    pb$tick()
    nmuts = length(dat[[x]])
    newseq = WT
    for(i in 1:nmuts){
      m = dat[[x]][i]
      mfrom=substr(gsub(m,pattern = "[0-9]",replacement = ""),1,1)
      minto=substr(gsub(m,pattern = "[0-9]",replacement = ""),2,2)
      pos = as.numeric(gsub(m,pattern = "[^0-9]",replacement = ""))
      if(newseq[pos+1]!=mfrom) stop("a residue in the reference sequence should match with the first letter in mutations")
      newseq[pos+1]=minto
    }
    newseq
  }
  rownames(seqdat)=1:length(dat)
  seq = sapply(1:nrow(seqdat), function(i) paste(seqdat[i,],collapse = ""))
  data.table(data.frame(V1=seq,stringsAsFactors = F))
}
