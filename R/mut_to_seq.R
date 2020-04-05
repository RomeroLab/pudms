#'@import data.table
#'@import pbapply
#'@import doParallel
mut_to_seq = function(txtdat,WT,nCores = 1){
  dat =strsplit(txtdat,split = ",")
  
  WT = sapply(X=1:nchar(WT),FUN = function(i) substr(WT,i,i))
  
  if(nCores>1){
    isParallel = TRUE
    nCores = min(nCores,detectCores())
    cl <- makeCluster(nCores)
    registerDoParallel(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())
    parallel::clusterExport(cl= cl,varlist = c("dat","WT"))
  }
  
  # for each row, we create a new sequence
  
  FUN = function(x){
    nmuts = length(dat[[x]]); 
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
  
  result = pblapply(X = 1:length(dat),FUN = FUN,cl = cl)
  seqdat = do.call("rbind",result)
  
  rownames(seqdat)=1:length(dat)
  seq = sapply(1:nrow(seqdat), function(i) paste(seqdat[i,],collapse = ""))
  data.table(data.frame(V1=seq,stringsAsFactors = F))
}
