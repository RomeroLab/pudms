#' Create a sparse model matrix in a parallel manner
#'@import doParallel
#'@import Matrix
#'@param order an integer; 1= main effects, 2= main effects + pairwise effects
#'@param data data to build a model matrix from
#'@param nCores number of threads for parallel computing
#'@export
sparse.model.matrix.parallel<-function(order = c(1,2), data, nCores =1, verbose = TRUE){
  
  data = lapply(1:ncol(data), function(x){
    newlevels = paste0(names(data[,x,drop=F]),levels(data[,x]))
    levels(data[,x]) = newlevels
    data[,x]
  }
  )
  data = data.frame(data)
  
  if(nCores>1){
    clustertype = ifelse(.Platform$OS.type=="windows", 'PSOCK', 'FORK') 
    cl <- makeCluster(nCores,type = clustertype, outfile = "")
    
    doParallel::registerDoParallel(cl)
    on.exit(stopCluster(cl))
  }else{
    registerDoSEQ()
  }
  
  X1=foreach(i=1L:ncol(data),.combine = "cbind")%dopar%{
    # print(i)
    t(Matrix::fac2sparse(data[,i],drop.unused.levels = TRUE)[-1,,drop=F])
    }
  
  if(order==2){
    
    uniqueLevels=foreach(i=1:(ncol(data)-1),.combine = "rbind")%:%foreach(j=(i+1):(ncol(data)),.combine = "rbind")%dopar%{
  
      a=interaction(data[,i],data[,j],sep = ":")
      baselevels = unique(c(paste(levels(data[,i])[1],levels(data[,j]),sep=":"),paste(levels(data[,i]),levels(data[,j])[1],sep=":")))
      levels(a)[levels(a)%in%baselevels] = "base"
      a=droplevels(a)
      c(i,j,length(levels(a)))
    }
    uniqueLevels = uniqueLevels[uniqueLevels[,3]>1,]
    
    X2=foreach(k=1:nrow(uniqueLevels),.combine = "cbind")%dopar%{

      i = uniqueLevels[k,1]; j = uniqueLevels[k,2]
      # cat(i,j,"\n")
      a=interaction(data[,i],data[,j],sep = ":")
      baselevels = unique(c(paste(levels(data[,i])[1],levels(data[,j]),sep=":"),paste(levels(data[,i]),levels(data[,j])[1],sep=":")))
      levels(a)[levels(a)%in%baselevels] = "base"
      # a=droplevels(a)
      # length(levels(a))
      t(Matrix::fac2sparse(a,drop.unused.levels = T)[-1,,drop=F])
    }
    
    # X2=foreach(i=1:(ncol(data)-1),.combine = "cbind")%:%foreach(j=(i+1):(ncol(data)),.combine = "cbind")%do%{
    #   cat(i,j,"\n")
    #   a=interaction(data[,i],data[,j],sep = ":")
    #   baselevels = unique(c(paste(levels(data[,i])[1],levels(data[,j]),sep=":"),paste(levels(data[,i]),levels(data[,j])[1],sep=":")))
    #   levels(a)[levels(a)%in%baselevels] = "base"
    #   
    #   t(Matrix::fac2sparse(a,drop.unused.levels = T)[-1,])
    # }
    X=cbind(X1,X2)
  }else{
    X = X1
  }
  
  return(X)
}
