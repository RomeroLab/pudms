
is_equal<-function(x,y){
  if(length(x)!=length(y)){
    return(FALSE)
  }else{
    return(all(x==y))
  }
}

check_duplicated_columns = function(X){
  
  # Check the type of X
  is.sparse = FALSE
  if (inherits(X, "sparseMatrix")) {
    is.sparse = TRUE
    X = as(X, "CsparseMatrix")
    X = as(X, "dgCMatrix")
  } else if (inherits(X, "dgeMatrix")) {
    X = as.matrix(X)
  }
  if (!(inherits(X,"matrix") || inherits(X, "dgCMatrix") )) {
    stop("X must be a matrix or a sparse matrix")
  }
  
  
  # return dup_indices, res
  
  if(!is.sparse){
    # Dense matrix  
    # obtain indices for duplicated columns
    dup_indices = which(duplicated(t(X)))
    
    # check which columns are identical to columns corr.to dup_indices
    if(length(dup_indices)>0){
      res = foreach(i=dup_indices)%do%{
        v = which(apply(X, 2, is_equal, y = X[,i])) # indices for the same column
        v[v!=i]
      }
    }else{
      res = NULL
    }
    
    
  }else{
    # Sparse matrix
    n <- nrow(X)
    p <- ncol(X)
    J <- rep(1:p, diff(X@p))
    I <- X@i + 1
    x <- X@x
    ColLst = split(I,J)
    names(ColLst) = colnames(X)
    dup_indices = which(duplicated(ColLst))
    
    if(length(dup_indices)>0){
      res = foreach(i=dup_indices)%do%{
        v = which(sapply(ColLst, is_equal, y = ColLst[[i]])) # indices for the same column
        v[v!=i]
      }
      names(res) = colnames(X)[dup_indices]
    }else{
      res = NULL
    }
    
  }
  
  if(!is.null(res)){
    res = lapply(res, names)
  }
  
  
  return(list(duplicated_indices = dup_indices,duplicated_columns = res))
}

# for (i in 1:length(res)){
#   cat(names(res[i]),":",res[[i]],"\n")
# }

