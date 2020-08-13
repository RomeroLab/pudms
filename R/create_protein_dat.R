#' Create a data set for protein
#'
#'@import magrittr
#'@import data.table
#'@importFrom utils read.table
#'@import foreach
#'@param path_l the file path for the labeled sequences. File type can be txt or txt.gzip
#'@param path_u the file path for the unlabeled sequences 
#'@param type type of the input files; "sequences" or "mutations" based on the type of the input; a default is "sequences"
#'@param WT a reference sequence if input type is "mutations"
#'@return a data table containing (sequence, labeled, unlabeled, seqId)
#'@export
create_protein_dat = function(path_l,path_u,type = "sequences",WT=NULL){
  
  type = match.arg(type,choices =  c("mutations","sequences") )
  if(type == "mutations"){
    if(is.null(WT)){stop("if type==mutations, WT needs to be provided")}
    txtdat_l = read.table(file = path_l,colClasses = "character")
    txtdat_u = read.table(file = path_u,colClasses = "character")
    cat("convert mutations into sequences for a labeled set\n")
    seq_l=mut_to_seq(txtdat = txtdat_l$V1,WT = WT)
    cat("convert mutations into sequences for an unlabeled set\n")
    seq_u=mut_to_seq(txtdat = txtdat_u$V1,WT = WT)
  }else{
    seq_l = fread(file = path_l,header = F)
    seq_u = fread(file = path_u,header = F)
  }
  
  # create an aggregated data with (sequence, number of sequence) for labeled and unlabeld
  g_l = seq_l[,.N,by=names(seq_l)] %>% setnames("V1","sequence") %>% setnames("N","labeled")
  g_u = seq_u[,.N,by=names(seq_u)] %>% setnames("V1","sequence") %>% setnames("N","unlabeled")
  # merge two aggregated datasets
  g_lu = merge(g_l,g_u,by="sequence",all = T)
  g_lu[is.na(g_lu$labeled),"labeled"] = 0
  g_lu[is.na(g_lu$unlabeled),"unlabeled"] =0 
  # check an error
  if(!sum(g_lu[,"labeled"])==length(seq_l$V1)){stop("error")}
  if(!sum(g_lu[,"unlabeled"])==length(seq_u$V1)){stop("error")}
  # assign a sequence ID: 1,...,number of unique sequences
  g_lu$seqId = 1:nrow(g_lu)
  g_lu
}

