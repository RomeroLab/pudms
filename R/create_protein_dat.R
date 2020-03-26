#'@import magrittr
#'@import data.table
#'@import foreach
#'@export
create_protein_dat = function(path_l,path_u,WT,type = c("mutations","sequences")){
  
  type = match.arg(type,choices =  c("mutations","sequences") )
  if(type == "mutations"){
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
  
  create_grouped_dat(seq_l,seq_u)
}

