#'@import magrittr
#'@import data.table
#'@export
create_grouped_dat <- function(seq_l,seq_u){
  # create an aggregated data with (sequence, number of sequence) for labeled and unlabeld
  g_l = seq_l[,.N,by=names(seq_l)] %>% setnames("V1","sequence") %>% setnames("N","labeled")
  g_u = seq_u[,.N,by=names(seq_u)] %>% setnames("V1","sequence") %>% setnames("N","unlabeled")
  # merge two aggregated datasets
  g_lu = merge(g_l,g_u,by="sequence",all = T)
  g_lu[is.na(labeled),"labeled"] = 0
  g_lu[is.na(unlabeled),"unlabeled"] =0 
  # check an error
  if(!sum(g_lu[,"labeled"])==length(seq_l$V1)){stop("error")}
  if(!sum(g_lu[,"unlabeled"])==length(seq_u$V1)){stop("error")}
  # assign a sequence ID: 1,...,number of unique sequences
  g_lu$seqId = 1:nrow(g_lu)
  g_lu
}
