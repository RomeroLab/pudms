extractIndices<-function(coef, position, excludeStates = "*", order = 1, position2 = NULL,verbose = TRUE){
  
  extract_pos = function(x){x[1]%>% gsub(pattern = "[^0-9]",replacement = "") %>% unlist}
  extract_aa = function(x){x[2]}
  allStates <- c("*", "A", "C", "D", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "E", "Y")
  includeStates = allStates[!allStates%in%excludeStates]
  colnames<-rownames(coef)[-1]
  
  if(order == 1){
    posm = colnames %>% strsplit(split="\\.") %>% sapply(extract_pos) %>% as.numeric()
    aa = colnames %>% strsplit(split="\\.") %>% sapply(extract_aa) 
    
    idxx= which((posm %in% position) & (aa %in% includeStates)) + 1 # 1 for intercept
    
  }else{
    
    pos1 = colnames %>%  strsplit(split=":") %>% lapply(function(x){x[1]}) %>% unlist
    pos2 = colnames %>%  strsplit(split=":") %>% lapply(function(x){x[2]}) %>% unlist
    midx = which(is.na(pos2))
    
    aam = pos1[midx] %>% strsplit(split="\\.") %>% sapply(extract_aa) 
    aa1 = pos1[-midx] %>% strsplit(split="\\.") %>% sapply(extract_aa) 
    aa2 = pos2[-midx] %>% strsplit(split="\\.") %>% sapply(extract_aa) 
    
    posm = pos1[midx] %>% strsplit(split="\\.") %>% sapply(extract_pos)  %>% as.numeric()
    pos1 = pos1[-midx] %>% strsplit(split="\\.") %>% sapply(extract_pos)  %>% as.numeric()
    pos2 = pos2[-midx] %>% strsplit(split="\\.") %>% sapply(extract_pos)  %>% as.numeric()
    
    if(is.null(position2)){
      filled_posm = c(posm, rep(-10, length(pos1))) 
      filled_aam  = c(aam,  rep(-10, length(aa1 )))
      idxx = which((filled_posm %in% position) & ( filled_aam %in% includeStates))+1
      
    }else{
      filled_pos1  = c(rep(-10,length(posm)),  pos1)
      filled_pos2  = c(rep(-10,length(posm)),  pos2)
      filled_aa1   = c(rep(-10,length(aam) ),  aa1 )
      filled_aa2   = c(rep(-10,length(aam) ),  aa2 )
      idxx= which((filled_pos1 %in% position) & (filled_pos2 %in% position2 ) & (filled_aa1 %in% includeStates) & (filled_aa2 %in% includeStates))+1
    }
    
  }
  if(verbose) cat("indices corresponding to:\n",paste(rownames(coef[idxx,,drop=F]),collapse = ","),"\n")
  return(idxx)
}
