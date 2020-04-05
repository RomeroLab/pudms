#'Return an adjusted auc value from a fitted model at py1 ==py1[py1idx]
#'@param vfit vpudms.fit object
#'@param py1idx an index of py1 for auc to be returned at
#'@export
#'
auc<-function(vfit,py1idx=NULL){
  
  #if py1idx is not null, check py1idx is valid
  if(!is.null(py1idx)){
    if(py1idx%in%1:length(vfit$py1)) stop("py1idx should be between 1 and length(input$py1)")
  }
    
  #if py1idx is null, py1.opt has to be not a null value
  if(is.null(py1idx) & is.null(vfit$py1.opt)){stop("py1.opt is NULL")}
  
  if(is.null(py1idx)){
    aucval = with(vfit, roc_curves[[which(py1 == py1.opt)]]$auc)
  }else{
    aucval = with(vfit, roc_curves[[py1idx]]$auc)
  }
  
  return(aucval)
  
}