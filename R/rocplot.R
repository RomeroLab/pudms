#' Draw an adjusted ROC curve
#'
#'@param roc_curve an roc curve from the output of the function adjusted_roc_curve()
#'@param py1 a numeric value
#'@import PRROC
#'@import ggplot2
#'@export
#'
rocplot = function(roc_curve,py1){
  
  x = seq(0,1,length.out = length(roc_curve$fpr_pu));
  ymax = (1/py1)*x*1*(x<py1)+1*(x>=py1)
  auc = unique(roc_curve$auc)
  
  dat= data.frame(fpr_pu = roc_curve$fpr_pu,
                  tpr_pu = roc_curve$tpr_pu,
                  fpr = roc_curve$fpr,
                  tpr = roc_curve$tpr)
  
  ggplot(dat,aes(.data$fpr_pu,.data$tpr_pu))+
    # geom_ribbon(aes(x = x, ymin = 0, ymax = ymax), alpha=0.3)+
    geom_path(size=1.5,color="grey40")+
    geom_path(aes(.data$fpr,.data$tpr),linetype="solid",color="black",size=1.5)+
    geom_abline(slope = 1, linetype="dashed",color="darkblue",size=1.5)+
    theme_classic(base_size = 15)+
    xlab("false-positive rate")+
    ylab("true-positive rate")+
    scale_x_continuous(expand=c(0,0)) +
    #custom breaks on y-axis
    scale_y_continuous(expand=c(0,0))+
    annotate("text",label=paste("auc (corrected): ",round(auc,3),sep=""), x=0.75,y=0.1,size=5)+
    theme(aspect.ratio = 1,plot.margin = unit(c(1,1,1,1), "cm"))+
    coord_fixed()
  
}
