#'@import PRROC
#'@import ggplot2
#'@export
#'
rocplot = function(roc_curve,py1, corrected = TRUE){
  
  x = seq(0,1,length.out = length(roc_curve$fpr_pu));
  ymax = (1/py1)*x*1*(x<py1)+1*(x>=py1)
  auc = unique(roc_curve$auc)
  
  ggplot(roc_curve,aes(fpr_pu,tpr_pu))+
    # geom_ribbon(aes(x = x, ymin = 0, ymax = ymax), alpha=0.3)+
    geom_path(size=1.5,color="grey40")+
    geom_path(aes(fpr,tpr),linetype="solid",color="black",size=1.5)+
    geom_abline(slope = 1, linetype="dashed",color="darkblue",size=1.5)+
    theme_classic(base_size = 20)+
    xlab("false-positive rate")+
    ylab("true-positive rate")+
    scale_x_continuous(expand=c(0,0)) +
    #custom breaks on y-axis
    scale_y_continuous(expand=c(0,0))+
    annotate("text",label=paste("auc (corrected): ",round(auc,3),sep=""), x=0.7,y=0.1,size=6)+
    theme(aspect.ratio = 1)+
    coord_fixed()
  
}
