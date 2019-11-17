#'@import PUlasso
#'@export
corrected_roc_curve = function(coef,
                               test_grouped_dat,
                               py1,
                               order=1) {
  
  rocdata = ROCdata(coef = coef,test_grouped_dat = test_grouped_dat,order = order)
  
  r = roc.curve(
    scores.class0 = rocdata$pr_test,
    weights.class0 = rocdata$labeled,
    scores.class1 = rocdata$pr_test,
    weights.class1 = rocdata$unlabeled,
    curve = T
  )
  roc_curve = data.frame(r$curve)
  colnames(roc_curve) = c("fpr_pu","tpr_pu","threshold")
  # modified ROC: Jain et al (2017)
  
  roc_curve$tpr = with(roc_curve, tpr_pu)
  roc_curve$fpr = with(roc_curve, (fpr_pu - py1*tpr)/(1-py1))
  
  roc_curve$auc_pu = r$auc
  roc_curve$auc = (r$auc-py1/2)/(1-py1)
  
  list(roc_curve= roc_curve, rocplot = rocplot(roc_curve,py1, corrected = TRUE))
}
