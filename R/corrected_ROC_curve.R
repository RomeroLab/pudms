#' Obtain a corrected ROC curve from rocdata
#' 
#' @import PRROC
#' @import PUlasso
#' @param coef a named vector for the fitted coefficient vector
#' @param test_grouped_dat a data table (sequence, labeled, unlabeled, seqId) containing validation examples
#' @param order a considered order of effects
#' @param basestate a character value
#' @param verbose a logical value
#' @param py1 the prevalence of positives in the unlabeled set (which will be used for a correction)
#' @return a list containing a data frame for generating a roc plot and the plot
#'@export
corrected_roc_curve = function(rocdata = NULL,
                               coef =NULL,
                               test_grouped_dat = NULL,
                               order=1,
                               basestate = NULL,
                               verbose=T,
                               py1
                               ) {
  
  # if rocdata is not provided, we create a new one
  if(is.null(rocdata)){
  rocdata = ROCdata(coef = coef,
                    test_grouped_dat = test_grouped_dat,
                    order = order,
                    basestate = basestate,
                    verbose = verbose)
  }
  
  # obtain an uncorrected ROC curve
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
  
  # return
  list(roc_curve= roc_curve, rocplot = rocplot(roc_curve,py1, corrected = TRUE))
}
