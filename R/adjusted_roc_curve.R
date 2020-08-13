#' Obtain an adjusted ROC curve from rocdata
#'
#' @import PRROC
#' @import PUlasso
#' @param rocdata NULL or an output from the rocdata function
#' @param coef a named vector for the fitted coefficient vector
#' @param test_grouped_dat a data table (sequence, labeled, unlabeled, seqId) containing validation examples
#' @param order a considered order of effects
#' @param refstate a character which will be used for the common reference state; the default is to use the most frequent amino acid as the reference state for each of the position.
#' @param verbose a logical value
#' @param py1 the prevalence of positives in the unlabeled set (which will be used for a correction)
#' @param plot a logical value. If TRUE, a plot will be returned
#' @return a list containing a data frame for generating a roc plot and the plot
#' @export
adjusted_roc_curve = function(rocdata = NULL,
                              coef = NULL,
                              test_grouped_dat = NULL,
                              Xprotein_test = NULL,
                              order = 1,
                              refstate = NULL,
                              verbose = T,
                              py1,
                              plot = TRUE) {
  # if rocdata is not provided, we create a new one
  if (is.null(rocdata)) {
    rocdata = rocdata(
      coef = coef,
      test_grouped_dat = test_grouped_dat,
      Xprotein_test = Xprotein_test,
      order = order,
      refstate = refstate,
      verbose = verbose
    )
  }
  
  # obtain an unadjusted ROC curve
  r = roc.curve(
    scores.class0 = rocdata$pr_test,
    weights.class0 = rocdata$labeled,
    scores.class1 = rocdata$pr_test,
    weights.class1 = rocdata$unlabeled,
    curve = T
  )
  roc_curve = data.frame(r$curve)
  colnames(roc_curve) = c("fpr_pu", "tpr_pu", "threshold")
  
  # modified ROC: Jain et al (2017)
  roc_curve$tpr = with(roc_curve, tpr_pu)
  roc_curve$fpr = with(roc_curve, (fpr_pu - py1 * tpr) / (1 - py1))
  
  roc_curve = as.list(roc_curve)
  roc_curve$auc_pu = r$auc
  roc_curve$auc = (r$auc - py1 / 2) / (1 - py1)
  
  if (plot) {
    rocplot = rocplot(roc_curve, py1)
  } else{
    rocplot = NULL
  }
  # return
  list(roc_data = rocdata,
       roc_curve = roc_curve,
       rocplot = rocplot)
}
