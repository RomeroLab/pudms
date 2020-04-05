#' @export
#' @method plot vpudms.fit
plot.vpudms.fit<-function(x,y,...){
  idx=which(x$py1 == x$py1.opt)
  rocplot(x$roc_curves[[idx]],x$py1.opt)
}