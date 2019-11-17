#' @export
#' @method coef pudms.fit
coef.pudms.fit<-function(object,...){
  coef(object$fitted)
}