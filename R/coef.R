#' @import PUlasso
#' @importFrom stats coef
#' @export
#' @method coef pudms.fit
coef.pudms.fit<-function(object,...){
  coef(object$fitted)
}