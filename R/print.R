#' @export
#' @method print pudms.fit
print.pudms.fit<-function(x,...){
  print(x$call)
}
#' @export
#' @method print vpudms.fit
print.vpudms.fit<-function(x,...){
  print(x$call)
}