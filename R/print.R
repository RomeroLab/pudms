#' @export
#' @method print pudms.fit
print.pudms.fit<-function(x,...){
  print(x$result_table[order(abs(x$result_table$p)),][1:10,])
}
