#' @export
#' @method print pudms.fit
print.pudms.fit<-function(x,...){
  cat("head(result table):\n ")
  print(x$result_table[order(abs(x$result_table$p)),][1:10,])
}
