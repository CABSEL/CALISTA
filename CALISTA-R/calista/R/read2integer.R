#' A calista Function
#'
#' @param
#' @keywords calista
#' @export
#' @examples
#' read2integer()



read2integer<-function(){
  writeLines('Specify the node pairs (e.g. 4 5):')
  Edge_1=strsplit(readLines(n=1), " ")
  #return(as.integer(n))
  Edge=c(as.integer(Edge_1[[1]][1]),as.integer(Edge_1[[1]][2]))
  return(Edge)
}




