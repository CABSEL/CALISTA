#' A calista Function
#'
#' @param i,MaxNumberOfEdges
#' @keywords calista
#' @export
#' @examples
#' CheckNumberOfEdge()


CheckNumberOfEdges<-function(i,MaxNumberOfEdges){
  if(i>MaxNumberOfEdges || i==0)
  {
    NoAdding=TRUE
  }else{
    NoAdding=FALSE
  }
}
