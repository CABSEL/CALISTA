#' A calista Function
#'
#' @param population
#' @keywords calista
#' @export
#' @examples
#'
#' get_consensus()

get_consensus<-function(population){
  if(is.null(population)){
    return(NULL)
  }
  runs=nrow(population)
  nvars=ncol(population)
  consensus=matrix(0,nrow = nvars,ncol = nvars)
  for (cycle in 1:runs) {
    for (j in 1:nvars) {
      for (i in 1:j) {
        if(population[cycle,j]==population[cycle,i] && i!=j){
          consensus[i,j]=consensus[i,j]+1
          consensus[j,i]=consensus[j,i]+1
        }
      }
    }
  }
  return(consensus+runs*diag(nvars))
}
