#' A calista Function
#'
#' @param v
#' @keywords calista
#' @export
#' @examples
#'
#' get_diff()




get_diff<-function(v){
  len=length(v)
  my_diff=matrix(0,1,len-1)
  for (i in 1:len-1) {
    my_diff[i]=v[i+1]-v[i]
  }
  my_diff
}
