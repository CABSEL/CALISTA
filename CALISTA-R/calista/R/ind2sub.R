#' A calista Function
#'
#' @param sz,ind
#' @keywords calista
#' @export
#' @examples
#' ind2sub()


ind2sub = function(sz,ind)
{
  ind = as.matrix(ind,ncol=1);
  sz = c(1,sz);
  den = 1;
  sub = c();
  for(i in 2:length(sz)){
    den = den * sz[i-1];
    num = den * sz[i];
    s = floor(((ind-1) %% num)/den) + 1;
    sub = cbind(sub,s);
  }
  return(sub);
}
