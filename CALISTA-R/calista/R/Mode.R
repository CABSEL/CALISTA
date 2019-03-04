#' A calista Function
#'
#' @param x
#' @keywords calista
#' @export
#' @examples
#' Mode()


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
