#' A calista Function
#'
#' @param DATA,INPUTS,Results
#' @keywords calista
#' @export
#' @examples
#' stopQuietly_calista()



stopQuietly_calista <- function(...) {
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "));
  stop(simpleError(blankMsg));
}
