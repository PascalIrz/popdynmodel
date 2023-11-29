#' Title
#'
#' @param fun
#' @param parameter
#'
#' @return
#' @export
#'
#' @examples
get_modparameters <- function(fun = NULL, parameter = NULL) {
  load("modfunparameters.rda")
  if (!is.null(fun)) {
    modfunparameters <- filter(modfunparameters, Function %in% fun)
  }
  if (!is.null(parameter)) {
    modfunparameters <- filter(modfunparameters, Para %in% parameter)
  }
  modfunparameters <-
    select(modfunparameters,-Para,-Function) %>% distinct() %>% print(right = FALSE, row.names = FALSE)
}
