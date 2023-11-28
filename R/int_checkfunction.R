#' Title
#'
#' @param df
#' @param vars_in_df
#' @param vars_na
#' @param vars_numeric
#' @param vars_duplicate
#' @param var_tmp
#' @param timestep
#' @param period
#' @param vars_pas
#' @param call
#'
#' @return
#' @export
#'
#' @examples
int_checkfunction <- function(df, vars_in_df, vars_na, vars_numeric, vars_duplicate, var_tmp, timestep, period, vars_pas=NULL, call = caller_env()) {
  if (FALSE %in% is.element(unlist(vars_in_df[!vars_in_df %in% "NULL"]),colnames(df))) {
    abort("must use existing variables", call = call)
  }
  df <- df[,which(colnames(df) %in% unlist(vars_in_df))] %>% distinct()
  #-----------------------------------------------------------------------------
  vars <- df[,which(colnames(df) %in% unlist(vars_na[!vars_na %in% NULL]))]
  if (TRUE %in% is.na(vars)) {
    vnames <- paste0(names(which(colSums(is.na(vars)) > 0)), collapse=", ")
    abort(paste0("NA or NaN are not allowed in '", vnames, "'", collapse=" "), call = call)
  }
  #-----------------------------------------------------------------------------
  vars <- lapply(df[,which(colnames(df) %in% unlist(vars_numeric[!vars_numeric %in% NULL]))], is.numeric)
  if (FALSE %in% vars) {
    vnames <- paste0(names(vars)[vars %in% TRUE], collapse=", ")
    abort(paste0("'", vnames, "' must be numeric", collapse = " "), call = call)
  }
  #-----------------------------------------------------------------------------
  if (TRUE %in% duplicated(df[,which(colnames(df) %in% unlist(vars_duplicate))])) {
    abort("must not have duplicate entries with multiple matches", call = call)
  }
  #-----------------------------------------------------------------------------
  if (max(diff(sort(unique(pull(df, !!var_tmp))))) > timestep) {
    abort("missing values are not allowed in 'var_tmp'", call = call)
  }
  #-----------------------------------------------------------------------------
  if (!is.null(period)) {
    if (FALSE %in% is.list(period)) {
      abort("'period' must be a list")
    }
    if (FALSE %in% (lapply(period, length) == 2) | FALSE %in% is.numeric(unlist(period))) {
      abort("'period' must be a numeric vector list of length 2", call = call)
    }
    if (FALSE %in% is.element(unlist(period), pull(df, !!var_tmp))) {
      abort("period values out of range", call = call)
    }
  }
  #-----------------------------------------------------------------------------
  if (!is.null(vars_pas)) {
    dpas <- reframe(df, across(unlist(vars_pas[1]), ~length(.x) == max(.x), .names = "pas"), .by = unlist(vars_pas[2:4]))
    if (FALSE %in% dpas$pas) {
      abord("missing passes in 'var_pas'", call = call)
    }
  }
  #-----------------------------------------------------------------------------
  return(df)
}
