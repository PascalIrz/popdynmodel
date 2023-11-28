#' Title
#'
#' @param df
#' @param var_id
#' @param var_tmp
#' @param vars
#' @param call
#'
#' @return
#' @export
#'
#' @examples
int_checkvarenv <- function(df, var_id, var_tmp, vars, call = caller_env()) {
  if (length(vars) > 0) {
    if (FALSE %in% is.element(c(unlist(vars[!vars %in% NULL]),quo_name(var_id),quo_name(var_tmp)),colnames(df))) {
      abort("must use existing variables",call = call)
    }
    vvars <- colnames(df)[which(colnames(df) %in% unlist(vars[!vars %in% NULL]))]
    tf <- df %>%
      summarise(across(all_of(vvars), ~ifelse(length(unique(na.omit(.x))) > 1, 1, 0)), .by = !!var_id) %>%
      summarise(across(all_of(vvars), ~sum(.x)))
    var_name <- names(tf)[which(tf == 0)]
    if (length(var_name) > 0) {
      df_var <- select(df, all_of(var_name), !!var_id) %>% na.omit() %>% distinct()
      df <- select(df, -var_name) %>% left_join(df_var)
    }
    tf <- df %>%
      summarise(across(all_of(vvars), ~ifelse(length(unique(na.omit(.x))) > 1, 1, 0)), .by = !!var_tmp) %>%
      summarise(across(all_of(vvars), ~sum(.x)))
    var_name <- names(tf)[which(tf == 0)]
    if (length(var_name) > 0) {
      df_var <- select(df, all_of(var_name), !!var_tmp) %>% na.omit() %>% distinct()
      df <- select(df, -var_name) %>% left_join(df_var)
    }
  }
  return(df)
}
