#' Title
#'
#' @param df
#' @param var_id
#' @param var_tmp
#' @param var_imp
#'
#' @return
#' @export
#'
#' @examples
mef_imputevalue <- function(df, var_id, var_tmp, var_imp) {
  check_required(var_id)
  check_required(var_tmp)
  check_required(var_imp)
  #-----------------------------------------------------------------------------
  var_id <- enquo(var_id)
  var_tmp <- enquo(var_tmp)
  var_imp <- enquo(var_imp)
  #-----------------------------------------------------------------------------
  df <- df %>% unite("id", quo_name(var_id), remove = FALSE) %>% unite("tmp", quo_name(var_tmp), remove = FALSE)
  vec_tmp <- sort(unique(df$tmp))
  vec_imp <- select(df, id, tmp, !!var_imp) %>%
    filter(!tmp == min(vec_tmp), is.na(df[,quo_name(var_imp)])) %>% distinct() %>% arrange(tmp)
  for (i in 1:nrow(vec_imp)) {
    df[(df$id %in% vec_imp[i,"id"] & df$tmp %in% vec_imp[i,"tmp"]),quo_name(var_imp)] <- unique(df[(df$id %in% vec_imp[i,"id"] & df$tmp %in% vec_tmp[which(vec_tmp == vec_imp[i,"tmp"]) - 1]),quo_name(var_imp)])
  }
  #-----------------------------------------------------------------------------
  vec_imp <- select(df, id, tmp, !!var_imp) %>%
    filter(!tmp == max(vec_tmp), is.na(df[,quo_name(var_imp)])) %>% distinct() %>% arrange(tmp)
  for (i in nrow(vec_imp):1) {
    df[(df$id %in% vec_imp[i,"id"] & df$tmp %in% vec_imp[i,"tmp"]),quo_name(var_imp)] <- unique(df[(df$id %in% vec_imp[i,"id"] & df$tmp %in% vec_tmp[which(vec_tmp == vec_imp[i,"tmp"]) + 1]),quo_name(var_imp)])
  }
  #-----------------------------------------------------------------------------
  df$id <- df$tmp <- NULL
  return(df)
}
