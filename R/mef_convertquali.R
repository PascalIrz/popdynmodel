#' Title
#'
#' @param df
#' @param var_quali
#'
#' @return
#' @export
#'
#' @examples
mef_convertquali <- function(df, var_quali) {
  var_quali <- enquo(var_quali)
  name_quali <- colnames(select(df, !!var_quali))
  for (i in name_quali) {
    modality <- unique(pull(df, i)[!is.na(pull(df, i))])
    for (j in 2:length(modality)) {
      df <- df %>% mutate(across(paste(i), ~ifelse(is.na(.x), NA, ifelse(.x %in% modality[j], 1, 0)), .names = paste(i,modality[j],sep="_")))
    }
  }
  return(df)
}
