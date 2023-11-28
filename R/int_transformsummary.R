#' Title
#'
#' @param mcmc_summary
#' @param datamodel
#' @param var_id
#' @param var_tmp
#' @param var_tax
#' @param period
#'
#' @return
#' @export
#'
#' @examples
int_transformsummary <- function(mcmc_summary, datamodel, var_id, var_tmp, var_tax, period) {
  #-----------------------------------------------------------------------------
  colnames(mcmc_summary) <- c("mean", "sd", "Q2.5", "Q25", "Q50", "Q75", "Q97.5", "Rhat", "n.eff")
  id <- datamodel$popdyn_int$id
  time <- datamodel$popdyn_int$time
  tax <- datamodel$popdyn_int$taxa
  region <- datamodel$popdyn_int$region
  guild <- datamodel$popdyn_int$guild
  varenv <- datamodel$popdyn_int$varenv
  varcol <- datamodel$popdyn_int$varcol
  vargrow <- datamodel$popdyn_int$vargrow
  vardet <- datamodel$popdyn_int$vardet
  start <- datamodel$popdyn_int$start
  end <- datamodel$popdyn_int$end
  protocol <- datamodel$popdyn_int$protocol
  #-----------------------------------------------------------------------------
  # get parameter names and position from row names
  mcmc_summary$para <- row.names(mcmc_summary)
  mcmc_summary$start <- str_locate(mcmc_summary$para, "\\[")[,1]
  mcmc_summary$end <- str_locate(mcmc_summary$para, "\\]")[,1]
  mcmc_summary$parameter <- substr(mcmc_summary$para, 1, mcmc_summary$start - 1)
  mcmc_summary$pos <- substr(mcmc_summary$para, mcmc_summary$start + 1, mcmc_summary$end - 1)
  position <- str_split_fixed(mcmc_summary$pos,",",Inf)
  mcmc_summary <- mcmc_summary %>%
    cbind(setNames(data.frame(position),paste("p",1:ncol(position),sep="")) %>%
            rowwise() %>% mutate(p2 = ifelse("p2" %in% names(.), p2, 0)) %>%
            rowwise() %>% mutate(p3 = ifelse("p3" %in% names(.), p3, 0))) %>%
    mutate(across(c("p1","p2","p3"), ~as.numeric(.x))) %>%
    relocate(parameter, .before = mean)
  #-----------------------------------------------------------------------------
  # get variable names from position
  mcmc_summary$paratax <- ifelse(mcmc_summary$parameter %in% c("z","p.per","p.col","p.ext","p.per_id","p.col_id","p.ext_id","alpha","beta","alpha_per","beta_per","alpha_col","beta_col","z_lambda","z_mulambda","turnover","OR","N","N_lambda_id","N_mulambda_id","N_lambda","N_mulambda","N_PGR","B","B_lambda_id","B_mulambda_id","B_lambda","B_mulambda","B_PGR","alpha_N","beta_N","alpha_B","beta_B","p.det","alpha_det","beta_det"), tax[mcmc_summary$p1], NA)
  mcmc_summary$paraid <- ifelse(mcmc_summary$parameter %in% c("z","p.per_id","p.ext_id","p.col_id","N","N_lambda_id","N_mulambda_id","B_lambda_id","B_mulambda_id"), id[mcmc_summary$p2], NA)
  mcmc_summary$paratime <- ifelse(mcmc_summary$parameter %in% c("z","z_lambda","turnover","z_lambda_gui","N","N_lambda_id","N_lambda","B","B_lambda_id","B_lambda","N_lambda_gui","B_lambda_gui"), time[mcmc_summary$p3], NA)
  #-----------------------------------------------------------------------------
  # get region names and subscripts
  sub_reg <- NULL
  if (!is.null(region)) {
    mcmc_summary$region <- ifelse(mcmc_summary$parameter %in% c("z_mulambda","z_lambda","turnover","OR","z_lambda_gui","z_mulambda_gui","GOR","N_lambda_gui","B_lambda_gui","N_mulambda_gui","B_mulambda_gui"), region[mcmc_summary$p2], NA)
    sub_reg <- select(mcmc_summary, region, p2) %>% na.omit() %>% distinct() %>% arrange(p2) %>% set_colnames(c("","subscript")) %>% set_rownames(NULL)
  }
  #-----------------------------------------------------------------------------
  # get guilds names and subscripts
  sub_gui <- NULL
  if (!is.null(guild)) {
    mcmc_summary$guild <- ifelse(mcmc_summary$parameter %in% c("z_lambda_gui","z_mulambda_gui","GOR","N_lambda_gui","N_mulambda_gui","B_lambda_gui","B_mulambda_gui"), guild[mcmc_summary$p1], NA)
    sub_gui <- select(mcmc_summary, guild, p1) %>% na.omit() %>% distinct() %>% arrange(p1) %>% set_colnames(c("","subscript")) %>% set_rownames(NULL)
  }
  #-----------------------------------------------------------------------------
  # get period from position and period subscripts
  sub_period <- NULL
  if (!is.null(period)) {
    mcmc_summary$period <- ifelse(mcmc_summary$parameter %in% c("z_mulambda","OR","z_mulambda_gui","GOR","N_mulambda_id","N_mulambda","N_PGR","N_mulambda_gui","B_mulambda_id","B_mulambda","B_PGR","B_mulambda_gui"), paste(time[start[mcmc_summary$p3] - 1],time[end[mcmc_summary$p3]],sep="-"), NA)
    sub_period <- select(mcmc_summary, period, p3) %>% na.omit() %>% distinct() %>% arrange(p3) %>% set_rownames(.$p3) %>% set_colnames(c("","subscript")) %>% set_rownames(NULL)
  }
  #-----------------------------------------------------------------------------
  # get names and subscripts of environmental covariates
  sub_var <- NULL
  if (!is.null(varenv)) {
    mcmc_summary$covariate <- ifelse(mcmc_summary$parameter %in% c("beta","beta_per","beta_col"), varenv[mcmc_summary$p2], NA)
    sub_var <- select(mcmc_summary, covariate, p2) %>% na.omit() %>% distinct() %>% arrange(p2) %>% set_colnames(c("","subscript")) %>% set_rownames(NULL)
  }
  if (!is.null(varcol)) {
    mcmc_summary$covariate <- ifelse(mcmc_summary$parameter %in% c("beta_col"), varcol[mcmc_summary$p2], NA)
    sub_var <- select(mcmc_summary, covariate, p2) %>% na.omit() %>% distinct() %>% arrange(p2) %>% set_colnames(c("","subscript")) %>% set_rownames(NULL)
  }
  if (!is.null(vargrow)) {
    mcmc_summary$covariate <- ifelse(mcmc_summary$parameter %in% c("beta_N","beta_B"), vargrow[mcmc_summary$p2], NA)
    sub_var <- select(mcmc_summary, covariate, p2) %>% na.omit() %>% distinct() %>% arrange(p2) %>% set_colnames(c("","subscript")) %>% set_rownames(NULL)
  }
  #-----------------------------------------------------------------------------
  # get names and subscripts of detection covariates
  sub_det <- NULL
  if (!is.null(vardet)) {
    mcmc_summary$detection_covariate <- ifelse(mcmc_summary$parameter %in% c("beta_det"), vardet[mcmc_summary$p2], NA)
    sub_det <- select(mcmc_summary, detection_covariate, p2) %>% na.omit() %>% distinct() %>% arrange(p2) %>% set_colnames(c("","subscript")) %>% set_rownames(NULL)
  }
  #-----------------------------------------------------------------------------
  # get names and subscripts of sampling protocol
  sub_pro <- NULL
  if (!is.null(protocol)) {
    mcmc_summary$protocol <- ifelse(mcmc_summary$parameter %in% c("p.det","alpha_det"), protocol[mcmc_summary$p2], NA)
    mcmc_summary$protocol <- ifelse(mcmc_summary$parameter %in% c("beta_det"), protocol[mcmc_summary$p3], mcmc_summary$protocol)
    sub_pro <- select(mcmc_summary, protocol, p2) %>% na.omit() %>% distinct() %>% arrange(p2) %>% set_colnames(c("","subscript")) %>% set_rownames(NULL)
  }
  #-----------------------------------------------------------------------------
  # set list of subscripts
  subscript <- list()
  subname <- vector("character",0)
  sub_tax <- select(mcmc_summary, paratax, p1) %>% na.omit() %>% distinct() %>% arrange(p1) %>% set_colnames(c("","subscript")) %>% set_rownames(NULL)
  if (FALSE %in% quo_is_null(var_tax) & nrow(sub_tax) > 0) {
    subscript$taxa <- sub_tax
    subname <- c(subname, quo_name(var_tax))
    mcmc_summary <- mutate(mcmc_summary, !!quo_name(var_tax) := paratax)
  }
  sub_id <- select(mcmc_summary, paraid, p2) %>% na.omit() %>% distinct() %>% arrange(p2) %>% set_colnames(c("","subscript")) %>% set_rownames(NULL)
  if (nrow(sub_id) > 0) {
    subscript$id <- sub_id
    subname <- c(subname, quo_name(var_id))
    mcmc_summary <- mutate(mcmc_summary, !!quo_name(var_id) := paraid)
  }
  sub_time <- select(mcmc_summary, paratime, p3) %>% na.omit() %>% distinct() %>% arrange(p3) %>% set_colnames(c("","subscript")) %>% set_rownames(NULL)
  if (nrow(sub_time) > 0) {
    subscript$time <- sub_time
    subname <- c(subname, quo_name(var_tmp))
    mcmc_summary <- mutate(mcmc_summary, !!quo_name(var_tmp) := paratime)
  }
  names(subscript) <- subname
  subscript$region <- sub_reg
  subscript$guild <- sub_gui
  subscript$period <- sub_period
  subscript$covariate <- sub_var
  subscript$protocol <- sub_pro
  subscript$detection_covariate <- sub_det
  #-----------------------------------------------------------------------------
  mcmc_summary <- select(mcmc_summary, -c(para, start, end, pos, p1, p2, p3, paratax, paraid, paratime))
  return(list(mcmc_summary = mcmc_summary, subscript = subscript))
}
