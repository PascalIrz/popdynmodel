#' Title
#'
#' @param df
#' @param var_id
#' @param var_tmp
#' @param var_env
#' @param var_tax
#' @param var_cnt
#' @param var_wei
#' @param var_surf
#' @param var_reg
#' @param var_guild
#' @param period
#' @param timestep
#' @param save_parameters
#' @param n_chain
#' @param n_iter
#' @param n_thin
#' @param n_burnin
#'
#' @return
#' @export
#'
#' @examples
modenv_popgrow <- function(df, var_id, var_tmp, var_env, var_tax=NULL, var_cnt=NULL, var_wei=NULL, var_surf=NULL, var_reg=NULL, var_guild=NULL, period=NULL, timestep=1, save_parameters = NULL, n_chain = 3, n_iter = 10000, n_thin = ceiling(n_iter/100), n_burnin = floor(n_iter/4)) {
  #-----------------------------------------------------------------------------
  # Check for missing required arguments
  check_required(var_id)
  check_required(var_tmp)
  check_required(var_env)
  if (quo_is_null(enquo(var_cnt)) & quo_is_null(enquo(var_wei))) {
    abort("'var_cnt' or 'var_wei' must be supplied")
  }
  #-----------------------------------------------------------------------------
  # Check for mistakes, if failure return an error message and stop
  df <- do.call(int_checkvarenv, list(df, var_id=enquo(var_id), var_tmp=enquo(var_tmp), vars=syms(var_env)))
  df <- do.call(int_checkfunction, list(df,
                                        vars_in_df=syms(c(var_id, var_tmp, var_tax, var_cnt, var_wei, var_env, var_surf, var_reg, var_guild)),
                                        vars_na=syms(c(var_id, var_tmp, var_tax, var_env, var_surf)),
                                        vars_numeric=syms(c(var_tmp, var_cnt, var_wei, var_env, var_surf)),
                                        vars_duplicate=syms(c(var_id, var_tmp, var_tax)),
                                        var_tmp=enquo(var_tmp), timestep, period,
                                        vars_pas=NULL))
  #-----------------------------------------------------------------------------
  var_id <- enquo(var_id)
  var_tmp <- enquo(var_tmp)
  var_tax <- enquo(var_tax)
  var_cnt <- enquo(var_cnt)
  var_wei <- enquo(var_wei)
  var_grow <- enquo(var_env)
  var_surf <- enquo(var_surf)
  var_reg <- enquo(var_reg)
  var_guild <- enquo(var_guild)
  if (FALSE %in% quo_is_null(var_cnt)) { var_pres <- var_cnt } else { var_pres <- var_wei }
  #-----------------------------------------------------------------------------
  # Write model and model data
  datamodel <- do.call(int_datamodel, list(df, occup=FALSE, grow=TRUE, modenv=NULL, modenvG=1, timestep, period, var_id, var_tmp, var_tax, var_pres, var_reg, var_guild, var_cnt, var_wei, var_surf, var_envO=NULL, var_envP=NULL, var_envC=NULL, var_grow))
  code <- do.call(int_popoccup, list(occup=FALSE,modenv=NULL,var_envO=NULL,var_envP=NULL,var_envC=NULL,var_guild=NULL))
  popdyn_code <- do.call(int_popgrow, list(code,modenvG=1,var_cnt,var_wei,var_guild))
  #-----------------------------------------------------------------------------
  # Define requested parameters
  if (is.null(save_parameters)) { save_parameters <- datamodel$popdyn_parameters } else {
    if (FALSE %in% is.element(save_parameters, datamodel$popdyn_parameters)) {
      para_name <- save_parameters[which(is.element(save_parameters, datamodel$popdyn_parameters) %in% FALSE)]
      abort(paste0("Some parameters are not in model: '", para_name, "'", collapse = " "))
    }
  }
  #-----------------------------------------------------------------------------
  # Fit model
  set.seed(123)
  popdyn <- nimbleModel(code = popdyn_code,
                        constants = datamodel$popdyn_const,
                        data = datamodel$popdyn_data,
                        inits = datamodel$popdyn_inits,
                        name = "popdyn", calculate = FALSE)
  popdynConf <- configureMCMC(popdyn, monitors = save_parameters)
  popdynMCMC <- buildMCMC(popdynConf)
  popdynComp <- compileNimble(popdyn)
  popdynModel <- compileNimble(popdynMCMC, project = popdyn, resetFunctions = TRUE)
  mcmc_chain <- runMCMC(popdynModel, nchains = n_chain, niter = n_iter, thin = n_thin, nburnin = n_burnin, setSeed = 123, samplesAsCodaMCMC = TRUE)
  if (n_chain == 1 & !is.null(dim(mcmc_chain))) {
    mcmc_na <- which(is.na(mcmc_chain[1,]))
    if (length(mcmc_na) > 0) { mcmc_chain <- mcmc_chain[,-mcmc_na] }
  }
  if (n_chain > 1 & !is.null(dim(mcmc_chain[[1]]))) { for (i in 1:n_chain) {
    mcmc_chain[[i]] <- mcmc_chain[[i]][,!colnames(mcmc_chain[[i]]) %in% names(which(is.na(mcmc_chain[[i]][1,])))]
  }}
  #-----------------------------------------------------------------------------
  # set summary data frame and list of subscripts
  try(mcmc_summary <- MCMCsummary(mcmc_chain, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), Rhat = TRUE, n.eff = TRUE))
  if (!exists("mcmc_summary")) {
    return(mcmc_chain)
    stop(expression(NULL))
  }
  list_summary <- do.call(int_transformsummary, list(mcmc_summary, datamodel, var_id, var_tmp, var_tax, period))
  output <- list(mcmc_summary = list_summary$mcmc_summary, mcmc_chain = mcmc_chain, subscript = list_summary$subscript)
  #-----------------------------------------------------------------------------
  return(output)
}
