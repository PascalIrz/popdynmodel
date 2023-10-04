get_ModParameters <- function(fun = NULL, parameter = NULL) {
  load("modfunparameters.rda")
  if (!is.null(fun)) {
    modfunparameters <- filter(modfunparameters, Function %in% fun)
  }
  if (!is.null(parameter)) {
    modfunparameters <- filter(modfunparameters, Para %in% parameter)
  }
  modfunparameters <- select(modfunparameters, -Para, - Function) %>% distinct() %>% print(right = FALSE, row.names = FALSE)
}
################################################################################
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
################################################################################
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
################################################################################
mod_popoccup <- function(df, var_id, var_tmp, var_tax, var_cnt, var_reg=NULL, var_guild=NULL, period=NULL, timestep=1, save_parameters=NULL, n_chain=3, n_iter=10000, n_thin=ceiling(n_iter/100), n_burnin=floor(n_iter/4)) {
  #-----------------------------------------------------------------------------
  # Check for missing required arguments
  check_required(var_id)
  check_required(var_tmp)
  check_required(var_tax)
  check_required(var_cnt)
  #-----------------------------------------------------------------------------
  # Check for mistakes, if failure return an error message and stop
  vars_in_df <- syms(c(var_id, var_tmp, var_tax, var_cnt, var_reg, var_guild))
  vars_na <- syms(c(var_id, var_tmp, var_tax))
  vars_numeric <- syms(c(var_tmp, var_cnt))
  vars_duplicate <- syms(c(var_id, var_tmp, var_tax))
  df <- do.call(int_checkfunction, list(df,vars_in_df,var_id=NULL,env_tmp=TRUE,vars_env=NULL,vars_na,vars_numeric,vars_duplicate,enquo(var_tmp),timestep,period,vars_pas=NULL))
  #-----------------------------------------------------------------------------
  var_id <- enquo(var_id)
  var_tmp <- enquo(var_tmp)
  var_tax <- enquo(var_tax)
  var_pres <- enquo(var_cnt)
  var_reg <- enquo(var_reg)
  var_guild <- enquo(var_guild)
  #-----------------------------------------------------------------------------
  # Write model and model data
  datamodel <- do.call(int_datamodel, list(df, occup=TRUE, grow=FALSE, modenv=NULL, modenvG=NULL, timestep, period, var_id, var_tmp, var_tax, var_pres, var_reg, var_guild, var_cnt=NULL, var_wei=NULL, var_surf=NULL, var_envO=NULL, var_envP=NULL, var_envC=NULL, var_grow=NULL))
  code <- do.call(int_popoccup, list(modenv=NULL,var_envO=NULL,var_envP=NULL,var_envC=NULL,var_guild))
  popdyn_code <- as.call(c(as.symbol("{"), code))
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
  # set summary data frame of taxa and guilds and list of subscripts
  mcmc_summary <- MCMCsummary(mcmc_chain, Rhat = TRUE, n.eff = TRUE)
  list_summary <- do.call(int_transformsummary, list(mcmc_summary, datamodel, var_id, var_tmp, var_tax, period)) 
  output <- list(mcmc_summary = list_summary$mcmc_summary, mcmc_chain = mcmc_chain, subscript = list_summary$subscript)
  #-----------------------------------------------------------------------------
  return(output)
}
################################################################################
modenv_popoccup <- function(df, var_id, var_tmp, var_tax, var_cnt, var_envO=NULL, var_envP=NULL, var_envC=NULL, env_tmp=TRUE, var_reg=NULL, var_guild=NULL, period=NULL, timestep=1, save_parameters=NULL, n_chain=3, n_iter=10000, n_thin=ceiling(n_iter/100), n_burnin=floor(n_iter/4)) {
  #-----------------------------------------------------------------------------
  # Check for missing required arguments
  check_required(var_id)
  check_required(var_tmp)
  check_required(var_tax)
  check_required(var_cnt)
  if (TRUE %in% quo_is_null(enquo(var_envO)) & TRUE %in% quo_is_null(enquo(var_envP)) & TRUE %in% quo_is_null(enquo(var_envC))) {
    abort("'var_envO', 'var_envP' or 'var_envC' must be supplied")
  }
  #-----------------------------------------------------------------------------
  # Check for mistakes, if failure return an error message and stop
  vars_in_df <- syms(c(var_id, var_tmp, var_tax, var_cnt, var_envO, var_envP, var_envC, var_reg, var_guild))
  vars_env <- syms(c(var_envO, var_envP, var_envC))
  vars_na <- syms(c(var_id, var_tmp, var_tax, var_envO, var_envP, var_envC))
  vars_numeric <- syms(c(var_tmp, var_cnt, var_envO, var_envP, var_envC))
  vars_duplicate <- syms(c(var_id, var_tmp, var_tax))
  df <- do.call(int_checkfunction, list(df,vars_in_df,enquo(var_id),env_tmp,vars_env,vars_na,vars_numeric,vars_duplicate,enquo(var_tmp),timestep,period,vars_pas=NULL))
  #-----------------------------------------------------------------------------
  var_id <- enquo(var_id)
  var_tmp <- enquo(var_tmp)
  var_tax <- enquo(var_tax)
  var_pres <- enquo(var_cnt)
  var_envO <- enquo(var_envO)
  var_envP <- enquo(var_envP)
  var_envC <- enquo(var_envC)
  var_reg <- enquo(var_reg)
  var_guild <- enquo(var_guild)
  #-----------------------------------------------------------------------------
  # Write model and model data
  datamodel <- do.call(int_datamodel, list(df, occup=TRUE, grow=FALSE, modenv=1, modenvG=NULL, timestep, period, var_id, var_tmp, var_tax, var_pres, var_reg, var_guild, var_cnt=NULL, var_wei=NULL, var_surf=NULL, var_envO, var_envP, var_envC, var_grow=NULL))
  code <- do.call(int_popoccup, list(modenv=1, var_envO, var_envP, var_envC, var_guild))
  popdyn_code <- as.call(c(as.symbol("{"), code))
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
  # set summary data frame of taxa and guilds and list of subscripts
  mcmc_summary <- MCMCsummary(mcmc_chain, Rhat = TRUE, n.eff = TRUE)
  list_summary <- do.call(int_transformsummary, list(mcmc_summary, datamodel, var_id, var_tmp, var_tax, period)) 
  output <- list(mcmc_summary = list_summary$mcmc_summary, mcmc_chain = mcmc_chain, subscript = list_summary$subscript)
  #-----------------------------------------------------------------------------
  return(output)
}
################################################################################
mod_popgrow <- function(df, var_id, var_tmp, var_tax, var_cnt=NULL, var_wei=NULL, var_surf=NULL, var_reg=NULL, var_guild=NULL, period=NULL, timestep=1, save_parameters = NULL, n_chain = 3, n_iter = 10000, n_thin = ceiling(n_iter/100), n_burnin = floor(n_iter/4)) {
  #-----------------------------------------------------------------------------
  # Check for missing required arguments
  check_required(var_id)
  check_required(var_tmp)
  check_required(var_tax)
  if (TRUE %in% quo_is_null(enquo(var_cnt)) & TRUE %in% quo_is_null(enquo(var_wei))) {
    abort("'var_cnt' or 'var_wei' must be supplied")
  }
  #-----------------------------------------------------------------------------
  # Check for mistakes, if failure return an error message and stop
  vars_in_df <- syms(c(var_id, var_tmp, var_tax, var_cnt, var_wei, var_surf, var_reg, var_guild))
  vars_na <- syms(c(var_id, var_tmp, var_tax, var_surf))
  vars_numeric <- syms(c(var_tmp, var_cnt, var_wei, var_surf))
  vars_duplicate <- syms(c(var_id, var_tmp, var_tax))
  df <- do.call(int_checkfunction, list(df,vars_in_df,var_id=NULL,env_tmp=TRUE,vars_env=NULL,vars_na,vars_numeric,vars_duplicate,enquo(var_tmp),timestep,period,vars_pas=NULL))
  #-----------------------------------------------------------------------------
  var_id <- enquo(var_id)
  var_tmp <- enquo(var_tmp)
  var_tax <- enquo(var_tax)
  var_cnt <- enquo(var_cnt)
  var_wei <- enquo(var_wei)
  var_surf <- enquo(var_surf)
  var_reg <- enquo(var_reg)
  var_guild <- enquo(var_guild)
  if (FALSE %in% quo_is_null(var_cnt)) { var_pres <- var_cnt } else { var_pres <- var_wei }
  #-----------------------------------------------------------------------------
  # Write model and model data
  datamodel <- do.call(int_datamodel, list(df, occup=FALSE, grow=TRUE, modenv=NULL, modenvG=NULL, timestep, period, var_id, var_tmp, var_tax, var_pres, var_reg, var_guild, var_cnt, var_wei, var_surf, var_envO=NULL, var_envP=NULL, var_envC=NULL, var_grow=NULL))
  popdyn_code <- do.call(int_popgrow, list(occup=FALSE,code=NULL,modenvG=NULL,var_cnt,var_wei,var_guild))
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
  mcmc_summary <- MCMCsummary(mcmc_chain, Rhat = TRUE, n.eff = TRUE)
  list_summary <- do.call(int_transformsummary, list(mcmc_summary, datamodel, var_id, var_tmp, var_tax, period)) 
  output <- list(mcmc_summary = list_summary$mcmc_summary, mcmc_chain = mcmc_chain, subscript = list_summary$subscript)
  #-----------------------------------------------------------------------------
  return(output)
} 
################################################################################
modenv_popgrow <- function(df, var_id, var_tmp, var_tax, var_env, var_cnt=NULL, var_wei=NULL, var_surf=NULL, var_reg=NULL, var_guild=NULL, period=NULL, timestep=1, save_parameters = NULL, n_chain = 3, n_iter = 10000, n_thin = ceiling(n_iter/100), n_burnin = floor(n_iter/4)) {
  #-----------------------------------------------------------------------------
  # Check for missing required arguments
  check_required(var_id)
  check_required(var_tmp)
  check_required(var_tax)
  check_required(var_env)
  if (TRUE %in% quo_is_null(enquo(var_cnt)) & TRUE %in% quo_is_null(enquo(var_wei))) {
    abort("'var_cnt' or 'var_wei' must be supplied")
  }
  #-----------------------------------------------------------------------------
  # Check for mistakes, if failure return an error message and stop
  vars_in_df <- syms(c(var_id, var_tmp, var_tax, var_cnt, var_wei, var_env, var_surf, var_reg, var_guild))
  vars_na <- syms(c(var_id, var_tmp, var_tax, var_env, var_surf))
  vars_numeric <- syms(c(var_tmp, var_cnt, var_wei, var_env, var_surf))
  vars_duplicate <- syms(c(var_id, var_tmp, var_tax))
  df <- do.call(int_checkfunction, list(df,vars_in_df,var_id=NULL,env_tmp=TRUE,vars_env=NULL,vars_na,vars_numeric,vars_duplicate,enquo(var_tmp),timestep,period,vars_pas=NULL))
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
  popdyn_code <- do.call(int_popgrow, list(occup=FALSE,code=NULL,modenvG=1,var_cnt,var_wei,var_guild))
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
  try(mcmc_summary <- MCMCsummary(mcmc_chain, Rhat = TRUE, n.eff = TRUE))
  if (!exists("mcmc_summary")) {
    return(mcmc_chain)  
    stop(expression(NULL))
  }
  list_summary <- do.call(int_transformsummary, list(mcmc_summary, datamodel, var_id, var_tmp, var_tax, period)) 
  output <- list(mcmc_summary = list_summary$mcmc_summary, mcmc_chain = mcmc_chain, subscript = list_summary$subscript)
  #-----------------------------------------------------------------------------
  return(output)
} 
################################################################################
mod_popdyn <- function(df, var_id, var_tmp, var_tax, var_cnt=NULL, var_wei=NULL, var_surf=NULL, var_reg=NULL, var_guild=NULL, period=NULL, timestep=1, save_parameters = NULL, n_chain = 3, n_iter = 10000, n_thin = ceiling(n_iter/100), n_burnin = floor(n_iter/4)) {
  #-----------------------------------------------------------------------------
  # Check for missing required arguments
  check_required(var_id)
  check_required(var_tmp)
  check_required(var_tax)
  if (TRUE %in% quo_is_null(enquo(var_cnt)) & TRUE %in% quo_is_null(enquo(var_wei))) {
    abort("'var_cnt' or 'var_wei' must be supplied")
  }
  #-----------------------------------------------------------------------------
  # Check for mistakes, if failure return an error message and stop
  vars_in_df <- syms(c(var_id, var_tmp, var_tax, var_cnt, var_wei, var_surf, var_reg, var_guild))
  vars_na <- syms(c(var_id, var_tmp, var_tax, var_surf))
  vars_numeric <- syms(c(var_tmp, var_cnt, var_wei, var_surf))
  vars_duplicate <- syms(c(var_id, var_tmp, var_tax))
  df <- do.call(int_checkfunction, list(df,vars_in_df,var_id=NULL,env_tmp=TRUE,vars_env=NULL,vars_na,vars_numeric,vars_duplicate,enquo(var_tmp),timestep,period,vars_pas=NULL))
  #-----------------------------------------------------------------------------
  var_id <- enquo(var_id)
  var_tmp <- enquo(var_tmp)
  var_tax <- enquo(var_tax)
  var_cnt <- enquo(var_cnt)
  var_wei <- enquo(var_wei)
  var_surf <- enquo(var_surf)
  var_reg <- enquo(var_reg)
  var_guild <- enquo(var_guild)
  if (FALSE %in% quo_is_null(var_cnt)) { var_pres <- var_cnt } else { var_pres <- var_wei }
  #-----------------------------------------------------------------------------
  # Write model and model data
  datamodel <- do.call(int_datamodel, list(df, occup=TRUE, grow=TRUE, modenv=NULL, modenvG=NULL, timestep, period, var_id, var_tmp, var_tax, var_pres, var_reg, var_guild, var_cnt, var_wei, var_surf, var_envO=NULL, var_envP=NULL, var_envC=NULL, var_grow=NULL))
  code <- do.call(int_popoccup, list(modenv=NULL,var_envO=NULL,var_envP=NULL,var_envC=NULL,var_guild))
  popdyn_code <- do.call(int_popgrow, list(occup=TRUE,code,modenvG=NULL,var_cnt,var_wei,var_guild))
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
  mcmc_summary <- MCMCsummary(mcmc_chain, Rhat = TRUE, n.eff = TRUE)
  list_summary <- do.call(int_transformsummary, list(mcmc_summary, datamodel, var_id, var_tmp, var_tax, period)) 
  output <- list(mcmc_summary = list_summary$mcmc_summary, mcmc_chain = mcmc_chain, subscript = list_summary$subscript)
  #-----------------------------------------------------------------------------
  return(output)
} 
################################################################################
modenv_popdyn <- function(df, var_id, var_tmp, var_tax, var_cnt=NULL, var_wei=NULL, var_env=NULL, var_envO=NULL, var_envP=NULL, var_envC=NULL, env_tmp=TRUE, var_surf=NULL, var_reg=NULL, var_guild=NULL, period=NULL, timestep=1, save_parameters = NULL, n_chain = 3, n_iter = 10000, n_thin = ceiling(n_iter/100), n_burnin = floor(n_iter/4)) {
  #-----------------------------------------------------------------------------
  # Check for missing required arguments
  check_required(var_id)
  check_required(var_tmp)
  check_required(var_tax)
  if (TRUE %in% quo_is_null(enquo(var_cnt)) & TRUE %in% quo_is_null(enquo(var_wei))) {
    abort("'var_cnt' or 'var_wei' must be supplied")
  }
  #-----------------------------------------------------------------------------
  # Check for mistakes, if failure return an error message and stop
  vars_in_df <- syms(c(var_id, var_tmp, var_tax, var_cnt, var_wei, var_env, var_envO, var_envP, var_envC, var_surf, var_reg, var_guild))
  vars_env <- syms(c(var_envO, var_envP, var_envC))
  vars_na <- syms(c(var_id, var_tmp, var_tax, var_env, var_envO, var_envP, var_envC, var_surf))
  vars_numeric <- syms(c(var_tmp, var_cnt, var_wei, var_env, var_envO, var_envP, var_envC, var_surf))
  vars_duplicate <- syms(c(var_id, var_tmp, var_tax))
  df <- do.call(int_checkfunction, list(df,vars_in_df,enquo(var_id),env_tmp,vars_env,vars_na,vars_numeric,vars_duplicate,enquo(var_tmp),timestep,period,vars_pas=NULL))
  #-----------------------------------------------------------------------------
  var_id <- enquo(var_id)
  var_tmp <- enquo(var_tmp)
  var_tax <- enquo(var_tax)
  var_cnt <- enquo(var_cnt)
  var_wei <- enquo(var_wei)
  var_grow <- enquo(var_env)
  var_envO <- enquo(var_envO)
  var_envP <- enquo(var_envP)
  var_envC <- enquo(var_envC)
  var_surf <- enquo(var_surf)
  var_reg <- enquo(var_reg)
  var_guild <- enquo(var_guild)
  if (FALSE %in% quo_is_null(var_cnt)) { var_pres <- var_cnt } else { var_pres <- var_wei }
  #-----------------------------------------------------------------------------
  # Write model and model data
  modenv <- modenvG <- NULL
  if (FALSE %in% quo_is_null(var_grow)) { modenvG <- 1 }
  if (FALSE %in% c(quo_is_null(var_envO), quo_is_null(var_envP), quo_is_null(var_envC))) { modenv <- 1 }
  datamodel <- do.call(int_datamodel, list(df, occup=TRUE, grow=TRUE, modenv, modenvG, timestep, period, var_id, var_tmp, var_tax, var_pres, var_reg, var_guild, var_cnt, var_wei, var_surf, var_envO, var_envP, var_envC, var_grow))
  code <- do.call(int_popoccup, list(modenv, var_envO, var_envP, var_envC, var_guild))
  popdyn_code <- do.call(int_popgrow, list(occup=TRUE,code,modenvG,var_cnt,var_wei,var_guild))
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
  try(mcmc_summary <- MCMCsummary(mcmc_chain, Rhat = TRUE, n.eff = TRUE))
  if (!exists("mcmc_summary")) {
    return(mcmc_chain)  
    stop(expression(NULL))
  }
  list_summary <- do.call(int_transformsummary, list(mcmc_summary, datamodel, var_id, var_tmp, var_tax, period)) 
  output <- list(mcmc_summary = list_summary$mcmc_summary, mcmc_chain = mcmc_chain, subscript = list_summary$subscript)
  #-----------------------------------------------------------------------------
  return(output)
}
################################################################################
################################################################################
################################################################################