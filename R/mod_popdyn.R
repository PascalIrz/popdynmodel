#' Fitting an abundance model
#'
#' In order to estimate population growth rates from the removal sampling-based model,
#'     use the popbay model specifying error=TRUE.
#'     Three MCMC chains (n_chain, default value) are considered to check for the convergence
#'     of the algorithm to a stationary posterior distribution.
#'
#' @param df Dataframe with the data and containing the variables passed as arguments.
#' @param var_id,var_tmp,var_tax,var_cnt,var_wei,var_surf,var_reg,var_guild Names of the variables
#'     of df respectively containing the site ID (compulsory), date/time of observations (compulsory),
#'     taxon (compulsory), count, weight, sampled area, region and functional guild of the taxon.
#'     Either `var_cnt` or `var_wei` must be given.
#' @param period
#' @param timestep
#' @param save_parameters Boolean. Should the parameters be saved?
#' @param n_chain,n_iter,n_thin,n_burnin Integers. Number of, respectively, MCMC chains (default = 3),
#'     iterations (default = 10000),
#'     thinning (eg n_thin = 50 indicates that the iterations are thinned to one every 50 iterations; default is n_iter/100),
#'     burn-in period (eg n_burnin = 1000 indicates that the first 1,000 iterations are discarded; default is n_iter/4).
#'
#' @importFrom rlang check_required quo_is_null enquo syms
#' @importFrom nimble nimbleModel configureMCMC buildMCMC compileNimble runMCMC
#' @importFrom MCMCvis MCMCsummary
#'
#' @return A list with the model output.
#' @export
#'
#' @examples
#' \dontrun{
#' # As fitting a Bayesian model can be time consuming, it is recommended to use few
#' #     iterations for testing models. In the example, the model is run with 6,000 iterations
#' #    (n_iter), thinned to one draw every 100th iteration (n_thin), after discarding a burn-in
#' #    period of 1,000 iterations (n_burnin).
#'
#'   my_model <- mod_popdyn(df = fauna,
#'                          var_id = site_id,
#'                          var_tmp = year,
#'                          var_cnt = abundance,
#'                          n_iter = 6000,
#'                          n_thin = 100,
#'                          n_burnin = 1000)
#' }
mod_popdyn <-
  function(df,
           var_id,
           var_tmp,
           var_tax = NULL,
           var_cnt = NULL,
           var_wei = NULL,
           var_surf = NULL,
           var_reg = NULL,
           var_guild = NULL,
           period = NULL,
           timestep = 1,
           save_parameters = NULL,
           n_chain = 3,
           n_iter = 10000,
           n_thin = ceiling(n_iter / 100),
           n_burnin = floor(n_iter / 4)) {
    #-----------------------------------------------------------------------------
    # Check for missing required arguments
    check_required(var_id)
    check_required(var_tmp)
    if (quo_is_null(enquo(var_cnt)) & quo_is_null(enquo(var_wei))) {
      abort("'var_cnt' or 'var_wei' must be supplied")
    }
    #-----------------------------------------------------------------------------
    # Check for mistakes, if failure return an error message and stop
    df <- do.call(
      int_checkfunction,
      list(
        df,
        vars_in_df = syms(
          c(
            var_id,
            var_tmp,
            var_tax,
            var_cnt,
            var_wei,
            var_surf,
            var_reg,
            var_guild
          )
        ),
        vars_na = syms(c(var_id, var_tmp, var_tax, var_surf)),
        vars_numeric = syms(c(var_tmp, var_cnt, var_wei, var_surf)),
        vars_duplicate = syms(c(var_id, var_tmp, var_tax)),
        var_tmp = enquo(var_tmp),
        timestep,
        period,
        vars_pas = NULL
      )
    )
    #-----------------------------------------------------------------------------
    var_id <- enquo(var_id)
    var_tmp <- enquo(var_tmp)
    var_tax <- enquo(var_tax)
    var_cnt <- enquo(var_cnt)
    var_wei <- enquo(var_wei)
    var_surf <- enquo(var_surf)
    var_reg <- enquo(var_reg)
    var_guild <- enquo(var_guild)
    if (quo_is_null(var_cnt)) {
      var_pres <- var_wei
    } else {
      var_pres <- var_cnt
    }
    #-----------------------------------------------------------------------------
    # Write model and model data
    datamodel <-
      do.call(
        int_datamodel,
        list(
          df,
          occup = TRUE,
          grow = TRUE,
          modenv = NULL,
          modenvG = NULL,
          timestep,
          period,
          var_id,
          var_tmp,
          var_tax,
          var_pres,
          var_reg,
          var_guild,
          var_cnt,
          var_wei,
          var_surf,
          var_envO = NULL,
          var_envP = NULL,
          var_envC = NULL,
          var_grow = NULL
        )
      )
    code <-
      do.call(
        int_popoccup,
        list(
          occup = TRUE,
          modenv = NULL,
          var_envO = NULL,
          var_envP = NULL,
          var_envC = NULL,
          var_guild
        )
      )
    popdyn_code <-
      do.call(int_popgrow,
              list(code, modenvG = NULL, var_cnt, var_wei, var_guild))
    #-----------------------------------------------------------------------------
    # Define requested parameters
    if (is.null(save_parameters)) {
      save_parameters <- datamodel$popdyn_parameters
    } else {
      if (FALSE %in% is.element(save_parameters, datamodel$popdyn_parameters)) {
        para_name <-
          save_parameters[which(is.element(save_parameters, datamodel$popdyn_parameters) %in% FALSE)]
        abort(paste0(
          "Some parameters are not in model: '",
          para_name,
          "'",
          collapse = " "
        ))
      }
    }
    #-----------------------------------------------------------------------------
    # Fit model
    set.seed(123)
    popdyn <- nimbleModel(
      code = popdyn_code,
      constants = datamodel$popdyn_const,
      data = datamodel$popdyn_data,
      inits = datamodel$popdyn_inits,
      name = "popdyn",
      calculate = FALSE
    )
    popdynConf <- configureMCMC(popdyn, monitors = save_parameters)
    popdynMCMC <- buildMCMC(popdynConf)
    popdynComp <- compileNimble(popdyn)
    popdynModel <-
      compileNimble(popdynMCMC, project = popdyn, resetFunctions = TRUE)
    mcmc_chain <-
      runMCMC(
        popdynModel,
        nchains = n_chain,
        niter = n_iter,
        thin = n_thin,
        nburnin = n_burnin,
        setSeed = 123,
        samplesAsCodaMCMC = TRUE
      )
    if (n_chain == 1 & !is.null(dim(mcmc_chain))) {
      mcmc_na <- which(is.na(mcmc_chain[1, ]))
      if (length(mcmc_na) > 0) {
        mcmc_chain <- mcmc_chain[, -mcmc_na]
      }
    }
    if (n_chain > 1 &
        !is.null(dim(mcmc_chain[[1]]))) {
      for (i in 1:n_chain) {
        mcmc_chain[[i]] <-
          mcmc_chain[[i]][, !colnames(mcmc_chain[[i]]) %in% names(which(is.na(mcmc_chain[[i]][1, ])))]
      }
    }
    #-----------------------------------------------------------------------------
    # set summary data frame and list of subscripts
    mcmc_summary <-
      MCMCvis::MCMCsummary(
        mcmc_chain,
        probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
        Rhat = TRUE,
        n.eff = TRUE
      )
    list_summary <-
      do.call(
        int_transformsummary,
        list(mcmc_summary, datamodel, var_id, var_tmp, var_tax, period)
      )
    output <-
      list(
        mcmc_summary = list_summary$mcmc_summary,
        mcmc_chain = mcmc_chain,
        subscript = list_summary$subscript
      )
    #-----------------------------------------------------------------------------
    return(output)
  }
