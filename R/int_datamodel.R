#' Title
#'
#' @param df
#' @param occup
#' @param grow
#' @param modenv
#' @param modenvG
#' @param timestep
#' @param period
#' @param var_id
#' @param var_tmp
#' @param var_tax
#' @param var_pres
#' @param var_reg
#' @param var_guild
#' @param var_cnt
#' @param var_wei
#' @param var_surf
#' @param var_envO
#' @param var_envP
#' @param var_envC
#' @param var_grow
#'
#' @return
#' @export
#'
#' @examples
int_datamodel <- function(df, occup, grow, modenv, modenvG, timestep, period, var_id, var_tmp, var_tax, var_pres, var_reg, var_guild, var_cnt, var_wei, var_surf, var_envO, var_envP, var_envC, var_grow) {
  guild <- varenv <- varcol <- vargrow <- NULL
  #-----------------------------------------------------------------------------
  if (TRUE %in% quo_is_null(var_tax)) {
    df$taxa <- 1
    var_tax <- "taxa"
    var_tax <- enquo(var_tax)
  }
  tax <- sort(unique(pull(df, !!var_tax)))
  ntax <- length(tax)
  #-----------------------------------------------------------------------------
  time <- seq(min(pull(df, !!var_tmp)),max(pull(df, !!var_tmp)),timestep)
  ntime <- length(time)
  #-----------------------------------------------------------------------------
  id <- sort(unique(pull(df, !!var_id)))
  nid <-  mapply(function(i) length(unique(pull(df[df[,quo_name(var_tax)] %in% i,], !!var_id))), tax, SIMPLIFY = "vector")
  idtax <- do.call("rbind",lapply(1:ntax, function(i) c(which(id %in% sort(unique(pull(df[df[,quo_name(var_tax)] %in% tax[i],], !!var_id)))), rep(NA, max(nid) - nid[i]))))
  #-----------------------------------------------------------------------------
  df_reg <- select(df, !!var_tax, !!var_id) %>% distinct() %>% mutate(reg = "global")
  if (FALSE %in% quo_is_null(var_reg)) {
    nr <- which(colnames(df) %in% colnames(select(df, !!var_reg)))
    df_reg <- df_reg %>%
      rbind(do.call("rbind",lapply(nr, function(i)  data.frame(filter(df, !is.na(df[,i])) %>%
                                                                 select(!!var_tax, !!var_id, all_of(i)) %>% distinct() %>%
                                                                 rename(reg = colnames(df[i])) %>%
                                                                 rowwise() %>% mutate(reg = paste(colnames(df[i]),reg,sep="_"))))))
  }
  region <- sort(unique(df_reg$reg))
  # regions per taxa
  nreg <- tapply(df_reg$reg, list(df_reg[,quo_name(var_tax)]), function(i) length(unique(i)))
  reg <- do.call("rbind",lapply(1:ntax, function(i) c(which(region %in% sort(unique(pull(df_reg[df_reg[,quo_name(var_tax)] %in% tax[i],], reg)))), rep(NA, max(nreg) - nreg[i]))))
  # id per taxa and regions
  nidreg <- do.call("rbind", lapply(1:ntax, function(i) c(do.call("cbind",lapply(1:nreg[i], function(j)
    length(unique(pull(df_reg[(df_reg[,quo_name(var_tax)] %in% tax[i] & df_reg$reg %in% region[reg[i,j]]),], !!var_id)))
  )), rep(NA, max(nreg) - nreg[i]))))
  idreg <- sapply(1:ntax, function(i) rbind(do.call("rbind",lapply(1:nreg[i], function(j)
    c(which(id %in% sort(unique(pull(df_reg[(df_reg[,quo_name(var_tax)] %in% tax[i] & df_reg$reg %in% region[reg[i,j]]),], !!var_id)))), rep(NA, max(nidreg, na.rm = T) - nidreg[i,j]))
  )), array(NA, dim=c(max(nreg) - nreg[i], max(nidreg, na.rm = T)))), simplify = "array")
  #-----------------------------------------------------------------------------
  popdyn_const <- list(step = timestep, ntax = ntax, nid = nid, ntime = ntime, idtax = idtax, nreg = nreg, reg = reg, nidreg = nidreg, idreg = idreg)
  #-----------------------------------------------------------------------------
  if (is.null(period)) {
    start <- c(2,NA)
    end <- c(ntime,NA)
    ndate <- c(ntime - 1,NA)
    nperiod <- 1
  } else {
    nperiod <- 1 + length(period)
    start <- c(2,mapply(function(i) which(time %in% (period[[i]][1] + timestep)), 1:length(period), SIMPLIFY = "vector"))
    end <- c(ntime,mapply(function(i) which(time %in% period[[i]][2]), 1:length(period), SIMPLIFY = "vector"))
    ndate <- mapply(function(i) end[i] - start[i], 1:nperiod,  SIMPLIFY = "vector")
  }
  popdyn_const <- c(popdyn_const, list(nperiod = nperiod, start = start, end = end))
  #-----------------------------------------------------------------------------
  if (FALSE %in% quo_is_null(var_guild)) {
    ng <- which(colnames(df) %in% colnames(select(df, !!var_guild)))
    df_guild <- rbind(do.call("rbind",lapply(ng, function(i) data.frame(filter(df, !is.na(df[,i])) %>%
                                                                          select(!!var_tax, !!var_id, all_of(i)) %>% distinct() %>%
                                                                          rename(guild = colnames(df[i])) %>%
                                                                          rowwise() %>% mutate(guild = paste(colnames(df[i]),guild,sep="_")))))) %>% merge(df_reg)
    # number of guilds
    guild <- sort(unique(df_guild$guild))
    ngui <- length(guild)
    # regions per guilds
    nregui <- tapply(df_guild$reg, list(df_guild$guild), function(i) length(unique(i)))
    regui <- do.call("rbind",lapply(1:ngui, function(i) c(which(region %in% sort(unique(df_guild$reg[df_guild$guild %in% guild[i]]))), rep(NA,max(nregui)-nregui[i]))))
    # taxa per guilds and region
    ngtax <- do.call("rbind",lapply(1:ngui, function(i) c(do.call("cbind",lapply(1:nregui[i], function(j)
      length(unique(pull(df_guild[(df_guild$guild %in% guild[i] & df_guild$reg %in% region[regui[i,j]]),], !!var_tax)))
    )), rep(NA, max(nregui) - nregui[i]))))
    gtax <- sapply(1:ngui, function(i) rbind(do.call("rbind",lapply(1:nregui[i], function(j)
      c(which(tax %in% sort(unique(pull(df_guild[(df_guild$guild %in% guild[i] & df_guild$reg %in% region[regui[i,j]]),], !!var_tax)))), rep(NA, max(ngtax, na.rm = T) - ngtax[i,j]))
    )), array(NA, dim=c(max(nregui) - nregui[i], max(ngtax, na.rm = T)))), simplify = "array")
    if (is.null(dim(gtax))) {
      gtax <- sapply(1:ngui, function(i) rbind(do.call("rbind",lapply(1:nregui[i], function(j)
        c(which(tax %in% sort(unique(pull(df_guild[(df_guild$guild %in% guild[i] & df_guild$reg %in% region[regui[i,j]]),], !!var_tax)))), rep(NA, max(ngtax, na.rm = T) - ngtax[i,j]))
      )), array(NA, dim=c(1, max(ngtax, na.rm = T)))), simplify = "array")
    }
    # add constants and parameters for guild
    popdyn_const <- c(popdyn_const, list(ngui = ngui, nregui = nregui, regui = regui, ngtax = ngtax, gtax = gtax))
  }
  #-----------------------------------------------------------------------------
  z <- tapply(df[,quo_name(var_pres)], list(df[,quo_name(var_tax)],df[,quo_name(var_id)],df[,quo_name(var_tmp)]), unique)
  z <- replace(z, which(z > 0), 1)
  popdyn_data <- list(z = z)
  popdyn_inits <- list(p.per = rep(0.5, ntax), p.col = rep(0.5, ntax), epsilon_per = nimMatrix(0, ntax, length(id)), epsilon_col = nimMatrix(0, ntax, length(id)))
  popdyn_parameters <- c("z","p.per","p.col","p.per_id","p.col_id")
  if (isTRUE(occup)) {
    popdyn_parameters <- c(popdyn_parameters, "p.ext_id","p.ext","z_mulambda","OR","turnover")
    popdyn_const <- c(popdyn_const, list(ndate = ndate))
    if (FALSE %in% quo_is_null(var_guild)) {
      popdyn_parameters <- c(popdyn_parameters,"z_lambda_gui","z_mulambda_gui","GOR")
    }
    if (!is.null(modenv)) {
      popdyn_parameters <- popdyn_parameters[!popdyn_parameters %in% c("p.col","p.per","p.ext","p.col_id","p.per_id","p.ext_id")]
      popdyn_inits <- popdyn_inits[names(popdyn_inits) %in% c("p.per","p.col","epsilon_per","epsilon_col") == FALSE]
      if (FALSE %in% quo_is_null(var_envP) & FALSE %in% quo_is_null(var_envC)) {
        varenv <- colnames(select(df, !!var_envP))
        nvar <- length(varenv)
        varcol <- c(varenv[varenv %in% colnames(select(df, !!var_envC))],colnames(select(df, !!var_envC))[!colnames(select(df, !!var_envC)) %in% varenv])
        nr <- which(colnames(df) %in% varcol)
        ncol <- length(nr)
        col <- sapply(nr, function(i) tapply(df[,i], list(df[,quo_name(var_id)],df[,quo_name(var_tmp)]), unique),simplify = "array")
        popdyn_inits <- list(alpha_per = rep(0, ntax), beta_per = nimMatrix(0, ntax, nvar),alpha_col =  rep(0, ntax), beta_col = nimMatrix(0, ntax, nvar))
        popdyn_parameters <- c(popdyn_parameters, "alpha_per", "beta_per","alpha_col", "beta_col")
        popdyn_const <- c(popdyn_const, list(ncol = ncol))
        popdyn_data <- c(popdyn_data, list(col = col))
      } else if (FALSE %in% quo_is_null(var_envP)) {
        varenv <- colnames(select(df, !!var_envP))
        nvar <- length(varenv)
        popdyn_inits <- list(alpha_per =  rep(0, ntax), beta_per = nimMatrix(0, ntax, nvar), p.col = rep(0.5, ntax), epsilon_col = nimMatrix(0, ntax, length(id)))
        popdyn_parameters <- c(popdyn_parameters, "alpha_per", "beta_per","p.col","p.col_id")
      } else if (FALSE %in% quo_is_null(var_envC)) {
        varenv <- colnames(select(df, !!var_envC))
        nvar <- length(varenv)
        popdyn_inits <- list(alpha_col =  rep(0, ntax), beta_col = nimMatrix(0, ntax, nvar), p.per = rep(0.5, ntax), epsilon_per = nimMatrix(0, ntax, length(id)))
        popdyn_parameters <- c(popdyn_parameters, "alpha_col", "beta_col","p.per","p.per_id")
      } else {
        varenv <- colnames(select(df, !!var_envO))
        nvar <- length(varenv)
        popdyn_inits <- list(alpha =  rep(0, ntax), beta = nimMatrix(0, ntax, nvar))
        popdyn_parameters <- c(popdyn_parameters, "alpha", "beta")
      }
      nr <- which(colnames(df) %in% varenv)
      var <- sapply(nr, function(i) tapply(df[,i], list(df[,quo_name(var_id)],df[,quo_name(var_tmp)]), unique),simplify = "array")
      popdyn_data <- c(popdyn_data,list(var = var))
      popdyn_const <- c(popdyn_const, list(nvar = nvar))
    }
  }
  #-----------------------------------------------------------------------------
  if (isTRUE(grow)) {
    if (FALSE %in% quo_is_null(var_surf)) {
      S <- tapply(df[,quo_name(var_surf)], list(df[,quo_name(var_tax)],df[,quo_name(var_id)],df[,quo_name(var_tmp)]), unique)
    } else {
      S <- nimArray(1, dim = c(ntax, length(id), ntime))
    }
    popdyn_data <- c(popdyn_data, list(S = S))
    #---------------------------------------------------------------------------
    if (!is.null(modenvG)) {
      varoccup <- c(varenv,varcol)
      vargrow <- c(varoccup[varoccup %in% colnames(select(df, !!var_grow))],colnames(select(df, !!var_grow))[!colnames(select(df, !!var_grow)) %in% varoccup])
      nr <- which(colnames(df) %in% vargrow)
      ngvar <- length(nr)
      gvar <- sapply(nr, function(i) tapply(df[,i], list(df[,quo_name(var_id)],df[,quo_name(var_tmp)]), unique),simplify = "array")
      popdyn_data <- c(popdyn_data,list(gvar = gvar))
      popdyn_const <- c(popdyn_const, list(ngvar = ngvar))
    }
    #---------------------------------------------------------------------------
    if (FALSE %in% quo_is_null(var_cnt)) {
      y <- tapply(df[,quo_name(var_cnt)], list(df[,quo_name(var_tax)],df[,quo_name(var_id)],df[,quo_name(var_tmp)]), unique)
      y <- replace(y, which(y == 0), 1)
      C <- log(replace(y, is.na(y), 1))
      popdyn_data <- c(popdyn_data, list(y = y))
      popdyn_parameters <- c(popdyn_parameters, "N_PGR","N_mulambda","N_lambda","N_mulambda_id","N_lambda_id","N")
      if (!is.null(modenvG)) {
        popdyn_inits <- c(popdyn_inits, list(C = C, alpha_N = rep(0, ntax), beta_N = nimMatrix(1, ntax, ngvar), tauN = nimMatrix(1, ntax, length(id))))
        popdyn_parameters <- c(popdyn_parameters, "alpha_N","beta_N")
      } else {
        popdyn_inits <- c(popdyn_inits, list(C = C, muN = nimMatrix(1, ntax, length(id)), tauN = nimMatrix(1, ntax, length(id))))
      }
      if (FALSE %in% quo_is_null(var_guild)) {
        popdyn_parameters <- c(popdyn_parameters, "N_lambda_gui","N_mulambda_gui","N_GGR")
      }
    }
    #-----------------------------------------------------------------------------
    if (FALSE %in% quo_is_null(var_wei)) {
      w <- tapply(df[,quo_name(var_wei)], list(df[,quo_name(var_tax)],df[,quo_name(var_id)],df[,quo_name(var_tmp)]), unique)
      w <- replace(w, which(w == 0), 1)
      W <- log(replace(w, is.na(w), 1))
      popdyn_data <- c(popdyn_data, list(w = w))
      popdyn_parameters <- c(popdyn_parameters, "B_PGR","B_mulambda","B_lambda","B_mulambda_id","B_lambda_id","B")
      if (!is.null(modenvG)) {
        popdyn_inits <- c(popdyn_inits, list(W = W, alpha_B = rep(0, ntax), beta_B = nimMatrix(1, ntax, ngvar), tauB = nimMatrix(1, ntax, length(id))))
        popdyn_parameters <- c(popdyn_parameters, "alpha_B","beta_B")
      } else {
        popdyn_inits <- c(popdyn_inits, list(W = W, muB = nimMatrix(1, ntax, length(id)), tauB = nimMatrix(1, ntax, length(id))))
      }
      if (FALSE %in% quo_is_null(var_guild)) {
        popdyn_parameters <- c(popdyn_parameters, "B_lambda_gui","B_mulambda_gui","B_GGR")
      }
    }
  }
  #-----------------------------------------------------------------------------
  popdyn_int <- list(taxa = tax, id = id, time = time, region = region, guild = guild, varenv = varenv, varcol = varcol, vargrow = vargrow, start = start, end = end, protocol = NULL, vardet = NULL)
  return(list(popdyn_data = popdyn_data, popdyn_const = popdyn_const, popdyn_inits = popdyn_inits, popdyn_parameters = popdyn_parameters, popdyn_int = popdyn_int))
}
################################################################################
