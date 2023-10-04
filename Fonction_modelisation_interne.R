int_popoccup <- function(modenv,var_envO,var_envP,var_envC,var_guild) { 
  if (is.null(modenv)) {
    #-----------------------------------------------------------------------------
    # Occupancy model without environmental effects
    occupancy <- nimbleCode(
      for(s in 1:ntax) {
        for(i in 1:nid[s]) {
          for (t in 2:ntime) {
            z[s,idtax[s,i],t] ~ dbern(p[s,idtax[s,i],t])
            p[s,idtax[s,i],t] <- z[s,idtax[s,i],t-1] * p.per_id[s,idtax[s,i]] + (1 - z[s,idtax[s,i],t-1]) * p.col_id[s,idtax[s,i]]
          }
          z[s,idtax[s,i],1] ~ dbern(p.per_id[s,idtax[s,i]])
          tauper[s,idtax[s,i]] ~ dgamma(0.1, 0.1)
          taucol[s,idtax[s,i]] ~ dgamma(0.1, 0.1)
          epsilon_per[s,idtax[s,i]] ~ dnorm(0, sd = tauper[s,idtax[s,i]])
          epsilon_col[s,idtax[s,i]] ~ dnorm(0, sd = taucol[s,idtax[s,i]])
          logit(p.per_id[s,idtax[s,i]]) <- p.per[s] + epsilon_per[s,idtax[s,i]]
          logit(p.col_id[s,idtax[s,i]]) <- p.col[s] + epsilon_col[s,idtax[s,i]]
          p.ext_id[s,idtax[s,i]] <- 1 - p.per_id[s,idtax[s,i]]
        }
        p.per[s] ~ dunif(0,1) 
        p.col[s] ~ dunif(0,1)
        p.ext[s] <- 1 - p.per[s]
      }
    )
  } else if (FALSE %in% quo_is_null(var_envP) & FALSE %in% quo_is_null(var_envC)) {
    #-----------------------------------------------------------------------------
    # Occupancy model with linear environmental effects on persistence and colonisation probabilities
    occupancy <- nimbleCode(
      for(s in 1:ntax) {
        for(i in 1:nid[s]) {
          for (t in 2:ntime) {
            z[s,idtax[s,i],t] ~ dbern(p[s,idtax[s,i],t])
            p[s,idtax[s,i],t] <- z[s,idtax[s,i],t-1] * p.per_id[s,idtax[s,i],t] + (1 - z[s,idtax[s,i],t-1]) * p.col_id[s,idtax[s,i],t]
            if (nvar > 1) {
              logit(p.per_id[s,idtax[s,i],t]) <- alpha_per[s] + sum(beta_per[s,1:nvar] * var[idtax[s,i],t,1:nvar]) 
            } else {
              logit(p.per_id[s,idtax[s,i],t]) <- alpha_per[s] + beta_per[s,1] * var[idtax[s,i],t,1] 
            } 
            if (ncol > 1) {
              logit(p.col_id[s,idtax[s,i],t]) <- alpha_col[s] + sum(beta_col[s,1:ncol] * col[idtax[s,i],t,1:ncol]) 
            } else {
              logit(p.col_id[s,idtax[s,i],t]) <- alpha_col[s] + beta_col[s,1] * col[idtax[s,i],t,1] 
            }  
          }
          z[s,idtax[s,i],1] ~ dbern(0.5)
        }
        taualpha_col[s] ~ dgamma(0.1, 0.1)
        alpha_col[s] ~ dnorm(0, sd = taualpha_col[s])
        taualpha_per[s] ~ dgamma(0.1, 0.1)
        alpha_per[s] ~ dnorm(0, sd = taualpha_per[s])
        for (k in 1:nvar) {
          taubeta_per[s,k] ~ dgamma(0.1, 0.1)
          beta_per[s,k] ~ dnorm(0, sd = taubeta_per[s,k])
        }
        for (k in 1:ncol) {
          taubeta_col[s,k] ~ dgamma(0.1, 0.1)
          beta_col[s,k] ~ dnorm(0, sd = taubeta_col[s,k])
        }
      }
    ) 
  } else if (FALSE %in% quo_is_null(var_envP)) {
    #-----------------------------------------------------------------------------
    # Occupancy model with linear environmental effects on persistence probability
    occupancy <- nimbleCode(
      for(s in 1:ntax) {
        for(i in 1:nid[s]) {
          for (t in 2:ntime) {
            z[s,idtax[s,i],t] ~ dbern(p[s,idtax[s,i],t])
            p[s,idtax[s,i],t] <- z[s,idtax[s,i],t-1] * p.per_id[s,idtax[s,i],t] + (1 - z[s,idtax[s,i],t-1]) * p.col_id[s,idtax[s,i]]
            if (nvar > 1) {
              logit(p.per_id[s,idtax[s,i],t]) <- alpha_per[s] + sum(beta_per[s,1:nvar] * var[idtax[s,i],t,1:nvar]) 
            } else {
              logit(p.per_id[s,idtax[s,i],t]) <- alpha_per[s] + beta_per[s,1] * var[idtax[s,i],t,1] 
            }  
          }
          z[s,idtax[s,i],1] ~ dbern(0.5)
          taucol[s,idtax[s,i]] ~ dgamma(0.1, 0.1)
          epsilon_col[s,idtax[s,i]] ~ dnorm(0, sd = taucol[s,idtax[s,i]])
          logit(p.col_id[s,idtax[s,i]]) <- p.col[s] + epsilon_col[s,idtax[s,i]]
        }
        p.col[s] ~ dunif(0,1)
        taualpha[s] ~ dgamma(0.1, 0.1)
        alpha_per[s] ~ dnorm(0, sd = taualpha[s])
        for (k in 1:nvar) {
          taubeta[s,k] ~ dgamma(0.1, 0.1)
          beta_per[s,k] ~ dnorm(0, sd = taubeta[s,k])
        }
      }
    ) 
  } else if (FALSE %in% quo_is_null(var_envC)) {
    #-----------------------------------------------------------------------------
    # Occupancy model with linear environmental effects on colonisation probability
    occupancy <- nimbleCode(
      for(s in 1:ntax) {
        for(i in 1:nid[s]) {
          for (t in 2:ntime) {
            z[s,idtax[s,i],t] ~ dbern(p[s,idtax[s,i],t])
            p[s,idtax[s,i],t] <- z[s,idtax[s,i],t-1] * p.per_id[s,idtax[s,i]] + (1 - z[s,idtax[s,i],t-1]) * p.col_id[s,idtax[s,i],t]
            if (nvar > 1) {
              logit(p.col_id[s,idtax[s,i],t]) <- alpha_col[s] + sum(beta_col[s,1:nvar] * var[idtax[s,i],t,1:nvar]) 
            } else {
              logit(p.col_id[s,idtax[s,i],t]) <- alpha_col[s] + beta_col[s,1] * var[idtax[s,i],t,1] 
            }  
          }
          z[s,idtax[s,i],1] ~ dbern(0.5)
          tauper[s,idtax[s,i]] ~ dgamma(0.1, 0.1)
          epsilon_per[s,idtax[s,i]] ~ dnorm(0, sd = tauper[s,idtax[s,i]])
          logit(p.per_id[s,idtax[s,i]]) <- p.per[s] + epsilon_per[s,idtax[s,i]]
        }
        p.per[s] ~ dunif(0,1)
        taualpha[s] ~ dgamma(0.1, 0.1)
        alpha_col[s] ~ dnorm(0, sd = taualpha[s])
        for (k in 1:nvar) {
          taubeta[s,k] ~ dgamma(0.1, 0.1)
          beta_col[s,k] ~ dnorm(0, sd = taubeta[s,k])
        }
      }
    )
  } else {
    #-----------------------------------------------------------------------------
    # Occupancy model with linear environmental effects on occupancy probability
    occupancy <- nimbleCode(
      for(s in 1:ntax) {
        for(i in 1:nid[s]) {
          for (t in 1:ntime) {
            z[s,idtax[s,i],t] ~ dbern(p[s,idtax[s,i],t])
            if (nvar > 1) {
              logit(p[s,idtax[s,i],t]) <- alpha[s] + sum(beta[s,1:nvar] * var[idtax[s,i],t,1:nvar]) 
            } else {
              logit(p[s,idtax[s,i],t]) <- alpha[s] + beta[s,1] * var[idtax[s,i],t,1] 
            }  
          }
        }
        taualpha[s] ~ dgamma(0.1, 0.1)
        alpha[s] ~ dnorm(0, sd = taualpha[s])
        for (k in 1:nvar) {
          taubeta[s,k] ~ dgamma(0.1, 0.1)
          beta[s,k] ~ dnorm(0, sd = taubeta[s,k])
        }
      }
    ) 
  }
  code <- list(occupancy)
  #-----------------------------------------------------------------------------
  # Occupancy rates at regional levels
  occupancy_regional <- nimbleCode(
    for(s in 1:ntax) {
      for (j in 1:nreg[s]) {
        for (t in 1:nperiod) {
          nz_reg[s,reg[s,j],t] <- max(z_reg[s,reg[s,j],start[t]:end[t]])
          z_mulambda[s,reg[s,j],t] <- nz_reg[s,reg[s,j],t] * prod(z_lambda[s,reg[s,j],start[t]:end[t]])^(1 / (step * ndate[t]))
          OR[s,reg[s,j],t] <- nz_reg[s,reg[s,j],t] * 100 * (z_mulambda[s,reg[s,j],t] - 1)
        }
        for (t in 2:ntime) {
          z_reg[s,reg[s,j],t] <- max(z_cal[s,reg[s,j],1:nidreg[s,j],t])
          nb_z[s,reg[s,j],t] <- sum(z_cal[s,reg[s,j],1:nidreg[s,j],t])
          z_lambda[s,reg[s,j],t] <- nb_z[s,reg[s,j],t] / max(1,nb_z[s,reg[s,j],t-1])
          cal.z_lambda[s,reg[s,j],t] <- z_reg[s,reg[s,j],t] * z_lambda[s,reg[s,j],t] + (1 - z_reg[s,reg[s,j],t])
          turnover[s,reg[s,j],t] <- sum(turn_cal[s,reg[s,j],1:nidreg[s,j],t]) / max(1,nb_z[s,reg[s,j],t])
          for (i in 1:nidreg[s,j]) {
            turn_cal[s,reg[s,j],i,t] <- z[s,idreg[j,i,s],t] * (1 - z[s,idreg[j,i,s],t-1])
            z_cal[s,reg[s,j],i,t] <- z[s,idreg[j,i,s],t]
          }
        }
        nb_z[s,reg[s,j],1] <- sum(z_cal[s,reg[s,j],1:nidreg[s,j],1])
        for (i in 1:nidreg[s,j]) {
          z_cal[s,reg[s,j],i,1] <- z[s,idreg[j,i,s],1]
        }
      }
    }
  )
  code <- c(code,list(occupancy_regional))
  #-----------------------------------------------------------------------------
  # Occupancy rates for guilds
  if (FALSE %in% quo_is_null(var_guild)) {
    occupancy_guild <- nimbleCode(
      for(g in 1:ngui) {
        for (j in 1:nregui[g]) {
          for (t in 1:nperiod) {
            z_guireg[g,regui[g,j],t] <- max(z_gui[g,regui[g,j],start[t]:end[t]])
            z_mulambda_gui[g,regui[g,j],t] <- z_guireg[g,regui[g,j],t] * prod(z_lambda_gui[g,regui[g,j],start[t]:end[t]])^(1 / (step * ndate[t]))
            GOR[g,regui[g,j],t] <- z_guireg[g,regui[g,j],t] * 100 * (z_mulambda_gui[g,regui[g,j],t] - 1)
          }
          for (t in 2:ntime) {
            z_lambda_gui[g,regui[g,j],t] <- prod(z_lambda_tax[g,1:ngtax[g,j],regui[g,j],t])^(1/ngtax[g,j])
            z_gui[g,regui[g,j],t] <- max(z_reg_tax[g,1:ngtax[g,j],regui[g,j],t])
            cal.z_lambda_gui[g,regui[g,j],t] <- z_gui[g,regui[g,j],t] * z_lambda_gui[g,regui[g,j],t] + (1 - z_gui[g,regui[g,j],t])
            
            for(s in 1:ngtax[g,j]) {
              z_lambda_tax[g,s,regui[g,j],t] <- z_lambda[gtax[j,s,g],regui[g,j],t]
              z_reg_tax[g,s,regui[g,j],t] <- z_reg[gtax[j,s,g],regui[g,j],t]
            }
          }
        }
      }
    )
    code <- c(code,list(occupancy_guild))
  }  
  return(code)
}
################################################################################
int_popgrow <- function(occup, code, modenvG, var_cnt, var_wei, var_guild) {
  ## Occupancy component
  if (isFALSE(occup)) {
    occupancy <- nimbleCode(
      for(s in 1:ntax) {
        for(i in 1:nid[s]) {
          for (t in 2:ntime) {
            z[s,idtax[s,i],t] ~ dbern(p[s,idtax[s,i],t])
            p[s,idtax[s,i],t] <- z[s,idtax[s,i],t-1] * p.per_id[s,idtax[s,i]] + (1 - z[s,idtax[s,i],t-1]) * p.col_id[s,idtax[s,i]]
          }
          z[s,idtax[s,i],1] ~ dbern(p.per_id[s,idtax[s,i]])
          p.per_id[s,idtax[s,i]] ~ dunif(0,1) 
          p.col_id[s,idtax[s,i]] ~ dunif(0,1)
        }  
      }
    )
    code <- list(occupancy)
  }
  #-----------------------------------------------------------------------------   
  ## Growth models from headcounts
  if (FALSE %in% quo_is_null(var_cnt)) {
    if (is.null(modenvG)) {
      #-------------------------------------------------------------------------
      # Growth model from headcounts without environmental effects
      counts <- nimbleCode(
        for (s in 1:ntax) {
          for (j in 1:nreg[s]) {
            for (n in 1:nperiod) {
              n_reg[s,reg[s,j],n] <- max(n_regt[s,reg[s,j],start[n]:end[n]])
              N_mulambda[s,reg[s,j],n] <- n_reg[s,reg[s,j],n] * prod(cal.N_lambda[s,reg[s,j],start[n]:end[n]])^(1/ max(1, step * sum(n_regt[s,reg[s,j],start[n]:end[n]])))
              N_PGR[s,reg[s,j],n] <- n_reg[s,reg[s,j],n] * 100 * (N_mulambda[s,reg[s,j],n] - 1)
            }
            for (t in 2:ntime) {
              n_regt[s,reg[s,j],t] <- max(n_zreg[s,reg[s,j],1:nidreg[s,j],t]) * max(n_zreg[s,reg[s,j],1:nidreg[s,j],t-1])
              N_lambda[s,reg[s,j],t] <- n_regt[s,reg[s,j],t] * sum(cal.N[s,reg[s,j],1:nidreg[s,j],t]) / max(1,sum(cal.N[s,reg[s,j],1:nidreg[s,j],t-1]))
              cal.N_lambda[s,reg[s,j],t] <- N_lambda[s,reg[s,j],t] + (1 - n_regt[s,reg[s,j],t])
            }
            for (t in 1:ntime) {
              for (i in 1:nidreg[s,j]) {
                n_zreg[s,reg[s,j],i,t] <- z[s,idreg[j,i,s],t]
                cal.N[s,reg[s,j],i,t] <- N[s,idreg[j,i,s],t] / S[s,idreg[j,i,s],t]
              }
            }
          }
          for (i in 1:nid[s]) {
            for (n in 1:nperiod) {
              n_id[s,idtax[s,i],n] <- max(n_idt[s,idtax[s,i],start[n]:end[n]])
              N_mulambda_id[s,idtax[s,i],n] <-  n_id[s,idtax[s,i],n] * prod(cal.Nlambda_id[s,idtax[s,i],start[n]:end[n]])^(1/(step * max(1,sum(n_idt[s,idtax[s,i],start[n]:end[n]]))))
            }
            for (t in 2:ntime) {
              C[s,idtax[s,i],t] ~ dgamma(1,1)
              log_N[s,idtax[s,i],t] <- z[s,idtax[s,i],t-1] * log(muN[s,idtax[s,i]] * y[s,idtax[s,i],t-1] * (S[s,idtax[s,i],t-1] / S[s,idtax[s,i],t])) + (1 - z[s,idtax[s,i],t-1]) * C[s,idtax[s,i],t]
              y[s,idtax[s,i],t] ~ dlnorm(log_N[s,idtax[s,i],t], tauN[s,idtax[s,i]])
              N[s,idtax[s,i],t] <- z[s,idtax[s,i],t] * y[s,idtax[s,i],t] 
              
              n_idt[s,idtax[s,i],t] <- z[s,idtax[s,i],t] * z[s,idtax[s,i],t-1]
              N_lambda_id[s,idtax[s,i],t] <- n_idt[s,idtax[s,i],t] * (y[s,idtax[s,i],t] / y[s,idtax[s,i],t-1]) * (S[s,idtax[s,i],t] / S[s,idtax[s,i],t-1])
              cal.Nlambda_id[s,idtax[s,i],t] <- N_lambda_id[s,idtax[s,i],t] + (1 - n_idt[s,idtax[s,i],t]) 
            }
            log_N[s,idtax[s,i],1] ~ dgamma(1,1)
            y[s,idtax[s,i],1] ~ dlnorm(log_N[s,idtax[s,i],1], tauN[s,idtax[s,i]])
            N[s,idtax[s,i],1] <- z[s,idtax[s,i],1] * y[s,idtax[s,i],1] 
            n_idt[s,idtax[s,i],1] <- 0
            muN[s,idtax[s,i]] ~ dgamma(0.01,0.01)
            tauN[s,idtax[s,i]] ~ dgamma(0.01,0.01)
          }
        }
      )
    } else {
      #---------------------------------------------------------------------------
      # Growth model from headcounts with linear environmental effects on rates
      counts <- nimbleCode(
        for (s in 1:ntax) {
          for (j in 1:nreg[s]) {
            for (n in 1:nperiod) {
              n_reg[s,reg[s,j],n] <- max(n_regt[s,reg[s,j],start[n]:end[n]])
              N_mulambda[s,reg[s,j],n] <- n_reg[s,reg[s,j],n] * prod(cal.N_lambda[s,reg[s,j],start[n]:end[n]])^(1/ max(1, step * sum(n_regt[s,reg[s,j],start[n]:end[n]])))
              N_PGR[s,reg[s,j],n] <- n_reg[s,reg[s,j],n] * 100 * (N_mulambda[s,reg[s,j],n] - 1)
            }
            for (t in 2:ntime) {
              n_regt[s,reg[s,j],t] <- max(n_zreg[s,reg[s,j],1:nidreg[s,j],t]) * max(n_zreg[s,reg[s,j],1:nidreg[s,j],t-1])
              N_lambda[s,reg[s,j],t] <- n_regt[s,reg[s,j],t] * sum(cal.N[s,reg[s,j],1:nidreg[s,j],t]) / max(1,sum(cal.N[s,reg[s,j],1:nidreg[s,j],t-1]))
              cal.N_lambda[s,reg[s,j],t] <- N_lambda[s,reg[s,j],t] + (1 - n_regt[s,reg[s,j],t])
            }
            for (t in 1:ntime) {
              for (i in 1:nidreg[s,j]) {
                n_zreg[s,reg[s,j],i,t] <- z[s,idreg[j,i,s],t]
                cal.N[s,reg[s,j],i,t] <- N[s,idreg[j,i,s],t] / S[s,idreg[j,i,s],t]
              }
            }
          }
          for (i in 1:nid[s]) {
            for (n in 1:nperiod) {
              n_id[s,idtax[s,i],n] <- max(n_idt[s,idtax[s,i],start[n]:end[n]])
              N_mulambda_id[s,idtax[s,i],n] <-  n_id[s,idtax[s,i],n] * prod(cal.Nlambda_id[s,idtax[s,i],start[n]:end[n]])^(1/(step * max(1,sum(n_idt[s,idtax[s,i],start[n]:end[n]]))))
            }
            for (t in 2:ntime) {
              if (ngvar > 1) {
                log(muN[s,idtax[s,i],t]) <- alpha_N[s] + sum(beta_N[s,1:ngvar] * gvar[idtax[s,i],t,1:ngvar])  
              } else {
                log(muN[s,idtax[s,i],t]) <- alpha_N[s] + beta_N[s,1] * gvar[idtax[s,i],t,1]  
              }
              C[s,idtax[s,i],t] ~ dgamma(1,1)
              log_N[s,idtax[s,i],t] <- z[s,idtax[s,i],t-1] * log(muN[s,idtax[s,i],t] * y[s,idtax[s,i],t-1] * (S[s,idtax[s,i],t-1] / S[s,idtax[s,i],t])) + (1 - z[s,idtax[s,i],t-1]) * C[s,idtax[s,i],t]
              y[s,idtax[s,i],t] ~ dlnorm(log_N[s,idtax[s,i],t], tauN[s,idtax[s,i]])
              N[s,idtax[s,i],t] <- z[s,idtax[s,i],t] * y[s,idtax[s,i],t] 
              
              n_idt[s,idtax[s,i],t] <- z[s,idtax[s,i],t] * z[s,idtax[s,i],t-1]
              N_lambda_id[s,idtax[s,i],t] <- n_idt[s,idtax[s,i],t] * (y[s,idtax[s,i],t] / y[s,idtax[s,i],t-1]) * (S[s,idtax[s,i],t] / S[s,idtax[s,i],t-1])
              cal.Nlambda_id[s,idtax[s,i],t] <- N_lambda_id[s,idtax[s,i],t] + (1 - n_idt[s,idtax[s,i],t]) 
            }
            log_N[s,idtax[s,i],1] ~ dgamma(1,1)
            y[s,idtax[s,i],1] ~ dlnorm(log_N[s,idtax[s,i],1], tauN[s,idtax[s,i]])
            N[s,idtax[s,i],1] <- z[s,idtax[s,i],1] * y[s,idtax[s,i],1] 
            n_idt[s,idtax[s,i],1] <- 0
            tauN[s,idtax[s,i]] ~ dgamma(0.01,0.01)
          }
          taualpha_N[s] ~ dgamma(0.1, 0.1)
          alpha_N[s] ~ dnorm(0, taualpha_N[s])
          for (k in 1:ngvar) {
            taubeta_N[s,k] ~ dgamma(0.1, 0.1)
            beta_N[s,k] ~ dnorm(0, taubeta_N[s,k])
          }
        }
      )
    }
    code <- c(code,list(counts))
    #---------------------------------------------------------------------------
    # Growth rates for guilds
    if (FALSE %in% quo_is_null(var_guild)) {
      counts_guild <- nimbleCode(
        for(g in 1:ngui) {
          for (j in 1:nregui[g]) {
            for (n in 1:nperiod) {
              n_gui[g,regui[g,j],n] <- max(n_guit[g,regui[g,j],start[n]:end[n]])
              N_mulambda_gui[g,regui[g,j],n] <- n_gui[g,regui[g,j],n] * prod(cal.N_lambda_gui[g,regui[g,j],start[n]:end[n]])^(1/max(1, step * sum(n_guit[g,regui[g,j],start[n]:end[n]])))
              N_GGR[g,regui[g,j],n] <- n_gui[g,regui[g,j],n] * 100 * (N_mulambda_gui[g,regui[g,j],n] - 1)
            }
            for (t in 2:ntime) {
              n_guit[g,regui[g,j],t] <- max(n_tax[g,regui[g,j],1:ngtax[g,j],t])
              N_lambda_gui[g,regui[g,j],t] <- n_guit[g,regui[g,j],t] * prod(cal.Nlambda_tax[g,regui[g,j],1:ngtax[g,j],t])^(1 / max(1,sum(n_tax[g,regui[g,j],1:ngtax[g,j],t])))
              cal.N_lambda_gui[g,regui[g,j],t] <- N_lambda_gui[g,regui[g,j],t] + (1 - n_guit[g,regui[g,j],t])
              
              for(s in 1:ngtax[g,j]) {
                n_tax[g,regui[g,j],s,t] <- n_regt[gtax[j,s,g],regui[g,j],t]
                cal.Nlambda_tax[g,regui[g,j],s,t] <- cal.N_lambda[gtax[j,s,g],regui[g,j],t]
              }
            }
          }
        }
      )
      code <- c(code,list(counts_guild))
    }
  }
  #---------------------------------------------------------------------------
  ## Growth models from biomass
  if (FALSE %in% quo_is_null(var_wei)) {
    if (is.null(modenvG)) {
      #-------------------------------------------------------------------------
      # Growth model from biomass without environmental effects
      biomass <- nimbleCode(
        for (s in 1:ntax) {
          for (j in 1:nreg[s]) {
            for (n in 1:nperiod) {
              b_reg[s,reg[s,j],n] <- max(b_regt[s,reg[s,j],start[n]:end[n]])
              B_mulambda[s,reg[s,j],n] <- b_reg[s,reg[s,j],n] * prod(cal.B_lambda[s,reg[s,j],start[n]:end[n]])^(1/ max(1, step * sum(b_regt[s,reg[s,j],start[n]:end[n]])))
              B_PGR[s,reg[s,j],n] <- b_reg[s,reg[s,j],n] * 100 * (B_mulambda[s,reg[s,j],n] - 1)
            }
            for (t in 2:ntime) {
              b_regt[s,reg[s,j],t] <- max(b_zreg[s,reg[s,j],1:nidreg[s,j],t]) * max(b_zreg[s,reg[s,j],1:nidreg[s,j],t-1])
              B_lambda[s,reg[s,j],t] <- b_regt[s,reg[s,j],t] * sum(cal.B[s,reg[s,j],1:nidreg[s,j],t]) / max(1,sum(cal.B[s,reg[s,j],1:nidreg[s,j],t-1]))
              cal.B_lambda[s,reg[s,j],t] <- B_lambda[s,reg[s,j],t] + (1 - b_regt[s,reg[s,j],t])
            }
            for (t in 1:ntime) {
              for (i in 1:nidreg[s,j]) {
                b_zreg[s,reg[s,j],i,t] <- z[s,idreg[j,i,s],t]
                cal.B[s,reg[s,j],i,t] <- B[s,idreg[j,i,s],t] / S[s,idreg[j,i,s],t]
              }
            }
          }
          for (i in 1:nid[s]) {
            for (n in 1:nperiod) {
              b_id[s,idtax[s,i],n] <- max(b_idt[s,idtax[s,i],start[n]:end[n]])
              B_mulambda_id[s,idtax[s,i],n] <-  b_id[s,idtax[s,i],n] * prod(cal.Blambda_id[s,idtax[s,i],start[n]:end[n]])^(1/(step * max(1,sum(b_idt[s,idtax[s,i],start[n]:end[n]]))))
            }
            for (t in 2:ntime) {
              W[s,idtax[s,i],t] ~ dgamma(1,1)
              log_B[s,idtax[s,i],t] <- z[s,idtax[s,i],t-1] * log(muB[s,idtax[s,i]] * w[s,idtax[s,i],t-1] * (S[s,idtax[s,i],t-1] / S[s,idtax[s,i],t])) + (1 - z[s,idtax[s,i],t-1]) * W[s,idtax[s,i],t]
              w[s,idtax[s,i],t] ~ dlnorm(log_B[s,idtax[s,i],t], tauB[s,idtax[s,i]])
              B[s,idtax[s,i],t] <- z[s,idtax[s,i],t] * w[s,idtax[s,i],t] 
              
              b_idt[s,idtax[s,i],t] <- z[s,idtax[s,i],t] * z[s,idtax[s,i],t-1]
              B_lambda_id[s,idtax[s,i],t] <- b_idt[s,idtax[s,i],t] * (w[s,idtax[s,i],t] / w[s,idtax[s,i],t-1]) * (S[s,idtax[s,i],t] / S[s,idtax[s,i],t-1])
              cal.Blambda_id[s,idtax[s,i],t] <- B_lambda_id[s,idtax[s,i],t] + (1 - b_idt[s,idtax[s,i],t]) 
            }
            log_B[s,idtax[s,i],1] ~ dgamma(1,1)
            w[s,idtax[s,i],1] ~ dlnorm(log_B[s,idtax[s,i],1], tauB[s,idtax[s,i]])
            B[s,idtax[s,i],1] <- z[s,idtax[s,i],1] * w[s,idtax[s,i],1] 
            b_idt[s,idtax[s,i],1] <- 0
            muB[s,idtax[s,i]] ~ dgamma(0.01,0.01)
            tauB[s,idtax[s,i]] ~ dgamma(0.01,0.01)
          }
        }
      )
    } else {
      #---------------------------------------------------------------------------
      # Growth model from biomass with linear environmental effects on rates
      biomass <- nimbleCode(
        for (s in 1:ntax) {
          for (j in 1:nreg[s]) {
            for (n in 1:nperiod) {
              b_reg[s,reg[s,j],n] <- max(b_regt[s,reg[s,j],start[n]:end[n]])
              B_mulambda[s,reg[s,j],n] <- b_reg[s,reg[s,j],n] * prod(cal.B_lambda[s,reg[s,j],start[n]:end[n]])^(1/ max(1, step * sum(b_regt[s,reg[s,j],start[n]:end[n]])))
              B_PGR[s,reg[s,j],n] <- b_reg[s,reg[s,j],n] * 100 * (B_mulambda[s,reg[s,j],n] - 1)
            }
            for (t in 2:ntime) {
              b_regt[s,reg[s,j],t] <- max(b_zreg[s,reg[s,j],1:nidreg[s,j],t]) * max(b_zreg[s,reg[s,j],1:nidreg[s,j],t-1])
              B_lambda[s,reg[s,j],t] <- b_regt[s,reg[s,j],t] * sum(cal.B[s,reg[s,j],1:nidreg[s,j],t]) / max(1,sum(cal.B[s,reg[s,j],1:nidreg[s,j],t-1]))
              cal.B_lambda[s,reg[s,j],t] <- B_lambda[s,reg[s,j],t] + (1 - b_regt[s,reg[s,j],t])
            }
            for (t in 1:ntime) {
              for (i in 1:nidreg[s,j]) {
                b_zreg[s,reg[s,j],i,t] <- z[s,idreg[j,i,s],t]
                cal.B[s,reg[s,j],i,t] <- B[s,idreg[j,i,s],t] / S[s,idreg[j,i,s],t]
              }
            }
          }
          for (i in 1:nid[s]) {
            for (n in 1:nperiod) {
              b_id[s,idtax[s,i],n] <- max(b_idt[s,idtax[s,i],start[n]:end[n]])
              B_mulambda_id[s,idtax[s,i],n] <-  b_id[s,idtax[s,i],n] * prod(cal.Blambda_id[s,idtax[s,i],start[n]:end[n]])^(1/(step * max(1,sum(b_idt[s,idtax[s,i],start[n]:end[n]]))))
            }
            for (t in 2:ntime) {
              if (ngvar > 1) {
                log(muB[s,idtax[s,i],t]) <- alpha_B[s] + sum(beta_B[s,1:ngvar] * gvar[idtax[s,i],t,1:ngvar])  
              } else {
                log(muB[s,idtax[s,i],t]) <- alpha_B[s] + beta_B[s,1] * gvar[idtax[s,i],t,1]  
              }
              W[s,idtax[s,i],t] ~ dgamma(1,1)
              log_B[s,idtax[s,i],t] <- z[s,idtax[s,i],t-1] * log(muB[s,idtax[s,i],t] * w[s,idtax[s,i],t-1] * (S[s,idtax[s,i],t-1] / S[s,idtax[s,i],t])) + (1 - z[s,idtax[s,i],t-1]) * W[s,idtax[s,i],t]
              w[s,idtax[s,i],t] ~ dlnorm(log_B[s,idtax[s,i],t], tauB[s,idtax[s,i]])
              B[s,idtax[s,i],t] <- z[s,idtax[s,i],t] * w[s,idtax[s,i],t] 
              
              b_idt[s,idtax[s,i],t] <- z[s,idtax[s,i],t] * z[s,idtax[s,i],t-1]
              B_lambda_id[s,idtax[s,i],t] <- b_idt[s,idtax[s,i],t] * (w[s,idtax[s,i],t] / w[s,idtax[s,i],t-1]) * (S[s,idtax[s,i],t] / S[s,idtax[s,i],t-1])
              cal.Blambda_id[s,idtax[s,i],t] <- B_lambda_id[s,idtax[s,i],t] + (1 - b_idt[s,idtax[s,i],t]) 
            }
            log_B[s,idtax[s,i],1] ~ dgamma(1,1)
            w[s,idtax[s,i],1] ~ dlnorm(log_B[s,idtax[s,i],1], tauB[s,idtax[s,i]])
            B[s,idtax[s,i],1] <- z[s,idtax[s,i],1] * w[s,idtax[s,i],1] 
            b_idt[s,idtax[s,i],1] <- 0
            tauB[s,idtax[s,i]] ~ dgamma(0.01,0.01)
          }
          taualpha_B[s] ~ dgamma(0.1, 0.1)
          alpha_B[s] ~ dnorm(0, taualpha_B[s])
          for (k in 1:ngvar) {
            taubeta_B[s,k] ~ dgamma(0.1, 0.1)
            beta_B[s,k] ~ dnorm(0, taubeta_B[s,k])
          }
        }
      )    
    }
    code <- c(code,list(biomass))
    #---------------------------------------------------------------------------
    # Growth rates for guilds
    if (FALSE %in% quo_is_null(var_guild)) {
      biomass_guild <- nimbleCode(
        for(g in 1:ngui) {
          for (j in 1:nregui[g]) {
            for (n in 1:nperiod) {
              b_gui[g,regui[g,j],n] <- max(b_guit[g,regui[g,j],start[n]:end[n]])
              B_mulambda_gui[g,regui[g,j],n] <- b_gui[g,regui[g,j],n] * prod(cal.B_lambda_gui[g,regui[g,j],start[n]:end[n]])^(1/max(1, step * sum(b_guit[g,regui[g,j],start[n]:end[n]])))
              B_GGR[g,regui[g,j],n] <- b_gui[g,regui[g,j],n] * 100 * (B_mulambda_gui[g,regui[g,j],n] - 1)
            }
            for (t in 2:ntime) {
              b_guit[g,regui[g,j],t] <- max(b_tax[g,regui[g,j],1:ngtax[g,j],t])
              B_lambda_gui[g,regui[g,j],t] <- b_guit[g,regui[g,j],t] * prod(cal.Blambda_tax[g,regui[g,j],1:ngtax[g,j],t])^(1 / max(1,sum(b_tax[g,regui[g,j],1:ngtax[g,j],t])))
              cal.B_lambda_gui[g,regui[g,j],t] <- B_lambda_gui[g,regui[g,j],t] + (1 - b_guit[g,regui[g,j],t])
              
              for(s in 1:ngtax[g,j]) {
                b_tax[g,regui[g,j],s,t] <- b_regt[gtax[j,s,g],regui[g,j],t]
                cal.Blambda_tax[g,regui[g,j],s,t] <- cal.B_lambda[gtax[j,s,g],regui[g,j],t]
              }
            }
          }
        }
      )
      code <- c(code,list(biomass_guild))
    }
  }
  return(as.call(c(as.symbol("{"), code)))
}
################################################################################
int_datamodel <- function(df, occup, grow, modenv, modenvG, timestep, period, var_id, var_tmp, var_tax, var_pres, var_reg, var_guild, var_cnt, var_wei, var_surf, var_envO, var_envP, var_envC, var_grow) {
  guild <- varenv <- varcol <- vargrow <- NULL
  #-----------------------------------------------------------------------------
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
  popdyn_parameters <- c("z","p.col_id","p.per_id")
  popdyn_inits <- list(p.per_id = nimMatrix(0.5, ntax, length(id)), p.col_id = nimMatrix(0.5, ntax, length(id))) 
  if (isTRUE(occup)) {                     
    popdyn_parameters <- c(popdyn_parameters, "p.ext_id","p.col","p.per","p.ext","z_mulambda","OR","turnover")
    popdyn_const <- c(popdyn_const, list(ndate = ndate))
    popdyn_inits <- list(p.per = rep(0.5, ntax), p.col = rep(0.5, ntax), epsilon_per = nimMatrix(0, ntax, length(id)), epsilon_col = nimMatrix(0, ntax, length(id)))
    if (FALSE %in% quo_is_null(var_guild)) {
      popdyn_parameters <- c(popdyn_parameters,"z_lambda_gui","z_mulambda_gui","GOR")
    }
    if (!is.null(modenv)) {
      popdyn_parameters <- popdyn_parameters[!popdyn_parameters %in% c("p.col","p.per","p.ext","p.col_id","p.per_id","p.ext_id")]
      if (FALSE %in% quo_is_null(var_envP) & FALSE %in% quo_is_null(var_envC)) {
        varenv <- colnames(select(df, !!var_envP))
        nvar <- length(varenv)
        varcol <- c(varenv[varenv %in% colnames(select(df, !!var_envC))],colnames(select(df, !!var_envC))[!colnames(select(df, !!var_envC)) %in% varenv])
        nr <- which(colnames(df) %in% varcol)
        ncol <- length(nr)
        col <- sapply(nr, function(i) tapply(df[,i], list(df[,quo_name(var_id)],df[,quo_name(var_tmp)]), unique),simplify = "array")
        popdyn_inits <- list(alpha_per =  rep(0, ntax), beta_per = nimMatrix(0, ntax, nvar),alpha_col =  rep(0, ntax), beta_col = nimMatrix(0, ntax, nvar)) 
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
int_checkfunction <- function(df, vars_in_df, var_id, env_tmp, vars_env, vars_na, vars_numeric, vars_duplicate, var_tmp, timestep, period, vars_pas=NULL, call = caller_env()) {
  if (FALSE %in% is.element(unlist(vars_in_df[!vars_in_df %in% "NULL"]),colnames(df))) {
    abort("must use existing variables", call = call)
  }
  df <- df[,which(colnames(df) %in% unlist(vars_in_df))] %>% distinct()
  #-----------------------------------------------------------------------------
  if (isFALSE(env_tmp)) {
    df_env <- df[,which(colnames(df) %in% c(quo_name(var_id), unlist(vars_env[!vars_env %in% NULL])))] %>% na.omit() %>% distinct()
    df[,which(colnames(df) %in% unlist(vars_env[!vars_env %in% NULL]))] <- NULL
    df <- left_join(df, df_env)
  }
  #-----------------------------------------------------------------------------
  if (TRUE %in% is.na(df[,which(colnames(df) %in% unlist(vars_na[!vars_na %in% NULL]))])) {
    vars_name <- paste0(names(which(colSums(is.na(df[,which(colnames(df) %in% unlist(vars_na[!vars_na %in% NULL]))])) > 0)), collapse=", ")
    abort(paste0("NA or NaN are not allowed in '", vars_name, "'", collapse=" "), call = call)
  } 
  #-----------------------------------------------------------------------------
  if (FALSE %in% lapply(df[,which(colnames(df) %in% unlist(vars_numeric[!vars_numeric %in% NULL]))], is.numeric)) { 
    vars_name <- vars_numeric[!vars_numeric %in% NULL][which(unlist(lapply(df[,which(colnames(df) %in% unlist(vars_numeric[!vars_numeric %in% NULL]))], is.numeric)) %in% FALSE)]
    abort(paste0("'", vars_name, "' must be numeric", collapse = " "), call = call)
  }
  #-----------------------------------------------------------------------------
  if (TRUE %in% duplicated(df[,which(colnames(df) %in% vars_duplicate)])) {
    abort("must not have duplicate entries with multiple matches", call = call)
  }
  #-----------------------------------------------------------------------------
  if (max(diff(sort(unique(pull(df, !!var_tmp))))) > timestep) {
    abort(paste0("missing values are not allowed in '", quo_name(var_tmp), "'", collapse = " "), call = call)
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
    dpas <- reframe(df, across(vars_pas[[1]], ~length(.x) == max(.x), .names = "pas"), .by = unlist(vars_pas[2:4]))
    if (FALSE %in% d$pas) {
      abord("missing passes in 'var_pas'", call = call)
    }
  }
  #-----------------------------------------------------------------------------
  return(df)
}
################################################################################
int_transformsummary <- function(mcmc_summary, datamodel, var_id, var_tmp, var_tax, period) {
  #-----------------------------------------------------------------------------
  colnames(mcmc_summary) <- c("mean", "sd", "Q2.5", "Q50", "Q97.5", "Rhat", "n.eff")
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
  mcmc_summary$paratax <- ifelse(mcmc_summary$parameter %in% c("z","p.per","p.col","p.ext","p.per_id","p.col_id","p.ext_id","alpha","beta","alpha_per","beta_per","alpha_col","beta_col","z_lambda","z_mulambda","turnover","OR","N","N_lambda_id","N_mulambda_id","N_lambda","N_mulambda","N_PGR","B","B_lambda_id","B_mulambda_id","B_lambda","B_mulambda","B_PGR","alpha_N","beta_N","alpha_B","beta_B"), tax[mcmc_summary$p1], NA)
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
  # get names and of environmental covariates 
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
  # set list of subscripts
  subscript <- list()
  subname <- vector("character",0)
  sub_tax <- select(mcmc_summary, paratax, p1) %>% na.omit() %>% distinct() %>% arrange(p1) %>% set_colnames(c("","subscript")) %>% set_rownames(NULL)
  if (nrow(sub_tax) > 0) {
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
  #-----------------------------------------------------------------------------
  mcmc_summary <- select(mcmc_summary, -c(para, start, end, pos, p1, p2, p3, paratax, paraid, paratime))
  return(list(mcmc_summary = mcmc_summary, subscript = subscript))
}
################################################################################
################################################################################
################################################################################
################################################################################
