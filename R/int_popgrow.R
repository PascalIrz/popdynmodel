#' Title
#'
#' @param code
#' @param modenvG
#' @param var_cnt
#' @param var_wei
#' @param var_guild
#'
#' @return
#' @export
#'
#' @examples
int_popgrow <- function(code, modenvG, var_cnt, var_wei, var_guild) {
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
