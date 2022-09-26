#' Writing an unbiased hierarchical Bayesian population dynamics models
#'
#' @param type Character
#' @param period Numeric
#' @param para.nat,para.bas,para.sta Logical
#'
#' @return gw.model Character vector
#' @export
#'
#' @examples
#' \donttest{
#' model <- writemodel_unbiased(type="A",
#' period=NA,
#' para.nat=TRUE,
#' para.bas=FALSE,
#' para.sta=FALSE)
#' }
#' \dontrun{
#' stat <- jags(data = varlist,
#' parameters.to.save=c("R.nat"),
#' model=textConnection(model))
#' write(model,"model.txt")
#' stat <- jags(data = varlist,
#' parameters.to.save=c("R.nat"),
#' model="model.txt")
#' }
writemodel_unbiased <- function(type,period,para.nat,para.bas,para.sta) {
	#-------------------------------------------------------------------------------------------------------------
	## Site-occupancy dynamics model
	#-------------------------------------------------------------------------------------------------------------
	occupancy <- "
	## Site occupancy dynamics
	for (s in 1:nsta) {
	p.per[s] ~ dunif(0,1)
	p.col[s] ~ dunif(0,1)
	for (t in start.1[s]:ntime) {
	z[t,s] ~ dbern(p[t,s])
	p[t,s] <- z[t-1,s] * p.per[s] + (1 - z[t-1,s]) * p.col[s]
	}}"
	#-------------------------------------------------------------------------------------------------------------
	## Modelling distribution changes rates at national scale
	#-------------------------------------------------------------------------------------------------------------
	occupancy.rate.nat <- "
	## Distribution change rates at national level
	zn.o[start.set-1] <- max(z[start.set-1,1:nsta])
	nb.sta[start.set-1] <- sum(z[start.set-1,1:nsta])
	for (t in start.set:ntime) {
	zn.o[t] <- max(z[t,1:nsta])
	nb.sta[t] <- sum(z[t,1:nsta])
	ro.nat[t] <- zn.o[t] * zn.o[t-1] * (nb.sta[t] / max(1,nb.sta[t-1])) + (1 - zn.o[t] * zn.o[t-1])
	no.nat[t] <- lap * zn.o[t] * zn.o[t-1]
	}
	Ro.nat <- pow(prod(ro.nat[start.set:ntime]),(1/max(1,sum(no.nat[start.set:ntime])))) "
	#-------------------------------------------------------------------------------------------------------------
	occupancy.rate.nat.period <- "
	## Intermediate distribution change rates at national level
	for (i in 1:nperiod) {
	Ro.nat.period[period.end[i]] <- pow(prod(ro.nat[period.start[i]:period.end[i]]),(1/max(1,sum(no.nat[period.start[i]:period.end[i]]))))
	}"
	#-------------------------------------------------------------------------------------------------------------
	## Modelling the distribution changes rates at regional scale
	#-------------------------------------------------------------------------------------------------------------
	occupancy.rate.bas <- "
	## Distribution change rates at regional level
	for (b in 1:nbas) {
	zb.o[start.set-1,b] <- max(z[start.set-1,sta.bas[b,1:nsta.bas[b]]])
	nb.sta.bas[start.set-1,b] <- sum(z[start.set-1,sta.bas[b,1:nsta.bas[b]]])
	for (t in start.set:ntime) {
	zb.o[t,b] <- max(z[t,sta.bas[b,1:nsta.bas[b]]])
	nb.sta.bas[t,b] <- sum(z[t,sta.bas[b,1:nsta.bas[b]]])
	ro.bas[t,b] <- zb.o[t,b] * zb.o[t-1,b] * (nb.sta.bas[t,b] / max(1,nb.sta.bas[t-1,b])) + (1 - zb.o[t,b] * zb.o[t-1,b])
	no.bas[t,b] <- lap * zb.o[t,b] * zb.o[t-1,b]
	}
	Ro.bas[b] <- pow(prod(ro.bas[start.set:ntime,b]),(1/max(1,sum(no.bas[start.set:ntime,b]))))
	}"
	#-------------------------------------------------------------------------------------------------------------
	occupancy.rate.bas.period <- "
	## Intermediate distribution change rates at regional level
	for (b in 1:nbas) {
	for (i in 1:nperiod) {
	Ro.bas.period[period.end[i],b] <- pow(prod(ro.bas[period.start[i]:period.end[i],b]),(1/max(1,sum(no.bas[period.start[i]:period.end[i],b]))))
	}}"
	#-------------------------------------------------------------------------------------------------------------
	## Abundance dynamics model
	#-------------------------------------------------------------------------------------------------------------
	a.dynamic <- "
	## Population dynamics from abundance time-series
	for (s in 1:nsta.1) {
	mu[sta.1[s]] ~ dunif(0,100)
	tau[sta.1[s]] ~ dgamma(0.01,0.01)
	y[start.1[sta.1[s]]-1,sta.1[s]] <- z[start.1[sta.1[s]]-1,sta.1[s]] * y.trans[start.1[sta.1[s]]-1,sta.1[s]]

	for (t in start.1[sta.1[s]]:ntime) {
	Y[t,sta.1[s]] <-  log(max(1, z[t,sta.1[s]] * (z[t-1,sta.1[s]] * mu[sta.1[s]] * y.trans[t-1,sta.1[s]]  * (S[t,sta.1[s]] / S[t-1,sta.1[s]]) + (1 - z[t-1,sta.1[s]]) * yc[t,sta.1[s]])))
	y.trans[t,sta.1[s]] ~ dlnorm(Y[t,sta.1[s]],tau[sta.1[s]])
	y[t,sta.1[s]] <- z[t,sta.1[s]] * y.trans[t,sta.1[s]]

	r.rel[t,sta.1[s]] <- (y[t,sta.1[s]] / max(1,y[t-1,sta.1[s]])) * (S[t-1,sta.1[s]] / S[t,sta.1[s]])
	r.sta[t,sta.1[s]] <- z[t,sta.1[s]] * z[t-1,sta.1[s]] * r.rel[t,sta.1[s]] + (1 - z[t,sta.1[s]] * z[t-1,sta.1[s]])
	nr.sta[t,sta.1[s]] <- lap * z[t,sta.1[s]] * z[t-1,sta.1[s]]
	r.cal[t,sta.1[s]] <- z[t,sta.1[s]] * z[t,sta.1[s]] * r.rel[t,sta.1[s]] * y[t-1,sta.1[s]]
	nr.cal[t,sta.1[s]] <- z[t,sta.1[s]] * z[t,sta.1[s]] * y[t-1,sta.1[s]]
	}
	R.sta[sta.1[s]] <- pow(prod(r.sta[start.set:ntime,sta.1[s]]),(1/max(1,sum(nr.sta[start.set:ntime,sta.1[s]]))))
	nb.r.sta[sta.1[s]] <- sum(nr.sta[start.set:ntime,sta.1[s]]) / lap
	}

	for (s in 1:nsta.2) {
	mu[sta.2[s]] ~ dunif(0,100)
	tau[sta.2[s]] ~ dgamma(0.01,0.01)
	y[start.1[sta.2[s]]-1,sta.2[s]] <- z[start.1[sta.2[s]]-1,sta.2[s]] * y.trans[start.1[sta.2[s]]-1,sta.2[s]]
	y[start.2[sta.2[s]]-1,sta.2[s]] <- z[start.2[sta.2[s]]-1,sta.2[s]] * y.trans[start.2[sta.2[s]]-1,sta.2[s]]
	r.sta[start.2[sta.2[s]]-1,sta.2[s]] <- 1
	r.cal[start.2[sta.2[s]]-1,sta.2[s]] <- 0
	nr.sta[start.2[sta.2[s]]-1,sta.2[s]] <- 0
	nr.cal[start.2[sta.2[s]]-1,sta.2[s]] <- 0

	for (t in start.1[sta.2[s]]:end.1[sta.2[s]]) {
	Y[t,sta.2[s]] <-  log(max(1, z[t,sta.2[s]] * (z[t-1,sta.2[s]] * mu[sta.2[s]] * y.trans[t-1,sta.2[s]]  * (S[t,sta.2[s]] / S[t-1,sta.2[s]]) + (1 - z[t-1,sta.2[s]]) * yc[t,sta.2[s]])))
	y.trans[t,sta.2[s]] ~ dlnorm(Y[t,sta.2[s]],tau[sta.2[s]])
	y[t,sta.2[s]] <- z[t,sta.2[s]] * y.trans[t,sta.2[s]]

	r.rel[t,sta.2[s]] <- (y[t,sta.2[s]] / max(1,y[t-1,sta.2[s]])) * (S[t-1,sta.2[s]] / S[t,sta.2[s]])
	r.sta[t,sta.2[s]] <- z[t,sta.2[s]] * z[t-1,sta.2[s]] * r.rel[t,sta.2[s]] + (1 - z[t,sta.2[s]] * z[t-1,sta.2[s]])
	nr.sta[t,sta.2[s]] <- lap * z[t,sta.2[s]] * z[t-1,sta.2[s]]
	r.cal[t,sta.2[s]] <- z[t,sta.2[s]] * z[t,sta.2[s]] * r.rel[t,sta.2[s]] * y[t-1,sta.2[s]]
	nr.cal[t,sta.2[s]] <- z[t,sta.2[s]] * z[t,sta.2[s]] * y[t-1,sta.2[s]]
	}

	for (t in start.2[sta.2[s]]:ntime) {
	Y[t,sta.2[s]] <-  log(max(1, z[t,sta.2[s]] * (z[t-1,sta.2[s]] * mu[sta.2[s]] * y.trans[t-1,sta.2[s]]  * (S[t,sta.2[s]] / S[t-1,sta.2[s]]) + (1 - z[t-1,sta.2[s]]) * yc[t,sta.2[s]])))
	y.trans[t,sta.2[s]] ~ dlnorm(Y[t,sta.2[s]],tau[sta.2[s]])
	y[t,sta.2[s]] <- z[t,sta.2[s]] * y.trans[t,sta.2[s]]

	r.rel[t,sta.2[s]] <- (y[t,sta.2[s]] / max(1,y[t-1,sta.2[s]])) * (S[t-1,sta.2[s]] / S[t,sta.2[s]])
	r.sta[t,sta.2[s]] <- z[t,sta.2[s]] * z[t-1,sta.2[s]] * r.rel[t,sta.2[s]] + (1 - z[t,sta.2[s]] * z[t-1,sta.2[s]])
	nr.sta[t,sta.2[s]] <- lap * z[t,sta.2[s]] * z[t-1,sta.2[s]]
	r.cal[t,sta.2[s]] <- z[t,sta.2[s]] * z[t,sta.2[s]] * r.rel[t,sta.2[s]] * y[t-1,sta.2[s]]
	nr.cal[t,sta.2[s]] <- z[t,sta.2[s]] * z[t,sta.2[s]] * y[t-1,sta.2[s]]
	}

	R.sta[sta.2[s]] <- pow(prod(r.sta[start.set:ntime,sta.2[s]]),(1/max(1,sum(nr.sta[start.set:ntime,sta.2[s]]))))
	nb.r.sta[sta.2[s]] <- sum(nr.sta[start.set:ntime,sta.2[s]]) / lap
	}"
	#-------------------------------------------------------------------------------------------------------------
	a.growth.sta.period <- "
	## Intermediate abundance growth rates at site level
	for (s in 1:nsta) {
	for (i in 1:nperiod) {
	R.sta.period[period.end[i],s] <- pow(prod(r.sta[period.start[i]:period.end[i],s]),(1/max(1,sum(nr.sta[period.start[i]:period.end[i],s]))))
	}}"
	#-------------------------------------------------------------------------------------------------------------
	## Modelling growth rates in abundance at national scale
	#-------------------------------------------------------------------------------------------------------------
	a.growth.nat <- "
	## Abundance growth rates at national level
	zn[start.set-1] <- max(z[start.set-1,1:nsta])
	for (t in start.set:ntime) {
	zn[t] <- max(z[t,1:nsta])
	r.nat[t] <- zn[t] * zn[t-1] * sum(r.cal[t,1:nsta]) / max(1,sum(nr.cal[t,1:nsta])) + ( 1 - zn[t] * zn[t-1])
	nr.nat[t] <- lap * zn[t] * zn[t-1]
	}
	R.nat <- pow(prod(r.nat[start.set:ntime]),(1/max(1,sum(nr.nat[start.set:ntime])))) "
	#-------------------------------------------------------------------------------------------------------------
	a.growth.nat.period <- "
	## Intermediate abundance growth rates at national level
	for (i in 1:nperiod) {
	R.nat.period[period.end[i]] <- pow(prod(r.nat[period.start[i]:period.end[i]]),(1/max(1,sum(nr.nat[period.start[i]:period.end[i]]))))
	}"
	#-------------------------------------------------------------------------------------------------------------
	## Modelling growth rates in abundance at regional scale
	#-------------------------------------------------------------------------------------------------------------
	a.growth.bas <- "
	## Abundance growth rates at regional level
	for (b in 1:nbas) {
	zb[start.set-1,b] <- max(z[start.set-1,sta.bas[b,1:nsta.bas[b]]])
	for (t in start.set:ntime) {
	zb[t,b] <- max(z[t,sta.bas[b,1:nsta.bas[b]]])
	r.bas[t,b] <- zb[t,b] * zb[t-1,b] * sum(r.cal[t,sta.bas[b,1:nsta.bas[b]]]) / max(1,sum(nr.cal[t,sta.bas[b,1:nsta.bas[b]]])) + (1 - zb[t,b] * zb[t-1,b])
	nr.bas[t,b] <- lap * zb[t,b] * zb[t-1,b]
	}
	R.bas[b] <- pow(prod(r.bas[start.set:ntime,b]),(1/max(1,sum(nr.bas[start.set:ntime,b]))))
	}"
	#-------------------------------------------------------------------------------------------------------------
	a.growth.bas.period <- "
	## Intermediate abundance growth rates at regional level
	for (b in 1:nbas) {
	for (i in 1:nperiod) {
	R.bas.period[period.end[i],b] <- pow(prod(r.bas[period.start[i]:period.end[i],b]),(1/max(1,sum(nr.bas[period.start[i]:period.end[i],b]))))
	}}"
	#-------------------------------------------------------------------------------------------------------------
	## Biomass dynamics model
	#-------------------------------------------------------------------------------------------------------------
	b.dynamic <- "
	## Population dynamics from biomass time-series
	for (s in 1:nsta.1) {
	muw[sta.1[s]] ~ dunif(0,1000)
	tauw[sta.1[s]] ~ dgamma(0.01,0.01)
	w[start.1[sta.1[s]]-1,sta.1[s]] <- z[start.1[sta.1[s]]-1,sta.1[s]] * w.trans[start.1[sta.1[s]]-1,sta.1[s]]

	for (t in start.1[sta.1[s]]:ntime) {
	W[t,sta.1[s]] <-  log(max(1, z[t,sta.1[s]] * (z[t-1,sta.1[s]] * muw[sta.1[s]] * w.trans[t-1,sta.1[s]]  * (S[t,sta.1[s]] / S[t-1,sta.1[s]]) + (1 - z[t-1,sta.1[s]]) * wc[t,sta.1[s]])))
	w.trans[t,sta.1[s]] ~ dlnorm(W[t,sta.1[s]],tauw[sta.1[s]])
	w[t,sta.1[s]] <- z[t,sta.1[s]] * w.trans[t,sta.1[s]]

	rw.rel[t,sta.1[s]] <- (w[t,sta.1[s]] / max(1,w[t-1,sta.1[s]])) * (S[t-1,sta.1[s]] / S[t,sta.1[s]])
	rw.sta[t,sta.1[s]] <- z[t,sta.1[s]] * z[t-1,sta.1[s]] * rw.rel[t,sta.1[s]] + (1 - z[t,sta.1[s]] * z[t-1,sta.1[s]])
	nrw.sta[t,sta.1[s]] <- lap * z[t,sta.1[s]] * z[t-1,sta.1[s]]
	rw.cal[t,sta.1[s]] <- z[t,sta.1[s]] * z[t,sta.1[s]] * rw.rel[t,sta.1[s]] * w[t-1,sta.1[s]]
	nrw.cal[t,sta.1[s]] <- z[t,sta.1[s]] * z[t,sta.1[s]] * w[t-1,sta.1[s]]
	}

	Rw.sta[sta.1[s]] <- pow(prod(rw.sta[start.set:ntime,sta.1[s]]),(1/max(1,sum(nrw.sta[start.set:ntime,sta.1[s]]))))
	nb.rw.sta[sta.1[s]] <- sum(nrw.sta[start.set:ntime,sta.1[s]]) / lap
	}

	for (s in 1:nsta.2) {
	muw[sta.2[s]] ~ dunif(0,1000)
	tauw[sta.2[s]] ~ dgamma(0.01,0.01)
	w[start.1[sta.2[s]]-1,sta.2[s]] <- z[start.1[sta.2[s]]-1,sta.2[s]] * w.trans[start.1[sta.2[s]]-1,sta.2[s]]
	w[start.2[sta.2[s]]-1,sta.2[s]] <- z[start.2[sta.2[s]]-1,sta.2[s]] * w.trans[start.2[sta.2[s]]-1,sta.2[s]]
	rw.sta[start.2[sta.2[s]]-1,sta.2[s]] <- 1
	rw.cal[start.2[sta.2[s]]-1,sta.2[s]] <- 0
	nrw.sta[start.2[sta.2[s]]-1,sta.2[s]] <- 0
	nrw.cal[start.2[sta.2[s]]-1,sta.2[s]] <- 0

	for (t in start.1[sta.2[s]]:end.1[sta.2[s]]) {
	W[t,sta.2[s]] <-  log(max(1, z[t,sta.2[s]] * (z[t-1,sta.2[s]] * muw[sta.2[s]] * w.trans[t-1,sta.2[s]]  * (S[t,sta.2[s]] / S[t-1,sta.2[s]]) + (1 - z[t-1,sta.2[s]]) * wc[t,sta.2[s]])))
	w.trans[t,sta.2[s]] ~ dlnorm(W[t,sta.2[s]],tauw[sta.2[s]])
	w[t,sta.2[s]] <- z[t,sta.2[s]] * w.trans[t,sta.2[s]]

	rw.rel[t,sta.2[s]] <- (w[t,sta.2[s]] / max(1,w[t-1,sta.2[s]])) * (S[t-1,sta.2[s]] / S[t,sta.2[s]])
	rw.sta[t,sta.2[s]] <- z[t,sta.2[s]] * z[t-1,sta.2[s]] * rw.rel[t,sta.2[s]] + (1 - z[t,sta.2[s]] * z[t-1,sta.2[s]])
	nrw.sta[t,sta.2[s]] <- lap * z[t,sta.2[s]] * z[t-1,sta.2[s]]
	rw.cal[t,sta.2[s]] <- z[t,sta.2[s]] * z[t,sta.2[s]] * rw.rel[t,sta.2[s]] * w[t-1,sta.2[s]]
	nrw.cal[t,sta.2[s]] <- z[t,sta.2[s]] * z[t,sta.2[s]] * w[t-1,sta.2[s]]
	}

	for (t in start.2[sta.2[s]]:ntime) {
	W[t,sta.2[s]] <-  log(max(1, z[t,sta.2[s]] * (z[t-1,sta.2[s]] * muw[sta.2[s]] * w.trans[t-1,sta.2[s]]  * (S[t,sta.2[s]] / S[t-1,sta.2[s]]) + (1 - z[t-1,sta.2[s]]) * wc[t,sta.2[s]])))
	w.trans[t,sta.2[s]] ~ dlnorm(W[t,sta.2[s]],tauw[sta.2[s]])
	w[t,sta.2[s]] <- z[t,sta.2[s]] * w.trans[t,sta.2[s]]

	rw.rel[t,sta.2[s]] <- (w[t,sta.2[s]] / max(1,w[t-1,sta.2[s]])) * (S[t-1,sta.2[s]] / S[t,sta.2[s]])
	rw.sta[t,sta.2[s]] <- z[t,sta.2[s]] * z[t-1,sta.2[s]] * rw.rel[t,sta.2[s]] + (1 - z[t,sta.2[s]] * z[t-1,sta.2[s]])
	nrw.sta[t,sta.2[s]] <- lap * z[t,sta.2[s]] * z[t-1,sta.2[s]]
	rw.cal[t,sta.2[s]] <- z[t,sta.2[s]] * z[t,sta.2[s]] * rw.rel[t,sta.2[s]] * w[t-1,sta.2[s]]
	nrw.cal[t,sta.2[s]] <- z[t,sta.2[s]] * z[t,sta.2[s]] * w[t-1,sta.2[s]]
	}

	Rw.sta[sta.2[s]] <- pow(prod(rw.sta[start.set:ntime,sta.2[s]]),(1/max(1,sum(nrw.sta[start.set:ntime,sta.2[s]]))))
	nb.rw.sta[sta.2[s]] <- sum(nrw.sta[start.set:ntime,sta.2[s]])
	}"
	#-------------------------------------------------------------------------------------------------------------
	b.growth.sta.period <- "
	## Intermediate abundance growth rates at site level
	for (s in 1:nsta) {
	for (i in 1:nperiod) {
	Rw.sta.period[period.end[i],s] <- pow(prod(rw.sta[period.start[i]:period.end[i],s]),(1/max(1,sum(nrw.sta[period.start[i]:period.end[i],s]))))
	}}"
	#-------------------------------------------------------------------------------------------------------------
	ab.growth.sta.period <- "
	## Intermediate abundance growth rates at site level
	for (s in 1:nsta) {
	for (i in 1:nperiod) {
	R.sta.period[period.end[i],s] <- pow(prod(r.sta[period.start[i]:period.end[i],s]),(1/max(1,sum(nr.sta[period.start[i]:period.end[i],s]))))
	Rw.sta.period[period.end[i],s] <- pow(prod(rw.sta[period.start[i]:period.end[i],s]),(1/max(1,sum(nrw.sta[period.start[i]:period.end[i],s]))))
	}}"
	#-------------------------------------------------------------------------------------------------------------
	## Modelling growth rates in biomass at national scale
	#-------------------------------------------------------------------------------------------------------------
	b.growth.nat <- "
	## Biomass growth rates at national level
	zn.w[start.set-1] <- max(z[start.set-1,1:nsta])
	for (t in start.set:ntime) {
	zn.w[t] <- max(z[t,1:nsta])
	rw.nat[t] <- zn.w[t] * zn.w[t-1] * sum(rw.cal[t,1:nsta]) / max(1,sum(nrw.cal[t,1:nsta])) + ( 1 - zn.w[t] * zn.w[t-1])
	nrw.nat[t] <- lap * zn.w[t] * zn.w[t-1]
	}
	Rw.nat <- pow(prod(rw.nat[start.set:ntime]),(1/max(1,sum(nrw.nat[start.set:ntime])))) "
	#-------------------------------------------------------------------------------------------------------------
	b.growth.nat.period <- "
	## Intermediate biomass growth rates at national level
	for (i in 1:nperiod) {
	Rw.nat.period[period.end[i]] <- pow(prod(rw.nat[period.start[i]:period.end[i]]),(1/max(1,sum(nrw.nat[period.start[i]:period.end[i]]))))
	}"
	#-------------------------------------------------------------------------------------------------------------
	## Modelling growth rates in biomass at regional scale
	#-------------------------------------------------------------------------------------------------------------
	b.growth.bas <- "
	## Biomass growth rates at regional level
	for (b in 1:nbas) {
	zb.w[start.set-1,b] <- max(z[start.set-1,sta.bas[b,1:nsta.bas[b]]])
	for (t in start.set:ntime) {
	zb.w[t,b] <- max(z[t,sta.bas[b,1:nsta.bas[b]]])
	rw.bas[t,b] <- zb.w[t,b] * zb.w[t-1,b] * sum(rw.cal[t,sta.bas[b,1:nsta.bas[b]]]) / max(1,sum(nrw.cal[t,sta.bas[b,1:nsta.bas[b]]])) + (1 - zb.w[t,b] * zb.w[t-1,b])
	nrw.bas[t,b] <- lap * zb.w[t,b] * zb.w[t-1,b]
	}
	Rw.bas[b] <- pow(prod(rw.bas[start.set:ntime,b]),(1/max(1,sum(nrw.bas[start.set:ntime,b]))))
	}"
	#-------------------------------------------------------------------------------------------------------------
	b.growth.bas.period <- "
	## Intermediate biomass growth rates at regional level
	for (b in 1:nbas) {
	for (i in 1:nperiod) {
	Rw.bas.period[period.end[i],b] <- pow(prod(rw.bas[period.start[i]:period.end[i],b]),(1/max(1,sum(nrw.bas[period.start[i]:period.end[i],b]))))
	}}"
	#-------------------------------------------------------------------------------------------------------------
	## Modelling growth rates in abundance and biomass at national scale
	#-------------------------------------------------------------------------------------------------------------
	ab.growth.nat <- "
	## Abundance growth rates at national level
	zn[start.set-1] <- max(z[start.set-1,1:nsta])
	for (t in start.set:ntime) {
	zn[t] <- max(z[t,1:nsta])
	r.nat[t] <- zn[t] * zn[t-1] * sum(r.cal[t,1:nsta]) / max(1,sum(nr.cal[t,1:nsta])) + ( 1 - zn[t] * zn[t-1])
	nr.nat[t] <- lap * zn[t] * zn[t-1]
	rw.nat[t] <- zn[t] * zn[t-1] * sum(rw.cal[t,1:nsta]) / max(1,sum(nrw.cal[t,1:nsta])) + ( 1 - zn[t] * zn[t-1])
	}
	R.nat <- pow(prod(r.nat[start.set:ntime]),(1/max(1,sum(nr.nat[start.set:ntime]))))
	Rw.nat <- pow(prod(rw.nat[start.set:ntime]),(1/max(1,sum(nr.nat[start.set:ntime])))) "
	#-------------------------------------------------------------------------------------------------------------
	ab.growth.nat.period <- "
	## Intermediate abundance growth rates at national level
	for (i in 1:nperiod) {
	R.nat.period[period.end[i]] <- pow(prod(r.nat[period.start[i]:period.end[i]]),(1/max(1,sum(nr.nat[period.start[i]:period.end[i]]))))
	Rw.nat.period[period.end[i]] <- pow(prod(rw.nat[period.start[i]:period.end[i]]),(1/max(1,sum(nr.nat[period.start[i]:period.end[i]]))))
	}"
	#-------------------------------------------------------------------------------------------------------------
	## Modelling growth rates in abundance and biomass at regional scale
	#-------------------------------------------------------------------------------------------------------------
	ab.growth.bas <- "
	## Abundance growth rates at regional level
	for (b in 1:nbas) {
	zb[start.set-1,b] <- max(z[start.set-1,sta.bas[b,1:nsta.bas[b]]])
	for (t in start.set:ntime) {
	zb[t,b] <- max(z[t,sta.bas[b,1:nsta.bas[b]]])
	r.bas[t,b] <- zb[t,b] * zb[t-1,b] * sum(r.cal[t,sta.bas[b,1:nsta.bas[b]]]) / max(1,sum(nr.cal[t,sta.bas[b,1:nsta.bas[b]]])) + (1 - zb[t,b] * zb[t-1,b])
	nr.bas[t,b] <- lap * zb[t,b] * zb[t-1,b]
	rw.bas[t,b] <- zb[t,b] * zb[t-1,b] * sum(rw.cal[t,sta.bas[b,1:nsta.bas[b]]]) / max(1,sum(nrw.cal[t,sta.bas[b,1:nsta.bas[b]]])) + (1 - zb[t,b] * zb[t-1,b])
	}
	R.bas[b] <- pow(prod(r.bas[start.set:ntime,b]),(1/max(1,sum(nr.bas[start.set:ntime,b]))))
	Rw.bas[b] <- pow(prod(rw.bas[start.set:ntime,b]),(1/max(1,sum(nr.bas[start.set:ntime,b]))))
	}"
	#-------------------------------------------------------------------------------------------------------------
	ab.growth.bas.period <- "
	## Intermediate abundance growth rates at regional level
	for (b in 1:nbas) {
	for (i in 1:nperiod) {
	R.bas.period[period.end[i],b] <- pow(prod(r.bas[period.start[i]:period.end[i],b]),(1/max(1,sum(nr.bas[period.start[i]:period.end[i],b]))))
	Rw.bas.period[period.end[i],b] <- pow(prod(rw.bas[period.start[i]:period.end[i],b]),(1/max(1,sum(nr.bas[period.start[i]:period.end[i],b]))))
	}}"
	#-------------------------------------------------------------------------------------------------------------
	## Selecting model components and writing model depending on "type", "para.nat" and "para.bas"
	#-------------------------------------------------------------------------------------------------------------
	if (type == "A") {
		mod <- c(occupancy,a.dynamic)
		if (para.nat == TRUE) { mod <- c(mod,a.growth.nat) }
		if (para.bas == TRUE) { mod <- c(mod,a.growth.bas) }
		if (!is.na(period)) {
			if (para.nat == TRUE) { mod <- c(mod,a.growth.nat.period) }
			if (para.bas == TRUE) { mod <- c(mod,a.growth.bas.period) }
			if (para.sta == TRUE) { mod <- c(mod,a.growth.sta.period) }
		}
		gw.model <- c("model {",mod,"}")
		#-------------------------------------------------------------------------------------------------------------
	}	else if (type == "B") {
		mod <- c(occupancy,b.dynamic)
		if (para.nat == TRUE) { mod <- c(mod,b.growth.nat) }
		if (para.bas == TRUE) { mod <- c(mod,b.growth.bas) }
		if (!is.na(period)) {
			if (para.nat == TRUE) { mod <- c(mod,b.growth.nat.period) }
			if (para.bas == TRUE) { mod <- c(mod,b.growth.bas.period) }
			if (para.sta == TRUE) { mod <- c(mod,b.growth.sta.period) }
		}
		gw.model <- c("model {",mod,"}")
		#-------------------------------------------------------------------------------------------------------------
	} else if (type == "AB") {
		mod <- c(occupancy,a.dynamic,b.dynamic)
		if (para.nat == TRUE) { mod <- c(mod,ab.growth.nat) }
		if (para.bas == TRUE) { mod <- c(mod,ab.growth.bas) }
		if (!is.na(period)) {
			if (para.nat == TRUE) { mod <- c(mod,ab.growth.nat.period) }
			if (para.bas == TRUE) { mod <- c(mod,ab.growth.bas.period) }
			if (para.sta == TRUE) { mod <- c(mod,ab.growth.sta.period) }
		}
		gw.model <- c("model {",mod,"}")
		#-------------------------------------------------------------------------------------------------------------
	} else if (type == "O") {
		mod <- occupancy
		if (para.nat == TRUE) { mod <- c(mod,occupancy.rate.nat) }
		if (para.bas == TRUE) { mod <- c(mod,occupancy.rate.bas) }
		if (!is.na(period)) {
			if (para.nat == TRUE) { mod <- c(mod,occupancy.rate.nat.period) }
			if (para.bas == TRUE) { mod <- c(mod,occupancy.rate.bas.period) }
		}
		gw.model <- c("model {",mod,"}")
		#-------------------------------------------------------------------------------------------------------------
	} else if (type == "OA") {
		mod <- c(occupancy,a.dynamic)
		if (para.nat == TRUE) { mod <- c(mod,occupancy.rate.nat,a.growth.nat) }
		if (para.bas == TRUE) { mod <- c(mod,occupancy.rate.bas,a.growth.bas) }
		if (!is.na(period)) {
			if (para.nat == TRUE) { mod <- c(mod,a.growth.nat.period,occupancy.rate.nat.period) }
			if (para.bas == TRUE) { mod <- c(mod,a.growth.bas.period,occupancy.rate.bas.period) }
			if (para.sta == TRUE) { mod <- c(mod,a.growth.sta.period) }
		}
		gw.model <- c("model {",mod,"}")
		#-------------------------------------------------------------------------------------------------------------
	}	else if (type == "OB") {
		mod <- c(occupancy,b.dynamic)
		if (para.nat == TRUE) { mod <- c(mod,occupancy.rate.nat,b.growth.nat) }
		if (para.bas == TRUE) { mod <- c(mod,occupancy.rate.bas,b.growth.bas) }
		if (!is.na(period)) {
			if (para.nat == TRUE) { mod <- c(mod,b.growth.nat.period,occupancy.rate.nat.period) }
			if (para.bas == TRUE) { mod <- c(mod,b.growth.bas.period,occupancy.rate.bas.period) }
			if (para.sta == TRUE) { mod <- c(mod,b.growth.sta.period) }
		}
		gw.model <- c("model {",mod,"}")
		#-------------------------------------------------------------------------------------------------------------
	} else {
		mod <- c(occupancy,a.dynamic,b.dynamic)
		if (para.nat == TRUE) { mod <- c(mod,occupancy.rate.nat,ab.growth.nat) }
		if (para.bas == TRUE) { mod <- c(mod,occupancy.rate.bas,ab.growth.bas) }
		if (!is.na(period)) {
			if (para.nat == TRUE) { mod <- c(mod,ab.growth.nat.period,occupancy.rate.nat.period) }
			if (para.bas == TRUE) { mod <- c(mod,ab.growth.bas.period,occupancy.rate.bas.period) }
			if (para.sta == TRUE) { mod <- c(mod,ab.growth.sta.period) }
		}
		gw.model <- c("model {",mod,"}")
	}
	#-------------------------------------------------------------------------------------------------------------
	return(gw.model) }
