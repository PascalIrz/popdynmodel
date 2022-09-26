#' Writing hierarchical Bayesian population dynamics models
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
#' model <- writemodel_removal(type="A",
#' period=NA,
#' para.nat=FALSE,
#' para.bas=TRUE,
#' para.sta=FALSE)
#' }
#' \dontrun{
#' stat <- jags(data = varlist,
#' parameters.to.save=c("R.bas"),
#' model=textConnection(model))
#' write(model,"model.txt")
#' stat <- jags(data = varlist,
#' parameters.to.save=c("R.bas"),
#' model="model.txt")
#' }
writemodel_removal <- function(type,period,para.nat,para.bas,para.sta) {
	#-------------------------------------------------------------------------------------------------------------
	## Site-occupancy dynamics model
	#-------------------------------------------------------------------------------------------------------------
	occupancy <- "
	## Site occupancy dynamics
	for (s in 1:nsta) {
		p.per[s] ~ dunif(0,1)
		p.col[s] ~ dunif(0,1)
		for (p in 1:npro[s]) {
			p.cap.sp[pro[s,p],s] ~ dunif(0,1)
		}
		for (t in start[s]:ntime) {
 			z[t,s] ~ dbern(p[t,s])
			p[t,s] <- z[t-1,s] * p.per[s] + (1 - z[t-1,s]) * p.col[s]
			for (n in 1:nsample[t,s]) {
				x[t,s,n] ~ dbern(z[t,s] * p.cap.sp[pro.time[t,s],s])
			}
		}
	}"
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
	## Modelling distribution changes rates at regional scale
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
	for (s in 1:nsta) {
		N[start[s]-1,s] ~ dpois(yini[s])
		mu.r[s] ~ dunif(0,10)
		tau.r[s] ~ dgamma(0.01,0.01)
		for (p in 1:npro[s]) {
			p.cap[pro[s,p],s] ~ dunif(0,1)
		}
		for (t in start[s]:ntime) {
			C[t,s] ~ dpois(yc[t,s])
			r[t,s] ~ dlnorm(mu.r[s],tau.r[s])
			N.rel[t,s] <- z[t,s] * (z[t-1,s] * r[t,s] * N[t-1,s] * (S[t,s] / S[t-1,s]) + (1 - z[t-1,s]) * C[t,s])
			N[t,s] <- round(N.rel[t,s])
			r.sta[t,s] <- z[t,s] * z[t-1,s] * (N[t,s] / max(1,N[t-1,s])) * (S[t-1,s] / S[t,s]) + (1 - z[t,s] * z[t-1,s])
			nr.sta[t,s] <- lap * z[t,s] * z[t-1,s]
		}
		for (t in (start[s]-1):ntime) {
			y[t,s,1] ~ dbin(p.cap[pro.time[t,s],s],N[t,s])
			Ny[t,s,1] <- N[t,s] - y[t,s,1]
		}
		for (t in 1:ntime.n[s]) {
		for (n in 2:nsample[time.n[s,t],s]) {
			y[time.n[s,t],s,n] ~ dbin(p.cap[pro.time[time.n[s,t],s],s],Ny[time.n[s,t],s,n-1])
			Ny[time.n[s,t],s,n] <- Ny[time.n[s,t],s,n-1] - y[time.n[s,t],s,n]
		}}
		R.sta[s] <- pow(prod(r.sta[start.set:ntime,s]),(1/max(1,sum(nr.sta[start.set:ntime,s]))))
		nb.r.sta[s] <- sum(nr.sta[start.set:ntime,s]) / lap
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
	Nn[start.set-1] <- sum(N[start.set-1,1:nsta])
	for (t in start.set:ntime) {
	zn[t] <- max(z[t,1:nsta])
	Nn[t] <- sum(N[t,1:nsta])
	r.nat[t] <- zn[t] * zn[t-1] * Nn[t] / max(1,Nn[t-1]) + ( 1 - zn[t] * zn[t-1])
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
	Nb[start.set-1,b] <- sum(N[start.set-1,sta.bas[b,1:nsta.bas[b]]])
	for (t in start.set:ntime) {
	zb[t,b] <- max(z[t,sta.bas[b,1:nsta.bas[b]]])
	Nb[t,b] <- sum(N[t,sta.bas[b,1:nsta.bas[b]]])
	r.bas[t,b] <- zb[t,b] * zb[t-1,b] * Nb[t,b] / max(1,Nb[t-1,b]) + (1 - zb[t,b] * zb[t-1,b])
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
	for (s in 1:nsta) {
	B[start[s]-1,s] ~ dpois(wini[s])
	mu.rw[s] ~ dunif(0,10)
	tau.rw[s] ~ dgamma(0.01,0.01)
	for (p in 1:npro[s]) {
		tau.cap[pro[s,p],s] ~ dgamma(0.01,0.01)
		var.cap[pro[s,p],s] <- 1 / tau.cap[pro[s,p],s]
	}
	for (t in start[s]:ntime) {
		Cb[t,s] ~ dpois(wc[t,s])
		rw[t,s] ~ dlnorm(mu.rw[s],tau.rw[s])
		B[t,s] <- z[t,s] * (z[t-1,s] * rw[t,s] * B[t-1,s] * (S[t,s] / S[t-1,s]) + (1 - z[t-1,s]) * Cb[t,s])
		rw.sta[t,s] <- z[t,s] * z[t-1,s] * (B[t,s] / max(1,B[t-1,s])) * (S[t-1,s] / S[t,s]) + (1 - z[t,s] * z[t-1,s])
		nrw.sta[t,s] <- lap * z[t,s] * z[t-1,s]
	}
	for (t in (start[s]-1):ntime) {
		w[t,s,1] ~ dnorm(B[t,s],tau.cap[pro.time[t,s],s]) T(0,)
		Bw[t,s,1] <- B[t,s] - w[t,s,1]
	}
	for (t in 1:ntime.n[s]) {
	for (n in 2:nsample[time.n[s,t],s]) {
		w[time.n[s,t],s,n] ~ dnorm(Bw[time.n[s,t],s,n-1],tau.cap[pro.time[time.n[s,t],s],s]) T(0,)
		Bw[time.n[s,t],s,n] <- Bw[time.n[s,t],s,n-1] - w[time.n[s,t],s,n]
	}}
	Rw.sta[s] <- pow(prod(rw.sta[start.set:ntime,s]),(1/max(1,sum(nrw.sta[start.set:ntime,s]))))
	nb.rw.sta[s] <- sum(nrw.sta[start.set:ntime,s]) / lap
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
	zn[start.set-1] <- max(z[start.set-1,1:nsta])
	Bn[start.set-1] <- sum(B[start.set-1,1:nsta])
	for (t in start.set:ntime) {
	zn[t] <- max(z[t,1:nsta])
	Bn[t] <- sum(B[t,1:nsta])
	rw.nat[t] <- zn[t] * zn[t-1] * Bn[t] / max(1,Bn[t-1]) + ( 1 - zn[t] * zn[t-1])
	nrw.nat[t] <- lap * zn[t] * zn[t-1]
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
	zb[start.set-1,b] <- max(z[start.set-1,sta.bas[b,1:nsta.bas[b]]])
	Bb[start.set-1,b] <- sum(B[start.set-1,sta.bas[b,1:nsta.bas[b]]])
	for (t in start.set:ntime) {
	zb[t,b] <- max(z[t,sta.bas[b,1:nsta.bas[b]]])
	Bb[t,b] <- sum(B[t,sta.bas[b,1:nsta.bas[b]]])
	rw.bas[t,b] <- zb[t,b] * zb[t-1,b] * Bb[t,b] / max(1,Bb[t-1,b]) + (1 - zb[t,b] * zb[t-1,b])
	nrw.bas[t,b] <- lap * zb[t,b] * zb[t-1,b]
	}
	Rw.bas[b] <- pow(prod(rw.bas[start.set:ntime,b]),(1/max(1,sum(nrw.bas[start.set:ntime,b]))))
	}"
	#-------------------------------------------------------------------------------------------------------------
	b.growth.bas.period <- "
	## Intermediate biomass growth rates at river-basins level
	for (b in 1:nbas) {
	for (i in 1:nperiod) {
	Rw.bas.period[period.end[i],b] <- pow(prod(rw.bas[period.start[i]:period.end[i],b]),(1/max(1,sum(nrw.bas[period.start[i]:period.end[i],b]))))
	}}"
	#-------------------------------------------------------------------------------------------------------------
	## Modelling growth rates in abundance and biomass at national scale
	#-------------------------------------------------------------------------------------------------------------
	ab.growth.nat <- "
	## Growth rates at national level
	zn[start.set-1] <- max(z[start.set-1,1:nsta])
	Nn[start.set-1] <- sum(N[start.set-1,1:nsta])
	Bn[start.set-1] <- sum(B[start.set-1,1:nsta])
	for (t in start.set:ntime) {
	zn[t] <- max(z[t,1:nsta])
	Nn[t] <- sum(N[t,1:nsta])
	r.nat[t] <- zn[t] * zn[t-1] * Nn[t] / max(1,Nn[t-1]) + ( 1 - zn[t] * zn[t-1])
	nr.nat[t] <- lap * zn[t] * zn[t-1]
	Bn[t] <- sum(B[t,1:nsta])
	rw.nat[t] <- zn[t] * zn[t-1] * Bn[t] / max(1,Bn[t-1]) + ( 1 - zn[t] * zn[t-1])
	}
	R.nat <- pow(prod(r.nat[start.set:ntime]),(1/max(1,sum(nr.nat[start.set:ntime]))))
	Rw.nat <- pow(prod(rw.nat[start.set:ntime]),(1/max(1,sum(nr.nat[start.set:ntime])))) "
	#-------------------------------------------------------------------------------------------------------------
	ab.growth.nat.period <- "
	## Intermediate growth rates at national level
	for (i in 1:nperiod) {
	R.nat.period[period.end[i]] <- pow(prod(r.nat[period.start[i]:period.end[i]]),(1/max(1,sum(nr.nat[period.start[i]:period.end[i]]))))
	Rw.nat.period[period.end[i]] <- pow(prod(rw.nat[period.start[i]:period.end[i]]),(1/max(1,sum(nr.nat[period.start[i]:period.end[i]]))))
	}"
	#-------------------------------------------------------------------------------------------------------------
	## Modelling the growth rates in abundance and biomass at regional scale
	#-------------------------------------------------------------------------------------------------------------
	ab.growth.bas <- "
	## Growth rates at regional level
	for (b in 1:nbas) {
	zb[start.set-1,b] <- max(z[start.set-1,sta.bas[b,1:nsta.bas[b]]])
	Nb[start.set-1,b] <- sum(N[start.set-1,sta.bas[b,1:nsta.bas[b]]])
	Bb[start.set-1,b] <- sum(B[start.set-1,sta.bas[b,1:nsta.bas[b]]])
	for (t in start.set:ntime) {
	zb[t,b] <- max(z[t,sta.bas[b,1:nsta.bas[b]]])
	Nb[t,b] <- sum(N[t,sta.bas[b,1:nsta.bas[b]]])
	r.bas[t,b] <- zb[t,b] * zb[t-1,b] * Nb[t,b] / max(1,Nb[t-1,b]) + (1 - zb[t,b] * zb[t-1,b])
	nr.bas[t,b] <- lap * zb[t,b] * zb[t-1,b]
	Bb[t,b] <- sum(B[t,sta.bas[b,1:nsta.bas[b]]])
	rw.bas[t,b] <- zb[t,b] * zb[t-1,b] * Bb[t,b] / max(1,Bb[t-1,b]) + (1 - zb[t,b] * zb[t-1,b])
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
	if (type == "O") {
		mod <- c(occupancy)
		if (para.nat == TRUE) { mod <- c(mod,occupancy.rate.nat) }
		if (para.bas == TRUE) { mod <- c(mod,occupancy.rate.bas) }
		if (!is.na(period)) {
			if (para.nat == TRUE) { mod <- c(mod,occupancy.rate.nat.period) }
			if (para.bas == TRUE) { mod <- c(mod,occupancy.rate.bas.period) }
		}
		gw.model <- c("model {",mod,"}")
		#-------------------------------------------------------------------------------------------------------------
	} else if (type == "A") {
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
	}	else if (type == "AB") {
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
	}	else if (type == "OA") {
		mod <- c(occupancy,a.dynamic)
		if (para.nat == TRUE) { mod <- c(mod,a.growth.nat,occupancy.rate.nat) }
		if (para.bas == TRUE) { mod <- c(mod,a.growth.bas,occupancy.rate.bas) }
		if (!is.na(period)) {
			if (para.nat == TRUE) { mod <- c(mod,a.growth.nat.period,occupancy.rate.nat.period) }
			if (para.bas == TRUE) { mod <- c(mod,a.growth.bas.period,occupancy.rate.bas.period) }
			if (para.sta == TRUE) { mod <- c(mod,a.growth.sta.period) }
		}
		gw.model <- c("model {",mod,"}")
		#-------------------------------------------------------------------------------------------------------------
	}	else if (type == "OB") {
		mod <- c(occupancy,b.dynamic)
		if (para.nat == TRUE) { mod <- c(mod,b.growth.nat,occupancy.rate.nat) }
		if (para.bas == TRUE) { mod <- c(mod,b.growth.bas,occupancy.rate.bas) }
		if (!is.na(period)) {
			if (para.nat == TRUE) { mod <- c(mod,b.growth.nat.period,occupancy.rate.nat.period) }
			if (para.bas == TRUE) { mod <- c(mod,b.growth.bas.period,occupancy.rate.bas.period) }
			if (para.sta == TRUE) { mod <- c(mod,b.growth.sta.period) }
		}
		gw.model <- c("model {",mod,"}")
		#-------------------------------------------------------------------------------------------------------------
	}	else {
		mod <- c(occupancy,a.dynamic,b.dynamic)
		if (para.nat == TRUE) { mod <- c(mod,ab.growth.nat,occupancy.rate.nat) }
		if (para.bas == TRUE) { mod <- c(mod,ab.growth.bas,occupancy.rate.bas) }
		if (!is.na(period)) {
			if (para.nat == TRUE) { mod <- c(mod,ab.growth.nat.period,occupancy.rate.nat.period) }
			if (para.bas == TRUE) { mod <- c(mod,ab.growth.bas.period,occupancy.rate.bas.period) }
			if (para.sta == TRUE) { mod <- c(mod,ab.growth.sta.period) }
		}
		gw.model <- c("model {",mod,"}")
	}
	#-------------------------------------------------------------------------------------------------------------
return(gw.model) }
