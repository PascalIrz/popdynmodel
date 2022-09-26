#' Creating input date files for Bayesian population dynamics models
#'
#' @import stats
#'
#' @param data A data frame
#' @param time,start,time.interval,period Numeric
#' @param type A character
#' @param para.bas A logical
#'
#' @return data.jags A list
#' @return sta Numeric or character vector
#' @return bas Numeric or character vector
#' @export
#'
#' @examples
#' \donttest{
#' data(riverfish)
#' varlist <- convertselect_removal(riverfish,
#' time=sort(unique(riverfish$date)),
#' start=min(riverfish$date),
#' time.interval=1,
#' type="A",
#' period=NA,
#' para.bas=TRUE)
#' }
convertselect_removal <- function(data,time,start,time.interval,type,period,para.bas) {
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Coding used functions
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# Replaces the missing initial values of a matrix mat with the replacement value
	initial_value <- function(mat,start.time,replacement) {
		if (is.na(mat[start.time])) { mat[start.time] <- 0 }
	return(mat) }
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Creating y, w and z variables
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (type %in% c("A","O","OA")) {
		data$y <- data$w <- data$abundance
	} else if (type %in% c("B","OB")) {
		data$y <- data$w <- data$biomass
	} else {
		data$y <- data$abundance
		data$w <- data$biomass
		data$wy <- sapply(data$sta_id, function(x) mean(data$biomass[data$sta_id == x]/data$abundance[data$sta_id == x],na.rm=T))
	}
	data$z <- ifelse(is.na(data$y), NA, ifelse(data$y > 0, 1, 0))
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Defining the site identifier vector (sta) and the number of sites (nsta)
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	sta <- sort(unique(data$sta_id))
	nsta <- length(sta)
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Defining the number of time step (ntime) and the position of the starting date in vector time (start.set)
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	ntime <- length(time)
	start.set <- which(time %in% (start + time.interval))
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Designing two-dimensional variable matrices (time, site)
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# true occupancy
	z <- tapply(replace(data$z,which(is.na(data$z)),0),list(data$date,data$sta_id), max)
	z <- replace(z,which(z==0),NA)
	# abundance of colonizers
	yc <- tapply(data$y, list(data$date,data$sta_id), sum)
	yc <- replace(yc,which(is.na(yc)),1)
	# biomass of colonizers
	wc <- tapply(data$w, list(data$date,data$sta_id), sum)
	wc <- replace(wc,which(is.na(wc)),1)
	# number of samples
	nsample <- tapply(data$sample, list(data$date,data$sta_id), max)
	# surface
	if (type == "O") {
		S <- NA
	} else {
		S <- tapply(data$surface, list(data$date,data$sta_id), mean)
	}
	# protocol
	protocol <- sort(unique(data$pro_id))
	pro.time <- tapply(data$pro_id, list(data$date,data$sta_id), function(i) which(protocol %in% i))
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Defining the start of time series (start.1), initial occupancy, initial abundances (yini), initial biomass (wini)
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	time.1 <- as.vector(sapply(sta, function(i) which(time %in% min(data$date[(data$sta_id == i & !is.na(data$y))]))))
	start.1 <- sapply(1+time.1, function(i) min(start.set,i))
	# assign initial value of yini (1 if missing)
	yini <- sapply(1:nsta, function(i) sum(data$y[(data$sta_id %in% sta[i] & data$date %in% time[start.1[i]-1])]))
	yini <- replace(yini,which(is.na(yini)),1)
	# assign initial value of wini (using the abundance-biomass relationship or 1 if missing)
	wini <- sapply(1:nsta, function(i) sum(data$w[(data$sta_id %in% sta[i] & data$date %in% time[start.1[i]-1])]))
	if (type %in% c("AB","OAB")) {
		wini <- sapply(1:nsta, function(i) ifelse(is.na(wini[i]), round(yini[i] * unique(data$wy[data$sta_id==sta[i]])), wini[i]))
	}
	# assign initial values of z if missing
	z <- sapply(seq_len(ncol(z)), function(i) initial_value(z[,i],start.1[i]-1,0))
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Defining the protocol identifer matrix (pro) and the number of protocols (npro) at sites-based scale
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	npro <- sapply(1:nsta, function(i) length(unique(data$pro_id[(data$sta_id %in% sta[i] & data$date >= time[start.1[i]-1])])))
	if (TRUE %in% all(npro == 1)) {
		pro <- matrix(sapply(sta, function(i) which(protocol %in% unique(data$pro_id[data$sta_id == i]))), nrow = nsta, ncol = 1)
	} else {
		pro <- t(sapply(1:nsta, function(i) c(which(protocol %in% sort(unique(data$pro_id[(data$sta_id %in% sta[i] & data$date >= time[start.1[i]-1])]))),rep(NA,max(npro)-npro[i]))))
	}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Defining times with multiple samples for each site (time.n) and the corresponding number (ntime.n)
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	sample.by.site <- aggregate(data$sample[!is.na(data$y)],list("sta_id"=data$sta_id[!is.na(data$y)],"date"=data$date[!is.na(data$y)]),max)
	ntime.n <- sapply(1:nsta, function(i) length(sample.by.site$date[(sample.by.site$sta_id %in% sta[i] & sample.by.site$date >= time[start.1[i]] & sample.by.site$x > 1)]))
	time.n <- t(sapply(1:nsta, function(i) c(which(time %in% sample.by.site$date[(sample.by.site$sta_id %in% sta[i] & sample.by.site$date >= time[start.1[i]] & sample.by.site$x > 1)]),rep(NA,max(ntime.n)-ntime.n[i]))))
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Designing three-dimensional variable matrices (time, site, sample)
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# observed occupancy
	x <- tapply(data$z, list(data$date,data$sta_id,data$sample), unique)
	# observed abundance
	y <- tapply(data$y, list(data$date,data$sta_id,data$sample), unique)
	# observed biomass
	w <- tapply(data$w, list(data$date,data$sta_id,data$sample), unique)
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Selecting and converting the required variables
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	list.O <- list(ntime=ntime,nsta=nsta,start.set=start.set,start=start.1,lap=time.interval,z=z,x=x)
	list.A <- list(y=y,yc=yc,yini=yini,S=S)
	list.B <- list(w=w,wc=wc,wini=wini,S=S)
	list.AB <- list(y=y,yc=yc,yini=yini,w=w,wc=wc,wini=wini,S=S)
	list.remove <- list(nsample=nsample,time.n=time.n,ntime.n=ntime.n,npro=npro,pro=pro,pro.time=pro.time)
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (type == "O") {
		data.jags <- c(list.O,list.remove)
	} else if (type == "A" | type == "OA") {
		data.jags <- c(list.O,list.A,list.remove)
	}	else if (type == "B" | type == "OB") {
		data.jags <- c(list.O,list.B,list.remove)
	} else {
		data.jags <- c(list.O,list.AB,list.remove)
	}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## If basins-wide parameters are requested : defining bas, nbas, sta.bas and nsta.bas
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (para.bas == TRUE) {
		# region identifier vector and number of regions
		bas <- sort(unique(data$bas_id))
		nbas <- length(bas)
		# number of sites per region
		nsta.bas <- sapply(bas, function(i) length(unique(data$sta_id[data$bas_id == i])))
		# two-dimensional matrix with the positions of site for each region
		sta.bas <- t(sapply(1:nbas, function(i) c(which(sta %in% sort(unique(data$sta_id[data$bas_id == bas[i]]))),rep(NA,max(nsta.bas)-nsta.bas[i]))))
		data.jags <- c(data.jags,list(nbas=nbas,nsta.bas=nsta.bas,sta.bas=sta.bas))
	} else { bas <- 0 }
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# If intermediate growth rates are requested : defining nperiod, period.start and period.end
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (!is.na(period)) {
		# start time of each period
		period.start <- 1 + seq((start.set-1),ntime-period,period)
		# end time of each period
		period.end <- period.start + period - 1
		# number of periods
		nperiod <- length(period.start)
		data.jags <- c(data.jags,list(nperiod=nperiod,period.start=period.start,period.end=period.end))
	}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
return(list(data.jags=data.jags,sta=sta,bas=bas,pro=protocol)) }

