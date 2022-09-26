#' Creating input date files for unbiased Bayesian population dynamics models
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
#' river.unbiased <- riverfish[riverfish$sample == 1,]
#' varlist <- convertselect_unbiased(river.unbiased,
#' time=sort(unique(riverfish$date)),
#' start=min(riverfish$date),
#' time.interval=1,
#' type="A",
#' period=NA,
#' para.bas=FALSE)
#' }
convertselect_unbiased <- function(data,time,start,time.interval,type,period,para.bas) {
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Coding used functions
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# Replaces the missing initial values of a matrix mat with the replacement value
	initial_value <- function(mat,min.time,start.time,replacement) {
		if (!is.na(start.time) & min.time >= start.time) { mat[start.time:min.time] <- replacement }
	return(mat) }
	# Replaces the missing initial values in biomass using the abundance-biomass relationship
	initial_biomass_value <- function(wmat,ymat,wy,start.time) {
		if (!is.na(wmat[start.time]) & !is.na(wmat[start.time])) { wmat[start.time] <- round(ymat[start.time] * wy) }
	return(wmat) }
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Creating y and w variables
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
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## If type is "O", surface and protocol are not required
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (type == "O") {
		data$surface <- 1
		data$pro_id <- 1
		nb.protocol <- data.frame("sta_id"=unique(data$sta_id),"x"=1)
	} else {
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Determining number of sampling protocols at sites-based scale and removing sites with more than 2 sampling protocols
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		data.period <- data[!data$date < start,]
		nb.protocol <- aggregate(data.period$pro_id,list("sta_id"=data.period$sta_id),function(x) length(unique(x)))
		data <- data[!data$sta_id %in% nb.protocol$sta_id[nb.protocol$x > 2],]
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Removing sites with more than 1 sampling protocol changes
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		sta.2.protocol <- nb.protocol$sta_id[nb.protocol$x == 2]
		nb.change <- sapply(sta.2.protocol, function (i) which(FALSE %in% (data.period$pro_id[data.period$sta_id == i] == lead(data.period$pro_id[data.period$sta_id == i]))))
		data <- data[!data$sta_id %in% sta.2.protocol[which(nb.change > 1)],]
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Removing dates with sampling protocol not used over the period defined by start
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		data$date.1 <- sapply(data$sta_id, function(i) min(data$date[(data$sta_id == i & data$pro_id %in% data$pro_id[(data$sta_id == i & data$date == start)])]))
		data$y <- ifelse(data$date >= start, data$y, ifelse(data$date < data$date.1, NA, data$y))
		data$w <- ifelse(data$date >= start, data$w, ifelse(data$date < data$date.1, NA, data$w))
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Removing sites with insufficient sampling dates for some protocols
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		date.by.protocol <- tapply(data$date[!is.na(data$y)], list(data$sta_id[!is.na(data$y)], data$pro_id[!is.na(data$y)]), function(i) length(unique(i)))
		remove.sta <- rownames(which(date.by.protocol < 2, arr.ind = TRUE))
		data <- data[!data$sta_id %in% remove.sta,]
	}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Defining the sites vectors (sta, sta.1, sta.2) and the corresponding number of sites (nsta, nsta.1, nsta.2)
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# site identifier vector
	sta <- sort(unique(data$sta_id))
	nsta <- length(sta)
	if (nsta == 0) {
		stop("no site meets the model conditions",call. = FALSE)
	}
	# positions of sites involving one sampling protocol in sta
	sta.1 <- which(sta %in% sort(nb.protocol$sta_id[nb.protocol$x == 1]))
	nsta.1 <- length(sta.1)
	# positions of sites involving one sampling protocol in sta
	sta.2 <- which(sta %in% sort(nb.protocol$sta_id[nb.protocol$x == 2]))
	nsta.2 <- length(sta.2)
	#-----------------------------------------------------------------------------------------------------------nsta--------------------------------------------------------------
	## Defining the number of time step (ntime) and the position of the starting date in vector time (start.set)
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	ntime <- length(time)
	start.set <- which(time %in% (start + time.interval))
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Designing two-dimensional variable matrices (time, sta)
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# site occupancy states
	z <- tapply(data$y, list(data$date,data$sta_id), function(i) ifelse(i > 0,1,0))
	# abundance
	y.trans <- tapply(data$y, list(data$date,data$sta_id), function(i) ifelse(i == 0,1,i))
	yc <- tapply(data$y, list(data$date,data$sta_id), function(i) ifelse(is.na(i),0,i))
	# biomass
	w.trans <- tapply(data$w, list(data$date,data$sta_id), function(i) ifelse(i == 0,1,i))
	wc <- tapply(data$w, list(data$date,data$sta_id), function(i) ifelse(is.na(i),0,i))
	# sampling surface
	S <- tapply(data$surface, list(data$date,data$sta_id), unique)
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Defining the start of time series at sites-based scale (start.1) and assigning missing initial values
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	time.1 <- as.vector(sapply(sta, function(i) which(time %in% min(data$date[(data$sta_id == i & data$pro_id %in% data$pro_id[(data$sta_id == i & data$date == start)] & !is.na(data$y))]))))
	start.1 <- sapply(1+time.1, function(i) min(start.set,i))
	# assign default initial values of z, y.trans and w.trans
	z <- sapply(seq_len(ncol(z)), function(i) initial_value(z[,i],time.1[i]-1,start.set-1,0))
	y.trans <- sapply(seq_len(ncol(y.trans)), function(i) initial_value(y.trans[,i],time.1[i]-1,start.set-1,1))
	w.trans <- sapply(seq_len(ncol(w.trans)), function(i) initial_value(w.trans[,i],time.1[i]-1,start.set-1,1))
	# assign initial biomass from the relationship between abundance and biomass if missing (for type "AB" or "OAB")
	if (type %in% c("AB","OAB")) {
		w.trans <- sapply(seq_len(ncol(w.trans)), function(i) initial_biomass_value(w.trans[,i],y.trans[,i],unique(data$wy[data$sta_id == sta[i]]),start.1[i]-1))
	}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Defining the end of times series for the first protocol (end.1), the start of time series for the second protocol (start.2) and assigning missing initial values
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	start.2 <- end.1 <- vector("numeric",nsta)
	if (nsta.2 > 0) {
		time.of.change <- sapply(1:nsta, function(i) ifelse (i %in% sta.2, which(time %in% min(data$date[(data$sta_id == sta[i] & !data$pro_id %in% data$pro_id[(data$sta_id == sta[i] & data$date == start)])])), NA))
		end.1 <- time.of.change - 1
		start.2 <- time.of.change + 1
		time.2 <- sapply(1:nsta, function(i) ifelse (i %in% sta.2, which(time %in% min(data$date[(data$sta_id == sta[i] & !data$pro_id %in% data$pro_id[(data$sta_id == sta[i] & data$date == start & !is.na(data$y))])])), NA))
		# assign default initial values of z, y.trans and w.trans of second times series if missing
		z <- sapply(seq_len(ncol(y.trans)), function(i) initial_value(z[,i],time.2[i],time.of.change[i],0))
		y.trans <- sapply(seq_len(ncol(y.trans)), function(i) initial_value(y.trans[,i],time.2[i],time.of.change[i],0))
		w.trans <- sapply(seq_len(ncol(w.trans)), function(i) initial_value(w.trans[,i],time.2[i],time.of.change[i],0))
		# assign initial biomass from the relationship between abundance and biomass of second times series if missing (for type "AB" or "OAB")
		if (type %in% c("AB","OAB")) {
			w.trans <- sapply(seq_len(ncol(w.trans)), function(i) initial_biomass_value(w.trans[,i],y.trans[,i],unique(data$wy[data$sta_id == sta[i]]),time.of.change[i]))
		}
	}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Selecting and converting the required variables
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	list.O <- list(ntime=ntime,nsta=nsta,start.set=start.set,start.1=start.1,lap=time.interval,z=z)
	list.A <- list(y.trans=y.trans,yc=yc,S=S,nsta.1=nsta.1,sta.1=sta.1,end.1=end.1,start.2=start.2,nsta.2=nsta.2,sta.2=sta.2)
	list.B <- list(w.trans=w.trans,wc=wc,S=S,nsta.1=nsta.1,sta.1=sta.1,end.1=end.1,start.2=start.2,nsta.2=nsta.2,sta.2=sta.2)
	list.AB <- list(w.trans=w.trans,wc=wc,y.trans=y.trans,yc=yc,S=S,nsta.1=nsta.1,sta.1=sta.1,end.1=end.1,start.2=start.2,nsta.2=nsta.2,sta.2=sta.2)
	#---------------------------------------------------------------------------------------------------------z----------------------------------------------------------------
	if (type == "A" | type == "OA") {
		data.jags <- c(list.O,list.A)
	}	else if (type == "B" | type == "OB") {
		data.jags <- c(list.O,list.B)
	} else if (type == "AB" | type == "OAB") {
		data.jags <- c(list.O,list.AB)
	} else  {
		data.jags <- list.O
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
return(list(data.jags=data.jags,sta=sta,bas=bas)) }
