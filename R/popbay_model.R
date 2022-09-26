#' Fitting hierarchical Bayesian population dynamics models
#'
#' @import stringr
#' @importFrom dplyr distinct relocate
#' @import coda
#' @import rjags
#' @importFrom R2jags attach.jags autojags jags jags2bugs
#'
#' @param data A data frame
#' @param type,file.model,file.mcmc Character
#' @param error,fit.model Logical
#' @param time.interval,start,end,period,n.chains,n.iter,n.burnin,n.thin Numeric
#' @param parameters.to.save Character vector
#'
#' @return statistics A data frame
#' @export
#' @examples
#' \donttest{
#' data <- data.frame("year"=rep(2010:2020,4),
#' "region"=rep(c("Occitanie","PACA"),each=22),
#' "site"=rep(c("S1","S2","S3","S4"),each=11),
#' "abundance"=sample(0:100,44))
#' new.data <- popbay_data(data,
#' ColNames=list(sta_id="site",date="year",bas_id="region"))
#' new.data$surface <- 1
#' new.data$pro_id <- 1
#' stat <- popbay_model(new.data,n.iter=1000,n.burnin=10,n.thin=10)
#' }
#' \dontrun{
#' data(riverfish)
#' stat <- popbay_model(riverfish,
#' type="AB",
#' n.iter=1000,n.burnin=10,n.thin=10)
#' stat <- popbay_model(riverfish,
#' type="AB",
#' period=6,
#' n.iter=1000,n.burnin=10,n.thin=10)
#' stat <- popbay_model(riverfish,
#' type="AB",
#' time.interval=2,
#' n.iter=1000,n.burnin=10,n.thin=10)
#' stat <- popbay_model(riverfish,
#' type="AB",
#' start=2010,
#' n.iter=1000,n.burnin=10,n.thin=10)
#' stat <- popbay_model(riverfish,
#' type="AB",
#' error=TRUE,
#' n.iter=1000,n.burnin=10,n.thin=10)
#' para <- popbay_model(riverfish,
#' type="OAB",
#' error=TRUE,
#' fit.model=FALSE,
#' file.model="complete")
#' load("complete_datamodel.rda")
#' stat <- jags(data = data.jags,
#' parameters.to.save=para,
#' model.file = "complete_model.txt")
#' }
popbay_model <- function(data,type="A",error=FALSE,time.interval=1,start,end,period=NA,
												 parameters.to.save,n.chains=3,n.iter=102000,n.burnin=2000,n.thin=1000,fit.model=TRUE,file.model,file.mcmc) {
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Checking the variables format, if the format is wrong the function returns an error message and stops
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (FALSE %in% test_var_is_numeric(data, "date")) {
		stop("'date' must be a numeric field",call. = FALSE)
	}
	if (FALSE %in% test_var_is_numeric(data, "abundance")) {
		stop("'abundance' must be a numeric field",call. = FALSE)
	}
	if (FALSE %in% test_var_is_numeric(data, "biomass")) {
		stop("'biomass' must be a numeric field",call. = FALSE)
	}
	if (FALSE %in% test_var_is_numeric(data, "surface")) {
		stop("'surface' must be a numeric field",call. = FALSE)
	}
	if (FALSE %in% test_var_is_numeric(data, "sample")) {
		stop("'sample' must be a numeric field",call. = FALSE)
	}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Checking for missing values in surface, pro_id and sample, if this occurs the function returns an error message and stops
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (TRUE %in% test_var_has_na(data, "surface")) {
		stop("'surface' must not have missing values",call. = FALSE)
	}
	if (TRUE %in% test_var_has_na(data, "pro_id")) {
		stop("'pro_id' must not have missing values",call. = FALSE)
	}
	if (TRUE %in% test_var_has_na(data, "sample")) {
		stop("'sample' must not have missing values",call. = FALSE)
	}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Checking for consistency between the sample positions and the number of samples, if this is not consistent the function returns an error message and stops
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (length(test_sample_position(data)) > 0) {
		stop("sample positions and number of samples are not consistent",call. = FALSE)
	}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Assigning default values to 'start' and 'end' if missing
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	# numeric: defines the starting date for estimating average rates
	if (missing(start)) {
		start <- min(data$date)
	}
	# numeric: defines the ending date for estimating average rates
	if (missing(end)) {
		end <- max(data$date)
	}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Defining the dates vector (time), filtering dates and removing duplicate rows
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	time <- c(rev(-1 * seq(from=-start,to=-min(data$date),by=time.interval)),seq(from=start+time.interval,to=end,by=time.interval))
	data <- data[(data$date %in% time & order(data$sta_id,data$date)),] %>%
		distinct()
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Checking for missing dates between start and end, if this occurs the function returns an error message and stops
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if(length(test_missing_date(data,time)) > 0) {
		stop("dataset must not have missing dates",call. = FALSE)
	}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Checking for duplicate dates with different values, if this occurs the function returns an error message
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (length(test_duplicate_date(data)) > 0) {
		stop("dataset must not have repeated dates with different values",call. = FALSE)
	}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Assigning default values to the character vector of names of the parameters to save if missing
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (missing(parameters.to.save)) {
		parameters.to.save <- vector("character",0)
		if (str_detect(type,"A") == "TRUE") {
			parameters.to.save <- c(parameters.to.save, c("R.nat","R.sta"))
			if (error == TRUE) {
				parameters.to.save <- c(parameters.to.save, c("p.cap"))
			}
			if (!is.na(period)) {
				parameters.to.save <- c(parameters.to.save, c("R.nat.period","R.sta.period"))
			}
			if (!is.null(data$bas_id)) {
				parameters.to.save <- c(parameters.to.save, c("R.bas"))
				if (!is.na(period)) { parameters.to.save <- c(parameters.to.save, c("R.bas.period")) }
			}
		}
		if (str_detect(type,"B") == "TRUE") {
			parameters.to.save <- c(parameters.to.save, c("Rw.nat","Rw.sta"))
			if (error == TRUE) {
				parameters.to.save <- c(parameters.to.save, c("var.cap"))
			}
			if (!is.na(period)) {
				parameters.to.save <- c(parameters.to.save, c("Rw.nat.period","Rw.sta.period"))
			}
			if (!is.null(data$bas_id)) {
				parameters.to.save <- c(parameters.to.save, c("Rw.bas"))
				if (!is.na(period)) { parameters.to.save <- c(parameters.to.save, c("Rw.bas.period")) }
			}
		}
		if (str_detect(type,"O") == "TRUE") {
			parameters.to.save <- c(parameters.to.save, c("Ro.nat"))
			if (error == TRUE) {
				parameters.to.save <- c(parameters.to.save, c("p.cap.sp"))
			}
			if (!is.na(period)) {
				parameters.to.save <- c(parameters.to.save, c("Ro.nat.period"))
			}
			if (!is.null(data$bas_id)) {
				parameters.to.save <- c(parameters.to.save, c("Ro.bas"))
				if (!is.na(period)) { parameters.to.save <- c(parameters.to.save, c("Ro.bas.period")) }
			}
		}}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Detecting if parameters are requested at site (para.sta) regional (para.bas) and national (para.nat) levels
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	para.sta <- ifelse(TRUE %in% str_detect(parameters.to.save, "sta|z|y|N|w|B|cap|p.per|p.col"), TRUE, FALSE)
	para.bas <- ifelse(TRUE %in% str_detect(parameters.to.save, "bas"), TRUE, FALSE)
	para.nat <- ifelse(TRUE %in% str_detect(parameters.to.save, "nat"), TRUE, FALSE)
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Selecting and converting the required variables for the unbiased growth model
	## 	and writing the model depending on the modelling type and the requested parameters
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (error == FALSE) {
		if (!is.null(data$sample)) { data <- data[data$sample == 1,] }
		data.model <- do.call(convertselect_unbiased,list(data,time,start,time.interval,type,period,para.bas))
		data.jags <- data.model$data.jags
		gw.model <- do.call(writemodel_unbiased,list(type,period,para.nat,para.bas,para.sta))
	#-------------------------------------------------------------------------------------------------------------
	## Selecting and converting to the jags format the required variables for the error growth model
	## 	and writing the model depending on the modelling type and the requested parameters
	#-------------------------------------------------------------------------------------------------------------
	} else {
		data.model <- do.call(convertselect_removal,list(data,time,start,time.interval,type,period,para.bas))
		data.jags <- data.model$data.jags
		gw.model <- do.call(writemodel_removal,list(type,period,para.nat,para.bas,para.sta))
	}
	#-------------------------------------------------------------------------------------------------------------
	## Saving the converted data and model files
	#-------------------------------------------------------------------------------------------------------------
	if(!missing(file.model)) {
		write(gw.model,paste(file.model,"_model.txt",sep=""))
		save(data.jags,file=paste(file.model,"_datamodel.rda",sep=""))
	}
	#-------------------------------------------------------------------------------------------------------------
	## Fitting the model (gw.model) to converted data (data.jags) and collecting the models outputs
	#-------------------------------------------------------------------------------------------------------------
	if (fit.model == TRUE) {
		set.seed(1989)
		run <- jags(data = data.jags,inits = NULL,parameters.to.save, model.file = textConnection(gw.model),n.chains,n.iter,n.burnin,n.thin)
		#------------------------------------------------------------------------------------------------------------
		# Saving the Markov chains file
		run.mcmc <- as.mcmc(run)
		if (!missing(file.mcmc)) {
			save(run.mcmc,file=paste(file.mcmc,"_mcmc.rda",sep=""))
		}
		#------------------------------------------------------------------------------------------------------------
		# Collecting the statistical outputs
		stat.run <- run$BUGSoutput
		stat <- stat.run$summary
		statistics <- as.data.frame(stat)
		colnames(statistics) <- c("mean","sd","Q2.5","Q25","Q50","Q75","Q97.5","Rhat","n.eff")
		# get parameter names from row names of statistics
		statistics$position <- as.vector(str_match(rownames(statistics), "\\d+,\\d+"))
		statistics$position <- ifelse(is.na(statistics$position),as.vector(str_match(rownames(statistics), "\\d+")),statistics$position)
		statistics$Para <- rownames(statistics)
		statistics$Parameters <- ifelse(is.na(statistics$position),statistics$Para,str_replace(sub("[[]","/",statistics$Para),paste("/",statistics$position,"]",sep=""),""))
		# get postions from row names of statistics
		statistics$p1 <- as.numeric(str_match(statistics$position,"\\d+"))
		statistics$p2 <- as.numeric(str_match(str_match(statistics$position,",\\d+"),"\\d+"))
		# get date, site identifer and/or sampling protocol for site-base scale parameters
		statistics$date <- NA
		if (para.sta == TRUE) {
			sta <- data.model$sta
			statistics$date <- ifelse(statistics$Parameters %in% c("R.sta.period","Rw.sta.period","r.sta","rw.sta","ro.sta","z","y","N","w","B"),time[statistics$p1],statistics$date)
			statistics$sta_id <- ifelse(statistics$Parameters %in% c("R.sta.period","Rw.sta.period","r.sta","rw.sta","ro.sta","z","y","N","w","B","p.cap.sp","p.cap","var.cap"),sta[statistics$p2],ifelse(statistics$Parameters %in% c("R.sta","Rw.sta","Ro.sta","p.per","p.col","nb.r.sta","nb.rw.sta"),sta[statistics$p1],NA))
			if (TRUE %in% str_detect(statistics$Parameters,"cap")) {
				pro <- data.model$pro
				statistics$pro_id <- ifelse(statistics$Parameters %in% c("p.cap.sp","p.cap","var.cap"),pro[statistics$p1],NA)
			}}
		# get date and region identifer for regional scale parameters
		if (para.bas == TRUE) {
			bas <- data.model$bas
			statistics$date <- ifelse(statistics$Parameters %in% c("R.bas.period","Rw.bas.period","Ro.bas.period","r.bas","rw.bas","ro.bas"),time[statistics$p1],statistics$date)
			statistics$bas_id <- ifelse(statistics$Parameters %in% c("R.bas.period","Rw.bas.period","Ro.bas.period","r.bas","rw.bas","ro.bas"),bas[statistics$p2],ifelse(statistics$Parameters %in% c("R.bas","Rw.bas","Ro.bas"),bas[statistics$p1],NA))
		}
		# get date for national scale parameters
		if (para.nat == TRUE) {
			statistics$date <- ifelse(statistics$Parameters %in% c("R.nat.period","Rw.nat.period","Ro.nat.period","r.nat"),time[statistics$p1],statistics$date)
		}
		statistics$position <- statistics$p1 <- statistics$p2 <- statistics$Para <- NULL
		if (TRUE %in% all(is.na(statistics$date))) { statistics$date <- NULL }
		return(statistics)
	}	else {
		return(parameters.to.save)
	}
}
