#' @import utils
if(getRversion() >= "2.15.1")  utils::globalVariables(c("sta_id"))

#' Checking and converting a dataset for popbay_model
#'
#' @import stats
#' @importFrom dplyr distinct arrange lead
#' @import magrittr
#'
#' @param data A data frame
#' @param ColNames A list
#' @param time.interval,start,end Numeric
#' @param na.replace A logical
#'
#' @return data A data frame
#' @export
#'
#' @examples
#' \donttest{
#' data <- data.frame("year"=rep(2010:2020,4),
#' "region"=rep(c("Occitanie","PACA"),each=22),
#' "site"=rep(c("S1","S2","S3","S4"),each=11),
#' "abundance"=sample(0:100,44))
#' new.data <- popbay_data(data,
#' ColNames=list(sta_id="site",date="year",bas_id="region"))
#' }
popbay_data <- function(data,ColNames,time.interval,start,end,na.replace=FALSE) {
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Standardising column name and checking if column names contains sta_id and date
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (!missing(ColNames)) {
		if (!is.null(ColNames$sta_id)) {
			colnames(data)[which(colnames(data) %in% ColNames$sta_id)] <- "sta_id"
		}
		if (!is.null(ColNames$date)) {
			colnames(data)[which(colnames(data) %in% ColNames$date)] <- "date"
		}
		if (!is.null(ColNames$sample)) {
			colnames(data)[which(colnames(data) %in% ColNames$sample)] <- "sample"
		}
		if (!is.null(ColNames$abundance)) {
			colnames(data)[which(colnames(data) %in% ColNames$abundance)] <- "abundance" }
		if (!is.null(ColNames$biomass)) { colnames(data)[which(colnames(data) %in% ColNames$biomass)] <- "biomass"
		}
		if (!is.null(ColNames$surface)) {
			colnames(data)[which(colnames(data) %in% ColNames$surface)] <- "surface"
		}
		if (!is.null(ColNames$pro_id)) {
			colnames(data)[which(colnames(data) %in% ColNames$pro_id)] <- "pro_id"
		}
		if (!is.null(ColNames$bas_id)) {
			colnames(data)[which(colnames(data) %in% ColNames$bas_id)] <- "bas_id"
		}
	}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Checking columns names, if column names do not contain sta_id and/or date the function returns an error message and stops
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if (FALSE %in% is.element(c("sta_id","date"),colnames(data))) {
		stop("data frame must contain the 'sta_id' and 'date' items",call. = FALSE)
	}
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
	## Assigning default values to missing function arguments
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		# numeric: defines the time interval between two consecutive dates
		if (missing(time.interval)) {
			time.interval <- 1
		}
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
		time <- seq(start,end,time.interval)
		data <- data[data$date %in% time,] %>%
			distinct() %>%
			arrange(sta_id, date)
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Checking for missing dates between start and end, if this occurs the function returns an error message
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		if(length(test_missing_date(data,time)) > 0) {
			message("Error: dataset includes missing dates (sta_id : ",
								paste0(test_missing_date(data,time),collapse=","),")")
			na.replace <- FALSE
		}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Checking for duplicate dates with different values, if this occurs the function returns an error message
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		if (length(test_duplicate_date(data)) > 0) {
			message("Error: dataset includes duplicate dates with different values (sta_id: ",
								paste0(test_duplicate_date(data),collapse=","),")")
			na.replace <- FALSE
		}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Checking for missing values in sample variable, if this occurs the function returns an error message
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		if (TRUE %in% test_var_has_na(data, "sample")) {
			message("Error: dataset includes missing 'sample' values")
			na.replace <- FALSE
		}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Checking for consistency between the sample positions and the number of samples, if this is not consistent the function returns an error message
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		if (length(test_sample_position(data)) > 0) {
			message("Error: sample positions and number of samples are not consistent (sta_id: ",
								paste0(test_sample_position(data),collapse=","),")")
			na.replace <- FALSE
		}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Checking for a single surface and protocol per site and date else the function returns an error message
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		if (length(test_var_multiple_values(data,"surface")) > 0) {
			message("Error : multiple 'surface'values per site and date (sta_id: ",
								paste0(test_var_multiple_values(data,"surface"),collapse=","),")")
			na.replace <- FALSE
		}
		if (length(test_var_multiple_values(data,"pro_id")) > 0) {
			message("Error : multiple 'pro_id'values per site and date (sta_id: ",
								paste0(test_var_multiple_values(data,"pro_id"),collapse=","),")")
			na.replace <- FALSE
		}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Checking for missing values in surface and pro_id variables, if this occurs and na.replace = TRUE missings values are imputed, else the function returns an error message
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		if (TRUE %in% test_var_has_na(data, "surface")) {
			if (na.replace == TRUE) {
				data <- imput_missing_value(data,"surface",time.interval,start,end)
			} else {
				message("Error: dataset includes missing 'surface' values")
			}
		}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		if (TRUE %in% test_var_has_na(data, "pro_id")) {
			if (na.replace == TRUE) {
				data <- imput_missing_value(data,"pro_id",time.interval,start,end)
			} else {
				message("Error: dataset includes missing 'pro_id' values")
			}
		}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Checking for sites without multiple-pass sampling
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		if (!is.null(data$sample)) {
			sta.id <- unique(data$sta_id)
			nb.sample <- sapply(sta.id, function(i) max(data$sample[data$sta_id == i]))
			if (TRUE %in% (nb.sample == 1)) {
				message(c("Warning (removal model): dataset includes sites without multiple-pass sampling (sta_id: ",
									paste0(sta.id[which(nb.sample == 1)],collapse=","),")"))
			}
		}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Checking for sites with more than 2 sampling protocols, if this occurs the function returns a warning message
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		if (!is.null(data$pro_id)) {
			data.na <- data[!is.na(data$pro_id),]
			sta.id <- unique(data.na$sta_id)
			nb.protocol <- sapply(sta.id, function(i) length(unique(data$pro_id[data.na$sta_id == i])))
			if (TRUE %in% (nb.protocol > 2)) {
				message("Warning (unbiased model): dataset includes sites with more than two sampling protocols (sta_id: ",
								paste0(sta.id[which(nb.protocol > 2)],collapse=","),")")
			}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	## Checking for sites with more than 1 sampling protocol changes, if this occurs the function returns a warning message
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
		nb.change <- sapply(sta.id, function (i) match(FALSE,data.na$pro_id[data.na$sta_id == i] == lead(data.na$pro_id[data.na$sta_id == i]),nomatch = 0))
			if (TRUE %in% (nb.change > 1)) {
				message("Warning (unbiased model): dataset includes sites with more than one sampling protocol changes (sta_id: ",
								paste0(sta.id[which(nb.change > 1)],collapse=","),")")
		}
	}
	#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	return(data)
}
