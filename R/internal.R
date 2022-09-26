#' Tests if a given variable in a dataframe is numeric
#'
#' @importFrom dplyr enquo pull
#' @import magrittr
#'
#' @param data Dataframe
#' @param variable Character
#'
#' @return data %>% pull(!!variable) %>% is.numeric() Logical
#' @keywords internal
#' @noRd
test_var_is_numeric <- function(data, variable) {
	enquo(variable)
	if (TRUE %in% is.element(variable, colnames(data))) {
		data %>% pull(!!variable) %>% is.numeric()
	}
}


#' Tests if a given variable in a dataframe has missing values
#'
#' @importFrom dplyr enquo pull
#' @import magrittr
#'
#' @param data Dataframe
#' @param variable Character
#'
#' @return NA %in% (data %>% pull(!!variable)) Logical
#' @keywords internal
#' @noRd
test_var_has_na <- function(data, variable) {
	enquo(variable)
	if (TRUE %in% is.element(variable, colnames(data))) {
		NA %in% (data %>% pull(!!variable))
	}
}


#' Tests if a dataframe includes missing dates
#'
#' @param data Dataframe
#' @param time Numeric vector
#'
#' @return sta.id[which(time.in.data == TRUE)] Vector
#' @keywords internal
#' @noRd
test_missing_date <- function(data, time) {
	sta.id <- unique(data$sta_id)
	time.in.data <- sapply(sta.id, function(i) length(unique(data$date[data$sta_id == i])) < length(time))
	return(sta.id[which(time.in.data == TRUE)])
}


#' Tests if a dataframe includes duplicate dates with different values
#'
#' @param data Dataframe
#'
#' @return unique(data[data$duplicate %in% duplicate[which(duplicate.date == TRUE)],"sta_id"]) Vector
#' @keywords internal
#' @noRd
test_duplicate_date <- function(data) {
	if (is.null(data$sample)) {
		data$duplicate <- data$sta_id
	} else {
		data$duplicate <- paste0(data$sta_id, data$sample)
	}
	duplicate <- unique(data$duplicate)
	duplicate.date <- unlist(lapply(duplicate, function(i) duplicated(data$date[data$duplicate == i])))
	return(unique(data[data$duplicate %in% duplicate[which(duplicate.date == TRUE)],"sta_id"]))
}


#' Tests if sample positions are consistent with the number of samples
#'
#' @importFrom stats aggregate
#'
#' @param data Dataframe
#'
#' @return sta Vector
#' @keywords internal
#' @noRd
test_sample_position <- function(data) {
	sta <- vector("numeric",0)
	if (TRUE %in% is.element("sample", colnames(data))) {
		number.of.sample <- aggregate(data$sample, list("sta_id"=data$sta_id, data$date), function(i) max(i) == length(i))
		sta <- unique(number.of.sample$sta_id[number.of.sample$x == FALSE])
	}
	return(sta)
}


#' Tests if a variables has multiple values per site and date
#'
#' @importFrom dplyr enquo
#'
#' @param data Dataframe
#' @param variable Character
#'
#' @return rownames(which(multiple.var == TRUE, arr.ind = TRUE)) Vector
#' @keywords internal
#' @noRd
test_var_multiple_values <- function(data,variable) {
	enquo(variable)
	sta.id <- vector("numeric",0)
	if (TRUE %in% is.element(variable, colnames(data))) {
		data$variable <- data[,which(colnames(data) %in% variable)]
		multiple.var <- tapply(data$var, list(data$sta_id,data$date), function(i) length(unique(i)) > 1)
		sta.id <- rownames(which(multiple.var == TRUE, arr.ind = TRUE))
	}
	return(sta.id)
}


#' Imputs missing values
#'
#' @param data Dataframe
#' @param variable Character
#' @param time.interval,start,end Numeric
#'
#' @return data Dataframe
#' @keywords internal
#' @noRd
imput_missing_value <- function(data,variable,time.interval,start,end) {
	n <- which(colnames(data) %in% variable)
	sta <- unique(data$sta_id[is.na(data[,n])])
	for (i in sta) {
		date <- unique(data$date[(data$sta_id == i & is.na(data[,n]))])
		for (t in date[!date == start]) {
			data[(data$sta_id == i & data$date == t),n] <- unique(data[(data$sta_id == i & data$date == (t - time.interval)),n])
		}
		for (t in rev(date[!date == end])) {
			data[(data$sta_id == i & data$date == t),n] <- unique(data[(data$sta_id == i & data$date == (t + time.interval)),n])
		}
	}
	return(data)
}
