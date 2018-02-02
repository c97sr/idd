#' extracts the incidence data for a given country and year
#' 
#' \code{extract_incidence} takes the flu incidence matrix and extracts
#' the incidence for a given flu season and a given country.
#' 
#' @param flu_data matrix with flu incidence data
#' @param country_code character vector of length 1: 3-letter country code
#' @param year numeric vector of length 1: year to extract
#' @return data frame with three columns:
#' t: numeric vector of length 52/53 (depending on the year): 
#' index running from 1 to 52/53
#' time_name: character vector of the same length as t:
#' week numbers of the data we've extracted
#' incidence: numeric vector of the same length as t:
#' incidence for those weeks
#' @export
#' 
extract_incidence <- function(flu_data,
                              country_code = "USA",
                              year = 2016) {
  flu_data <- as.data.frame(flu_data)
  year_names <- rownames(flu_data)
  row_name_start <- paste0(year, "-27")
  row_name_end <- paste0(year + 1, "-26")
  row_index_start <- which(rownames(flu_data) == row_name_start)
  row_index_end <- which(rownames(flu_data) == row_name_end)
  incidence <- flu_data[seq(row_index_start, row_index_end), 
                        colnames(flu_data) == country_code]
  time_name_vec <- year_names[seq(row_index_start, row_index_end)]
  
  incidence_data <- data.frame(t = seq_along(time_name_vec), 
                        time_name = time_name_vec, 
                        incidence = incidence)
  return(incidence_data)
}

#' plots the incidence
#' 
#' \code{plot_incidence} plots the incidence extracted by \code{extract_incidence}
#' 
#' @param incidence_data data frame extracted by \code{extract_incidence}
#' @param log_scale logical vector of length 1: if TRUE, plot on log scale, 
#' otherwise plot on linear scale
#' @return ggplot object
#' @import ggplot2
#' @export
plot_incidence <- function(incidence_data, log_scale = FALSE) {

  label_x_axis_every <- 5
  label_index <- seq(1, nrow(incidence_data), by = label_x_axis_every)
  g <- ggplot(incidence_data, aes(x = t)) + geom_point(aes(y = incidence)) +
    scale_x_continuous("Week", breaks = incidence_data[label_index, "t"],
                       labels = incidence_data[label_index, "time_name"]) +
    theme(axis.text.x= element_text(angle = 90)) 
  if(log_scale) {
    g <- g + scale_y_log10("incidence")
  }
  if("forecast" %in% colnames(incidence_data)) {
    g <- g + geom_point(aes(y = forecast), color = "blue")
  }
  return(g)
}

#' performs linear regression for a given country and year
#' 
#' \code{linear_regression} takes the flu incidence matrix and performs
#' linear regression for a given flu season and a given country.
#' We also specify the "current week" and how many previous weeks' worth of data
#' to use.
#' 
#' @param flu_data matrix with flu incidence data
#' @param country_code character vector of length 1: 3-letter country code
#' @param year numeric vector of length 1: year to plot
#' @param current_week numeric vector of length 1: week number of the current week
#' @param n_weeks_prior numeric vector of length 1: number of weeks' previous data to use
#' @param log_transform logical vector of length 1: if TRUE, perform linear 
#' regression on log transformed data, otherwise perform linear regression on 
#' original data
#' @return linear regression object (see documentation for \code{lm})
#' @export
#' 
linear_regression <- function(flu_data,
                              country_code = "USA",
                              year = 2016,
                              current_week = 47,
                              n_weeks_prior = 4,
                              log_transform = TRUE) {
  incidence_data <- extract_incidence(flu_data, country_code, year)
  week_string <- as.character(current_week)
  if(current_week < 10) {
    week_string <- paste0("0", week_string)
  }
  end_row <- grep(paste0("-", week_string), incidence_data$time_name)
  start_row <- end_row - n_weeks_prior
  incidence_data <- incidence_data[seq(start_row, end_row),]
  if(log_transform) {
    incidence_data$incidence <- log(incidence_data$incidence)
  }
  
  linear_regression_output <- lm(incidence ~ t, data = incidence_data)
  
  return(linear_regression_output)
}

#' extract forecasted points from linear regression object
#' 
#' \code{extract_forecasted_points} takes the linear regression object
#' returned by \code{linear_regression} and extracts model predictions
#' 
#' @param lm_output linear regression object returned by \code{linear_regression}
#' @param incidence_data incidence data extracted by \code{extract_incidence}
#' @param log_transform logical vector of length 1: if TRUE, linear 
#' regression was performed on log transformed data, otherwise linear regression
#' was performed on original data
#' @return data frame with predicted points and original data
#' @export
#' 
extract_forecasted_points <- function(lm_output, incidence_data, log_transform) {
  
  intercept <- lm_output$coefficients[["(Intercept)"]]
  gradient <- lm_output$coefficients[["t"]]
  last_timepoint_used <- lm_output$model$t[length(lm_output$model$t)]
  times_to_forecast <- seq(last_timepoint_used + 1, nrow(incidence_data))
  forecasted_points <- intercept + gradient * times_to_forecast
  if(log_transform) {
    forecasted_points <- exp(forecasted_points)
  }
  forecasted_points <- c(rep(NA, last_timepoint_used), forecasted_points)
  incidence_data$forecast <- forecasted_points
  incidence_data
}