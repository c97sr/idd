# Copyright Steven Riley (sr@stevenriley.net)
#
# This file is part of the library idd.
#
# idd is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This work is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this idsource.  If not, see <http://www.gnu.org/licenses/>.

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
  g <- ggplot(incidence_data, aes(x = t))
  if("forecast" %in% colnames(incidence_data)) {
    g <- g + geom_point(aes(y = incidence, color = time_used_to_forecast)) +
      scale_color_manual(breaks = c(F, T), values = c("black", "red"))
  } else {
    g <- g + geom_point(aes(y = incidence))
  }
  g <- g + scale_x_continuous("Week", breaks = incidence_data[label_index, "t"],
                       labels = incidence_data[label_index, "time_name"]) +
    theme(axis.text.x= element_text(angle = 90),
          legend.position = "none") 
  if(log_scale) {
    g <- g + scale_y_log10("incidence")
  }
  if("forecast" %in% colnames(incidence_data)) {
    g <- g + geom_point(aes(y = forecast), color = "blue")
  }
  return(g)
}

#' subset incidence data for weeks of interest
#' 
#' \code{subset_incidence_data} takes the incidence data and subsets it for the weeks of interest.
#' 
#' @param incidence_data data frame extracted by \code{extract_incidence}
#' @param current_week numeric vector of length 1: week number of the current week
#' @param n_weeks_prior numeric vector of length 1: number of weeks' previous data to use
#' @return data frame for weeks of interest
#' @export
#' 
subset_incidence_data <- function(incidence_data, current_week, n_weeks_prior) {
  week_string <- as.character(current_week)
  if(current_week < 10) {
    week_string <- paste0("0", week_string)
  }
  end_row <- grep(paste0("-", week_string), incidence_data$time_name)
  start_row <- end_row - n_weeks_prior
  incidence_data <- incidence_data[seq(start_row, end_row),]
  return(incidence_data)
}

#' performs linear regression for a given country and year
#' 
#' \code{linear_regression} takes the incidence data (already subsetted for 
#' country and year) and performs linear regression.
#' We also specify the "current week" and how many previous weeks' worth of data
#' to use.
#' 
#' @param incidence_data data frame extracted by \code{extract_incidence}
#' @param current_week numeric vector of length 1: week number of the current week
#' @param n_weeks_prior numeric vector of length 1: number of weeks' previous data to use
#' @param log_transform logical vector of length 1: if TRUE, perform linear 
#' regression on log transformed data, otherwise perform linear regression on 
#' original data
#' @return linear regression object (see documentation for \code{lm})
#' @export
#' 
linear_regression <- function(incidence_data,
                              current_week = 47,
                              n_weeks_prior = 4,
                              log_transform = TRUE) {

  incidence_data <- subset_incidence_data(incidence_data, current_week,
                                          n_weeks_prior)
  
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
#' @param weeks_ahead numeric vector of length 1: number of weeks to forecast ahead
#' @param log_transform logical vector of length 1: if TRUE, linear 
#' regression was performed on log transformed data, otherwise linear regression
#' was performed on original data
#' @return data frame with predicted points and original data
#' @export
#' 
extract_forecasted_points <- function(lm_output, incidence_data, weeks_ahead,
                                      log_transform) {
  
  intercept <- lm_output$coefficients[["(Intercept)"]]
  gradient <- lm_output$coefficients[["t"]]
  last_timepoint_used <- lm_output$model$t[length(lm_output$model$t)]
  if(missing(weeks_ahead)) {
    weeks_ahead <- nrow(incidence_data) - last_timepoint_used
  }
  times_to_forecast <- seq_len(weeks_ahead) + last_timepoint_used
  forecasted_points <- intercept + gradient * times_to_forecast
  if(log_transform) {
    forecasted_points <- exp(forecasted_points)
  }
  
  time_used_to_forecast <- rep(FALSE, nrow(incidence_data))
  time_used_to_forecast[lm_output$model$t] <- TRUE
  incidence_data$time_used_to_forecast <- time_used_to_forecast
  incidence_data$forecast <- rep(NA, nrow(incidence_data))
  incidence_data$forecast[last_timepoint_used + seq_len(weeks_ahead)] <- forecasted_points
  return(incidence_data)
}

#' calculate log likelihood of data given model prediction (assuming Poisson intensity)
#' 
#' \code{calc_log_likelihood} calculates the log likelihood of the incidence data 
#' given model prediction (assuming Poisson intensity)
#' 
#' @param data numeric vector: incidence data
#' @param model_prediction numeric vector of same length as \code{data}: 
#' incidence according to model 
#' @return log likelihood: numeric vector of length 1
#' @export
#' 
calc_log_likelihood <- function(data, model_prediction) {
  stopifnot(length(data) == length(model_prediction) &&
              identical(data, round(data)))
  log_likelihood <- sum(dpois(data, model_prediction, log = TRUE))
  return(log_likelihood)
}

#' simplified version of SEIR model with default inputs
#' 
#' \code{simple_comp_seir} is a simplified version of \code{comp_seir}
#' with default inputs
#' 
#' @param R_0 numeric vector of length 1: basic reproduction number
#' @return A list of two elements. The first is the incidence of infection
#' and the second is the timepoints to which the incidence refers.
#' @export
simple_comp_seir <- function(R_0) {
  latent_period <- 1.6 # Cori et al. (2012) 10.1016/j.epidem.2012.06.001
  infectious_period <- 1 # Cori et al. (2012) 10.1016/j.epidem.2012.06.001
  pop_size <- 323.1e6 # USA population size 2016
  model_output <- comp_seir(De=latent_period,
            Tg=infectious_period + latent_period,
            R0=R_0,
            N=pop_size,
            I0=10,
            trickle=0,
            trickleStart=0,
            reprate=1,
            dt=1,
            R1=R_0,
            t1=999,
            R2=R_0,
            t2=9999,
            noTimeSteps=360,
            noReals=1,
            A=0,
            deterministic=TRUE
  )
  return(model_output[c(1,2)])
}

#' #' fit the SEIR model to data
#' #' 
#' #' \code{fit_seir} fits an SEIR model to incidence data
#' #' 
#' #' @param incidence_data data frame extracted by \code{extract_incidence}
#' #' @param current_week numeric vector of length 1: week number of the current week
#' #' @param n_weeks_prior numeric vector of length 1: number of weeks' previous data to use
#' #' @return A list of two elements. The first is the incidence of infection
#' #' and the second is the timepoints to which the incidence refers.
#' #' @export
#' fit_seir <- function(incidence_data, current_week, n_weeks_prior) {
#'   
#'   evaluate_model_and_calc_log_likelihood <- function(R_0) {
#'     model_output <- simple_comp_seir(R_0)
#'     model_prediction <- model_output$inf_inc[seq(1, n_weeks_prior * 7 + 1, by = 7)]
#'     log_likelihood <- calc_log_likelihood(data, model_prediction)
#'   }
#'   
#'   incidence_data <- subset_incidence_data(incidence_data, current_week,
#'                                           n_weeks_prior)
#'   browser()
#'   return(NULL)
#' }
