# Copyright Steven Riley (sr@stevenriley.net), 
# Caroline Walters (c.e.walters6@gmail.com),
# Ada Yan (ada.w.yan@gmail.com)
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
extract_incidence <- function(flu_data,
                              country_code = "ISR",
                              year = 2016) {
  flu_data <- as.data.frame(flu_data)
  year_names <- rownames(flu_data)
  # start plotting at week 27 of the current year
  row_name_start <- paste0(year, "-27") 
  # stop plotting at week 26 of the next year
  row_name_end <- paste0(year + 1, "-26")
  # find the corresponding weeks in the data
  row_index_start <- which(rownames(flu_data) == row_name_start)
  row_index_end <- which(rownames(flu_data) == row_name_end)
  # extrac the week number and incidence for those weeks
  incidence <- flu_data[seq(row_index_start, row_index_end), 
                        colnames(flu_data) == country_code]
  time_name_vec <- year_names[seq(row_index_start, row_index_end)]
  
  incidence_data <- data.frame(t = seq_along(time_name_vec), 
                               time_name = time_name_vec, 
                               incidence = incidence)
  return(incidence_data)
}

#' plots the incidence for a given country (all years)
#'
#' \code{plot_incidence_all} plots the incidence for a given country (all years)
#'
#' @param flu_data matrix with flu incidence data
#' @param country_code character vector of length 1: 3-letter country code
#' @return ggplot object
#' @import ggplot2
#' @export
plot_incidence_all <- function(flu_data, country_code = "ISR") {

  years <- seq(2010, 2016)
  all_incidence <- lapply(years, function(x) extract_incidence(flu_data,
                                            country_code = country_code,
                                                    year = x))
  label_with_season <- function(df, label) {
    df$season <- label
    df
  }
  all_incidence <- Map(label_with_season, all_incidence, years)
  all_incidence <- do.call(rbind, all_incidence)
  all_incidence$season <- factor(all_incidence$season, levels = years)

  g <- ggplot(all_incidence, aes(x = t)) +
    geom_line(aes(y = incidence, color = season), na.rm = TRUE) +
    scale_x_continuous("Weeks since week 26 of calendar year", 
                       breaks = seq(0, max(all_incidence$t), by = 5))
  return(g)
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
  # if the data frame contains a prediction, plot the data points used to
  # predict in red, and the rest in black
  if("prediction" %in% colnames(incidence_data)) {
    g <- g + geom_point(aes(y = incidence, color = time_used_to_predict),
                        na.rm = TRUE) +
      scale_color_manual(breaks = c(F, T), values = c("black", "red"))
  } else {
    # otherwise plot all data points in black
    g <- g + geom_point(aes(y = incidence), na.rm = TRUE)
  }
  g <- g + scale_x_continuous("Week", breaks = incidence_data[label_index, "t"],
                              labels = incidence_data[label_index, "time_name"]) +
    theme(axis.text.x= element_text(angle = 90),
          legend.position = "none") 
  if(log_scale) {
    g <- g + scale_y_log10("incidence")
  }
  # if the data frame contains a prediction, plot it in blue
  if("prediction" %in% colnames(incidence_data)) {
    g <- g + geom_point(aes(y = prediction), color = "blue", na.rm = TRUE)
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
subset_incidence_data <- function(incidence_data, current_week, n_weeks_prior) {
  end_row <- extract_week_index(current_week, incidence_data$time_name)
  start_row <- end_row - n_weeks_prior
  incidence_data <- incidence_data[seq(start_row, end_row),]
  return(incidence_data)
}

#' extract index of week number from week names
#' 
#' \code{extract_week_index} extract index of week number from week names
#' 
#' @param week_no numeric vector of length 1: week number of the week of interest
#' @param time_name character vector with names of weeks
#' @return numeric vector of length 1: index
#' @export
extract_week_index <- function(week_no, time_name) {
  week_string <- as.character(week_no)
  if(week_no < 10) {
    week_string <- paste0("0", week_string)
  }
  week_index <- grep(paste0("-", week_string), time_name)
  return(week_index)
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
linear_regression <- function(incidence_data,
                              current_week = 47,
                              n_weeks_prior = 4,
                              log_transform = TRUE) {
  
  # extract the incidence data from n weeks prior to now
  incidence_data <- subset_incidence_data(incidence_data, current_week,
                                          n_weeks_prior)
  
  # log transform the incidence if necessary
  if(log_transform) {
    incidence_data$incidence <- log(incidence_data$incidence)
  }
  
  # perform linear regression
  linear_regression_output <- lm(incidence ~ t, data = incidence_data)
  
  # returns intercept and gradient, plus other statistics
  return(linear_regression_output)
}

#' extract predicted points from linear regression object
#' 
#' \code{extract_predicted_points} takes the linear regression object
#' returned by \code{linear_regression} and extracts model predictions
#' 
#' @param lm_output linear regression object returned by \code{linear_regression}
#' @param incidence_data incidence data extracted by \code{extract_incidence}
#' @param weeks_ahead numeric vector of length 1: number of weeks to predict ahead
#' @param log_transform logical vector of length 1: if TRUE, linear 
#' regression was performed on log transformed data, otherwise linear regression
#' was performed on original data
#' @return data frame with predicted points and original data
#' @export
extract_predicted_points <- function(lm_output, incidence_data, weeks_ahead,
                                      log_transform) {
  
  # extract the intercept and gradient from the linear regression output
  intercept <- lm_output$coefficients[["(Intercept)"]]
  gradient <- lm_output$coefficients[["t"]]
  
  # find the current week (last timepoint used for predicting)
  last_timepoint_used <- lm_output$model$t[length(lm_output$model$t)]
  if(missing(weeks_ahead)) {
    weeks_ahead <- nrow(incidence_data) - last_timepoint_used
  }
  # find the times for which we did the prediction
  times_to_predict <- seq_len(weeks_ahead) + last_timepoint_used
  predicted_points <- intercept + gradient * times_to_predict
  if(log_transform) {
    predicted_points <- exp(predicted_points)
  }
  
  # make a column in the data frame recording which time points were used to predict
  time_used_to_predict <- rep(FALSE, nrow(incidence_data))
  time_used_to_predict[lm_output$model$t] <- TRUE
  incidence_data$time_used_to_predict <- time_used_to_predict
  
  # paste the predicted data points at the right times into the data frame
  incidence_data$prediction <- rep(NA, nrow(incidence_data))
  incidence_data$prediction[last_timepoint_used + seq_len(weeks_ahead)] <- predicted_points
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
calc_log_likelihood <- function(data, model_prediction) {
  # check that the number of data points is the same as the number of 
  # prediction points
  stopifnot(length(data) == length(model_prediction) &&
              identical(data, round(data)))
  
  # exclude na data points
  is_na_data <- which(is.na(data))
  if(length(is_na_data) > 0) {
    data <- data[-is_na_data]
    model_prediction <- model_prediction[-is_na_data]
  }
  
  # to avoid extremly low log likelihoods, set predicted incidence to 1 if it
  # is below 1
  model_prediction[model_prediction < 1] <- 1
  # actual counts = reported counts / reporting rate -- correct for this
  # when calculating likelihood
  reporting_rate <- 0.006
  # log likelihood assuming Poisson-distributed errors
  log_likelihood <- sum(dpois(round(data / reporting_rate), model_prediction/ reporting_rate, log = TRUE))
  return(log_likelihood)
}

#' construct a likelihood profile for R_0 given data
#' 
#' \code{likelihood_profile_seir} calculates the log likelihood of the model 
#' given the data, for a range of R_0. The other parameters are fixed
#' 
#' @param incidence_data incidence data extracted by \code{extract_incidence}
#' @param current_week numeric vector of length 1: week number of the current week
#' @param starting_week numeric vector of length 1: 
#' guess for week number when the epidemic started. Use data from starting week
#' to current week to predict
#' @param R_0_min numeric vector of length 1: lower bound of R_0 values over
#' which to construct likelihood profile
#' @param R_0_max numeric vector of length 1: upper bound of R_0 values over
#' which to construct likelihood profile
#' @return a data frame with the following columns:
#' R_0_vec: numeric vector of R_0 values which we're scanning over
#' log_likelihood_vec: log likelihood for those values of R_0
#' @export
likelihood_profile_seir <- function(incidence_data, current_week, starting_week, R_0_min, R_0_max) {
  
  # extract incidence data for weeks used to predict
  starting_week_index <- extract_week_index(starting_week, incidence_data$time_name)
  current_week_index <- extract_week_index(current_week, incidence_data$time_name)
  
  n_weeks_prior <- current_week_index - starting_week_index
  
  incidence_data <- subset_incidence_data(incidence_data, current_week,
                                          n_weeks_prior)
  
  # a function which solves the SEIR model and calculates the log likelihood
  # for a given value of R_0
  evaluate_model_and_calc_log_likelihood <- function(R_0) {
    model_prediction <- solve_seir_wrapper(R_0, n_weeks_prior)
    log_likelihood <- calc_log_likelihood(incidence_data$incidence, model_prediction)
    return(log_likelihood)
  }
  
  # solve the SEIR model and calculate the log likelihood for 
  # 100 evenly spaced values of R_0
  R_0_vec <- seq(R_0_min, R_0_max, length.out = 100)
  log_likelihood_vec <- vapply(R_0_vec, evaluate_model_and_calc_log_likelihood, numeric(1))
  likelihood_profile_output <- data.frame(R_0_vec = R_0_vec, 
                                          log_likelihood_vec = log_likelihood_vec)
  return(likelihood_profile_output)
}

#' plot a likelihood profile for R_0 given data
#' 
#' \code{plot_likelihood_profile} plots the log likelihood of the model 
#' given the data, for a range of R_0. The other parameters are fixed
#' 
#' @param likelihood_profile_output data frame returned by 
#' \code{likelihood_profile_seir}
#' @return ggplot object
#' @export
plot_likelihood_profile <- function(likelihood_profile_output) {
  ggplot(likelihood_profile_output, aes(x = R_0_vec, y = log_likelihood_vec)) +
    geom_line() + xlab ("R_0") + ylab("log likelihood")
}

#' fit the SEIR model to data
#'
#' \code{fit_seir} fits an SEIR model to incidence data
#'
#' @param incidence_data data frame extracted by \code{extract_incidence}
#' @param current_week numeric vector of length 1: week number of the current week
#' @param starting_week numeric vector of length 1: 
#' guess for week number when the epidemic started. Use data from starting week
#' to current week to predict
#' @param R_0_min numeric vector of length 1: lower bound of R_0 values over
#' which to search
#' @param R_0_max numeric vector of length 1: upper bound of R_0 values over
#' which to search
#' @return the value of R_0 with the maximum likelihood
#' @export
fit_seir <- function(incidence_data, current_week, starting_week, R_0_min, R_0_max) {
  
  # extract incidence data for weeks used to predict
  starting_week_index <- extract_week_index(starting_week, incidence_data$time_name)
  current_week_index <- extract_week_index(current_week, incidence_data$time_name)
  
  n_weeks_prior <- current_week_index - starting_week_index
  
  incidence_data <- subset_incidence_data(incidence_data, current_week,
                                          n_weeks_prior)
  
  # a function which solves the SEIR model and calculates the log likelihood
  # for a given value of R_0
  evaluate_model_and_calc_log_likelihood <- function(R_0) {
    model_prediction <- solve_seir_wrapper(R_0, n_weeks_prior)
    log_likelihood <- calc_log_likelihood(incidence_data$incidence, model_prediction)
    return(log_likelihood)
  }
  
  # find the value of R_0 which maximises the log likelihood.
  # NOte that the optimise function tries to minimise, so we need the minus sign.
  max_likelihood_output <- optimise(function(R_0) -evaluate_model_and_calc_log_likelihood(R_0), c(R_0_min, R_0_max))
  R_0 <- max_likelihood_output$minimum
  return(R_0)
}

#' extract predicted points for SEIR model
#' 
#' \code{extract_predicted_points_seir} takes the R_0
#' returned by \code{fit_seir} and extracts model predictions
#' 
#' @param R_0 numeric vector of length 1: basic reproduction number
#' @param incidence_data data frame extracted by \code{extract_incidence}
#' @param starting_week numeric vector of length 1: 
#' guess for week number when the epidemic started. Use data from starting week
#' to current week to predict
#' @param current_week numeric vector of length 1: week number of the current week
#' @param weeks_ahead numeric vector of length 1: number of weeks to predict ahead
#' @return data frame with predicted points and original data
#' @export
#' 
extract_predicted_points_seir <- function(R_0, incidence_data, 
                                           current_week,
                                           starting_week, 
                                           weeks_ahead) {
  # find the times for which we did the prediction
  starting_week_index <- extract_week_index(starting_week, incidence_data$time_name)
  current_week_index <- extract_week_index(current_week, incidence_data$time_name)
  
  # make a column in the data frame recording which time points were used to predict
  time_used_to_predict <- rep(FALSE, nrow(incidence_data))

  time_used_to_predict[seq(starting_week_index, current_week_index)] <- TRUE
  incidence_data$time_used_to_predict <- time_used_to_predict
  
  # solve the SEIR model from the starting week
  n_weeks_prior <- current_week_index - starting_week_index
  predicted_points <- solve_seir_wrapper(R_0, weeks_ahead + n_weeks_prior)
  # exclude points before current week from prediction
  predicted_points <- predicted_points[-seq_len(n_weeks_prior + 1)]
  
  # paste the predicted data points at the right times into the data frame
  incidence_data$prediction <- rep(NA, nrow(incidence_data))
  incidence_data$prediction[current_week_index + seq_along(predicted_points)] <- predicted_points
  return(incidence_data)
}


#' solve SEIR model for incidence
#' 
#' \code{solve_seir_model} solves the SEIR model for the weekly incidence
#' 
#' @param R_0 numeric vector of length 1: basic reproduction number
#' @param latent_period numeric vector of length 1: latent period
#' @param infectious_period numeric vector of length 1: infectious period
#' @param N numeric vector of length 1: population size
#' @param I_0 numeric vector of length 1: initial number of infectious individuals
#' @param n_weeks numeric vector of length 1: number of weeks for which to solve the model
#' @return data frame with predicted points and original data
#' @importFrom deSolve ode
#' @export
#' 
solve_seir_model <- function(R_0, latent_period, infectious_period, N, I_0, n_weeks) {
  # a function which returns the rate of change of each compartment, and the
  # rate of change of the cumulative incidence
  seir_model <- function(time, x, params) {
    with(as.list(c(x, params)), {
      dS <- -beta * S * I / N
      dE <- beta* S * I / N - E / latent_period
      dI <- E / latent_period - I / infectious_period
      dcumincidence <- E / latent_period
      list(c(dS, dE, dI, dcumincidence))
    })
  }
  
  # initial values
  iv <- c(S = N - I_0, E = 0, I = I_0, cumincidence = 0)
  # times for which to solve the differential equations
  times <- (seq_len(n_weeks + 1) - 1) * 7
  # specify model parameters
  params <- c(beta = R_0 / infectious_period, 
              latent_period = latent_period, 
              infectious_period = infectious_period, 
              N = N)
  out <- deSolve::ode(iv, times, seir_model,params)
  # find the incidence from the cumulative incidence
  incidence <- c(0, diff(out[, "cumincidence"]))
  return(incidence)
}

#' solve SEIR model for incidence, specifying default parameters except for 
#' \code{R_0} and \code{n_weeks}
#' 
#' \code{solve_seir_wrapper} solves the SEIR model for the weekly incidence,
#' specifying default parameters except for \code{R_0} and \code{n_weeks}
#' 
#' @param R_0 numeric vector of length 1: basic reproduction number
#' @param n_weeks numeric vector of length 1: number of weeks for which to solve the model
#' @return data frame with predicted points and original data
#' @export
#' 
solve_seir_wrapper <- function(R_0, n_weeks) {
  # specify default parameters
  latent_period <- 1.6
  infectious_period <- 1
  N <- 8.5e6
  I_0 <- 100
  reporting_rate <- 0.006
  # solve ODEs
  solve_seir_model(R_0, latent_period, infectious_period, N, I_0, n_weeks) * reporting_rate
}

#' calculate proportion of prediction points which are within threshold
#' 
#' \code{calc_prediction_accuracy} calculates the proportion of prediction points 
#' which are within a given percentage threshold of the observed incidence
#' 
#' @param prediction_df prediction data frame returned by 
#' \code{extract_predicted_points} or \code{extract_predicted_points_seir}
#' @param tolerance numeric vector of length 1: parameter specifying range
#' within which the predicted incidence should fall.  For example, if 
#' tolerance = 0.25, a predicted point is deemed to be accurate if the 
#' predicted incidence is between 75\% and 125\% of the observed incidence.
#' @return numeric vector of length 1: the proportion of prediction points 
#' which are within a given percentage threshold of the observed incidence
#' @export
#' 
calc_prediction_accuracy <- function(prediction_df, tolerance = 0.25) {
    ## number of predicted points within threshold
    n_accurate_points <- sum(prediction_df$prediction < prediction_df$incidence * (1 + tolerance) &
                           prediction_df$prediction > prediction_df$incidence * (1 - tolerance), 
                           na.rm = TRUE)
    ## number of predicted points
    n_points <- sum(!is.na(prediction_df$prediction) & !is.na(prediction_df$incidence))
    proportion <- n_accurate_points / n_points
    return(proportion)
}

#' extract subset of prediction data frame containing predictions
#' #' 
#' \code{calc_prediction_accuracy} calculates the proportion of prediction points 
#' which are within a given percentage threshold of the observed incidence
#' 
#' @param prediction_df prediction data frame returned by 
#' \code{extract_predicted_points} or \code{extract_predicted_points_seir}
#' @return data frame: incidence and prediction for time points where data is 
#' available
#' @export
#' 
extract_subset_prediction <- function(prediction_df) {
  prediction_df[!is.na(prediction_df$incidence) & !is.na(prediction_df$prediction),
              -which((colnames(prediction_df) == "time_used_to_predict"))]
}
