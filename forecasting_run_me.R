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

rm(list = ls())
devtools::load_all()

# load the incidence data
data("fluIliCountryData")
incidence_data <- extract_incidence(fluIliCountryData, 
                    country = "ISR", 
                    year = 2016)

# plot the incidence data
plot_incidence(incidence_data)

# perform linear regression starting from calendar week 47,
# using data from 4 weeks prior
linear_regression_output <- linear_regression(incidence_data,
                              current_week = 47,
                              n_weeks_prior = 4,
                              log_transform = FALSE)

# extract forecasted points from linear regression -- forecast 5 weeks ahead
forecast_df <- extract_forecasted_points(linear_regression_output, 
                                         incidence_data,
                                         weeks_ahead = 5,
                                         log_transform = FALSE)

# plot forecast
linear_regression_4_weeks_figure <- plot_incidence(forecast_df)
print(calc_forecast_accuracy(forecast_df))

# perform linear regression starting from calendar week 47,
# using data from 8 weeks prior
linear_regression_output <- linear_regression(incidence_data,
                                              current_week = 47,
                                              n_weeks_prior = 8,
                                              log_transform = FALSE)

# extract forecasted points from linear regression, forecasting 5 weeks ahead
forecast_df <- extract_forecasted_points(linear_regression_output, 
                                         incidence_data,
                                         weeks_ahead = 5,
                                         log_transform = FALSE)

# plot forecast
plot_incidence(forecast_df)
print(calc_forecast_accuracy(forecast_df))

# Question: does the forecast capture the qualitative trend in the data?
# Question: how does changing the number of previous data points to include 
# change the forecast accuracy?
# Question: how does moving the current week forward change the forecast accuracy?
# (keeping the number of previous data points to include constant)

# perform linear regression starting from calendar week 47,
# using data from 4 weeks prior, log transforming data
linear_regression_output <- linear_regression(incidence_data,
                                              current_week = 47,
                                              n_weeks_prior = 4,
                                              log_transform = TRUE)

# extract forecasted points from linear regression, forecasting 5 weeks ahead
forecast_df <- extract_forecasted_points(linear_regression_output, 
                                         incidence_data,
                                         weeks_ahead = 5,
                                         log_transform = TRUE)

# plot forecast
plot_incidence(forecast_df)
print(calc_forecast_accuracy(forecast_df))

# Question: what are the qualitative changes in the forecast when we log transform the data?
# Question: does the forecast capture the qualitative trend in the data?
# Question: how does changing the number of previous data points to include 
# change the forecast accuracy?
# Question: how does moving the current week forward change the forecast accuracy?
# (keeping the number of previous data points to include constant)
# Question: what happens if we try to forecast many more weeks forward?

# fit seir model
# first, plot likelihood profile to determine reasonable bounds for maximum likelihood algorithm
# the likelihood profile plots the log likelihood over a range of R_0
current_week <- 47
starting_week <- 37 # our guess for the start time of the epidemic. We assume that
# the epidemic starts with 100 infectious individuals at this time.
R_0_min <- 1
R_0_max <- 5
likelihood_profile_output <- likelihood_profile_seir(incidence_data, current_week, starting_week, R_0_min, R_0_max)
plot_likelihood_profile(likelihood_profile_output)

# Question: from looking at the likelihood profile, what is the most likely 
# value of R_0?
# Is this a reasonable value of R_0?

# when we fit the SEIR model to the data, we will need to choose reasonable 
# bounds of R_0 to search over.
# plotting the likelihood profile helps us choose these bounds.
# the bounds should include the peak of the log likelihood which we have seen
# visually.  Ideally it should not include very flat regions of parameter space
# where the fitting algorithm will have trouble finding the maximum 
# (in this case, very high values of R_0)

# choose reasonable bounds for maximum likelihood algorithm, then fit R_0
R_0_min <- 1
R_0_max <- 2
R_0 <- fit_seir(incidence_data, current_week, starting_week, R_0_min, R_0_max)

# forecast using found value of R_0
forecast_df <- extract_forecasted_points_seir(R_0, incidence_data, 
                                           current_week,
                                           starting_week, 
                                           weeks_ahead = 10)
# plot forecast
plot_incidence(forecast_df)
print(calc_forecast_accuracy(forecast_df))

# Question: what features of the incidence curve are captured by the SEIR model
# forecast which were not captured by the linear regression forecast?
# Question: what happens as we advance the current week?
# Question: how sensitive is the forecast to our guess of when the epidemic started?
# Question: what might be causes of error in the forecast?