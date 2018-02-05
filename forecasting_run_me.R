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

# extract forecasted points from linear regression -- forecast 10 weeks ahead
forecast_df <- extract_forecasted_points(linear_regression_output, 
                                         incidence_data,
                                         weeks_ahead = 10,
                                         log_transform = FALSE)

# plot forecast
linear_regression_4_weeks_figure <- plot_incidence(forecast_df)

# perform linear regression starting from calendar week 47,
# using data from 8 weeks prior
linear_regression_output <- linear_regression(incidence_data,
                                              current_week = 47,
                                              n_weeks_prior = 8,
                                              log_transform = FALSE)

# extract forecasted points from linear regression, forecasting 10 weeks ahead
forecast_df <- extract_forecasted_points(linear_regression_output, 
                                         incidence_data,
                                         weeks_ahead = 10,
                                         log_transform = FALSE)

# plot forecast
plot_incidence(forecast_df)

# perform linear regression starting from calendar week 47,
# using data from 8 weeks prior, log transforming data
linear_regression_output <- linear_regression(incidence_data,
                                              current_week = 47,
                                              n_weeks_prior = 8,
                                              log_transform = TRUE)

# extract forecasted points from linear regression, forecasting 10 weeks ahead
forecast_df <- extract_forecasted_points(linear_regression_output, 
                                         incidence_data,
                                         weeks_ahead = 10,
                                         log_transform = TRUE)

# plot forecast
plot_incidence(forecast_df)

# fit seir model
# first, plot likelihood profile to determine reasonable bounds for maximum likelihood algorithm
current_week <- 48
starting_week <- 37
R_0_min <- 1
R_0_max <- 5
likelihood_profile_output <- likelihood_profile_seir(incidence_data, current_week, starting_week, R_0_min, R_0_max)
plot_likelihood_profile(likelihood_profile_output)

# choose reasonable bounds for maximum likelihood algorithm, then fit R_0
R_0_min <- 1
R_0_max <- 2
R_0 <- fit_seir(incidence_data, current_week, starting_week, R_0_min, R_0_max)

# forecast using found value of R_0
forecast_df <- extract_forecasted_points_seir(R_0, incidence_data, 
                                           current_week,
                                           starting_week, 
                                           weeks_ahead = 20)
# plot forecast
plot_incidence(forecast_df)