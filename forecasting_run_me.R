rm(list = ls())
devtools::load_all()

# load the incidence data
data("fluIliCountryData")
incidence_data <- extract_incidence(fluIliCountryData, 
                    country = "USA", 
                    year = 2016)

# plot the incidence data
data_figure <- plot_incidence(incidence_data, log_scale = FALSE)

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
linear_regression_8_weeks_figure <- plot_incidence(forecast_df)

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
linear_regression_log_transform_8_weeks_figure <- plot_incidence(forecast_df)