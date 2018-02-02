rm(list = ls())
devtools::load_all()
data("fluIliCountryData")
incidence_data <- extract_incidence(fluIliCountryData, 
                    country = "USA", 
                    year = 2016)
figure1 <- plot_incidence(incidence_data, log_scale = FALSE)

linear_regression_output <- linear_regression(fluIliCountryData,
                              country_code = "USA",
                              year = 2016,
                              current_week = 47,
                              n_weeks_prior = 4,
                              log_transform = FALSE)

forecast_df <- extract_forecasted_points(linear_regression_output, 
                                         incidence_data,
                                         log_transform = FALSE)
figure2 <- plot_incidence(forecast_df, log_scale = FALSE)