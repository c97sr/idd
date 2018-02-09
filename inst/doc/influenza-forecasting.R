## ---- echo = FALSE-------------------------------------------------------
library(idd)

## ------------------------------------------------------------------------
data("fluIliCountryData")
plot_incidence_all(fluIliCountryData, "ISR")

## ----plot incidence------------------------------------------------------
incidence_data <- extract_incidence(fluIliCountryData, 
                    country = "ISR", 
                    year = 2016)

plot_incidence(incidence_data, log_scale = FALSE)

## ------------------------------------------------------------------------
linear_regression_output <- linear_regression(incidence_data,
                              current_week = 49,
                              n_weeks_prior = 4,
                              log_transform = FALSE)

## ----eval = FALSE--------------------------------------------------------
#  linear_regression_output

## ------------------------------------------------------------------------
prediction_df <- extract_predicted_points(linear_regression_output, 
                                         incidence_data,
                                         weeks_ahead = 8,
                                         log_transform = FALSE)

## ----eval = FALSE--------------------------------------------------------
#  plot_incidence(prediction_df)

## ----eval = FALSE--------------------------------------------------------
#  extract_subset_prediction(prediction_df)

## ------------------------------------------------------------------------
linear_regression_output <- linear_regression(incidence_data,
                              current_week = 49,
                              n_weeks_prior = 4,
                              log_transform = TRUE)

prediction_df <- extract_predicted_points(linear_regression_output,
                                         incidence_data,
                                         weeks_ahead = 8,
                                         log_transform = TRUE)

## ----eval = FALSE--------------------------------------------------------
#  plot_incidence(prediction_df)

## ----eval = FALSE--------------------------------------------------------
#  calc_prediction_accuracy(prediction_df)

## ---- echo = FALSE-------------------------------------------------------
params_df <- data.frame(Parameter = c("$R_0 = \\beta \\tau_I$",
                                      "$\\tau_E$",
                                      "$\\tau_I$",
                                      "$N$",
                                      "$I_0$",
                                      "$r$"),
                        Description = c("Basic reproduction number",
                                        "Latent period",
                                        "Infectious period",
                                        "Population size",
                                        "Initial number of infectious individuals",
                                        "Proportion of cases which are reported"),
                                      Value = c("varied",
                                                "1.6 day$^{-1}$",
                                                "1 day$^{-1}$",
                                                "$8.5 \\times 10^6$",
                                                "100",
                                                "0.006"))
knitr::kable(params_df)

## ------------------------------------------------------------------------
current_week <- 49
starting_week <- 37 # our guess for the start time of the epidemic. We assume that
# the epidemic starts with 100 infectious individuals at this time.
R_0_min <- 1
R_0_max <- 5
likelihood_profile_output <- likelihood_profile_seir(incidence_data, current_week, starting_week, R_0_min, R_0_max)

## ----eval = FALSE--------------------------------------------------------
#  plot_likelihood_profile(likelihood_profile_output)

## ------------------------------------------------------------------------
R_0_min <- 1
R_0_max <- 2
R_0 <- fit_seir(incidence_data, current_week, starting_week, R_0_min, R_0_max)

## ----eval = FALSE--------------------------------------------------------
#  R_0

## ------------------------------------------------------------------------
# predict using found value of R_0
prediction_df <- extract_predicted_points_seir(R_0, incidence_data, 
                                           current_week,
                                           starting_week, 
                                           weeks_ahead = 8)

## ----eval = FALSE--------------------------------------------------------
#  # plot prediction
#  plot_incidence(prediction_df)
#  calc_prediction_accuracy(prediction_df)

