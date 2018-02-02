#' plots the incidence for a given flu season and a given country
#' 
#' \code{plot_incidence} takes the flu incidence matrix and makes a plot of
#' the incidence for a given flu season and a given country.
#' 
#' @param flu_data matrix with flu incidence data
#' @param country_code character vector of length 1: 3-letter country code
#' @param year numeric vector of length 1: year to plot
#' @param log_scale logical vector of length 1: if TRUE, plot on log scale, 
#' otherwise plot on linear scale
#' @return ggplot object
#' @import ggplot2
#' @export
plot_incidence <- function(flu_data, 
                           country_code = "USA", 
                           year = 2016, 
                           log_scale = TRUE){
  flu_data <- as.data.frame(flu_data)
  row_name_start <- paste0(year, "-27")
  row_name_end <- paste0(year + 1, "-26")
  incidence <- flu_data[seq(which(rownames(flu_data) == row_name_start),
                            which(rownames(flu_data) == row_name_end)), 
                        colnames(flu_data) == country_code]
  time_vec <- seq_along(incidence)
  plot_df <- data.frame(t = time_vec, incidence = incidence)
  g <- ggplot(plot_df, aes(x = t, y = incidence)) + geom_point()
  if(log_scale) {
    g <- g + scale_y_log10("incidence")
  }
  g
}
