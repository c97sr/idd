#' Simulation model example data frame
#'
#' This object contains synthetic data on the location of farms in Kent, UK.
#' These are not real locations but they have a similar density distribution.
#' These data are intended to be used in a teaching example about
#' simulation models of infectious disease.
#'
#' @name smeDf
#' 
#' @docType data
#' 
#' @author Steven Riley \email{sr@stevenriley.net}
#'
#' @references These synthetic data have been constructed to have a similar
#' distribution of spatial density as the actual distribution in this part of
#' the UK.
#' \url{https://www.gov.uk/government/statistical-data-sets/structure-of-the-livestock-industry-in-england-at-december}
NULL

#' Forecast model example data frame
#'
#' These are influenza-like-illness data gathered from many countries across the 
#' world and curated by the WHO as FluID.
#'
#' @name fluIliCountryData
#' 
#' @docType data
#' 
#' @author Steven Riley \email{sr@stevenriley.net}
#'
#' @references These incidence data are from publically available WHO syndromic data
#' \url{https://www.who.org}
NULL

#' Economist weekly COVID-19 excess mortality data for 88 weeks from week 1 2020
#'
#' Downloaded from our world in data
#'
#' @format
#' A data 88 rows and 238 columns
#' \describe{
#'   \item{country}{Country names are row labels}
#'   \item{week}{Weeks starting from week 1 2020 are row labels}
#' }
#' @source <https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/excess_mortality/excess_mortality_economist_estimates.csv>
"country_excess_econ"