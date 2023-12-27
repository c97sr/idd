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

#' World Health Organization TB data
#'
#' A subset of data from the World Health Organization Global Tuberculosis
#' Report ...
#'
#' @format ## `who`
#' A data frame with 7,240 rows and 60 columns:
#' \describe{
#'   \item{country}{Country name}
#'   \item{iso2, iso3}{2 & 3 letter ISO country codes}
#'   \item{year}{Year}
#'   ...
#' }
#' @source <https://www.who.int/teams/global-tuberculosis-programme/data>
"country_excess_death"