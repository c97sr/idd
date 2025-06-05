## code to prepare `DATASET` dataset goes here

library("dplyr")
library("tidyr")

## Download the lookup table for three letter country codes
fn_country_lookup <- "https://raw.githubusercontent.com/lukes/ISO-3166-Countries-with-Regional-Codes/master/all/all.csv"
country_lookup <- read.csv(fn_country_lookup)

#' Work to load up excess mortality data. This will go futher up in the final
#' version of thius
## localfn_w <- "/mnt/c/Users/sriley/Documents/tmpdata/excess_death_econ.csv"
## wwebfn <- "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/excess_mortality/excess_mortality_economist_estimates.csv"
z <- read.csv(wwebfn)
## write.csv(w,file=localfn_w,row.names=FALSE)

## z <- read.csv(localfn_w)
names(z)
dim(z)

head(table(z$country))

names(country_lookup)
head(table(country_lookup$name))

length(table(z$country))

country_excess_econ <- z %>%
  dplyr::filter(date < as.Date("2021-08-31")) %>%
  dplyr::select("date","country",
               "cumulative_estimated_daily_excess_deaths_per_100k") %>%
  dplyr::rename("cumexcess" =
                  "cumulative_estimated_daily_excess_deaths_per_100k") %>%
  tidyr::spread(country,cumexcess)

dim(country_excess_econ)
colnames <- names(country_excess_econ)
windex <- which(colnames=="World")
dindex <- which(colnames=="date")
countries <- colnames[-c(windex,dindex)]
ncountries <- length(countries)

## Download quarterly GDP data
## file available here https://data.oecd.org/gdp/quarterly-gdp.htm
## but not available via and API
fn_gdp_data <- "~/gdrive/Briefcase/DP_LIVE_27122023181343752.csv"
gdp_data <- read.csv(fn_gdp_data)
dim(gdp_data)



usethis::use_data(country_excess_econ, overwrite = TRUE)
