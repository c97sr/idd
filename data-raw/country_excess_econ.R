## code to prepare `DATASET` dataset goes here

library("dplyr")
library("tidyr")

#' Work to load up excess mortality data. This will go futher up in the final
#' version of thius
localfn_w <- "/mnt/c/Users/sriley/Documents/tmpdata/excess_death_econ.csv"
wwebfn <- "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/excess_mortality/excess_mortality_economist_estimates.csv"
w <- read.csv(wwebfn)
write.csv(w,file=localfn_w,row.names=FALSE)

z <- read.csv(localfn_w)
names(z)
dim(z)

length(table(z$country))

country_excess_death <- z %>%
  dplyr::filter(date < as.Date("2021-08-31")) %>%
  dplyr::select("date","country",
               "cumulative_estimated_daily_excess_deaths_per_100k") %>%
  dplyr::rename("cumexcess" =
                  "cumulative_estimated_daily_excess_deaths_per_100k") %>%
  tidyr::spread(country,cumexcess)

dim(country_excess_death)
colnames <- names(country_excess_death)
windex <- which(colnames=="World")
dindex <- which(colnames=="date")
countries <- colnames[-c(windex,dindex)]
ncountries <- length(countries)

usethis::use_data(country_excess_econ, overwrite = TRUE)
