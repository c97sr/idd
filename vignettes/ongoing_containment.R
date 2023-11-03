#' # Ongoing containment as a policy objective
#'
#' Script designed to be run with rmarkdown::render(), then run to
#'
#' ## Background
#'
#' This analysis was developed between mid-February and early morning of March
#' 10th 2020. Initially, it was used to give early estimates of the possible
#' impact of school closures.
#'
#' Next steps
#'
#' * add the alternative excess mortality data
#' * add the processed excess mortality as a data object to the package
#' * fix the cumulative excess mortality so it works even when the ICU limit
#' hit
#'
#' ## Load functions and population data
#'
#' First clear memory.
rm(list=ls(all=TRUE))

#' Load functions to run the hybrid model from local packages, from github or directly from local copies
#' if that is how the code has been distributed.
## source("~/git/idd/R/cov_hybrid.r")
library(devtools)
library(dplyr)

#' Install SR's idd package either from github or locally
## install_github("c97sr/idd")
## install("~/gdrive/git/idd", dependencies=FALSE)
## library("idd")

#' Make an initial run for 2 years with an R0 of 2.0 in a UK-like
#' population, with a trickle seed of 10 cases per week and nop
#' reactive behaviour change
dt <- 0.1
pop <- 66000000*c(0.25,0.25,0.25,0.25)
seed <- rep(6000/length(pop),length(pop))
xvals <- seq(0,540,dt)
y1 <- idd::cov_hybrid(
    vecN=pop,
    vecI0=seed,
    R0=2.0,
    Rp=1.4,
    Rover=0.95,
    trig_pres=99999999,
    icu_cap=4357/53000000,
    sevBchange=FALSE,
    vecPS=rep(0.05,length(pop)),
    trickle=10,
    vecTcalc=xvals)
#' Then run with a reactive drop when the virus is known to be present
devtools::load_all()
y2 <- cov_hybrid(
    vecN=pop,
    vecI0=seed,
    R0=2.0,
    Rp=1.4,
    Rover=0.95,
    trig_pres=5,
    icu_cap=4357/53000000,
    sevBchange=FALSE,
    vecPS=rep(0.05,length(pop)),
    trickle=10,
    vecTcalc=xvals)
#' Then run with a reactive drop when the virus is known to be present
y3 <- cov_hybrid(
    vecN=pop,
    vecI0=seed,
    R0=2.0,
    Rp=0.95,
    Rover=0.95,
    trig_pres=5,
    icu_cap=4357/53000000,
    sevBchange=FALSE,
    vecPS=rep(0.05,length(pop)),
    trickle=10,
    vecTcalc=xvals)
#' Then run with a reactive drop when the virus is known to be present
y4 <- cov_hybrid(
    R0=2.0,
    vecN=pop,
    vecI0=seed,
    Rp=1.4,
    Rover=0.99,
    trig_pres=5,
    icu_cap=4357/53000000,
    sevBchange=TRUE,
    vecPS=rep(0.05,length(pop)),
    trickle=10,
    vecTcalc=xvals)
#' Then run with a reactive drop when the virus is known to be present
y4A <- cov_hybrid(
    R0=2.0,
    vecN=pop,
    vecI0=seed,
    Rp=1.4,
    Rover=0.99,
    trig_pres=5,
    icu_cap=2*4357/53000000,
    sevBchange=TRUE,
    vecPS=rep(0.05,length(pop)),
    trickle=10,
    vecTcalc=xvals)
y4B <- cov_hybrid(
    R0=2.0,
    vecN=pop,
    vecI0=seed,
    Rp=1.4,
    Rover=1.1,
    trig_pres=5,
    icu_cap=4357/53000000,
    sevBchange=TRUE,
    vecPS=rep(0.05,length(pop)),
    trickle=10,
    vecTcalc=xvals)

#' Plot the baselines
par(mfrow=c(1,2))
for (s in c("","y")) {
    plot(y1$t/30,rowSums(y1$inf[,,1]/dt+1),
         log=s,col="red",type="l",lwd=2,xlim=c(0,18),
         xlab="Month",ylab="Daily incidence (+1)")
    points(y2$t/30,rowSums(y2$inf[,,1]/dt+1),col="blue",type="l",lwd=2)
    points(y3$t/30,rowSums(y3$inf[,,1]/dt+1),col="green",type="l",lwd=2)
    points(y4$t/30,rowSums(y4$inf[,,1]/dt+1),col="cyan",type="l",lwd=2)
}

#' Plot the sensitvities
par(mfrow=c(1,2))
for (s in c("","y")) {
    plot(y1$t/30,rowSums(y1$inf[,,1]/dt+1),
         log=s,col="gray",type="l",lwd=2,xlim=c(0,18),
         xlab="Month",ylab="Daily incidence (+1)")
    points(y2$t/30,rowSums(y2$inf[,,1]/dt+1),col="gray",type="l",lwd=2)
    points(y3$t/30,rowSums(y3$inf[,,1]/dt+1),col="gray",type="l",lwd=2)
    points(y4$t/30,rowSums(y4$inf[,,1]/dt+1),col="gray",type="l",lwd=2)
    points(y4A$t/30,rowSums(y4A$inf[,,1]/dt+1),col="orange",type="l",lwd=2)
    points(y4B$t/30,rowSums(y4B$inf[,,1]/dt+1),col="magenta",type="l",lwd=2)
}

#' Define the format for the table
tabdisp <- function(x)(formatC(signif(x,3), format="f", big.mark=",", digits=0))

#' Output the total amount of infection
tabdisp(sum(y1$inf[,,1]))
tabdisp(sum(y2$inf[,,1]))
tabdisp(sum(y3$inf[,,1]))
tabdisp(sum(y4$inf[,,1]))
tabdisp(sum(y4A$inf[,,1]))
tabdisp(sum(y4B$inf[,,1]))

#' Output the total number of deaths
tabdisp(y1$sevIn*0.2 + y1$sevOut)
tabdisp(y2$sevIn*0.2 + y2$sevOut)
tabdisp(y3$sevIn*0.2 + y3$sevOut)
tabdisp(y4$sevIn*0.2 + y4$sevOut)
tabdisp(y4A$sevIn*0.2 + y4A$sevOut)
tabdisp(y4B$sevIn*0.2 + y4B$sevOut)

#' Output the proportion of herd immunity acheived at 18 months
round(sum(y1$inf[,,1])/(sum(pop)/2)*100)
round(sum(y2$inf[,,1])/(sum(pop)/2)*100)
round(sum(y3$inf[,,1])/(sum(pop)/2)*100)
round(sum(y4$inf[,,1])/(sum(pop)/2)*100)
round(sum(y4A$inf[,,1])/(sum(pop)/2)*100)
round(sum(y4B$inf[,,1])/(sum(pop)/2)*100)

#' Include some model calibration runs as evidence of correctness because
#' there is a little bit of model here and its import that the NGM calc
#' is correct


#' Work to load up excess mortality data. This will go futehr up in the final
#' version of thius
localfn <- "C:/Users/sr/Documents/tmpdata/excess_death.csv"
localfn_w <- "C:/Users/sr/Documents/tmpdata/excess_death_econ.csv"
## webfn <- "https://raw.githubusercontent.com/owid/covid-19-data/4132d2c726fa406d5ff934b4f6411f463992a07d/public/data/excess_mortality/excess_mortality.csv"
## x <- read.csv(webfn)
## write.csv(x,file=localfn,row.names = FALSE)
wwebfn <- "https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/excess_mortality/excess_mortality_economist_estimates.csv"
w <- read.csv(wwebfn)
write.csv(w,file=localfn_w,row.names=FALSE)

y <- read.csv(localfn)
z <- read.csv(localfn_w)
names(z)
dim(z)

length(table(z$country))

## all(dim(x) == dim(y))
## all(names(x) == names(y))
df <- y %>%
  dplyr::filter(time_unit=="weekly", date < as.Date("2021-08-31")) %>%
  dplyr::select("location", "cum_excess_per_million_proj_all_ages","date") %>%
  tidyr::spread(location,cum_excess_per_million_proj_all_ages)

dim(df)
names(df)

plot(y2$t,rowSums(y2$csevtreat[,,1]),type="l",ylim=c(0,10000))

plot(df$`United Kingdom`,type="l",col="blue",lwd=3)
points(df$Australia,type="l",col="red",lwd=3)

#' Now make some simple plots of cumulative covid deaths from the model
