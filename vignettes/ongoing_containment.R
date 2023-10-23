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
#' Needs
#' - local install of packages for idd
#' - local install of dplyr
#'
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
devtools::load_all()

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

#' Load a database of social mixing patterns
library(socialmixr)
data(polymod)
age_bounds <- c(0,5,11,18,30)
nACs <- length(age_bounds)
mixr_polyuk <- contact_matrix(
    polymod,
    countries="Great Britain",
    age.limits = age_bounds,
    symmetric=TRUE)

#' Extract basic contact matrices
matpolyuk <- mixr_polyuk$matrix
poppolyuk <- mixr_polyuk$demography$population

#' ## Define susceptibiity profiles and school closure mixing
#'
#' Establish some baseline values for susceptibility
susBase <- rep(1,nACs)
susLowKids <- susBase
susLowKids[1:4] <- 0.333

#' Run two test cases for the R0 parameterization
yA <- cov_hybrid(
    vecN=poppolyuk,
    R0=1.01,
    matCt=matpolyuk,
    vecTcal=seq(0,20,0.001),
    vecInfNess=rep(1,nACs),
    vecSusTy=susLowKids,
    deterministic=TRUE)
yB <- cov_hybrid(
    vecN=poppolyuk,
    R0=0.99,
    matCt=matpolyuk,
    vecTcal=seq(0,20,0.001),
    vecInfNess=rep(1,nACs),
    vecSusTy=susLowKids,
    deterministic=TRUE)

plot(rowSums(yA$inf[,,1]),log="y",col="red",ylim=c(0.1,5))
points(rowSums(yB$inf[,,1]),col="blue")

#' Check to see if attack rates are about right with uniform mixing and
#' uniform susceptibility and a perfectly balanced population. Seems to be OK
#' to within half a percent.
yC <- cov_hybrid(
    R0=1.8)
a <- sum(yC$inf[,,1])/sum(yC$pop)
epsilon <- (a + exp(-1.8*a) - 1)/a
epsilon

#' Work to load up excess mortality data. This will go futehr up in the final
#' version of thius
## webfn <- "https://raw.githubusercontent.com/owid/covid-19-data/4132d2c726fa406d5ff934b4f6411f463992a07d/public/data/excess_mortality/excess_mortality.csv"
## x <- read.csv(webfn)
## write.csv(x,file=localfn,row.names = FALSE)
localfn <- "C:/Users/sr/Documents/tmpdata/excess_death.csv"
y <- read.csv(localfn)
dim(y)
names(y)

## all(dim(x) == dim(y))
## all(names(x) == names(y))

df <- y %>%
  dplyr::filter(time_unit=="weekly", date < as.Date("2021-02-28")) %>%
  dplyr::select("location", "cum_excess_per_million_proj_all_ages","date") %>%
  tidyr::spread(location,cum_excess_per_million_proj_all_ages)

plot(df$`United Kingdom`,type="l",col="blue",lwd=3)
points(df$Australia,type="l",col="red",lwd=3)

#' Now make some simple plots of cumulative covid deaths from the model
