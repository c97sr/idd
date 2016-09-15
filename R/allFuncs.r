# Copyright Steven Riley (sr@stevenriley.net)
#
# This file is part of the library idd.
#
# idd is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This work is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this idsource.  If not, see <http://www.gnu.org/licenses/>.

#' Solves the SEIR model
#'
#' Creates both deterministic and stochastic solutions to a common
#' epidemic model.
#'
#' @param De duration latent
#' @param Tg generation time
#' @param R0 Basic reproductive number
#' @param N size of the population
#' @param I0 number of initial infectious individuals
#' @param dt time step
#' @param R1 second value for the reproductive nuber
#' @param t1 time of change to the second R value
#' @param R2 third value for the reproductive number
#' @param t2 time of change to the third reproductive number
#' @param noTimeSteps number of timesteps
#' @param trickle expected number of trickle infections per day
#' @param trickleStart start time for trickle infection
#' @param reprate reporting rate
#' @param A coefficient of seasonaliity
#' @param noReals number of realizations, NOT YET IMPLEMENTED
#' @param deterministic mode of solutions
#'
#' @return A list of two elements. The first is the incidence of infection
#' and the second is the timepoints to which the incidence refers.
#'
#' @examples comp_seir()
#' @examples comp_seir(R0=0.999)
comp_seir <- function(
  De=7,
  Tg=15,
  R0=1.8,
  N=1000000,
  I0=10,
  trickle=0,
  trickleStart=0,
  reprate=1,
  dt=1,
  R1=1.0,
  t1=999,
  R2=1.0,
  t2=9999,
  noTimeSteps=360,
  noReals=1,
  A=0,
  deterministic=TRUE
) {

  # Define the variable to be returned
  rtn_inf_inc <- vector(mode="numeric",length=noTimeSteps+1)
  rtn_beta <- vector(mode="numeric",length=noTimeSteps+1)

  # Define derived parameters
  Di <- Tg - De

  # Set up state variables
  state_S <- N - I0
  state_E <- 0
  state_I <- I0
  state_R <- 0

  # Assign the seed infection to time 0
  rtn_inf_inc[1] <- I0

  # Start the main loop
  for (i in 1:noTimeSteps) {

    # Trickle indicator
    if (i*dt >= trickleStart) {
      indTrick <- 1
    } else {
      indTrick <- 0
    }

    # Set Beta
    if (i < t1/dt) {
      beta = R0 / Di * (1+A*sin(i*dt*2*pi/360))
    } else if (i < t2/dt) {
      beta = R1 / Di * (1+A*sin(i*dt*2*pi/360))
    } else {
      beta = R2 / Di * (1+A*sin(i*dt*2*pi/360))
    }

    rtn_beta[i] <- beta

    # Calculate hazard rates
    # browser()
    hazInf     <- beta * state_I / N
    hazBecInf  <- 1 / De
    hazRec     <- 1 / Di
    expTrickle <- indTrick * trickle

    # Calculate probabilities
    pInf    <- 1 - exp(-dt * hazInf)
    pBecInf <- 1 - exp(-dt * hazBecInf)
    pRec    <- 1 - exp(-dt * hazRec)

    # Calculate the expected or simulated number of events
    if (deterministic) {
      nInf 	<- state_S * pInf
      nBecInf   <- state_E * pBecInf
      nRec 	<- state_I * pRec
      nTrickle  <- expTrickle
      rtn_inf_inc[i+1] <- nInf*reprate
    } else {
      nInf <- rbinom(1,state_S,pInf)
      nBecInf <- rbinom(1,state_E,pBecInf)
      nRec <- rbinom(1,state_I,pRec)
      nTrickle <- rpois(1,expTrickle)
      rtn_inf_inc[i+1] <- rbinom(1,nInf,reprate)
    }

    # Update the state variables
    state_S <- state_S - nInf - nTrickle
    state_E <- state_E + nInf - nBecInf + nTrickle
    state_I <- state_I + nBecInf - nRec
    state_R <- state_R - nRec

    # Update the output variables


  } # End main loop

  # Return incidence
  list(
      inf_inc = rtn_inf_inc,
      time=(0:noTimeSteps)*dt,
      beta=rtn_beta)

}

#' Illustration of possible zika invasion dynamics
#'
#' Simulates a stochastic compartmental model to reproduce an
#' illustration of how two very similar communities might have
#' very different initial epidemics.
#'
#' @param its number of iterations for the second population
#' @param file name of file to write pdf or "screen"
#' @param N size of population
#'
#' @example zika_invasion(its=100,file="~/tmp/zika_example.pdf")
#' @example zika_invasion(its=100,file="screen",N=1000000)
zika_invasion <- function(
    its=10,
    file="screen",
    N=1000000,
    seed=19283719,
    highlights=c(1,2,3),
    highcolors=c("red","magenta","cyan"),
    maxts=360*2,
    R0axes=0:4,
    incaxes=seq(0,3.5,0.5)) {
  set.seed(seed)
  nosens <- length(highlights)
  senstab <- array(dim=c(its,maxts+1))
  transval <- vector(mode="numeric",length=maxts+1)
  De <- 11
  Tg <- 16
  Di <- Tg - De
  r1 <- comp_seir(
    De=De, Tg=Tg, A=0.9, I0 = 10, trickle=0/7,
    noTimeSteps=maxts,deterministic=FALSE,R0=2,N=N
  )
  if (file != "screen") pdf_fig_proof(file=file,pw=9)
  plot(r1$beta*Di,type="n",ylim=c(min(R0axes),max(R0axes)),axes=FALSE,xlab="",ylab="")
  axis(4,at=R0axes,col=grey.colors(10)[5],las=1)
  polygon(
      c(0:maxts,maxts:0),c(r1$beta*Di,rep(0,maxts+1)),
      col=grey.colors(10)[5],border=NA)
  par(new=TRUE)
  plot(r1$inf_inc/N*100,type="l",col="black",
       xlab="Day", ylab="Incidence",
       lwd=2,axes=FALSE,ylim=c(min(incaxes),max(incaxes)))
  axis(1)
  axis(2,las=1,at=incaxes)
  for (i in 1:its) {
    r2 <- comp_seir(
      De=De, Tg=Tg, A=0.9, I0 = 0, trickle=3/7,
      noTimeSteps=maxts,deterministic=FALSE,
      trickleStart=70,R0=2,N=N)
      points(r2$inf_inc/N*100,type="l",col=grey.colors(10)[9])
      senstab[i,] <- r2$inf_inc/N*100
  }
  for (i in 1:length(highlights)) {
      points(senstab[highlights[i],],type="l",col=highcolors[i],lwd=2)
  }
  if (file != "screen") dev.off()
}

#' Solves the basic individual-based serial interval model
#'
#' Simulates a sequence of individual infection events.
#' Infections are drawn from simple offspring and serial
#' interval distributions.
#'
#' @param Tg generation time
#' @param R0 the basic reproductive number
#' @param Nmax the maximum number if infections in a cluster
#'
#' @return Returns a list of times of infection and indices.
#' The position in the list of infection times corresponds to
#' the index given in the second argument.
#'
#' @examples ind.tau()
ind_tau <- function(
	Tg=5,
	R0=2.6,
	Nmax=30
) {

  # Define the outputs to be returned
  rtn_times <- numeric(mode="numeric",length=Nmax)
  rtn_infectors <- numeric(mode="numeric",length=Nmax)

  # Initiate the
  rtn_times[1] <- 0
  rtn_infectors[1] <- 9999

  completed_infector <- 0
  current_infector <- 1

  # Main loop to go through every infector
  while (completed_infector <= Nmax) {

    # Not yet completed
    completed_infector <- completed_infector + 1
  }

  # Return list of infectors and times of infection
  list(t=rtn_times, i = rtn_infectors )

}

#' Translates a number from a log or linear scale to a unit scale
#'
#' @param y is the number on the non-unit scale
#' @param min is the assumed minimum value of the non-unit scale
#' @param max is the assumed maximum value of the non-unit scale
#' @param logbase is the base of the log scale if used
#' @param log is a boolean for whether the scale is log or linear
#' Returns the value on the unit scale
to_unit_scale <- function(y,min=1,max=100,logbase=10,logflag=FALSE) {
  if (logflag) {
    rtn <- (log(y,logbase)-log(min,logbase))/(log(max,logbase)-log(min,logbase))
  } else {
    rtn <- (y-min)/(max-min)
  }
  rtn
}


#' Proposes parameter updates for an MCMC algorithm
#'
#' Uses either linear or log-linear random walks to propose
#' multivariate jumps in parameter space.
#'
#' @param pt a table of parameter names, values and min and maxes.
#' See Details.
#' @param fmask a vector of strings of the names of parameters within the table
#' to be fitted.
#' The default value is to assume all the parameters are to be fitted
#'
#' @details The
#'
#' @return A vector of proposed parameter values. The vector is the
#' same length as ptab so it copies in the values of the parameters
#' that are not being updated.
#'
#' @examples sirTab <- data.frame(val=c(1.8,2.6), max=c(5.0,7.0),
#' min=c(0.9,5.0), step=c(0.1,0.1), log=c(FALSE,FALSE),
#' row.names = c("R0","Tg"))
#' mh_update(sirTab)
#'
mh_update <- function(pt,fmask=1:(dim(pt)[1])) {

	bps <- pt[,"val"]
	rtn <- bps
	nofit <- length(rtn)

	for (i in 1:nofit) {

		# Set up and transform to unit scale
		rv <- runif(1)
		rv <- (rv-0.5)* pt[i,"step"]
		x <- bps[fmask[i]]
		x <- to_unit_scale(x,min=pt[i,"min"],max=pt[i,"max"],logflag=pt[i,"log"])
		x <- x + rv

		# Bouncing boundary conditons
		# if (x < 0) x <- -x
		# if (x > 1) x <- 2 - x

		# Cyclical boundary conditions
		if (x < 0) x <- 1 + x
		if (x > 1) x <- x - 1

		# Test for errors and return to originl scales
		if (x < 0 || x > 1) stop("problem here")
		rtn[fmask[i]] <- from_unit_scale(x,min=pt[i,"min"],max=pt[i,"max"],logflag=pt[i,"log"])

	}

	rtn
}

#' Solves the classic relationship between R0 and cumulative attack rate
#'
#' @param R0 basic reproductive number
#'
#' @example attack_rate(1.8)
#'
#' @example attack_rate(4.0)
#' @example attack_rate(3.0)
attack_rate <- function(R0=1.8) {
  f <- function(a){1-exp(-R0*a)-a}
  tmp <- uniroot(f,c(0.000001,0.9999999),tol=0.001)
  tmp$root
}
