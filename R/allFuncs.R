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
#' @param noReals number of realizations, NOT YET IMPLEMENTED
#' @param deterministic mode of solution, NOT YET IMPLEMENTED
#'
#' @return A list of two elements. The first is the incidence of infection and the
#' second is the timepoints to which the incidence refers.
#'
#' @examples
#' comp.seir()
#' comp.seir(R0=0.999)
comp.seir <- function(
  De=1.48,
  Tg=2.6,
  R0=1.8,
  N=6800000,
  I0=10,
  dt=1,
  R1=1.0,
  t1=999,
  R2=1.0,
  t2=9999,
  noTimeSteps=10,
  noReals=1,
  deterministic=TRUE
) {

  # Define the variable to be returned
  rtn_inf_inc <- vector(mode="numeric",length=noTimeSteps+1)

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

    # Set Beta
    if (i < t1*dt) {
      beta = R0 / Di
    } else if (i < t2*dt) {
      beta = R1 / Di
    } else {
      beta = R2 / Di
    }

    # Calculate hazard rates
    hazInf     <- beta * state_I / N
    hazBecInf  <- 1 / De
    hazRec     <- 1 / Di

    # Calculate probabilities
    pInf    <- 1 - exp(-dt * hazInf)
    pBecInf <- 1 - exp(-dt * hazBecInf)
    pRec    <- 1 - exp(-dt * hazRec)

    # Calculate the expected or simulated number of events
    if (deterministic) {
      nInf 	<- state_S * pInf
      nBecInf <- state_E * pBecInf
      nRec 	<- state_I * pRec
    } else {
      nInf <- rbinom(1,state_S,pInf)
      nBecInf <- rbinom(1,state_E,pBecInf)
      nRec <- rbinom(1,state_I,pRec)
    }

    # Update the state variables
    state_S <- state_S - nInf
    state_E <- state_E + nInf - nBecInf
    state_I <- state_I + nBecInf - nRec
    state_R <- state_R - nRec

    # Update the output variables
    rtn_inf_inc[i+1] <- nInf

  } # End main loop

  # Return incidence
  list(inf_inc = rtn_inf_inc, time=(0:noTimeSteps)*dt)

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
#' @examples
#' ind.tau()
ind.tau <- function(
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

    #

    # Not yet completed
    completed_infector <- completed_infector + 1
  }

  # Return list of infectors and times of infection
  list(t=rtn_times, i = rtn_infectors )

}

#' Proposes parameter updates for an MCMC algorithm
#'
#' Uses either linear or log-linear random walks to propose 
#' multivariate jumps in parameter space.
#'
#' @param pt a table of parameter names, values and min and maxes. 
#' See Details.
#' @param fm a vector of strings of the parameters within the table 
#' to be fitted.
#' The default value is to assume all the parameters are to be fitted
#'
#' @details The 
#'
#' @return A vector of proposed parameter values. The vector is the 
#' same length as ptab so it copies in the values of the parameters 
#' that are not being updated.
mh_update <- function(pt,fmask=1:(dim(pt)[1])) {

	bps <- pt[,"val"]
	rtn <- bps
	nofit <- length(rtn)

	for (i in 1:nofit) {

		# Set up and transform to unit scale
		rv <- runif(1)
		rv <- (rv-0.5)* pt[i,"step"]
		x <- bps[fmask[i]]
		x <- SR_to_unit(x,min=pt[i,"min"],max=pt[i,"max"],logflag=pt[i,"log"])
		x <- x + rv

		# Bouncing boundary conditons
		# if (x < 0) x <- -x
		# if (x > 1) x <- 2 - x

		# Cyclical boundary conditions
		if (x < 0) x <- 1 + x
		if (x > 1) x <- x - 1

		# Test for errors and return to originl scales
		if (x < 0 || x > 1) stop("problem here")
		rtn[fmask[i]] <- SR_from_unit(x,min=pt[i,"min"],max=pt[i,"max"],logflag=pt[i,"log"])

	}

	rtn
}

fnProposeParamUpdatesSingle <- function(ptab,index) {

	# Takes biological parameters, a mask which is a list of the parameters
	# being fitted, a table of the fitted parameters and their ranges and
	# the number of parameters to be fitted.
	# Returns a proposed vector of parameters.
	# Comment in and out bouncing and cyclical boundary conditions
	# for the random walk

	bps <- ptab[,"val"]
	rtn <- bps
	if (index < 0 || index > dim(ptab)[1]) {
		stop("index must be less than or equal to number of parameters")
	}

	# Set up and transform to unit scale
	rv <- runif(1)
	rv <- (rv-0.5)* ptab[index,"step"]
	x <- bps[index]
	x <- SR_to_unit(x,min=ptab[index,"min"],max=ptab[index,"max"],logflag=ptab[index,"log"])
	x <- x + rv

	# Cyclical boundary conditions
	if (x < 0) x <- 1 + x
	if (x > 1) x <- x - 1

	# Test for errors and return to originl scales
	if (x < 0 || x > 1) stop("problem here")
	rtn[index] <- SR_from_unit(x,min=ptab[index,"min"],max=ptab[index,"max"],logflag=ptab[index,"log"])

	rtn
}


