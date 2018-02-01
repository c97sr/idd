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
#' @details Details go here XXX new edit YYYY
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

if (FALSE) {
    library(devtools)
    detach("package:idd", unload=TRUE)
    install_github("c97sr/idd",force=TRUE)
    install.packages("../../idd",repos=NULL,force=TRUE)
    library(idd)
}
