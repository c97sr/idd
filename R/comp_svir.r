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

#' Solves the SVIR model
#'
#' Creates both deterministic and stochastic solutions to a common
#' epidemic model.
#'
#' @param Tg generation time
#' @param R0 Basic reproductive number
#' @param N size of the population
#' @param I0 number of initial infectious individuals
#' @param dt time step
#' @param A coefficient of seasonaliity
#' @param noReals number of realizations, NOT YET IMPLEMENTED
#' @param deterministic mode of solutions
#'
#' @details Details go here XXX new edit YYYY
#'
#' @return A list of two elements. The first is the incidence of infection
#' and the second is the timepoints to which the incidence refers.
#'
#' @examples comp_svir()
#' @examples comp_svir(R0=0.999)
comp_svir <- function(
  Tg=2.8,
  R0=1.8,
  N=100000,
  I0=1,
  V0=0,
  dt=0.1,
  startVac=9999,
  vacRate=0,
  noTimeSteps=10*180,
  noReals=1,
  deterministic=TRUE,
  doPlot=FALSE
) {

    ## Define the variable to be returned
    rtn_inf_inc <- vector(mode="numeric",length=noTimeSteps+1)

    ## Set up state variables
    state_S <- N - I0 - V0
    state_I <- I0
    state_R <- 0
    state_V <- V0

    ## Assign the seed infection to time 0
    rtn_inf_inc[1] <- I0

    ## Start the main loop
    for (i in 1:noTimeSteps) {

        ## Set variables
        beta = R0 / Tg
        if (i*dt > startVac) {
            tmpVacRate <- vacRate
        } else {
            tmpVacRate <- 0
        }

        ## Calculate hazard rates
        hazInf <- beta * state_I / N
        hazRec <- 1 / Tg
        hazVac <- tmpVacRate / N
        hazLeaveSus <- hazInf + hazVac
        hazLeaveInf <- hazRec + hazVac
        
        ## Calculate probabilities
        pLeaveSus <- 1 - exp(-dt * hazLeaveSus)
        pLeaveInf <- 1 - exp(-dt * hazLeaveInf)
        pVac <- 1 - exp(-dt * hazVac)

        ## Calculate the expected or simulated number of events
        if (deterministic) {
            nInf <- state_S * pLeaveSus * hazInf / hazLeaveSus
            nVacSus <- state_S * pLeaveSus * hazVac / hazLeaveSus
            nRec <- state_I * pLeaveInf * hazRec / hazLeaveInf
            nVacInf <- state_I * pLeaveInf * hazVac / hazLeaveInf
            nVacRec <- state_R * pVac
            rtn_inf_inc[i+1] <- nInf
        } else {
            error("comp_svir: this case not yet implemented 91837")
        }

        ## Update the state variables
        state_S <- state_S - nInf - nVacSus
        state_I <- state_I + nInf - nVacInf - nRec
        state_R <- state_R + nRec - nVacRec
        state_V <- state_V + nVacSus + nVacInf + nVacRec

    } # End main time loop

    ## Plot if required
    if (doPlot) {
        plot(
          (0:noTimeSteps)*dt,
          rtn_inf_inc/dt,type="l",
          xlab="Time",ylab="Incidence (per day)"
        )
    }
    
    ## Return incidence
    list(
        inf_inc = rtn_inf_inc,
        time=(0:noTimeSteps)*dt
    )

}
