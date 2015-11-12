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
