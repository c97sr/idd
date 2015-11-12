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
