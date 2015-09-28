# Function to give a complete tree for a simple branching process
# assumes poisson distributed offspring distribution and serial interval
# distribution

ind.tau <- function(
	Tg=5,		# Generation time
	R0=2.6,		# Basic reproductive number
	Nmax=30		# Maximum number of infections in the cluster
) {
	
	# Define the outputs to be returned
	rtn_times <- numeric(mode="numeric",length=Nmax)
	rtn_infectors <- numeric(mode="numeric",length=Nmax)
	
	# Initiate the 
	rtn_times[1] <- 0
	rtn_infectors[1] <- 1
	
	completed_infector <- 0
	current_infector <- 1
	
	# Main loop to go through every infector 
	while (completed_infector <= Nmax) {

	# Not yet completed
		
	}

    # Return list of infectors and times of infection
    list(t=rtn_times, i = rtn_infectors )

}
