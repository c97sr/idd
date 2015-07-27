# comp.seir()
comp.seir <- function(	
		De=1.48,			# Duration latent
		Tg=2.6,				# Generation time
		R0=1.8,				# Basic reproductive number
		N=6800000,			# Population size
		I0=10,				# Initial number infective
		dt=1,				# Timestep
		R1=1.0,         	# Second R value
		t1=999,         	# time point of change
		R2=1.0,         	# Third R value
		t2=9999,        	# time point of change
		noTimeSteps=10, 	# Number of timesteps
		noReals=1,			# Number of realizations, ignored if solution is deterministic NOT YET IMPLEMENTED
		deterministic=TRUE	# Mode of solution
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
		hazInf 		<- beta * state_I / N
		hazBecInf 	<- 1 / De
		hazRec 		<- 1 / Di 
		
		# Calculate probabilities
		pInf 	<- 1 - exp(-dt * hazInf)
		pBecInf <- 1 - exp(-dt * hazBecInf)
		pRec  <- 1 - exp(-dt * hazRec)
		
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
