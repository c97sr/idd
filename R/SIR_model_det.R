####
#Deterministic SIR model
####
#Author: Harriet L Mills
#Date: 03/01/2014
####

#### example use:

# library(deSolve)

#initial_conditions = c(S=0.99, I=0.01, R=0) #, d_I=0) #from lsoda: the initial (state) values for the ODE system. If initial_conditions has a name attribute, the names will be used to label the output matrix.

# times=0:100 #from lsoda: times at which explicit estimates for y are desired. The first value in times must be the initial time.

# p_list = c(N=1, beta = 0.2, rho = 0.08, mu_S = 0.001, mu_I = 0.002, mu_R = 0.001)

# solution <- lsoda(y=initial_conditions, times=times, func=SIR_model_det, parms=p_list, atol=1e-80) 


#### function
SIR_model_det <- function(t, y, p_list){
  #t is the current time point
  #y is vector of the states
  #p_list is the parameter list
  
  rtn = vector() #creates a vector to store the derivatives in
  
  ## model
  lambda  <- p_list["beta"] * y["I"]/p_list["N"] 
  
  birth  <- p_list["mu_S"] * y["S"] + p_list["mu_I"] * y["I"]  + p_list["mu_R"] * y["R"]        
  
  #Susceptible class
  rtn["S"]   <- birth - lambda * y["S"] - p_list["mu_S"] * y["S"]                   
  
  #Infected class
  rtn["I"]   <- lambda * y["S"] - p_list["rho"] * y["I"] - p_list["mu_I"] * y["I"]   
  
  #Recovered class
  rtn["R"]   <- p_list["rho"] * y["I"] - p_list["mu_R"] * y["R"]                     
  
#   #records the incidence of infection
#   rtn["d_I"] <- - y["d_I"] + lambda * y["S"]   
  
  
  list(rtn)  	# must return a list
}
