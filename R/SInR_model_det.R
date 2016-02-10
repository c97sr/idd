####
#Deterministic SI^nR model
####
#Author: Harriet L Mills
#Date: 03/01/2014
####

#### Please note: if you know the number of groups required, it is better (more efficient) to write this function explicitly for that number of groups.

#### example use:
# library(deSolve)

# for a four stage infection
# initial_conditions=c(S=0.99, I1=0.01, I2=0, I3=0, I4=0, R=0) #, d_I=0) #from lsoda: the initial (state) values for the ODE system.
# 
# times=0:100 #from lsoda: times at which explicit estimates for y are desired. The first value in times must be the initial time.
# 
# p_list = c(N=1, beta = 0.2, number_infclasses = 4, 
#            alpha1 = 1, alpha2 = 0.2, alpha3 = 1.2, alpha4=1.5, #infectivity parameters
#            rho1 = 0.0, rho2 = 0.0, rho3 = 0.0, rho4 = 0.0, #recovery parameters
#            sigma1 = 0.1, sigma2 = 0.1, sigma3 = 0.1,  # progression parameters (notice no progression from final stage)
#            mu_S = 0.001, mu_I1 = 0.002, mu_I2 = 0.002, mu_I3 = 0.004, mu_I4 = 0.004, mu_R=0.00) # death rate parameters
# 
# solution <- lsoda(y=initial_conditions, times=times, func=SInR_model_det, parms=c(p_list), atol=1e-80) 



#### function
SInR_model_det <- function(t, y, p_list){
  #t is the current time point in the integration (using lsoda to solve)
  #y is vector of the states
  #p_list is the parameter list
  
  rtn = vector() #creates a vector to store the derivatives in
  
  ## model
  
  # to find force of infection we must sum over all infectious classes, scaling by alpha
  lambda = p_list["beta"] * sum( p_list[paste("alpha", 1:p_list["number_infclasses"], sep="")] * y[paste("I", 1:p_list["number_infclasses"], sep="")]/p_list["N"] )
  
  birth = p_list["mu_S"] * y["S"] + sum( p_list[paste("mu_I", 1:p_list["number_infclasses"], sep="")] * y[paste("I", 1:p_list["number_infclasses"], sep="")] ) + p_list["mu_R"] * y["R"] 

  #Susceptible class
  rtn["S"]   <- birth - lambda * y["S"] - p_list["mu_S"] * y["S"]

  #infectious classes
  for (i in 1:p_list["number_infclasses"]){
    
    infclass=paste("I", i, sep="") #determine what infecious class we are in
    
    if (i==1){ # if in first stage we include the new infections
      
      rtn[infclass]  <- lambda * y["S"] - (p_list["rho1"] + p_list["sigma1"] + p_list["mu_I1"]) * y[infclass]   
      
    } else if (i==p_list["number_infclasses"]){ #if in final stage, no progression
      
      rtn[infclass]  <- p_list[paste("sigma", i - 1, sep="")] * y[paste("I", i-1, sep="")] - (p_list[paste("rho", i, sep="")] + p_list[paste("mu_I", i, sep="")]) * y[infclass] 
    
    } else { #in any other infectious stage individuals progress in and out, recover and die
      
      rtn[infclass]  <- p_list[paste("sigma", i - 1, sep="")] * y[paste("I", i-1, sep="")] - (p_list[paste("rho", i, sep="")] + p_list[paste("sigma", i, sep="")] + p_list[paste("mu_I", i, sep="")]) * y[infclass]   #Infected class
    
    }
    
  }
  #recovered class
  rtn["R"] <- sum( p_list[paste("rho", 1:p_list["number_infclasses"], sep="")] * y[paste("I", 1:p_list["number_infclasses"], sep="")] ) - p_list["mu_R"] * y["R"] 
  
#   #records the incidence of infection (those passing into the first infectious stage)
#   rtn["d_I"] <- - y["d_I"] + lambda * y["S"]     
  
  
  list(rtn)  # must return a list
}






