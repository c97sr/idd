cov.clusters <- function(inchains,sim) {

    ## Not sure this one is needed
    ## Going to try to add this functionality to the function above first
    
    ## Setup aux variables
    chains <- inchains$snd
    noreals <- length(chains)
    inf <- sim$inf
    tmax <- max(sim$inf)

    ## Declare the return types for the function
    rtn_nl_1s <- vector(length=tmax,mode="numeric")
    rtn_nl_2s <- vector(length=tmax,mode="numeric")
    rtn_nl_3s <- vector(length=tmax,mode="numeric")
    
    for (i in 1:noreals) {
        real <- chains[[i]]
        secs <- length(real)
        for (j in 1:nsecs) {
            ## Coding here if needed
        }
    }

    rtn
    
}
