#' @export
ncov.hybrid <- function(
                        R0 = 2.1,
                        vecN=c(5000000,5000000),
                        vecI0=rep(15,length(vecN)),
                        vecInfNess=rep(1,length(vecN)),
                        vecPs=rep(0.1,length(vecN)),
                        vecPa=rep(0.1,length(vecN)),
                        matCt=matrix(1/length(vecN),
                                     nrow=length(vecN),
                                     ncol=length(vecN)),
                        vecTcalc=seq(0,360,0.1),
                        vecRtrel=rep(1,length(vecTcalc)),
                        D_E=2,
                        D_I1=5,
                        D_I2=3,
                        D_HR=20,
                        D_HD=23,
                        deterministic=FALSE,
                        nReals=2
                        ) {
    
    ## Define housekeeping variables
    nAges <- length(vecN)
    nTimes <- length(vecTcalc)
    dt <- vecTcalc[2] - vecTcalc[1]

    ## Test for consistency in the data
    ## TODO add in the tests here

    ## Declare items to be returned
    rtn_inf <- array(dim=c(nTimes,nAges,nReals))

    ## Reparameterise from R0 to beta
    beta <- R0 / D_I1
    
    ## Initiate the realisation loop
    for (i in 1:nReals) {

        ## Initiate the state variables
        S <- vecN - vecI0
        E <- rep(0,nAges)
        I_1M <- vecI0
        I_2M <- rep(0,nAges)
        R <- rep(0,nAges)

        ## Initiate constant hazards
        vecHazExE = rep(1/D_E,nAges)
        vecHazExI_1M = rep(1/D_I1,nAges)
        vecHazExI_2M = rep(1/D_I2,nAges)

        ## Initiate constant probs
        vecProbExE = 1 - exp(-dt * vecHazExE)
        vecProbExI = 1 - exp(-dt * vecHazExI)
        vecProbExI_1M = 1 - exp(-dt * vecHazExI_1M)
        vecProbExI_2M = 1 - exp(-dt * vecHazExI_2M)
        
        ## Initiate the return data structures
        rtn_inf[1,,i] <- vecI0
        
        ## Initiate the time loop
        j <- 2
        while (j <= nTimes) {

            ## Set current beta
            beta_cur <- vecRtrel[j] * beta
            
            ## Calculate variable hazards
            vecFoi <-  beta_cur * vecInfNess * (
                (I_1M %*% matCt) / vecN +
                (I_2M %*% matCt) / vecN
                )
            
            ## Calculate variable probabilites
            pVecFoi <- 1 - exp(-dt*vecFoi)
            
            ## Either draw random numbers or calculate averages
            if (deterministic) {
                noInf <- S * pVecFoi
                noExE <- E * vecProbExE
                noExI_1M <- I * vecProbExI_1M
                noExI_2M <- I * vecProbExI_2M
            } else {
                noInf <- rbinom(nAges,S,pVecFoi)
                noExE <- rbinom(nAges,E,vecProbExE)
                noExI_1M <- rbinom(nAges,I,vecProbExI_1M)
                noExI_2M <- rbinom(nAges,I,vecProbExI_2M)
            }

            ## Update the state variables
            S <- S - noInf
            E <- E + noInf - noExE
            I_1M <- I + noExE - noExI_1M
            I_2M <- I + noExI_2M - noExI_2M
            R <- R + noExI_2M

            ## Latest revision up to here
            
            ## Record the other outpur variables
            rtn_inf[j,,i] <- noInf

            ## Close off the loop
            j <- j+1
            
        }
    }

    ## Make derived outputs
    rtn_inf_cum <- array(dim=c(nTimes,nAges,nReals))
    rtn_inf_cum[1,,] <- rtn_inf[1,,]
    for (i in 2:nTimes) {
        rtn_inf_cum[i,,] <- rtn_inf_cum[i-1,,] + rtn_inf[i,,]
    }
    
    ## Return the outputs
    list(inf=rtn_inf,
         cinf=rtn_inf_cum,
         t=vecTcalc)

}
