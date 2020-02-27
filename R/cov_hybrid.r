#' @export
cov_hybrid <- function(
                       R0 = 2.6,
                       vecN=60000000*c(0.25,0.25,0.25,0.25),
                       vecI0=rep(1,length(vecN)),
                       vecInfNess=rep(1,length(vecN)),
                       vecSusTy=rep(1,length(vecN)),
                       vecPS=rep(0.1,length(vecN)),
                       vecPM=rep(0.1,length(vecN)),
                       vecPA=rep(0.1,length(vecN)),
                       matCt=matrix(1/length(vecN),
                                     nrow=length(vecN),
                                     ncol=length(vecN)),
                       matCtClosure=matrix(1/length(vecN),
                                     nrow=length(vecN),
                                     ncol=length(vecN)),
                       scLim=c(99999,99999),
                       vecTcalc=seq(0,720,0.1),
                       vecRtrel=rep(1,length(vecTcalc)),
                       D_E=2,
                       D_I1=5,
                       D_I2=3,
                       D_HR=20,
                       D_HD=23,
                       deterministic=TRUE,
                       plot=FALSE,
                       nReals=2,
                       dummy=NULL
                       ) { 
    
    ## Define housekeeping variables
    nAges <- length(vecN)
    nTimes <- length(vecTcalc)
    totalN <- sum(vecN)
    dt <- vecTcalc[2] - vecTcalc[1]
    vecOne <- rep(1,nAges)
    vecPF <- vecOne - vecPS - vecPM - vecPA
    vecPMcond <- vecPM / (vecOne - vecPS)
    vecPAcond <- vecPA / (vecOne - vecPS - vecPM)
    sevClassChars <- c("S","M","F","A")
    nSevClasses <- length(sevClassChars)
    matPsev <- data.frame(S=vecPS,M=vecPM,F=vecPF,A=vecPA)
    
    ## Test for consistency in the input parameters
    if (min(vecPF) < 0) {stop("Age specific severity parameters are not consistent")}
    if (deterministic) {nReals <- 1}

    ## Declare items to be returned
    rtn_inf <- array(dim=c(nTimes,nAges,nReals))
    rtn_pop <- array(dim=c(nTimes,nAges,nReals))
    
    ## Create next generation matrix, assuming each severity level is an
    ## infectious type. Needs to be nested to 4 levels for infector and infectee
    ngm = matrix(ncol=nAges*nSevClasses,nrow=nAges*nSevClasses)
    for (i in 1:nSevClasses) {
        tor <- sevClassChars[i]
        for (j in 1:nAges) {
            for (k in 1:nSevClasses) {
                tee <- sevClassChars[k] 
                for (l in 1:nAges) {
                    ngm_col <- (i-1)*nAges + j
                    ngm_row <- (k-1)*nAges + l

                    ## Not entirely sure why this doesn't have a term for the
                    ## relative size of each population, but will worry about that
                    ## later! vecN[l] / totalN. Seems to pass the tests.
                    ngm[ngm_row,ngm_col] <-
                        matCt[j,l] * (D_I1+D_I2) * matPsev[l,k] *
                        vecInfNess[j] * vecSusTy[l]

                    ## Close all the loops
                 }
            }
        }
    }

    ## Set the value of beta using the NGM matrix
    R0_ngm <- Re((eigen(ngm))$values[1])
    beta <- R0 / R0_ngm
    beta_check <- R0 / (D_I1 + D_I2)

    ## Initiate the realisation loop
    for (i in 1:nReals) {

        ## Initiate the state variables
        S <- vecN - vecI0
        E <- rep(0,nAges)
        I_1M <- vecI0
        I_2M <- rep(0,nAges)
        I_1F <- rep(0,nAges)
        I_2F <- rep(0,nAges)
        I_1A <- rep(0,nAges)
        I_2A <- rep(0,nAges)
        I_1S <- rep(0,nAges)
        I_2S <- rep(0,nAges)
        R <- rep(0,nAges)

        ## Initiate constant hazards
        vecHazExE = rep(1/D_E,nAges)
        vecHazExI_1M = rep(1/D_I1,nAges)
        vecHazExI_2M = rep(1/D_I2,nAges)
        vecHazExI_1F = rep(1/D_I1,nAges)
        vecHazExI_2F = rep(1/D_I2,nAges)
        vecHazExI_1A = rep(1/D_I1,nAges)
        vecHazExI_2A = rep(1/D_I2,nAges)
        vecHazExI_1S = rep(1/D_I1,nAges)
        vecHazExI_2S = rep(1/D_I2,nAges)

        ## Initiate constant probs
        vecProbExE = 1 - exp(-dt * vecHazExE)
        vecProbExI_1M = 1 - exp(-dt * vecHazExI_1M)
        vecProbExI_2M = 1 - exp(-dt * vecHazExI_2M)
        vecProbExI_1F = 1 - exp(-dt * vecHazExI_1F)
        vecProbExI_2F = 1 - exp(-dt * vecHazExI_2F)
        vecProbExI_1S = 1 - exp(-dt * vecHazExI_1S)
        vecProbExI_2S = 1 - exp(-dt * vecHazExI_2S)
        vecProbExI_1A = 1 - exp(-dt * vecHazExI_1A)
        vecProbExI_2A = 1 - exp(-dt * vecHazExI_2A)
        
        ## Initiate the return data structures
        rtn_inf[1,,i] <- vecI0
        rtn_pop[1,,i] <- sum(
            S,E,
            I_1M,I_2M,I_1F,I_2F,
            I_1S,I_2S,I_1A,I_2A,
            R)
        
        ## Initiate the time loop
        j <- 2
        while (j <= nTimes) {

            ## Set current beta
            beta_cur <- vecRtrel[j] * beta

            ## Set mixing matrix
            t_cur <- vecTcalc[j]
            if ((t_cur < scLim[1]) || (t_cur >= scLim[2])) {
                matCtTmp <- matCt
            } else {
                matCtTmp <- matCtClosure
            }
            
            ## Calculate variable hazards
            vecFoi <-  beta_cur * vecSusTy * (
                ((I_1M * vecInfNess) %*% matCtTmp) / vecN +
                ((I_2M * vecInfNess) %*% matCtTmp) / vecN +
                ((I_1F * vecInfNess) %*% matCtTmp) / vecN +
                ((I_2F * vecInfNess) %*% matCtTmp) / vecN +
                ((I_1S * vecInfNess) %*% matCtTmp) / vecN +
                ((I_2S * vecInfNess) %*% matCtTmp) / vecN +
                ((I_1A * vecInfNess) %*% matCtTmp) / vecN +
                ((I_2A * vecInfNess) %*% matCtTmp) / vecN
            )
            
            ## Calculate variable probabilites
            pVecFoi <- 1 - exp(-dt*vecFoi)
            
            ## Either draw random numbers or calculate averages
            ## Currently just below here
            if (deterministic) {
                noInf <- S * pVecFoi
                noExE <- E * vecProbExE
                noEntI1S <- noExE * vecPS
                noEntI1M <- noExE * vecPM
                noEntI1A <- noExE * vecPA
                noEntI1F <- noExE * (vecOne - vecPM - vecPS - vecPA)
                noExI_1M <- I_1M * vecProbExI_1M
                noExI_2M <- I_2M * vecProbExI_2M
                noExI_1F <- I_1F * vecProbExI_1F
                noExI_2F <- I_2F * vecProbExI_2F
                noExI_1S <- I_1S * vecProbExI_1S
                noExI_2S <- I_2S * vecProbExI_2S
                noExI_1A <- I_1A * vecProbExI_1A
                noExI_2A <- I_2A * vecProbExI_2A
            } else {
                noInf <- rbinom(nAges,S,pVecFoi)
                noExE <- rbinom(nAges,E,vecProbExE)
                noEntI1S <- rbinom(nAges, noExE, vecPS)
                noEntI1M <- rbinom(nAges, noExE - noEntI1S, vecPMcond)
                noEntI1A <- rbinom(nAges, noExE - noEntI1S - noEntI1M, vecPAcond)
                noEntI1F <- noExE - noEntI1S - noEntI1M - noEntI1A                        
                noExI_1M <- rbinom(nAges,I_1M,vecProbExI_1M)
                noExI_2M <- rbinom(nAges,I_2M,vecProbExI_2M)
                noExI_1F <- rbinom(nAges,I_1F,vecProbExI_1F)
                noExI_2F <- rbinom(nAges,I_2F,vecProbExI_2F)
                noExI_1S <- rbinom(nAges,I_1S,vecProbExI_1S)
                noExI_2S <- rbinom(nAges,I_2S,vecProbExI_2S)
                noExI_1A <- rbinom(nAges,I_1A,vecProbExI_1A)
                noExI_2A <- rbinom(nAges,I_2A,vecProbExI_2A)
            }

            ## Update the state variables
            ## Problems are probably here
            S <- S - noInf
            E <- E + noInf - noExE
            I_1M <- I_1M + noEntI1M - noExI_1M
            I_2M <- I_2M + noExI_1M - noExI_2M
            I_1F <- I_1F + noEntI1F - noExI_1F
            I_2F <- I_2F + noExI_1F - noExI_2F
            I_1S <- I_1S + noEntI1S - noExI_1S
            I_2S <- I_2S + noExI_1S - noExI_2S
            I_1A <- I_1A + noEntI1A - noExI_1A
            I_2A <- I_2A + noExI_1A - noExI_2A
            R <- R + noExI_2M + noExI_2F + noExI_2S + noExI_2A

            ## Record the other outpur variables
            rtn_inf[j,,i] <- noInf
            rtn_pop[j,,i] <- sum(S,E,
                                 I_1M,I_2M,I_1F,I_2F,
                                 I_1S,I_2S,I_1A,I_2A,
                                 R)
            
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

    ## Plot a quick diagnostic, mainly for dbugging the flow equations
    if (plot) {
        plot(
            vecTcalc,rtn_pop[,1,1]+1,
            type="l",log="y",col="red",
            ylim=c(1,max(rtn_pop[,1,1]))
        )
        points(vecTcalc/30,rowSums(rtn_inf[,,1])+1)
    }
    
    ## Return the outputs, expand as necessary
    list(inf=rtn_inf,
         cinf=rtn_inf_cum,
         t=vecTcalc,
         pop=vecN)

}
