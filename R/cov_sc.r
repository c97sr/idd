#' Explores school closure options for a viral respiratory pandemic
#' 
#' Used in Feb 2020.
#' @export
cov_sc <- function(vecR0=seq(1.5,3,0.5),
                   vecTStart=seq(21,21+2*7,7),
                   vecDuration=seq(14,14+3*28,28),
                   vecSusChild=c(1,0.333)
                   ) {

    ## Make the dataframe for the output. Define a signle row and then 
    ## make
    ## the return object by zeroing off the row
    one_row <- data.frame(R0=-1,sus=-1,tstart=-1,dur=-1,bpeak=-1,scpeak=-1,deltacar=-1)
    rtn_df <- one_row[-1,]
    
    # Make a toy example to test the R0 parameterisation
    matopen <- matrix(c(4,2,2,1),nrow=2,ncol=2)
    matclose <- matrix(c(1,2,2,1),nrow=2,ncol=2)

    ## Setup the loops
    for (R0 in vecR0) {
        for (suschild  in vecSusChild) {

            ## Run the baseline scenario and record the key outputs
            y_base <- cov_hybrid(
                R0=R0,
                matCt=matopen,
                matCtClosure=matclose,
                scLim=c(999,1000),
                vecSusTy=c(suschild,1)
            )
            
            for (ts in vecTStart) {
                for (dur in vecDuration) {

                    ## Run the intervention scenarios and record the key values
                    y_int <- cov_hybrid(
                        R0=R0,
                        matCt=matopen,
                        matCtClosure=matclose,
                        scLim=c(ts,ts+dur),
                        vecSusTy=c(suschild,1)
                    )

                    ## Add the output to a tmp row
                    one_row[1,"R0"] <- R0
                    one_row[1,"sus"] <- suschild
                    one_row[1,"tstart"] <- ts
                    one_row[1,"dur"] <- dur
                    one_row[1,"deltacar"] <-
                        (sum(y_base$inf[,,1]) - sum(y_int$inf[,,1])) / sum(y_int$pop) 
                    rtn_df <- rbind(rtn_df,one_row)
                    
                }
            }
        }
    }

    ## Return the dataframe of results
    rtn_df
    
}
