cov.chains  <- function(
                        simres,
                        max_ons = 49,
                        max_trv = 60) {

    require(data.tree)
    
    inf <- simres$inf
    ons <- simres$ons
    trv <- simres$trv
    wiw <- simres$wiw
    cty <- simres$cty

    nreals <- dim(inf)[2]
    noinf <- dim(inf)[1]
    vec_all_inf <- 1:noinf

    ## move this to the sim function    
    rtn_exp <- vector(length=nreals, mode="numeric")
    rtn_ons <- vector(length=nreals, mode="numeric")
    rtn_trv <- vector(length=nreals, mode="numeric")
    rtn_2nd <- NULL 

    ## Need to edit below to extract linelist of exported chains
    ## This needs to be a database of case trees from successfully
    ## exported infected individuals
    ## Not currently working for the number of secondary cases

    for (i in 1:nreals) {
        
        mask_exp  <- (ons[,i] <= max_ons) & (trv[,i] <= max_trv)
        mask_ons  <- ons[,i] <= max_ons
        
        rtn_exp[i] <- sum(mask_exp)
        rtn_ons[i] <- sum(mask_ons)
        rtn_trv[i] <- mean(trv[mask_exp,i]-ons[mask_exp,i])
        index_exp <- (1:noinf)[mask_exp]
        noexp <- length(index_exp)
        tmp_rtn_2nd <- NULL
        if (noexp > 1) {
            ## Need to add in the tree structure here
            ## Make a table of the time of the number of severe cases per cluster
            ## and the time of the onset of each of the severe cases
            ## If this is done on a per cluster basis, it should allow for a
            ## a pretty decent analysis of likleihood of detection and the distribution
            ## of cluster sizes.
            for (j in index_exp) {
                ## TODO Start cluster here and add to the  cluster one layer at a
                ## time until there are no more cases left
                ## while (more)
                tmp_2nd <- as.vector(c(j,vec_all_inf[wiw[,i]==j]))
                ## edits here to build the cluster
                tmp_rtn_2nd <- c(tmp_rtn_2nd,list(tmp_2nd))
            }
        }
        rtn_2nd <- c(rtn_2nd,list(tmp_rtn_2nd))

    }

    ## Assign infections to be one of three severiy levels according to the inputted
    ## weightings.
    
    list(
        exp=rtn_exp,
        ons=rtn_ons,
        trv=rtn_trv,
        snd=rtn_2nd
    )

}
