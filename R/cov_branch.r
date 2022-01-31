#' Function to simulate an ncov-like brtanching process
#' This might be useful for the border crossing work
#' Need to go through this and see if I cam make it more 
#' relevant to omicron
#' @export
cov_branch <- function(
                     nreal=2,
                     ninf=100,
                     tinfmax=999,
                     seed=10,
                     popN=-1,
                     gtd=dpois(0:20,8.1),
                     off=dpois(0:10,2.6),
                     ipd=dpois(0:10,5.8),
                     pass_per_day=3301,
                     R0 = -1,
                     k = 1,
                     popsize=19000000,
                     casetypes=c(1,1,1)
                     ) {

    ## Check for inconsistencies in the arguments
    if ((is.null(off) && R0 < 0)) {
        stop("conflict in vector off and R0 values")
    } 

    ## Set broad options
    if (R0 > 0) {
        disoff = TRUE
        ## Here set k vlaue
    }
    
    ## Define ranges for the waiting time ditrsibutions
    dom_gtd  <- 0:(length(gtd)-1)
    dom_off  <- 0:(length(off)-1)
    dom_ipd  <- 0:(length(ipd)-1)

    ## Setup output objects
    rtn_inf <- matrix(data=-1,nrow=ninf,ncol=nreal)
    rtn_wiw <- matrix(data=-1,nrow=ninf,ncol=nreal)
    rtn_wiw_v2 <- array(dim=c(ninf,length(off),nreal))
    rtn_ons <- matrix(data=9999,nrow=ninf,ncol=nreal)
    rtn_trv <- matrix(data=9999,nrow=ninf,ncol=nreal)
    rtn_cty <- matrix(data=sample(1:length(casetypes),ninf*nreal,
                              replace=TRUE,prob=casetypes)
                 ,nrow=ninf,ncol=nreal)
    
    ## Define possion rate travel
    lambda_travel <- pass_per_day / popsize

    ## Initiate loop for realizations
    for (i in 1:nreal) {

        ## Initiaate the seed
        for (j in 1:seed) {
            rtn_inf[j,i] <- 1
            rtn_wiw[j,i] <- 0
            rtn_trv[j,i] <- rtn_inf[j,i] + round(rexp(1,lambda_travel))
            rtn_ons[j,i] <- rtn_inf[j,i] + sample(dom_ipd,1,prob=ipd)
        }
        cur_inf  <- 1
        cur_sus  <- seed+1
        still_run  <- TRUE
        while (still_run) {

            ## TODO code in finite populations size here
            
            inf_to_assign <- sample(dom_off,1,prob=off)
            while (still_run && (inf_to_assign > 1)) {
                rtn_inf[cur_sus,i] <- rtn_inf[cur_inf,i] + sample(dom_gtd,1,prob=gtd)
                if (rtn_inf[cur_sus,i] > tinfmax) {still_run <- FALSE} 
                rtn_wiw[cur_sus,i] <- cur_inf
                rtn_wiw_v2[cur_inf,inf_to_assign,i] <- cur_sus
                rtn_ons[cur_sus,i] <- rtn_inf[cur_sus,i] + sample(dom_ipd,1,prob=ipd)
                rtn_trv[cur_sus,i] <- rtn_inf[cur_sus,i] + round(rexp(1,lambda_travel))
                inf_to_assign  <- inf_to_assign - 1
                cur_sus  <-  cur_sus + 1
                if (cur_sus > ninf) {still_run  <- FALSE}
            }
            cur_inf  <- cur_inf + 1
            if (cur_inf > ninf) {
                stop("this loop should be deleted, probably not needed")
                still_run  <- FALSE
            }
        }
    }

    ## Return the different matrices
    list(
        inf=rtn_inf,
        wiw=rtn_wiw,
        wiw2 = rtn_wiw_v2,
        ons=rtn_ons,
        trv=rtn_trv,
        cty=rtn_cty
    )    
}

if (FALSE) {
    cov_branch()    
}
