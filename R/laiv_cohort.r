#' This function is supposed to return simulated cohort studies
#' of young children in the UKL that we can use to evaluate the
#' possible impact of repeated LAIV vaccination on the population
#' The default arguments are currently a little clumsy
laiv_cohort <- function(np = 100,
                        ny = 6,
                        inf_sus = 0.3333,
                        an_attack = rep(0.15,ny),
                        an_vacc = c(rep(0,2),rep(0.5,4))) {

  ## Setup the variables to be returned
  rtninf <- matrix(0,nrow = np,ncol=ny)
  rtnvac <- matrix(0,nrow = np,ncol=ny)

  ## Cycle through the years simulating events
  ## Note that in the description of this we use 0
  ## as the firts index
  for (y in 1:ny) {
    this_ar <- an_attack[y]
    this_vacc <- an_vacc[y]
    for (p in 1:np) {

      ## Look at vaccination first
      if (runif(1) < this_vacc) {rtnvac[p,y] <- 1}

      ## Now look at infection uniform for prior years
      if (y < ny) {
        if (runif(1) < this_ar) {rtninf[p,y] <- 1}

        ## But not uniform for the last year
      } else {
        if (rtninf[p,ny-1] > 0) {
          if (runif(1) < (this_ar*inf_sus)) {rtninf[p,y] <- 1}
        } else {
          if (runif(1) < this_ar) {rtninf[p,y] <- 1}
        }
      }
    }
  }

  ## Pick the natural histories we want to compare
  mask_inf_novacc <- (rtninf[,ny-1] > 0) & (rtnvac[,ny-1] < 1)
  mask_noinf_vacc <- (rtninf[,ny-1] < 1) & (rtnvac[,ny-1] > 0)

  ## Extract the two populations and outcomes based on natural
  ## history
  n_inf_novacc <- sum(mask_inf_novacc)
  n_noinf_vacc <- sum(mask_noinf_vacc)
  inf_inf_novacc <- sum(rtninf[mask_inf_novacc,ny])
  inf_noinf_vacc <- sum(rtninf[mask_noinf_vacc,ny])

  vecx = c(inf_inf_novacc,inf_noinf_vacc)
  vecn = c(n_inf_novacc,n_noinf_vacc)

  ## Then look at the power to distinguish between them
  power_calc <- prop.test(
    x = vecx,
    n = vecn,
    alternative = "l",
    conf.level=0.95)

  ## Return the simulated events and the summary statistics
  list(
    ## inf = rtninf,
    ## vac = rtnvac,
    pval = power_calc$p.val,
    raw = data.frame(x=vecx,n=vecn)
  )

}
