#' Define a function to take a table and add in binomial estimates and
#' CIs
#' @export
add.CIs <- function(tab,method="wilson") {
    rtn <- tab
    tmp <- dim(tab)
    nrows <- tmp[1]
    ncols <- tmp[2]
    rowP <- vector(mode="numeric", length=ncols)
    rowUB <- vector(mode="numeric", length=ncols)
    rowLB <- vector(mode="numeric", length=ncols)
    for (i in 1:ncols) {
        tmpbin <- propCI(x = tab["Detected",i],
                         n = tab["Detected",i] + tab["Not Detected",i],
                         method=method)
        rowP[i] <- tmpbin$p
        rowUB[i] <- tmpbin$upper
        rowLB[i] <- tmpbin$lower
    }
    rtn <- rbind(rtn, p=rowP, ub=rowUB, lb=rowLB)
    rtn
}

sim_las <- function(n_la=315,n_sim, max_per = 10) {
  
  mat_sim <- matrix()
  
}
