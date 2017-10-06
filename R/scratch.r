if (FALSE) {
    
    ## Clear out memory
    rm(list=ls(all=TRUE))

    ## First install the idd package either from local or from
    ## github. Assumes that pwd is the top directory of the IDD
    ## package
    strIddPackage <- "local"
    library("devtools")
    if (strIddPackage=="local") {
        install("./")
    } else {
        install_github("c97sr/idd")
    }

    # reload(idd)
    detach("package:idd",unload=TRUE)
    library(idd)
    library(gtools)

    ## Then separate bits of analysis below here

    # Run On's code
    # Solve x for y=Mx 
    M <- matrix(c(3.2,2.1,5.7,1.0,2.5,3.0,6.7,2.3,3.8),ncol=3,nrow=3)
    y <- c(1.3,2.7,0)
    solve(M,y)
    ## However, x includes some negative -1.1902549  1.5976893  0.5240487
    ## We try different combination of x which follow multinomial
    ## distribution
    ## using rdirichlet. We run all sample space by setting a very
    ## large n in rdiric## hlet. ## function and run
    ## all possible Mx and compare this predicted y with the true y 
    ## and see the function summation of square of (true y - predicted y)
    ## can be minimized
    min_error <- function(y,n){
        
        temp <- rep(0,n)        
        counter <- rep(0,n)        
    
        ## Below doesn't seem to be needed
        ## difference <- 5
        ## counter <- 0
        rtn <- rdirichlet(n,rep(1,3))     
        
        for(counter in 1:n) {
            temp[counter] <- sum(
                (y-rtn[counter,]%*%M)*(y-rtn[counter,]%*%M)
            ) 
        }
      
        list(a=rtn[which(temp == min(abs(temp))),],min = min(abs(temp)))
  
    }

    res <- min_error(y,1000000)
    ## a is the estimate ######

    res$a
    res$min
    ## TRUE y ###
    ## Need to chat from here XXX

    y
    ## Predicted y ###

    a%*%M
    
    M

    y

}
