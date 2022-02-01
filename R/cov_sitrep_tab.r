#' export
cov_sitrep_tab <- function(df,
                           fldname="new_dx_total",
                           long=TRUE,
                           propNAs=TRUE) {

    ## Make it work for wide and long data
    if (long) {
        metcol <- "metric_value"
        dftmp <- df[df$metric_name==fldname,]
        trust_field <- "org_code"
    } else {
        metcol <- fldname
        dftmp <- df
        trust_field <- "trust_code"
    }

    ## 
    allDates <- as.Date(
        seq(from=as.numeric(min(df$date)),to=as.numeric(max(df$date))),
        origin="1970-01-01"
    )
    allTCodes <- names(table(df[,trust_field]))
    nd <- length(allDates)
    nt <- length(allTCodes)

    #' Define a function to make tables for different fields in the sitrep
    #' data
    tabRtn <- matrix(nrow=nt,ncol=nd)
    tabNoNAs <- matrix(0,nrow=nt,ncol=nd)
    tabNoVals <- matrix(0,nrow=nt,ncol=nd)

    for (i in 1:nObs) {
        curCol <- match(dftmp$date[i],allDates)
        curRow <- match(dftmp[i,trust_field],allTCodes)
        curVal <- dftmp[i,metcol]
        curNA <- is.na(curVal)
        if (curNA) {
            tabNoNAs[curRow,curCol] <- tabNoNAs[curRow,curCol] + 1
            ## Consider a conditional loop on how to propagate NAs
            if (propNAs) {
            }
        } else {
            tabNoVals[curRow,curCol] <- tabNoVals[curRow,curCol] + 1
            if (is.na(tabRtn[curRow,curCol])) {
                tabRtn[curRow,curCol] <- curVal
            } else {
                tabRtn[curRow,curCol] <- tabRtn[curRow,curCol] + curVal
            }
        }
    }
    
    ## Return the different elements
    list(
        tab=tabRtn,
        tabVals=tabNoVals,
        tabNAs=tabNoNAs,
        dates=allDates,
        trusts=allTCodes)

}

