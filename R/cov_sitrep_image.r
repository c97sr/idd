#' @export
cov_sitrep_image <- function(tab,
                             dates,
                             trusts,
                             title="",
                             crchar= "",
                             torder=1:length(trusts),
                             titlewidth=50) {

    ## Manage long titles
    newtitle <- substr(title,1,titlewidth)
    title <- substr(title,titlewidth+1,999)
    while (nchar(title) > 0) {
        newtitle <- paste(newtitle,substr(title,1,titlewidth),sep=crchar)
        title <- substr(title,titlewidth+1,999) 
    }

    ## Translate one time
    tab <- t(tab)

    ## Get the max values
    maxvals <- colSums(tab,na.rm=TRUE)

    image(tab[,torder],axes=FALSE)
    mtext(newtitle, cex=0.5, line=1)
    axis(1,at=seq(0,1,1/(length(dates)-1)),labels=dates,las=2,cex.axis=0.5)
    axis(2,at=seq(0,1,1/(length(trusts)-1)),labels=trusts[torder],las=2,cex.axis=0.5)
    axis(4,at=seq(0,1,1/(length(trusts)-1)),labels=maxvals[torder],las=2,cex.axis=0.5)
    mtext("Trust",side=2,line=2)
    mtext("Max row value",side=4,line=2)
    ## legend needed
}
