#' @export
load_uk_cov_data <- function(
                             at="region",
                             met1="newCasesBySpecimenDate",
                             met2="",
                             met3="",
                             met4="") {
    AreaType1 <- at
    Metric1 <- ifelse(met1=="","",paste0("&metric=",met1))
    Metric2 <- ifelse(met2=="","",paste0("&metric=",met2))
    Metric3 <- ifelse(met3=="","",paste0("&metric=",met3))
    Metric4 <- ifelse(met4=="","",paste0("&metric=",met4))
    DynUrl1 <- paste0(
        "https://api.coronavirus.data.gov.uk/v2/data?areaType=",
        AreaType1,Metric1,Metric2,Metric3,Metric4,"&format=csv")
    rtn <- read.csv(DynUrl1)
    rtn
}

## Script to test the functions developed here
if (FALSE) {
    x <- load_uk_cov_data(at="ltla")
    dim(x)    
}
