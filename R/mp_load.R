
if (FALSE) {

  ## install.packages("readr")
  library("readr")
  
  fn_remote <- "https://raw.githubusercontent.com/globaldothealth/monkeypox/main/latest.csv"
  fn_local <- "C:/localfiles/tmp/dat_mpx.csv"
  
  # dat_mpx_rem <- read_csv(fn_remote)
  # write_csv(dat_mpx_rem,fn_local)
  dat_mpx <- read_csv(fn_local)
  # all.equal(dat_mpx_rem,dat_mpx)
  
  dim(dat_mpx)
  names(dat_mpx)
  table(dat_mpx$Country,useNA = "always")
  table(dat_mpx$Date_entry,useNA = "always")
  table(dat_mpx$Country,is.na(dat_mpx$Date_onset),useNA = "always")  
  
}
