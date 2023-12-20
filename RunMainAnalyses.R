##rm(list=ls())

# key configuration options:
# rerun table generation from the beginning or not?
data_already_extracted <- FALSE
# if rerunning, start from existing filtered files by LTLA or not?
skip_create_filtered <- FALSE
# save files as gzip or not?
save_as_gz <- TRUE

###############  PREPARE SCRIPT

using<-function(...) {
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,require,character.only=TRUE))
  need<-libs[req==FALSE]
  if(length(need)>0){ 
    install.packages(need)
    suppressWarnings(lapply(need,require,character.only=TRUE))
  }
}

using("dplyr", "reticulate", "readr", "tibble", "tidyr", "plotly", "viridis", "glue", "here", "data.table", "stringr", 
      "tidyverse","GGally","lubridate","R.utils")


############### MAIN CODE


# remember to load private configuration before this point!
##source("./configuration.R")

# set extension
extension_csv <- ifelse(save_as_gz,"csv.gz","csv")

# set base working directory
setwd(dir_main)

# load list of LTLAs
ltlas <- sort(as.vector(read.table(ltla_file)[,1]))
ltlas <- ltlas[ltlas!="E06000053"] # exclude Scilly

# set max delay between exposure and notification
maxdelay_EN <- 14
# set max delay between notifications and positive tests
maxdelay_NT <- 14
maxdelay_ET <- maxdelay_EN + maxdelay_NT
# set max number of exposures per day
max_exposure <- 48

# set cutoff date for the analysis
start_date_notifications <- as.Date("2021-01-01")
start_date_notifications_iOS <- as.Date("2021-01-01")
cutoff_date_notifications <- as.Date("2023-03-14")
cutoff_date_positivetests <- cutoff_date_notifications + maxdelay_NT

# months used for the analysis
months_iOS <- c(paste(c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),"21",sep=""),
                paste(c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),"22",sep=""),
                paste(c("Jan","Feb","Mar"),"23",sep=""))
months <- months_iOS
# high quality data in this period, Apr/May considered together
monthshq <- c(paste(c("Apr/May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),"21",sep=""),
                paste(c("Jan","Feb"),"22",sep="")) 

# get initial time in order to compute runtime
begins<-Sys.time()

# load utils/functions
source(paste0(dir_code,"/functions_utils.R"))

# extract statistics 
source(paste0(dir_code,"/load_exposure_data.R"))

# prepare contact patterns per individual 
source(paste0(dir_code,"/extract_contact_patterns.R"))

# prepare other data (geography, case counts, daily analytics...)
source(paste0(dir_code,"/prepare_external_data.R"))

#...and merge
source(paste0(dir_code,"/merge_external_data.R"))

#...and postprocess into tables
source(paste0(dir_code,"/postprocess_individual_data.R"))

# get final time and print runtime
ends<-Sys.time()
print(ends-begins)

#################
# STATISTICAL ANALYSES
#################

# generate tables/results for the analyses
source(paste0(dir_code,"/generate_statistics_for_paper_risk.R"))

# generate figures from the tables
source(paste0(dir_code,"/plots_for_paper_risk.R"))

# analyse and generate figures about ROC curves
source(paste0(dir_code,"/plot_predictors_ROC_curves.R"))



