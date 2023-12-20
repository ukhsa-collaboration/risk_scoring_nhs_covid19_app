####################
# DEFINE PARAMETERS
####################

# rescale risk scores to new scale for the labels (old risk score 90 -> calibrated risk score 1), with 1 decimal digit 
rescaledRiskScore <- function(x){round(x/90*10)/10}

# boundaries of risk buckets
bucket_boundaries <- c(-1,100,120,140,160,200,300,400,500,Inf)
# number of risk buckets actually used (lower one includes only non-risky windows)
n_buckets <- length(bucket_boundaries)-2
# labels of buckets
bucket_labels <- paste(rescaledRiskScore(bucket_boundaries[-length(bucket_boundaries)]),rescaledRiskScore(bucket_boundaries[-1]),sep="-")
names_bucket_labels <- paste("[",bucket_boundaries[-length(bucket_boundaries)],",",bucket_boundaries[-1],")",sep="")
names(bucket_labels) <- names_bucket_labels 
# names of buckets
bucket_names <- paste("riskScore_",rescaledRiskScore(bucket_boundaries[1+c(1:n_buckets)]),sep="")

####################
# DEFINE FUNCTIONS
####################

# save table
SaveTable <- function(x,name=deparse(substitute(x))){
  write_csv(x,file=paste0(dir_out,"/OutputTable_",name,".csv"))
}
# read table saved by SaveTable
ReadTable <- function(x){
  assign(x,read_csv(file=paste0(dir_out,"/OutputTable_",x,".csv")), envir = .GlobalEnv)
}

# save plot in png and pdf
SaveFigure <- function(x,name=deparse(substitute(x)),...){
  ggsave(x,file=paste0(dir_out,"/OutputFigure_",name,".png"),units="cm",...)
  ggsave(x,file=paste0(dir_out,"/OutputFigure_",name,".pdf"),units="cm",...)
}

# distanceScore to distance in m (assume 1m for distances of <1m)
score2distance <- function(x){
  1/sqrt(x)
}

# compute delay in days
diff_days <- function(x,y){
  round(as.numeric(difftime(as.Date(x),as.Date(y),units='d')))
}

# extract month or month/year
getmonth <- function(x){
  z<-paste(month(x,label=TRUE),substr(year(x),3,4),sep="")
  z
}

# extract month or month/year for high quality period, merge Apr/May
getmonthhq <- function(x){
  z<-paste(month(x,label=TRUE),substr(year(x),3,4),sep="")
  ifelse(z %in% c("Apr21","May21"),"Apr/May21",z)
}

# extract OS
getOS <-function(x){
  if_else(str_starts(x,"iPhone"),"iOS","Android")
}

# extract day of the week
getdayweek <- function(x){
  wday(as.Date(x), week_start = 1)
}

# extract week of the year
getweek <-function(x){
  as.Date("2021-01-04") + 7*( diff_days(as.Date(x),"2021-01-04") %/% 7 )
}

# convert day of the week to name
dayweek2text <- function(x){
  vdays <- wday(ISOdate(1, 1, 1:7), label = TRUE, abbr = FALSE, week_start = 1)
  vdays[x]
}
# convert day of the week to short name
dayweek2short <- function(x){
  vdays <- wday(ISOdate(1, 1, 1:7), label = TRUE, abbr = TRUE, week_start = 1)
  vdays[x]
}

# function extracting scanInstances summaries
summarise_scanInstances <- function(scanInstances) {
  # extract vectors, restrict s to [0,30min]
  s_vec <- lapply(str_match_all(scanInstances, "\"s\":([-]?[0-9]+)"),function(x){pmax(0,pmin(1800,as.integer(x[,2])))})
  m_vec <- lapply(str_match_all(scanInstances, "\"m\":([-]?[0-9]+)"),function(x){as.integer(x[,2])})
  t_vec <- lapply(str_match_all(scanInstances, "\"t\":([-]?[0-9]+)"),function(x){as.integer(x[,2])}) 
  num_scans <- sapply(s_vec,length)
  stopifnot(all(num_scans > 0L))
  stopifnot(num_scans == sapply(m_vec,length))
  stopifnot(num_scans == sapply(t_vec,length))
  list(number_scans=num_scans,
       duration_scans=(sapply(s_vec,function(x){sum(x[-1])})),
       avg_min_attenuation=(sapply(m_vec,mean)),
       avg_typical_attenuation=(sapply(t_vec,mean)),
       min_min_attenuation=(sapply(m_vec,min)),
       min_typical_attenuation=(sapply(t_vec,min))
  )
}


# turn infectiousness into numeric score (normalization different from paper by 2.5x)
infectiousness_value <- function(x){
  c(high=1,standard=0.4,none=0)[x]
}

# define positivity (0 or 1) for individuals
tests_positive <- function(positive,number_exposures){
  0+(positive*number_exposures>2-0.01 | positive>1-0.01)
}

# filter exposures by LTLA from full data
extract_filtered_exposures_full <- function(){
  
  initial_time <- Sys.time()
  
  # initialise
  stats_exposures<-c()
  counter<-0
  
  # load data for full exposures by LTLA
  for(ltla in ltlas){
    counter<-counter+1
    if(!file.exists(paste0(dir_ltla,"/filtered_events_data_all_",ltla,".",extension_csv))){
        print(paste0("processing LTLA #",counter," ",ltla))
        # read file all exposures
        data <- read_csv(paste0(dir_full,"/events_data_all_",ltla,extension_csv))
        setDT(data)
        # filter null scanInstances
        data <- data[ grepl("s",scanInstances), ]
        # edit scanInstances
        data <- data[ , scanInstances := str_replace_all(scanInstances, fixed('""'), '"')]
        # extract summary scanInstances
        data <- data[ , c("number_scans", "duration_scans", "avg_min_attenuation",
                          "avg_typical_attenuation","min_min_attenuation",
                          "min_typical_attenuation") :=
                        summarise_scanInstances(scanInstances)]
        #extract risk buckets
        data <- data[ , risk_bucket := cut(riskScore,bucket_boundaries,right=F) ]
        # select only relevant columns
        filtered_data <- as_tibble(data) %>% 
          mutate(
            ltla=ltla,
            exposure_date=as.Date(payload_date),
            received_date=as.Date(BlobLastModifiedUtcTime),
            distanceScore=riskScore/duration_scans/infectiousness_value(infectiousness)
          ) %>%
          filter(
            diff_days(received_date,exposure_date)>=0 &
              diff_days(received_date,exposure_date)<=if_else(type=="exposureWindow",maxdelay_EN,maxdelay_ET) &
              received_date<=if_else(type=="exposureWindow",cutoff_date_notifications,cutoff_date_positivetests)
          )
        write_csv(filtered_data %>% select(-scanInstances), paste0(dir_ltla,"/filtered_events_data_all_",ltla,".",extension_csv))
        print(difftime(Sys.time(),initial_time))
    }
  }
}


# filter exposures by LTLA from matched data
extract_filtered_exposures_matched <- function(){
  
  initial_time <- Sys.time()
  
  # initialise
  stats_exposures<-c()
  counter<-0
  
  # load data for full exposures by LTLA
  for(ltla in ltlas){
    counter<-counter+1
    if(!file.exists(paste0(dir_ltla,"/filtered_unique_individuals_in_events_",ltla,".",extension_csv))){
      print(paste0("processing LTLA #",counter," ",ltla))
      # read file all exposures
      data <- read_csv(paste0(dir_matched,"/unique_individuals_in_events_",ltla,extension_csv))
      setDT(data)
      # filter null scanInstances
      data <- data[ scanInstances!="", ]
      # edit scanInstances
      data <- data[ , scanInstances := str_replace_all(scanInstances, fixed('""'), '"')]
      # extract summary scanInstances
      data <- data[ , c("number_scans", "duration_scans", "avg_min_attenuation",
                        "avg_typical_attenuation","median_min_attenuation",
                        "median_typical_attenuation") :=
                      summarise_scanInstances(scanInstances)]
      # extract risk buckets
      data <- data[ , risk_bucket := cut(riskScore,bucket_boundaries,right=F)]
      # select only relevant columns
      filtered_data <- as_tibble(data) %>% 
        mutate(
          ltla=ltla,
          exposure_date=as.Date(payload_date),
          notification_date=as.Date(BlobLastModifiedUtcTime_2),
          positive_date=as.Date(BlobLastModifiedUtcTime_3),
          distanceScore=riskScore/duration_scans/infectiousness_value(infectiousness),
          positive_test=(!is.na(as.Date(BlobLastModifiedUtcTime_3)))
        ) %>%
        filter(
          diff_days(notification_date,exposure_date)>=0 &
            diff_days(notification_date,exposure_date)<=maxdelay_EN &
            ifelse(positive_test,
                   diff_days(positive_date,notification_date)>=0 &
                     diff_days(positive_date,notification_date)<=maxdelay_NT,
                   TRUE) &
            notification_date<=cutoff_date_notifications                                      
        ) %>%
        group_by(individual_id,exposure_date) %>%
        # get number of exposure windows per individual and day
        mutate(daily_notifications=n(),) %>%
        ungroup() %>%
        # filter no more than 48 exposure windows per day
        filter(daily_notifications<=48) %>%
        select(-daily_notifications)
      fwrite(filtered_data %>% select(-scanInstances), paste0(dir_ltla,"/filtered_unique_individuals_in_events_",ltla,".",extension_csv))
    }
    print(difftime(Sys.time(),initial_time))
  }
}















# extract statistics by LTLA from full data
# (function mystat working on tibble)
extract_statistics_full <- function(mystat, save.file="", save.files.ltla=""){
  initial_time <- Sys.time()
  # initialise
  stats_exposures<-c()
  counter<-0
  # load data for full exposures by LTLA
  for(ltla in ltlas){
    counter<-counter+1
    print(paste0("processing LTLA #",counter," ",ltla))
    filtered_data <- as_tibble(fread(paste0(dir_ltla,"/filtered_events_data_all_",ltla,".",extension_csv)))
    # extract summary statistics
    temp_stats_exposures <- filtered_data %>% mystat()
    if(save.files.ltla!=""){
      write_csv(temp_stats_exposures,file=paste0(dir_ltla,"/",save.files.ltla,"_",ltla,".",extension_csv))
    }
    if(save.file!=""){
      stats_exposures <- bind_rows(stats_exposures, temp_stats_exposures)
      write_csv(stats_exposures,save.file)
    }
    print(difftime(Sys.time(),initial_time))
  }
  if(save.file=="") {
    return(stats_exposures)
  }
}


# extract statistics by LTLA from matched data
# (function mystat working on tibble)
extract_statistics_matched <- function(mystat, save.file="", save.files.ltla=""){
  initial_time <- Sys.time()
  # initialise
  stats_exposures<-c()
  counter<-0
  # load data for full exposures by LTLA
  for(ltla in ltlas){
    counter<-counter+1
    print(paste0("processing LTLA #",counter," ",ltla))
    filtered_data <- as_tibble(fread(paste0(dir_ltla,"/filtered_unique_individuals_in_events_",ltla,".",extension_csv)))
    # extract summary statistics
    temp_stats_exposures <- filtered_data %>% mystat()
    if(save.files.ltla!=""){
      write_csv(temp_stats_exposures,file=paste0(dir_ltla,"/",save.files.ltla,"_",ltla,".",extension_csv))
    }
    if(save.file!=""){
      stats_exposures <- bind_rows(stats_exposures, temp_stats_exposures)
      write_csv(stats_exposures,save.file)
    }
    print(difftime(Sys.time(),initial_time))
  }
  if(save.file=="") {
    return(stats_exposures)
  }
}

# NOT USED
## pre-filter functions to load exposure data
# pre-filter iOS only
pre_filter_iOS <- function(data){
  return(data %>% filter(str_starts(deviceModel,"iPhone")))
}
# pre-filter Android only
pre_filter_Android <- function(data){
  return(data %>% filter(!str_starts(deviceModel,"iPhone")))
}
# pre-filter only dates relevant for the analysis
pre_filter_dates_full <- function(data){
  return(data %>% filter(received_date>=start_date_notifications))
  # & received_date<=as.Date("2022-01-10")))
}
# pre-filter only dates relevant for the analysis
pre_filter_dates_matched <- function(data){
  return(data %>% filter(notification_date>=start_date_notifications))
  # & notification_date<=as.Date("2022-01-10")))
}
# pre-filter only OS and dates relevant for the analysis
pre_filter_iOS_full <- function(data){
  return(data %>% filter(
    str_starts(deviceModel,"iPhone") & received_date>=start_date_notifications_iOS
  ))
  # & received_date<=as.Date("2022-01-10")))
}
# pre-filter only OS and dates relevant for the analysis
pre_filter_iOS_matched <- function(data){
  return(data %>% filter(
    str_starts(deviceModel,"iPhone") & notification_date>=start_date_notifications_iOS
  ))
  # & notification_date<=cutoff_date_notifications))
}


##########################################
# define functions to process data by LTLA
##########################################

apply_by_ltla <- function(filename, fun_extract, fun_rejoin, gz=FALSE){
  data<-c()
  if(gz){
    extension <- ".csv.gz"
  } else {
    extension <- ".csv"
  }
  for(ltla in ltlas){
    data<-bind_rows(data,
                    fun_extract(read_csv(paste0(filename,"_",ltla,extension),progress = F,col_types = cols()))
    )
    cat("*")
  }
  cat(" done.")
  return(fun_rejoin(data))
}
