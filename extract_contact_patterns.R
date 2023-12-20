##################################
# code to extract contact patterns
##################################

# peak day = day with most exposure windows
# get number of exposure windows on the peak day
get_peak_duration <- function(biweekly_counts){
  maxbc<-max(biweekly_counts)
}

# get the peak day of the week
get_peak_day <- function(biweekly_counts){
  (which.max(biweekly_counts)-1) %% 7 + 1
}

# get how many days of exposure
get_num_days <- function(biweekly_counts){
  sum(biweekly_counts>0)
}

# define contact pattern classification
classify_pattern <- function(maxbc,ndays){
  if(maxbc>=16){ return("household") } else {
    if(ndays>1) { return("recurring")} else {
      if(maxbc==1){ return("fleeting") } else {
        return("one-day")
      }
    }
  }
}
names_patterns <- c("household","recurring","one-day","fleeting")

# function to extract patterns
mystats_contact_patterns_from_matched_exposures <- function(data){
  data %>%
    filter(getmonth(notification_date) %in% months) %>%
    # Jan 3 2021 was a Sunday, so days are coded Mon:Sun=1:7
    mutate(relative_date=diff_days(exposure_date,"2021-01-03")) %>% 
    group_by(individual_id,ltla) %>%
    # extract first, last, random day of exposure, as well as number of days of exposure
    mutate(first_expdate=min(relative_date),last_expdate=max(relative_date),random_expdate=sample(rep(relative_date,each=2),size=1),ndays=length(unique(relative_date))) %>%
    # get the relative day in a Mon-Sun,Mon-Sun 2 week period coded as dw101-dw114
    mutate(relative_day=paste("dw",100+relative_date-7*(ceiling(min(relative_date)/7)-1),"",sep="")) %>%
    group_by(individual_id,ltla,first_expdate,last_expdate,random_expdate,ndays,relative_day) %>%
    summarise(n=n()) %>%
    ungroup() %>%
    select(individual_id,ltla,first_expdate,last_expdate,random_expdate,ndays,relative_day,n) %>%
    pivot_wider(names_from = relative_day, names_sort=TRUE, values_from = n, values_fill = 0)  %>% 
    select(all_of(c("individual_id","ltla","first_expdate","last_expdate","random_expdate","ndays",paste("dw",101:114,sep="")))) %>%
    rowwise() %>%
    # extract peak day and number of exposures at peak day
    mutate(peak_day=get_peak_day(c_across(starts_with("dw"))),peak_duration=get_peak_duration(c_across(starts_with("dw")))) %>%
    # extract classification
    mutate(classification=classify_pattern(peak_duration,ndays)) %>%
    ungroup()
  }

# build contact pattern tables
if(!data_already_extracted){
  extract_statistics_matched(
    function(x){mystats_contact_patterns_from_matched_exposures(pre_filter_dates_matched(x))},
    save.files.ltla="table_contactpatterns"
    )
}

