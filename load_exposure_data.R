##########################################
#define functions to extract data
##########################################


# define basic summary statistics for full exposures
mystats_full_exposures <- function(data){
  data_all <- data %>% filter(type=="exposureWindow")
  data_pos <- data %>% filter(type!="exposureWindow")
  data_all_stats <- data_all %>% 
    group_by(exposure_date) %>%
    summarise(
      ltla=ltla[1],
      # mean delay EN by day of exposure
      avg_delay_EN=mean(diff_days(received_date,exposure_date)),
      # number of exposures by day of exposure
      daily_exposures=n(),
      # mean riskScore by day of exposure
      avg_riskScore=mean(riskScore),
      # median riskScore by day of exposure
      median_riskScore=median(riskScore),
      # mean attenuation by day of exposure
      avg_attenuation=mean(avg_typical_attenuation),
      # mean min attenuation by day of exposure
      avg_minattenuation=median(avg_min_attenuation),
      # mean distanceScore by day of exposure
      avg_distanceScore=mean(distanceScore),
      # median distanceScore by day of exposure
      median_distanceScore=median(distanceScore),
      # mean duration by day of exposure
      avg_duration_scans=mean(duration_scans),
      # median duration by day of exposure
      median_duration=median(duration_scans),
      # mean infectiousness by day of exposure
      avg_infectiousness=mean(infectiousness_value(infectiousness))
    )
  data_pos_stats<-data_pos %>% 
    group_by(exposure_date) %>%
    summarise(
      # mean delay ET by day of exposure
      avg_delay_ET_positive=mean(diff_days(received_date,exposure_date)),
      # number of exposures by day of exposure
      daily_exposures_positive=n(),
      # mean riskScore by day of exposure
      avg_riskScore_positive=mean(riskScore),
      # median riskScore by day of exposure
      median_riskScore_positive=median(riskScore),
      # mean attenuation by day of exposure
      avg_attenuation_positive=mean(avg_typical_attenuation),
      # mean min attenuation by day of exposure
      avg_minattenuation_positive=median(avg_min_attenuation),
      # mean distanceScore by day of exposure
      avg_distanceScore_positive=mean(distanceScore),
      # median distanceScore by day of exposure
      median_distanceScore_positive=median(distanceScore),
      # mean duration by day of exposure
      avg_duration_positive=mean(duration_scans),
      # median duration by day of exposure
      median_duration_positive=median(duration_scans),
      # mean infectiousness by day of exposure
      avg_infectiousness_positive=mean(infectiousness_value(infectiousness))
    )
  return(full_join(data_all_stats,data_pos_stats,by=(exposure_date="exposure_date"),suffix=c(".all",".pos")))
}


# define basic summary statistics for matched exposures
mystats_matched_exposures <- function(data){
  data_allpos <- data %>%
      group_by(exposure_date) %>%
      summarise(
        ltla=ltla[1],
        # mean delay EN by day of exposure
        avg_delay_EN=mean(diff_days(notification_date,exposure_date)),
        # mean delay ET by day of exposure
        avg_delay_ET=mean(diff_days(positive_date,exposure_date),na.rm=TRUE),
        # number of exposures by day of exposure
        daily_exposures=n(),
        # mean riskScore by day of exposure
        avg_riskScore=mean(riskScore),
        # median riskScore by day of exposure
        median_riskScore=median(riskScore),
        # mean attenuation by day of exposure
        avg_attenuation=mean(avg_typical_attenuation),
        # mean min attenuation by day of exposure
        avg_minattenuation=median(avg_min_attenuation),
        # mean distanceScore by day of exposure
        avg_distanceScore=mean(distanceScore),
        # median distanceScore by day of exposure
        median_distanceScore=median(distanceScore),
        # mean duration by day of exposure
        avg_duration=mean(duration_scans),
        # median duration by day of exposure
        median_duration=median(duration_scans),
        # mean infectiousness by day of exposure
        avg_infectiousness=mean(infectiousness_value(infectiousness)),
        # mean delay EN by day of exposure
        avg_delay_EN_positive=mean((diff_days(notification_date,exposure_date))[positive_test]),
        # mean delay ET by day of exposure
        avg_delay_ET_positive=mean((diff_days(positive_date,exposure_date)[positive_test]),na.rm=TRUE),
        # number of exposures by day of exposure
        daily_exposures_positive=sum(positive_test),
        # mean riskScore by day of exposure
        avg_riskScore_positive=mean(riskScore[positive_test]),
        # median riskScore by day of exposure
        median_riskScore_positive=median(riskScore[positive_test]),
        # mean attenuation by day of exposure
        avg_attenuation_positive=mean(avg_typical_attenuation[positive_test]),
        # mean min attenuation by day of exposure
        avg_minattenuation_positive=median(avg_min_attenuation[positive_test]),
        # mean distanceScore by day of exposure
        avg_distanceScore_positive=mean(distanceScore[positive_test]),
        # median distanceScore by day of exposure
        median_distanceScore_positive=median(distanceScore[positive_test]),
        # mean duration by day of exposure
        avg_duration_positive=mean(duration_scans[positive_test]),
        # median duration by day of exposure
        median_duration_positive=median(duration_scans[positive_test]),
        # mean infectiousness by day of exposure
        avg_infectiousness_positive=mean(infectiousness_value(infectiousness[positive_test]))
      )
  return(data_allpos)
}


# define basic summary statistics for matched exposures by individual
mystats_matched_exposures_individual <- function(data){
  data %>%
    group_by(individual_id,ltla) %>%
    summarise(
      # OS
      OS=ifelse(mean(str_starts(deviceModel,"iPhone"))>0.5,"iOS","Android"),
      # positivity
      positive=mean(positive_test),
      # notification date
      notification_date=mean(as.Date(notification_date)),
      # EN delay
      delay_EN=mean(diff_days(notification_date,exposure_date)),
      # of exposures
      number_exposures=n(),
      # means
      meanRiskScore = mean(riskScore),
      meanDuration = mean(duration_scans),
      meanInfectiousness = mean(infectiousness_value(infectiousness)),
      meanDistanceScore = mean(distanceScore),
      # distribution of riskScore buckets
      distrRiskScore=data.frame(as.list(setNames(table(c(names_bucket_labels,risk_bucket))[names_bucket_labels[-1]]-1,bucket_names))),
      # distribution of infectiousness
      distrInfectiousness=data.frame(as.list(table(c("high","standard",infectiousness))[c("high","standard")]-1)),
      # distribution of duration
      distrDuration=data.frame(as.list(setNames(hist(duration_scans,breaks=c(0,900,Inf),plot=FALSE)$counts,c("half_duration","full_duration")))),
      # distribution of distanceScore
      distrDistanceScore=data.frame(as.list(setNames(hist(distanceScore,breaks=1/rev(c(0,1.5,2,2.5,3,3.5,Inf))^2,plot=FALSE)$counts,paste0("distance_",rev(c(1,1.5,2,2.5,3,3.5)))))),
      # max risk and distance score
      maxRiskScore = max(riskScore),
      maxDistanceScore = max(distanceScore)
    ) %>% 
    unpack(cols = c(distrRiskScore,distrInfectiousness,distrDuration,distrDistanceScore))
}


##########################################
# extract data
##########################################
if (!data_already_extracted){

  if(!skip_create_filtered){
    extract_filtered_exposures_full()
    extract_filtered_exposures_matched()
  }
    
  stats_exposures_full <- extract_statistics_full(mystats_full_exposures,
    save.file=paste0(dir_ltla,"/simplestats_full_table.",extension_csv))
  stats_exposures_matched <- extract_statistics_matched(mystats_matched_exposures,
    save.file=paste0(dir_ltla,"/simplestats_matched_table.",extension_csv))
  extract_statistics_matched(
    function(x){mystats_matched_exposures_individual(pre_filter_dates_matched(x))},
    save.files.ltla="table_individual"
    )

} else {
  
  stats_exposures_full <- as_tibble(fread("simplestats_full_table.csv"))
  stats_exposures_matched <- as_tibble(fread("simplestats_matched_table.csv"))
  
}




