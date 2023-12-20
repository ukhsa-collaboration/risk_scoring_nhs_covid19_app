##########################################
# extract exposure risk tables
##########################################
setwd(dir_ltla)

if(!data_already_extracted){
  table_outcome_riskScores <-
    apply_by_ltla("table_individual",
                  function(x){x %>%
                      filter(getmonth(notification_date) %in% months) %>% 
                      # assign positivity to contacts
                      mutate(positive=tests_positive(positive,number_exposures)) %>%
                      # add some extra stats
                      mutate(iOS=0+(OS=="iOS"),month=getmonth(notification_date),duration=number_exposures*meanDuration,region=region_from_ltla[ltla],ruc=ruc_from_ltla_simple[ltla],cumRiskScore=meanRiskScore*number_exposures) %>%
                      rename(date=notification_date) %>%
                      select(date,ltla,iOS,positive,month,number_exposures,duration,maxRiskScore,cumRiskScore,starts_with("riskScore_"),peak_day,peak_duration,first_expdate,last_expdate,random_expdate,ndays,classification,bg_rate_cases_app,
                             meanRiskScore,maxDistanceScore,meanDistanceScore,region,ruc)
                  }, function(y){
                    y %>% 
                      rename(avgRiskScore=meanRiskScore,maxProximityScore=maxDistanceScore,avgProximityScore=meanDistanceScore,RuralUrbanScore=ruc,firstExpDate=first_expdate,lastExpDate=last_expdate,randomExpDate=random_expdate)
                  },
              gz=save_as_gz)
  write_csv(table_outcome_riskScores,file=paste0("table_riskScores.",extension_csv))
  rm(table_outcome_riskScores)
  
  table_outcome_distances <-
    apply_by_ltla("table_individual",
                          function(x){x %>%
                              filter(getmonth(notification_date) %in% months) %>% 
                              # assign positivity to contacts
                              mutate(positive=tests_positive(positive,number_exposures)) %>%
                              # add some extra stats
                              mutate(month=getmonth(notification_date),duration=number_exposures*meanDuration,minDistance=1/sqrt(maxDistanceScore)) %>%
                              rename(date=notification_date) %>%
                              select(date,ltla,positive,month,number_exposures,duration,minDistance,starts_with("distance_"),peak_day,peak_duration,classification,bg_rate_cases_app)
                          }, function(y){
                            y
                          },
                      gz=save_as_gz)
  write_csv(table_outcome_distances,file=paste0("table_distances.",extension_csv))
  rm(table_outcome_distances)
}