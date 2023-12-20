if(!data_already_extracted){ 
  
  # merge data with matched exposures
  for(ltla in ltlas){
    print(ltla)
    myfilename <- paste0(dir_ltla,"/filtered_unique_individuals_in_events_",ltla,".",extension_csv)
    if(file.exists(myfilename)){
      # load original file
      oldfile_ltla <- as_tibble(fread(myfilename))
      # load contact pattern file
      contactfile <- read_csv(paste0(dir_ltla,"/table_contactpatterns_",ltla,".",extension_csv),progress = F,col_types = cols())
      # add information on classification and contact patterns 
      newfile <- left_join(oldfile_ltla,contactfile %>% select(individual_id,ltla,classification,peak_day,peak_duration,first_expdate,last_expdate,random_expdate,ndays))
      # add background cases
      newfile <- left_join(newfile %>% mutate(BlobLastModifiedUtcDate_2=as.Date(BlobLastModifiedUtcDate_2)),background_cases %>% mutate(date=as.Date(date)),by=c(localAuthority="ltla",BlobLastModifiedUtcDate_2="date"))
      # add OS and day of the week of the exposure
      newfile <- newfile %>% mutate(OS=getOS(deviceModel),dayweek=getdayweek(payload_date))
      write_csv(newfile, file=paste0(dir_ltla,"/filtered_unique_individuals_in_events_",ltla,".",extension_csv))
    }
  }
  
  # merge data with individual summaries
  for(ltla in ltlas){
    myfilename <- paste0(dir_ltla,"/table_individual_",ltla,".",extension_csv)
    if(file.exists(myfilename)){
      print(ltla)
      # load original file
      oldfile_ltla <- as_tibble(fread(myfilename))
      # load contact pattern file
      contactfile <- read_csv(paste0(dir_ltla,"/table_contactpatterns_",ltla,".",extension_csv),progress = F,col_types = cols())
      # add information on classification and contact patterns 
      newfile <- left_join(oldfile_ltla,contactfile %>% select(individual_id,ltla,classification,peak_day,peak_duration,first_expdate,last_expdate,random_expdate,ndays))
      # add background cases
      newfile <- left_join(newfile %>% mutate(notification_date=as.Date(notification_date)),background_cases %>% mutate(date=as.Date(date)),by=c(ltla="ltla",notification_date="date"))
      write_csv(newfile, file=paste0(dir_ltla,"/table_individual_",ltla,".",extension_csv))
    }
  }
  
}