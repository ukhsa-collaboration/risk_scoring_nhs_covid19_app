setwd(dir_main)

# read other files with daily analytics or external data
AAE <- read_csv(AAEfile)

public_cases <- read_csv(public_cases_file)
geo <- read_csv(geo_file)
# create matching tables between LTLAs and regions
table_ltla_region <- geo %>% 
  select(lad20cd,region) %>% 
  distinct() %>% 
  mutate(region=str_replace_all(region," ","_"))
region_from_ltla <- setNames(
  pull(table_ltla_region,region),
  pull(table_ltla_region,lad20cd)
)
regions <- sort(unique(region_from_ltla[ltlas]))

# rural-urban scores in many different flavours
ruc <- read.table(ruc_file,sep="\t",header=T,quote="|")
rownames(ruc) <- ruc$LAD11CD
ruc_from_ltla <- as.character(ruc[ltlas,"RUC11CD"])
# _main is urban,rural,conurbation
ruc_from_ltla_main <- as.character(ruc[ltlas,"RUC11CD"])
ruc_from_ltla_main[ruc_from_ltla_main %in% c('1','2')]<-"rural"
ruc_from_ltla_main[ruc_from_ltla_main %in% c('3','4')]<-"urban"
ruc_from_ltla_main[ruc_from_ltla_main %in% c('5','6')]<-"conurbation"
# _simple version is purely numeric
ruc_from_ltla_simple <- as.numeric(ruc_from_ltla)
ruc_from_ltla_simple[region_from_ltla[ltlas]=="London"]<-7
ruc_from_ltla_simple[region_from_ltla[ltlas]=="Wales"]<-0
ruc_from_ltla_simple[is.na(ruc_from_ltla_simple)]<-3.5
# basic version includes Wales and London as separate entities, but merges the remaining conurbations
ruc_from_ltla[region_from_ltla[ltlas]=="Wales"]<-"Wales"
ruc_from_ltla[region_from_ltla[ltlas]=="London"]<-"London"
ruc_from_ltla[ruc_from_ltla %in% c('5','6')]<-"5/6"
names(ruc_from_ltla)<-ltlas
names(ruc_from_ltla_simple)<-ltlas
names(ruc_from_ltla_main)<-ltlas
list_locations <- sort(unique(ruc_from_ltla[ruc_from_ltla!='NA']))

# weekly smoothed cases by ltla/date
daily_rates_cases <- full_join(AAE %>% select(date,ltla,population,uptake,test_positive),public_cases) %>%
  arrange(ltla,date) %>%
  mutate(date=as.Date(date),
         daily_rate_cases_app=frollmean(test_positive/population/uptake,7,align="center"),
         daily_rate_cases_under15=frollmean(casesUnder15/population,7,align="center"),
         daily_rate_cases_over15=frollmean(casesOver15/population,7,align="center")
  ) %>%
  filter(getmonth(date) %in% c(months,months_iOS)) %>%
  select(-uptake,-population,-test_positive,-casesUnder15,-casesOver15)

# background probability of testing positive among non-notified app users in a 2-week window after the relevant date
background_cases <- AAE %>% 
  select(date,ltla,population,uptake,test_positive,positive_after_EN,notifications) %>%
  arrange(ltla,date) %>%
  mutate(date=as.Date(date),
         cum_notifications=frollsum(notifications,14,align="right"),
         bg_rate_cases_app=frollsum((test_positive-positive_after_EN)/(population*uptake-cum_notifications),14,align="left"),
  ) %>%
  filter(getmonth(date) %in% c(months,months_iOS)) %>%
  select(-uptake,-population,-test_positive,-positive_after_EN,-notifications,-cum_notifications)


