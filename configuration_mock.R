############################## 
# MOCK CONFIGURATION INFO
############################## 

# FILES AND DIRECTORIES
# set base working directory
dir_main <- "YOUR/MAIN/DIRECTORY/FOR/THIS/ANALYSIS"
# set directory containing original data
dir_data <- dir_main
# set code directory
dir_code <- "YOUR/MAIN/CODE/DIRECTORY/FOR/THIS/ANALYSIS"
# set directory containing metadata
dir_meta <- paste0(dir_code,"/metadata")
# set subdirectories with exposure data
dir_matched <- dir_data
dir_full <- dir_data
dir_ltla <- dir_main
# set output directory
dir_out <- paste0(dir_main,"/FiguresTables") 
# set file AAE data and public cases data
AAEfile <- paste0(dir_main,"/clean_timeseries_24Sep20_9May23.csv")
# set file for public cases, obtained from coronavirus.data.gov.uk using 
# > write_csv(data %>% filter(!(age %in% c("unassigned","00_59","60+"))) %>% mutate(under15=ifelse(age %in% c("00_04","05_09","10_14"),"casesUnder15","casesOver15")) %>% group_by(areaCode,date,under15) %>% summarise(cases=sum(cases)) %>% pivot_wider(names_from=under15,values_from=cases) %>% rename(ltla=areaCode),"cases_by_ltla_under15.csv")
public_cases_file <- paste0(dir_meta,"/cases_by_ltla_under15.csv")
# set file used for geographical mappings
geo_file <- paste0(dir_meta,"/Merged_PD_Demo_Geo_2021-01-21_v2.csv")
# set file with rural/urban score, see https://www.ons.gov.uk/methodology/geography/geographicalproducts/ruralurbanclassifications/2011ruralurbanclassification 
ruc_file <- paste0(dir_meta,"/RUC11_LAD11_ENv2.txt")
# set file containing list of LTLAs
ltla_file <- paste0(dir_meta,"/list_ltlas.csv")


