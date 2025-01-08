#####Adding CTD data (not 0-10m mounted on bongo)1#######
library(tidyverse)
credentials::set_github_pat() # for privaate repository
#devtools::install_github("clehoux/zooenv")
year_to_report=2024
year_of_assess=2025
egg <- read_rds(paste0("data/",year_of_assess,"/PL_Bongo_Capelin_larvae_Counts_L2_", year_to_report, ".RDS")) %>% filter(year > 2021) %>% 
  filter(!is.na(day)) %>%
  mutate(
    date = paste(year, str_pad(month,width=2, pad="0",side="left") , str_pad(day, width=2, pad="0", side="left"), sep = "/"),
    time = ifelse(is.na(start_time), "12:00", str_pad(start_time, width=4, pad="0")),
    time=ifelse(str_detect(time, pattern=":"),time, paste(substring(time, 1,2),substring(time, 3,4), sep=":")),
    depth.max = ifelse(is.na(sounding), start_depth + 5, sounding)
  ) %>%
  filter(!is.na(latitude), !is.na(longitude), !is.na(date), !is.na(depth.max)) %>%
  filter(!is.na(date))

library(zooenv)

#check if uynique ID are indeed unique
#set_sgdo_pass()
length(unique(egg$sample_id))==nrow(egg)
eggenv<- SGDO_sql(latitude=egg$latitude, longitude=egg$longitude, date=egg$date, time=egg$time, timezone=rep("UTC", nrow(egg)),depth.max=egg$depth.max, ID=egg$sample_id, geotol=0.02, timetol=24, station=egg$consecutive, range_ctd=c(0,50), range_bot=c(0,100))

egg_to_save<- left_join(egg, eggenv %>%  dplyr::rename(sample_id =ID)) %>%dplyr::select(-depth.max)
library(janitor)
get_dupes(egg_to_save %>%  dplyr::select(T_0_50, S_0_50))  
# 5 c'est pas si mal

t010 <- read_rds(paste0("data/",year_of_assess,"/PL_Bongo_Scomber_eggs_larvae_Counts_L2_", year_to_report, ".RDS")) %>% filter(year > 2021) %>%  dplyr::select(sample_id, temperature0_10, salinity0_10)
  
egg_to_save<- left_join(egg_to_save,t010)

write.table(egg_to_save %>%  arrange(year, station), paste0("data/",year_of_assess,"/PL_Bongo_Capelin_larvae_Counts_L2_",year_to_report,"_CTD.tsv"), row.names=F, sep="\t",dec=".")  
write_rds(egg_to_save %>%  arrange(year, station), paste0("data/",year_of_assess,"/PL_Bongo_Capelin_larvae_Counts_L2_",year_to_report,"_CTD.RDS"))  
