library(tidyverse)
library(lubridate)
library(magrittr)
library(stringr)
library(sf)
library(nngeo)
extract_T0_10 <-  function(year){
#format changed. stations are not used to merge with ctd, instead, positions, date. and nearest join are used  
  
  
#the format changed starting in 2021 so 2 file lists because date is missing 
file_list <- paste0("data/CTD_casts_0_10m/", list.files(path = "data/CTD_casts_0_10m/"))
file_list <-  file_list[!grepl(file_list, pattern="Gregoire")]
egg <- read_rds(paste0("data/",year,"/PL_Bongo_Scomber_eggs_larvae_Counts_L2_", year, ".RDS"))
dat_out <- data.frame()

#no data in 95,97 and 2020

#read and get uniform format
for(y in c(1987:1994,1996, 1998:2019, 2021:year)){
  if(y <= 2000){
    filey<- file_list[grep(file_list, pattern=y)]
    daty <-  read.delim(filey, header = FALSE, sep=",",skip = 1, col.names = c("file","station","datetime","lat","lon","temperature0_10","salinity","sigma_t")) %>% 
      mutate(station=as.character(station), 
             date=substring(datetime, 1,11), year=y) %>% dplyr::select(-datetime)
  }
  if(y >2000 & y < 2021){
    filey<- file_list[grep(file_list, pattern=y)]
    daty <-  read.delim(filey, header = FALSE, sep=",",skip = 1, col.names = c("file","station","date","time","lat","lon","temperature0_10","salinity","sigma_t")) %>%       
      mutate(station=as.character(station), year=y) %>%  dplyr::select(-time)
  }
  if(y %in% 2021:2022){
    filey<- file_list[grep(file_list, pattern=y)]
    daty <-  read.delim(filey, header = FALSE, sep=",",skip = 1, col.names = c("file","station","temperature0_10","salinity","sigma_t")) %>%  
      mutate(station=as.character(station), date=NA, lat=NA, lon=NA, year=y)
  }
  if(y  ==2023){
    filey<- file_list[grep(file_list, pattern=y)]
    daty <-  read_excel(filey, sheet=2) 
    colnames(daty) = c("file","station","datetime","lat","lon","prof_max","temperature0_10","salinity","sigma_t")  
      daty<- daty %>% mutate(station=as.character(station), date=NA, lat=NA, lon=NA, year=y)
  }
  if(y >=2024){
    filey<- file_list[grep(file_list, pattern=y)]
    daty <-  read.delim(filey, header = FALSE, sep=",",skip = 1, col.names = c("file","station","datetime","lat","lon","prof_max","temperature0_10","salinity","sigma_t")) %>%       
      mutate(station=as.character(station), year=y, date=substring(datetime, 10)) %>%  dplyr::select(-datetime)
  }  
  
  dat_out<- bind_rows(dat_out, daty) 
}

dat_out <- dat_out %>% mutate(#year=as.numeric(str_sub(date, 8,11)),
  month=str_sub(date,4,6),
  month = as.numeric(ifelse(month %in% c("Jun", "JUN"),as.numeric(6),
                            ifelse(month %in% c("JUL", "Jul"),as.numeric(7), NA))),
  day=str_sub(date,1,2),
  date = ymd(ifelse(!is.na(date),as.character(paste(year,month,day,sep = "-")), NA)))

#dat_out %<>%  filter(!grepl(file, pattern="IML2022024-4r"))#dplyr::filter(Station<13) %>% 

#multiple ctds at a given day/station in some cases...take mean


#theorical station position for years in which we don't have coordinates in CTD files.
lookup<- read.delim("data/lookup_station_egg.txt") %>%  mutate(station=as.character(station))
dat_out <- left_join(dat_out, lookup %>%dplyr::select(station, latitude,longitude)) %>%
  mutate(
    lat = coalesce(lat, latitude),
    lon = coalesce(lon, longitude)) %>%
  dplyr::select(-longitude, -latitude) %>% 
  mutate(date=as.character(date))  %>%  filter(!is.na(lat)) %>%  dplyr::rename(temp_station=station)# for correct export from nn_join


lookup_f = st_transform(st_as_sf(lookup, coords = c("longitude", "latitude"),crs=projdeg), crs=lcc)

dat_out_f = st_transform(st_as_sf(dat_out, coords = c("lon", "lat"),crs=projdeg), crs=lcc)

dat_out<- dat_out %>%  mutate(newstation=st_join(dat_out_f,lookup_f, join = st_nn, maxdist=40000)$station,
                         stratum=st_join(dat_out_f,lookup_f, join = st_nn, maxdist=40000)$stratum,
                         stratum_area=st_join(dat_out_f,lookup_f, join = st_nn, maxdist=40000)$stratum_area
                    )

#nearest join. can be the next neighbour station.
egg_joined <-  data.frame()
for(r in 1:nrow(egg)){
  eggr<- egg[r,]
  
  dat_y<- dat_out %>%  filter(year==eggr$year)
  
  if(any(!is.na(dat_y$date))){
    nearest_date <- as.character(c(seq.Date(ymd(eggr$date), by="-1 day", length.out = 2), seq.Date(ymd(eggr$date), by="1 day", length.out = 2)[2]))
    dat_y<- dat_y %>% filter(date %in% nearest_date)
  }
  
  egg_f = st_transform(st_as_sf(eggr, coords = c("longitude", "latitude"),crs=projdeg), crs=lcc)
  if(nrow(dat_y) > 0) temp_f = st_transform(st_as_sf(dat_y, coords = c("lon", "lat"),crs=projdeg), crs=lcc)
  if(nrow(dat_y) == 0) temp_f = data.frame()
  
  joined <- eggr %>%   dplyr::mutate(temperature0_10 =ifelse(nrow(temp_f) >0,st_join(egg_f,temp_f, join = st_nn, maxdist=40000)$temperature0_10, NA),
                                     salinity0_10 =ifelse(nrow(temp_f) >0,st_join(egg_f,temp_f, join = st_nn, maxdist=40000)$salinity, NA)) #%>% 
  #if_else make this function not work in case nrow==0
   # dplyr::mutate(date_ctd =as.character(ifelse(nrow(temp_f) >0,st_join(egg_f,temp_f, join = st_nn, maxdist=40000)$date.y, NA)))#à des fins de valdiation seulement
  
  egg_joined <-  bind_rows(egg_joined,joined)
 }

#adding 79-86 from Gregoire et al. 2013 tables. #trajet wsas verified with dates in Gregoire et al figures 
t79_86<- read.delim("data/CTD_casts_0_10m/temperature0_10_1979_1986_Gregoire_et_al_2013.txt") %>%  filter(!is.na(temperature0_10))

egg_joined79_86<- full_join(egg_joined %>%  filter(year < 1987) %>%  dplyr::select(-temperature0_10),t79_86 %>%  dplyr::mutate(station=as.character(station),
                                                                                                                               trajet=as.character(trajet))) 
egg_joined<- bind_rows(egg_joined79_86, egg_joined %>%  filter(year>1986)) %>%  filter(!is.na(sample_id))#to remove station that were in the temperature file but were remove from the egg dataset for quality control reasons


#using the stratum average if no nearest neighbour.
egg_joined<- left_join(egg_joined,egg_joined %>% group_by(year, stratum, trajet) %>%  summarize(ctemp = mean(temperature0_10, na.rm=T),
                                                                                                csal = mean(salinity0_10, na.rm=T))) %>% 
mutate(
  temperature0_10 = coalesce(temperature0_10, ctemp),
  salinity0_10 = coalesce(salinity0_10, csal)) %>%  dplyr::select(-csal, -ctemp)
  
write_rds(egg_joined, paste0("data/",year,"/PL_Bongo_Scomber_eggs_larvae_Counts_L2_", year, ".RDS"))

write.table(egg_joined, paste0("data/",year,"/PL_Bongo_Scomber_eggs_larvae_Counts_L2_", year, ".tsv"), sep="\t", dec=".", row.names=F)

}
