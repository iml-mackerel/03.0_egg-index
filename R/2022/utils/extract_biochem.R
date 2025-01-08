library(tidyverse)
library(lubridate)
library(magrittr)
library(stringr)
source("bdOracle.Rprofile") #containing password and username for data extraction
source("R/biochem/PL_Get_SampleID_Batch.R")
source('R/biochem/PL_Get_Counts_Batch.R')
source("R/biochem/PL_Read_Filter.R")
source("R/biochem/PL_Taxonomic_Grouping.R")
#######Biochem extraction##########

#path_to_file for the excel file of the last year of sampling in case it is not available in biochem.if NULL , skipped and taken only from biochem

extract_biochem <-  function(year,path_to_file=NULL){

#get metadata  
PL_Get_SampleID_Batch(sql_file=paste0("sql/",year,"/BOCHEM_Bongo.sql"), output_file=paste0("data/",year,"/PL_Bongo_SampleID_", year, ".tsv"))
#metadata extracted
bongocheck<- read.delim(paste0("data/",year,"/PL_Bongo_SampleID_", year, ".tsv"))

#get counts
PL_Get_Counts_Batch(ifile=paste0("data/",year,"/PL_Bongo_SampleID_", year, ".tsv"), output_file=paste0("data/",year,"/PL_Bongo_Counts_L0_", year, ".tsv"),
                    sql_file= paste0("sql/",year,"/PL_Counts_from_SampleID.sql"))
#counts extracted
bongo<- read.delim(paste0("data/",year,"/PL_Bongo_Counts_L0_", year, ".tsv"))

##Cleaning up, largely based on A.Smith script
sgsl_ichthyo<- bongocheck %>%  
  mutate(date = as.Date(paste(day, month, year, sep="/"), format = "%d/%m/%Y"),
         doy=yday(date),start_time = str_pad(start_time, width=4, side="left", pad="0"),
         end_time = str_pad(end_time, width=4, side="left", pad="0"),
         hour_start =as.numeric(substring(start_time, 1,2)),
         minute_start=as.numeric(substring(start_time, 3,4)),
         hour_end =as.numeric(substring(end_time, 1,2)),
         minute_end=as.numeric(substring(end_time, 3,4)),
         hour_start = ifelse(hour_start == 0,24, hour_start),
         hour_end = ifelse(hour_end == 0,24,
                           ifelse(hour_end==1 & hour_start==24, 25,hour_end)), 
         set_duration= ((hour_end*3600) + (minute_end*60)) - ((hour_start*3600) + (minute_start*60))) %>% 
  dplyr::select(-end_time, -hour_start, -minute_start, -hour_end, -minute_end)# calculate length of time bongos were sampling


# simplify consecutiveutive variable (consecutive = mission specific order the sets were done in, letters at the end indicate a second set or pass at the same station or using both tribord and babord bongo nets etc)
sgsl_ichthyo %<>% mutate(consecutive = str_remove(consecutive,"consec"))
# check for double occurrences of consecutives (I checked this and it corresponds to a second pass generally)
sgsl_ichthyo  %<>% 
  mutate(
    extra_consecutive = case_when(
      stringi::stri_detect_regex(consecutive, "[abcdefghijklmnopqrstuvwxyz]") ~ "0",
      TRUE ~ "1"
    )
  )

# In some years, station names were appended with various suffixes 
# '' | 'a' indicates first pass
# 'b' & 'c' indicate 2nd and 3rd pass respectively (only one case of the latter)
# '_B' & '_T' indicate babord and tribord respectively (i.e. which bongo net was used. Default is babord)
# 'COR' | 'cor' indicate a correction (not sure if a second pass at the same station or just a new line of data)

# Remove all lowercase characters (keep upper case for now because they will be used later for filtering)
sgsl_ichthyo  %<>% 
  mutate(
    station = str_remove_all(collector_station_name, "[abcdefghijklmnopqrstuvwxyz]")
  )

# Essentially anything alphanumeric is a station from the second pass or a special experiment so create a new variable to ID them
sgsl_ichthyo  %<>% 
  mutate(
    extra_station = case_when(
      stringi::stri_detect_regex(collector_station_name, "[bcdefghijklmnopqrstuvwxyz]") ~ "0", # notice that 'a' is removed here
      TRUE ~ "1"
    )
  )

# Determine passes via suffix, defaults to = 1, i.e. when there is no suffix then it is first pass
sgsl_ichthyo  %<>% 
  mutate(
    pass1 = case_when(
      stringi::stri_endswith_fixed(collector_station_name, "a",case_insensitive=F) ~ "1",
      stringi::stri_endswith_fixed(collector_station_name, "b",case_insensitive=F) ~ "2",
      stringi::stri_detect_regex(collector_station_name, "b_") ~ "2",
      stringi::stri_endswith_fixed(collector_station_name, "c",case_insensitive=F) ~ "3",
      TRUE ~ "1"
    )
  )

# Determine passes via order of occurrence (i.e. order by date then assign 1 to first occurrence of a station and 2 to last occurrence of the same station). 
# Note, I checked the counts of stations and all counts are either 1 or 2 except for one invertebrate taxon in 1989 where n=4. df <- sgsl_ichthyo  %>%  arrange(date,taxons) %>% group_by(year, event_collector_stn_name, taxons) %>% dplyr::summarise(n=n())
sgsl_ichthyo %<>% arrange(date) %>% 
  group_by(year, station) %>% 
  mutate(pass2 = case_when(
    doy == min(doy)~ "1",
    doy == max(doy)~ "2")) %>% ungroup()

# Combine both sources of pass information and create new variable 
sgsl_ichthyo %<>% mutate(pass = case_when(pass1=="1"&pass2=="1"~"1", # if both are 1 then 1
                                          pass1=="2"|pass2=="2"~"2", # if either indicate 2 then 2
                                          is.na(pass2)~pass1))       # if NA then whatever pass1 says 



# actual variable for trajet 1 vs 2
# function to extract last n characters
substrRight <- function(x, n){
  substr(x, nchar(x) - n + 1, nchar(x))
}

sgsl_ichthyo %<>%
  mutate(trajet = substrRight(collector_deployment_id,2))
unique(sgsl_ichthyo$trajet)

# simplify
sgsl_ichthyo %<>% mutate(trajet = 
                           fct_collapse(trajet, 
                                        "1" = c("RT", "T1"),
                                        "2" = c("T2", "T3")))

# Is there a station correction?
sgsl_ichthyo  %<>% 
  mutate(
    correction = case_when(
      stringi::stri_detect_regex(collector_station_name, "cor",case_insensitive=T) ~ "1",
      TRUE ~ "0"
    )
  )

# Side of Bongo net used for measure (i.e. babord or tribord) based on comments. Default is babord
sgsl_ichthyo  %<>% 
  mutate(
    bongo1 = case_when(
      stringi::stri_detect_regex(data_manager_comment, "babord",case_insensitive=F) ~ "babord",
      stringi::stri_detect_regex(data_manager_comment, "tribord",case_insensitive=F) ~ "tribord",
      TRUE ~ "babord"
    )
  )

# Side of Bongo net used for measure (i.e. babord or tribord) based on station name 
sgsl_ichthyo  %<>% 
  mutate(
    bongo2 = case_when(
      stringi::stri_detect_regex(collector_station_name, "_B",case_insensitive=F) ~ "babord",
      stringi::stri_detect_regex(collector_station_name, "_T",case_insensitive=F) ~ "tribord"
    )
  )

sgsl_ichthyo %<>% mutate(bongo = ifelse(is.na(bongo1),bongo2,bongo1)) # defaults to data manager comment variable

# Make sure everything checks out (look for NAs and cases where variables don't match)
df <- distinct(sgsl_ichthyo %>%  dplyr::select(date, consecutive, extra_station, station, collector_station_name, pass1, pass2, pass, data_manager_comment))

# Remove clutter
rm(df)
sgsl_ichthyo %<>% dplyr::select(-pass1,-pass2,-bongo2, -bongo1)


# There are a few additional missions and/or experiments that are found in this dataset 
# Make easy way to find them 1 = yes, 0 = no
sgsl_ichthyo %<>% 
  mutate(extra_1991a = ifelse(collector_station_name %in%
                                c("7.5.00", "7.5.01", "7.5.02", "7.5.03", "7.5.04", "7.5.05", "7.5.06", "7.5.07", "7.5.08", "7.5.09", "7.5.10",
                                  "7.5.11", "7.5.12", "7.5.13", "7.5.14", "7.5.15", "7.5.16", "7.5.17", "7.5.18", "7.5.19", "7.5.20", "7.5.21",
                                  "7.5.22", "7.5.23", "7.5.24", "7.5.25", "7.5.26", "7.5.27", "7.5.28", "7.5.29", "7.5.30", "7.5.31", "7.5.32",
                                  "7.5.33", "7.5.34", "7.5.35"), 1,0),
         extra_1991b = ifelse(year == 1991 & collector_station_name %in% 
                                c("7.4A", "7.4B", "7.4C", "7.4D", "7.4E", "7.4F", "7.4G", "7.4H","7.4I","7.4J", "8.4A", "8.4B", "8.4C", "8.4D","8.4E","8.4F","8.4G","8.4H","8.4I"),1,0),
         extra_1986 = ifelse(year == 1986 & collector_station_name %in% c("13.1","14.1","14.2","14.3","14.4","15.1","15.2","15.3","15.4","15.5","16.1","16.2","16.3","16.4"),1,0), 
         extra_northumberland = ifelse(collector_station_name %in% c("5.0","6.0","7.0","8.0","9.0"),1,0), 
         extra_ns_1998 = ifelse(str_starts(collector_station_name,"SMB"),1,0), 
         extra_st_georges_bay_1994 = ifelse(collector_station_name %in% c("3.11","3.12","3.13","3.14","3.15","3.16","3.17","3.18","3.19","3.20","3.21","3.22","3.23","3.24","3.25","3.26","3.27","3.28","3.29"),1,0), 
         extra_2000 = ifelse(str_starts(collector_station_name,"Supp"),1,0),
         extra_scotian_shelf_2009 = ifelse(grepl(mission_name , pattern="international"),1,0)) %>% 
  
  mutate(extra= rowSums(across(extra_1991a:extra_scotian_shelf_2009))) %>%  filter(extra < 1) %>%  dplyr::select(-c(extra_1991a:extra_scotian_shelf_2009))


# Define new station variable to be able to merge with other data (e.g. temperature from ctd casts)
sgsl_ichthyo  %<>% 
  mutate(
    station_name = str_remove_all(collector_station_name, "_B")) %>% 
  mutate(
    station_name = str_remove_all(station_name, "_T"))

sgsl_ichthyo  %<>% 
  mutate(
    station = str_remove_all(station, "_B")) %>% 
  mutate(
    station = str_remove_all(station, "_T")) %>% 
  mutate(
    station = str_remove_all(station, " COR")) %>% 
  mutate(
    station = str_remove_all(station, "A")) %>% 
  mutate(
    station = str_remove_all(station, "b"))

sgsl_ichthyo %<>% mutate(station = 
                           fct_collapse(station, 
                                        "7.7" = c("7.7", "7.71")))




###adjust extra_consecutive that were not caught 
sgsl_ichthyo %<>% mutate(extra_consecutive = ifelse((year==1983 & station_name==7.71) |
                                                      (year==1989 & station_name==5.4 & consecutive==22) |
                                                      (year==1993 & station_name=="6.5c")|
                                                      (year==1999 & collector_station_name=="1.2_T") |
                                                      (year==2000 & collector_station_name=="6.2_T") |
                                                      (year==2000 & collector_station_name=="7.1_T") |
                                                      (year==2000 & collector_station_name=="2.5b_T"), 0, extra_consecutive),
                         trajet = ifelse((year==1989 & consecutive ==51)|
                                            (year==2002 & station %in% c(2.1,2.2)), "1", trajet))


#check if pass and T correspond
test<- sgsl_ichthyo %>%  mutate(pass = as.numeric(str_remove_all(pass,"[[:alpha:]]")),
                         trajet=as.numeric(str_remove_all(trajet, "[[:alpha:]]")),
                                                diff=pass-trajet,) %>%  filter(diff!=0)

test2 <-  inner_join(test %>%  dplyr::select(year, station), sgsl_ichthyo)
 #most of the time trajet is OK except 89 5.4 and 2002 changed in line above 


sgsl_ichthyo %<>% dplyr::select(-pass)


###adding 1979###has already been processed change format
library(readxl)
X1979 <- read_excel("data/1979.xlsx", sheet = "data") %>%  filter(!(station==3.6 & total==0))

X1979 %<>% pivot_longer(total:n_15, names_to="molt_number", values_to="counts") %>%  
   rename(latitude=lat, longitude=lon, start_depth=sample_depth) %>% 
   mutate(counts=counts/start_depth*volume,#these were real counts/m2 but will be applied later
          split_fraction = 1,
         trajet = "1",
         consecutive=as.character(consec),
         sample_id=paste(year,consecutive),
         station=as.character(station),
         taxonomic_name="Scomber scombrus",
         stage="egg", sex="UNASSIGNED",
         data_manager_comment="Data not verified by data mangement",
         extra_consecutive="1")# not an extra
 
 
#######adding last sampled year if not available in biochem############
if(!is.null(path_to_file)){  
meta <- read_excel(path_to_file, 
                         sheet = "Data")
meta  %<>% ungroup %>%  transmute(year = as.numeric(Année), 
                            month = as.numeric(str_sub(Date,6,7)),
                            day = as.numeric(str_sub(Date,9,10)),
                            date = as.Date(paste(day,month,year, sep = "/" ), format = "%d/%m/%Y"),
                            doy = as.numeric(yday(date)),
                            start_time=format(Heure, format="%H:%M") ,
                            station = as.character(Station),
                            consecutive = as.character(`Conséc.`),
                            bongo = recode(`Bab./Trib.`, `B`="babord", `T`="tribord"),
                            set_duration = `Durée totale (sec)`,
                            start_depth = `Prof. bongo (m)`,
                            volume = `Volume (m3) Bionet`,
                            sounding = `Prof. station (m)`,
                            comments=Comentaires,
                      sample_id= paste(year, consecutive, bongo),) #vérifier si le nom de colonne est le même à chaque année.


  maqytr <- read_excel(path_to_file, 
                               sheet = "Plancton")

  
  p1 <- maqytr %>% mutate(year = as.numeric(Année), 
                                     month = as.numeric(str_sub(Date,6,7)),
                                     day = as.numeric(str_sub(Date,9,10)),
                                     date = as.Date(paste(day,month,year, sep = "/" ), format = "%d/%m/%Y"),
                                     doy = as.numeric(yday(date)),
                                     station = as.character(Station),
                                     consecutive = as.character(`Conséc.`),
                                     bongo = recode(`Bâb./Trib.`, `B`="babord", `T`="tribord"), # notice the â
                                     split_fraction = `SOUS-ECHAN.  Œufs Maq.`,
                                     #total = `Total Oeufs    Maquereau`,
                                     taxonomic_name="Scomber scombrus",
                                     stage="egg",
                          sex="UNASSIGNED",
                          collector_comment=ifelse(`Arrêt du tri (+ de 6h)`=="O", paste("Arrêt du tri car plus de 6hres.", Remarques), Remarques),
                          sample_id= paste(year, consecutive, bongo),) %>% 
    pivot_longer(cols = c(`STADE 1`: `STADE 5`), names_to = "molt_number", values_to = "counts") %>% 
    mutate(start_time=format(Heure, format="%H:%M"))#%
  
  #split différetn pour alrves et oeufs
  
  
  p2 <- maqytr %>% transmute(year = as.numeric(Année), 
                                     month = as.numeric(str_sub(Date,6,7)),
                                     day = as.numeric(str_sub(Date,9,10)),
                                     date = as.Date(paste(day,month,year, sep = "/" ), format = "%d/%m/%Y"),
                                     doy = as.numeric(yday(date)),
                                     station = as.character(Station),
                                     consecutive = as.character(`Conséc.`),
                                     bongo = recode(`Bâb./Trib.`, `B`="babord", `T`="tribord"), # notice the â
                                     split_fraction = `SOUS-ECHAN.      Larves`,
                                     counts = `S. scombrus`,
                                     taxonomic_name="Scomber scombrus",
                                     stage="larvae",molt_number="", sex="UNASSIGNED",
                                     sample_id= paste(year, consecutive, bongo))
                             
                             
  newyear <- full_join(left_join(meta, p1), left_join(meta, p2)) %>%
  mutate(collector_comment = paste(comments, collector_comment)) %>%
  dplyr::select(year:sounding, split_fraction:counts) %>%
  mutate(
    gear_model = "Bongo 61cm",
    mesh_size = 333,
    extra_station = "1", extra_consecutive = "1", trajet = "1", station_name = station, sample_id = paste(year, consecutive),
    molt_number = recode(molt_number, `STADE 1` = "STAGE I",`STADE 2` = "STAGE II",`STADE 3` = "STAGE III",`STADE 4` = "STAGE IV",`STADE 5` = "STAGE V")) %>% 
    filter(!(is.na(taxonomic_name) & bongo=="tribord"))####whatch out for 8.6!!!!!!!!


  coord <- sgsl_ichthyo %>% group_by(station) %>% dplyr::summarise(longitude=mean(longitude,na.rm=T),latitude = mean(latitude,na.rm=T))
 newyear <- left_join(newyear,coord)
  
  

  all_sgsl<- full_join(full_join(left_join(sgsl_ichthyo,bongo),X1979),newyear)
}  

if(is.null(path_to_file)) all_sgsl<- full_join(left_join(sgsl_ichthyo,bongo),X1979)
  
# Strata - defined by P. Ouellet, ResDoc 87/62 # modifié par CL 2022, il y avait des erreurs
  all_sgsl %<>% mutate(stratum = ifelse(station %in%
                                            c("1.1", "1.2", "1.3", "1.4", "1.5", "2.1", "2.2", "2.3", "2.4", "2.5", "2.6", "3.1", "3.2", "3.3", "3.4", "3.5", "3.8","3.9", "4.1", "4.2", "4.9", "5.7","6.8","7.7", "8.1","8.7","9.1", "9.2"), "1", 
                                          ifelse(station %in%
                                                   c("3.6", "3.7","4.3", "4.4", "4.5", "4.6","4.8","5.1", "5.2","5.6","6.1", "6.2","6.3","6.7", "7.1","7.2","7.5", "7.6","8.2","8.6","9.3", "9.4","10.1"), '2',
                                                 ifelse(station %in%
                                                          c("4.7","5.3", "5.4", "5.5","6.4", "6.5", "6.6", "7.3", "7.4","8.3", "8.4", "8.5", "9.5","11.1","12.1"),"3", NA)))) 

  all_sgsl %<>% mutate(stratum_area = ifelse(stratum == "1",29.61e+9,
                                               ifelse(stratum == "2", 21.91e+9,
                                                      ifelse(stratum == "3", 17.93e+9, NA)))) # from ouellet area of stratum modifié il y avait une errreur
#for validation
#  distinct( all_sgsl %>%  filter(year==1987) %>% dplyr::select(latitude, longitude, stratum, year)) %>% 
#  ggplot(aes(x=longitude, y=latitude)) +geom_point(aes(col=stratum))+ facet_wrap(~year, ncol=5)
  

  #fill in volume with NA a partir des volumes notés
  years_to_verify<- distinct(all_sgsl %>%  filter(is.na(volume) & !is.na(counts)) %>%  dplyr::select(year, station, trajet))
  years_to_verify
  x1982<- left_join(years_to_verify, read_excel("S:/Pélagiques/Plancton/Relevés/Relevé 1982/Métadonnées_1982.xls", sheet="Bongo82 P273 et P275") %>% 
          dplyr::mutate(station=as.character(STATION), newvolume=VOLUME, trajet=as.character(TRAJET), comment="Volume not vertified by data mangement") %>%  
         dplyr::select(station, newvolume, trajet, comment))
  
all_sgsl <- full_join(all_sgsl, x1982) %>%  mutate(volume= coalesce(volume, replace=newvolume), data_manager_comment=paste(data_manager_comment, comment))  %>%  dplyr::select(-newvolume, -comment) 
  
  
#####################  Calculate Abundance  ####################

# Calculate concentration and density by dividing the counts by the fraction of the sample treated, dividing by the filtered volume sampled by the gear, and multiplying by the depth sampled
sgsl_ichthyo_L1 <-  all_sgsl %>%  mutate(
  n = counts/split_fraction,
  n_m3 = ifelse(is.na(volume) & n==0, 0,n/volume),#if no eggs found but volume is missing, set Nm2 to zero
  n_m2 = n_m3*(start_depth),
) %>% arrange(date, consecutive) %>%  filter(!is.na(counts))

sgsl_ichthyo_L1<- sgsl_ichthyo_L1[!duplicated(sgsl_ichthyo_L1 %>%  dplyr::select(sample_id, taxonomic_name, stage, molt_number)),] #patching up for 2005 for now. update in biochem pending


#getmetadata for joining filering    
sgsl_ichthyo_meta <-  distinct(sgsl_ichthyo_L1 %>%  filter(!(year==year & stage=="larvae")) %>% #problems with merging with data not in biochem. remove comments for larvae
                                 dplyr::select(-c(split_fraction:consec), -c(n:n_m2)))

#filter file with eggs and larvae
df.filter <- PL_Read_Filter(filter_file = paste0("sql/",year,"/Biochem_maq_eggs_filter_file.txt")) #two last rows of filter file are for 1979

df.data.L2 <- left_join(PL_Taxonomic_Grouping(
  df.data = sgsl_ichthyo_L1 %>% dplyr::select(., sample_id, taxonomic_name, stage, molt_number, sex, n_m2),
  df.filter
), sgsl_ichthyo_meta) %>% filter(extra_consecutive==1) %>% 
  dplyr::select(sample_id, mission_name, consecutive,station,bongo,trajet, date,year:day, doy, start_time, set_duration, latitude:start_depth, sounding:volume, set_duration, mesh_size,stratum, stratum_area, collector_name:mission_comment2, maq_eggs:maq_larvae)


write.table(df.data.L2, paste0("data/",year,"/PL_Bongo_Scomber_eggs_larvae_Counts_L2_", year, ".tsv"), sep="\t", dec=".", row.names=F)
}
