#password and username for Oracle database
source("bdOracle.Rprofile")


year_to_report <- 2022 # maximum year to extract data
year_of_assess <- 2023 # to link to file. year of the assessment

#####################  LOAD PACKAGES  ####################
library(boot) # for inv.logit
library(tidyverse)
library(magrittr)
library(ggpmisc)
library(ggpubr)
source(paste0("R/",year_of_assess,"/utils/basemap.R")) #if error "plot new has not been call", restart R, package compatibility issues
source(paste0("R/",year_of_assess,"/utils/mackerel_fun_incubation.R")) # Mackerel incubation .
source(paste0("R/",year_of_assess,"/utils/spatial_projections.R"))#



####### Data extraction######
source(paste0("R/",year_of_assess,"/utils/extract_biochem.R"))
#for early analysis set current year excel file to read data
#extract_biochem(year_to_report = year_to_report, path_to_file = "S:/Pélagiques/Plancton/Relevés/Relevé 2022/Métadonnées_IML2022-024_2022_11_07.xlsx")
#when data are loaded in Biochem set path_to_file to NULL
extract_biochem(year_to_report = year_to_report, path_to_file = NULL)


## create a lookup file
lookup <- read.delim(paste0("data/",year_of_assess,"/PL_Bongo_Scomber_eggs_larvae_Counts_L2_", year_to_report, ".tsv")) %>%
  group_by(station) %>%
  summarize(
    stratum = first(stratum),
    stratum_area = first(stratum_area),
    depth = mean(sounding, na.rm = T),
    latitude = mean(latitude, na.rm = T),
    longitude = mean(longitude, na.rm = T)
  )

write.table(lookup, "data/lookup_station_egg.txt", row.names = F, sep = "\t", dec = ".")

source(paste0("R/",year_of_assess,"/utils/extract_T0_10.R"))
extract_T0_10(year_to_report = year_to_report, year_of_assess=year_of_assess)
# the resulting file is in data PL_Bongo_Scomber_eggs_larvae_Counts_L2_....

# Stage 1 and 5 mackerel eggs. Response variable of interest is n_m2, i.e. the number of eggs per square meter after standardizing for station depth and volume
egg <- read.delim(paste0("data/",year_of_assess,"/PL_Bongo_Scomber_eggs_larvae_Counts_L2_", year_to_report, ".tsv"))# %>%

View(egg %>% mutate(NAtemp = if_else(is.na(temperature0_10), "absent", "present")) %>% group_by(year, trajet, NAtemp) %>%  tally())

egg <- egg %>%  filter(!year %in% c(1982), !(trajet == 2 & year %in% c(1983,1984, 1985, 1987,1993,2000, 2007))) # 2001, 2006, trajet 2 2000 les données sont là mais on m'a pas la température
#87, 1993 small table size 38, 29

####volumes exploration
ggplot(data=egg %>%  filter(trajet==1), aes(x=set_duration, y=volume))+geom_point(aes(col=start_depth)) +facet_wrap(~year)

ggplot(data=egg %>%  filter(trajet==1), aes(y=volume, x=as.factor(year)))+geom_boxplot()


#incubation time.
egg %<>%
  mutate(
    I = I_lockwood(temperature0_10),
    IM = I_mendiola(temperature0_10),
    DEP = maq_eggs_stage1_5 / I * 24,
    DEPM = maq_eggs_stage1_5 / IM * 24
  ) # Calculate daily egg production (DEP) by station

#find out potential outliers
egg_d <- egg %>% filter(set_duration < 2000) # gives a weird value but egg=0 so no need to remove
lm_out <- lm(volume ~ set_duration, data = egg_d)

simu <- data.frame(set_duration = 0:2000)

simu$fit <- predict(lm_out, newdata = simu, se = T)$fit
simu$se <- predict(lm_out, newdata = simu, se = T)$se

egg_d <- egg_d %>% mutate(
  PA = ifelse(DEP == 0, 0, 1),
  colmatage = ifelse(grepl(collector_comment, pattern = "colm"), "Y", "N")
)

my.formula <- y ~ x # for smooth in ggplot

ggplot(egg_d, aes(y = volume, x = set_duration)) +
  geom_point(aes(col = maq_eggs_stage1_5)) +
  #geom_jitter(aes(col = as.factor(colmatage)), width = 20) +
  #geom_label(aes(label = as.factor(station)), width = 20) +
  theme_few() +
  geom_line(data = simu, aes(x = set_duration, y = fit), col = "red") +
  geom_line(data = simu, aes(y = fit + (3 * se), x = set_duration), col = "red", lty = 2) +
  geom_line(data = simu, aes(y = fit - (3 * se), x = set_duration), col = "red", lty = 2) +
  geom_abline(intercept = lm_out$coefficients[1] - 180, slope = lm_out$coefficients[2], col = "purple") +
  geom_abline(intercept = lm_out$coefficients[1] - 250, slope = lm_out$coefficients[2], col = "green")+
geom_hline(yintercept = 80, lty=2, col="orange")


egg<- egg %>%  mutate(outliers=ifelse(volume < 80, "Y", 
                                ifelse(volume > 800, "Y",
                                       ifelse(set_duration < 500 & volume > 500, "Y",
                                              ifelse(set_duration > 1000 & volume < 200 , "Y", "N")))),
                      outliersEN=ifelse(outliers=="Y",paste0("outliers pass ", trajet), "keep"),
                      outliersFR=ifelse(outliers=="Y",paste0("valeurs aberrantes trajet ", trajet), "conservées"),
                      outliersBI=ifelse(outliers=="Y",paste0("valeurs aberrantes / outliers trajet/pass ", trajet), "conservées / keep"))

po<- ggplot(egg %>%  filter(set_duration < 2000), aes(y = volume, x = set_duration)) +
    #geom_jitter(aes(col = as.factor(colmatage)), width = 20) +
  #geom_label(aes(label = as.factor(station)), width = 20) +
  theme_few() +
  stat_poly_eq(
    formula = my.formula,
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
    parse = TRUE, vstep = 0.08
  )+
  geom_line(data = simu, aes(x = set_duration, y = fit), col = "red") +
  geom_line(data = simu, aes(y = fit + (3 * se), x = set_duration), col = "red", lty = 2) +
  geom_line(data = simu, aes(y = fit - (3 * se), x = set_duration), col = "red", lty = 2) +
  geom_hline(yintercept = 80, lty=2, col="blue")+
  geom_hline(yintercept = 250, lty=5, col="black", lwd=0.2)+
  geom_vline(xintercept = 600, lty=5, col="black", lwd=0.2)+
  theme(legend.position = "bottom",
        legend.background = element_blank())
 

poEN<- po+xlab("Duration (s)") + ylab("Volume (m³)") +  geom_point(aes(col = outliersEN, size=outliersEN)) +
  scale_color_manual(values=c("black", "red", "blue"), name="")+
  scale_size_manual(values=c(0.2,1,1), name="")
ggsave(paste0("figures/",year_of_assess,"/volume_outliersEN.png"), width=6, height=5, dpi=600)
poFR<- po+xlab("Dur\u00e9e (s)") + ylab("Volume (m³)")+  geom_point(aes(col = outliersFR, size=outliersFR)) +
  scale_color_manual(values=c("black", "red", "blue"), name="")+
  scale_size_manual(values=c(0.2,1,1), name="")
ggsave(paste0("figures/",year_of_assess,"/volume_outliersFR.png"), width=6, height=5, dpi=600)
poBI<- po+xlab("Dur\u00e9e / Duration (s)") + ylab("Volume (m³)")+  geom_point(aes(col = outliersBI, size=outliersBI)) +
  scale_color_manual(values=c("black", "red", "blue"), name="")+
  scale_size_manual(values=c(0.2,1,1), name="")
ggsave(paste0("figures/",year_of_assess,"/volume_outliersBI.png"), width=6, height=5, dpi=600)



#  filtering
egg %<>% dplyr::filter(!is.na(DEP), DEP < 1500) # there was one really large value, volume is suspiciously high considering the depth, that could not be predicted with ZAG, removed. 

sample_id_to_remove <- egg %>%  filter((set_duration < 500 & volume >550) |(volume  >800 & set_duration < 1000))

egg <-  egg %>%  filter(!sample_id %in% sample_id_to_remove$sample_id)

egg <- left_join(egg, egg %>% group_by(station) %>% summarize(etopo = mean(sounding, na.rm = T))) %>% mutate(depth = coalesce(sounding, etopo), DEPbackup = DEP)


##calculate 95% of max depth per stations
filtermax<- egg %>%  filter(year > 2001) %>%  group_by(station) %>%  summarize(meandepth= mean(start_depth, na.rm=T),
                                                                               sddepth= sd(start_depth, na.rm=T),
                                                                               max99depth= meandepth + (3.291*sddepth),
                                                                               min99depth= meandepth - (3.291*sddepth)) %>%  
                                                                     dplyr::select(station, max99depth, min99depth)
 

# Data for cross validation for INLA####
year_to_cv <- full_join(egg, expand_grid(station = unique(egg$station), unique(egg %>% dplyr::select(year, trajet)))) %>%
  group_by(year, trajet) %>%
  filter(is.na(DEP), year != year_to_report) %>%
  tally() %>%
  filter(n <= 1)
# set to NA the value you want to fit
eggcv <- egg %>%
  filter(year %in% year_to_cv$year) %>%
  group_by(year) %>%
  sample_n(size = 10) %>%
  mutate(DEP = NA)
eggcv_to_fit <- egg %>% filter(!sample_id %in% eggcv$sample_id)
eggcv_all <- as.data.frame(bind_rows(eggcv, eggcv_to_fit))



########### DEP############
source(paste0("R/",year_of_assess,"/utils/basemap.R"))
source(paste0("R/",year_of_assess,"/INLA/Mesh.R"))
trajet2filter <- data.frame(year = sort(unique(egg$year)), trajet = c(1, 1, 2, 2, 1, 1, 2, 2, 2, 1, 2, 1, 2, 2, 2, 1, 1, rep(1, length(sort(unique(egg$year))) - 17)), trajet2_filter = "keep")

eggt1 <- egg %>% filter(trajet == 1)
#series that include the 2nd pass when avilable. 
eggt2 <- full_join(egg, trajet2filter) %>% filter(trajet2_filter == "keep", !is.na(station))

## figure for area
lookup <- read.delim("data/lookup_station_egg.txt")

basemap2 +
  geom_path(data = fortify(bat50), aes(x = long, y = lat, group = group), col = "grey45", linewidth = 0.4) +
  geom_path(data = fortify(bat200n), aes(x = long, y = lat, group = group), col = "grey85", linewidth = 0.4) +
  geom_polygon(data = fortify(nec), aes(x = long, y = lat, group = group), fill = "bisque2", col = "bisque4") +
  geom_label(data = lookup, aes(x = longitude, y = latitude, label = station), size = 3, label.padding = unit(0.1, "lines"))

ggsave(paste0("figures/",year_of_assess,"/area_of_interest.png"), width = 6, height = 6, units = "in", dpi = 600)
# full_join(eggt1, eggt2 %>%  filter(trajet==2)) %>% ggplot()+ geom_histogram(aes(x=year, group=as.factor(trajet),fill=as.factor(trajet)), position="dodge",stat="count")
tableN <- full_join(eggt1, eggt2 %>% filter(trajet == 2)) %>%
  group_by(year, trajet) %>%
  tally() %>%
  pivot_wider(names_from = trajet, values_from = n)

write.table(tableN, paste0("results/",year_of_assess,"/tableN_trajets.txt"), sep = "\t", dec = ".", row.names = F)

###INLA mesh###
loc_egg1 <- make_boundary(dat = eggt1, year_to_report, trajet = 1)
make_mesh(range = 150, Loc = loc_egg1, year_to_report, trajet = 1)

loc_egg2 <- make_boundary(dat = eggt2, year_to_report, trajet = 2)
make_mesh(range = 150, Loc = loc_egg2, year_to_report, trajet = 2)

# prediction grid for INLA
predstation <- expand_grid(distinct(egg %>% dplyr::select(year)), egg %>% group_by(station) %>% dplyr::summarise(x = mean(longitude, na.rm = T), y = mean(latitude, na.rm = T)) %>% mutate(station = as.factor(station)))
coord <- create_coord(predstation[, c("x", "y")], crsout = lcc)
predstation <- predstation %>% mutate(X.m = coord$x, Y.m = coord$y)
saveRDS(predstation, paste0("data/",year_of_assess,"/INLA/prediction_grid_station", year_to_report, ".rds"), version = 2)

#Test with AR in INLA, there is no difference not to be done in  later stock assessment
#source("R/explorations_archived/INLA_ZAG_covar_Stations_testAR.R")
#miss<- full_join(expand_grid(year=c(1995, 1997, 2020), station=unique(egg$station)), lookup)
#INLA_ZAGAR(dat=full_join(eggt1, miss) , res_mesh=150, res_pred="station", Rvar="DEP", year_to_report, trajet=1, serie="")
  

#######execute INLA model##########
source(paste0("R/",year_of_assess,"/INLA/INLA_ZAG_covar_Stations.R"))
INLA_ZAG(dat = eggt1, res_mesh = 150, res_pred = "station", year_to_report,year_of_assess, Rvar = "DEP", trajet = 1)
#estimating density instead of DEP
INLA_ZAG(dat = eggt1, res_mesh = 150, res_pred = "station", year_to_report,year_of_assess, Rvar = "maq_eggs_stage1_5", trajet = 1) #no effect
INLA_ZAG(dat = eggt1 %>%  filter(outliers=="N"), res_mesh = 150, res_pred = "station", year_to_report, year_of_assess,Rvar = "DEP", trajet = 1, serie="outliers_rm")
INLA_ZAG(dat = eggt1, res_mesh = 150, res_pred = "station", year_to_report,year_of_assess, Rvar = "DEPM", trajet = 1)
INLA_ZAG(dat = eggt2, res_mesh = 150, res_pred = "station", year_to_report,year_of_assess, Rvar = "DEP", trajet = 2)
INLA_ZAG(dat = eggt2, res_mesh = 150, res_pred = "station", year_to_report, year_of_assess,Rvar = "DEPM", trajet = 2)
#outliers of start depth removed and estimated. 0 stays 0
eggt1filter<- left_join(eggt1, filtermax) %>%  mutate(outliers_maxdepth= if_else(start_depth <= max99depth, "N",
                                                                              if_else(start_depth > max99depth & DEP ==0, "N","Y")),
                                                      outliers_mindepth= if_else(start_depth >= min99depth, "N",
                                                                              if_else(start_depth < min99depth  & depth-start_depth <10, "N","Y")),  # remainining NA are all > 50m so second criteria does not apply
                                                      outliers_depth= if_else(outliers_maxdepth=="Y", "Depth > CI",
                                                                             if_else(outliers_mindepth =="Y", "Depth < CI",
                                                                                     "Depth within CI"))
                                                      )

eggt1filter<-eggt1filter %>%  mutate(outliers_depth = if_else(DEP==0, "No egg", outliers_depth),
                                      outliers_depth= recode_factor(outliers_depth, `Depth within CI`="Depth within CI",`Depth < CI`="Depth < CI", `Depth > CI`="Depth > CI", `No egg`="No egg"),
                                      outliers_depthFR= recode_factor(outliers_depth, `Depth within CI` = "Profondeur à l'intérieur des IC", 
                                      `Depth < CI` = "Profondeur < IC",
                                      `Depth > CI` = "Profondeur > IC", `Sans oeuf`="Sans oeuf"))

#map of outliers values:
basemap2 + geom_point(data=eggt1filter, aes(col=outliers_depth,shape=outliers_depth,x=longitude, y=latitude))+facet_wrap(~year) + scale_shape_manual(values=c(1,16, 16,3), name="") +scale_color_manual(values=c("black","dodgerblue2", "darkred", "black"), name="") +xlab("") + ylab("") +theme(legend.position = c(0.5,0.07), legend.text = element_text(size=15))
ggsave(paste0("figures/",year_of_assess,"/position_outliers_depths.png"), dpi=600, units="in",width=12, height=13)

basemap2 + geom_point(data=eggt1filter, aes(col=outliers_depthFR,shape=outliers_depthFR,x=longitude, y=latitude))+facet_wrap(~year) + scale_shape_manual(values=c(1,16, 16,3), name="") +scale_color_manual(values=c("black","dodgerblue2", "darkred", "black"), name="") +xlab("") + ylab("") +theme(legend.position = c(0.5,0.07), legend.text = element_text(size=15))
ggsave("figures/",year_of_assess,"position_outliers_depthsFR.png", dpi=600, units="in",width=12, height=13)



##histogram of outliers values
poEN <- eggt1filter %>%  
            ggplot()+geom_histogram(aes(x=start_depth, y = ifelse(after_stat(count) > 0, after_stat(count), NA),fill=outliers_depth, col=outliers_depth), breaks=seq(1,201, 10)) + geom_vline(xintercept=51, lty=2, lwd=0.8, col="grey")+
            facet_wrap(~year) +scale_fill_manual(name="",values=c("white","dodgerblue2", "darkred", "black"))+
  scale_color_manual(name="",values=c("black","dodgerblue2", "darkred", "black"))+
basetheme + ylab("N stations")+ xlab("Maximum depth sampled (m)")  +theme(legend.position = c(0.5,0.07), legend.text = element_text(size=12))
            ggsave("figures/",year_of_assess,"outliers_depth_histoEN.png", dpi=600, units="in",width=10, height=10)
                         
 poFR <- eggt1filter %>%  
              ggplot()+geom_histogram(aes(x=start_depth, y = ifelse(after_stat(count) > 0, after_stat(count), NA),fill=outliers_depthFR, col=outliers_depthFR), breaks=seq(1,201, 10)) + geom_vline(xintercept=51, lty=2, lwd=0.8, col="grey")+
              facet_wrap(~year) +scale_fill_manual(name="",values=c("white","dodgerblue2", "darkred","black"))+
              scale_color_manual(name="",values=c("black","dodgerblue2", "darkred","black"))+
              basetheme + ylab("N stations")+ xlab("Profondeur maximale échantillonnée (m)")  +theme(legend.position = c(0.5,0.07), legend.text = element_text(size=12))
            ggsave("figures/",year_of_assess,"outliers_depth_histoFR.png", dpi=600, units="in",width=10, height=10)
      
INLA_ZAG(dat = eggt1filter %>%  filter(outliers_depth=="Depth within CI"), res_mesh = 150, res_pred = "station", year_to_report, Rvar = "DEP", trajet = 1, serie="outliers_depth_rm")


############# INLA cross-validation###########
#do not need to do every year.
source(paste0("R/",year_of_assess,"INLA/INLA_ZAG_covar_Stations_CV.R"))
INLA_ZAG_CV(dat = eggcv_all %>% filter(trajet == 1), res_mesh = 150, res_pred = "station", trajet = 1, year_of_assess)
# to get result use summary fitted values
load(paste0("results/",year_of_assess,"/INLA/DEP_INLA_mesh150km_pred_station_ZAG_",year_to_report,"_CV_trajet1.RData"))
# Extract mean from fitted values
idb.prd <- inla.stack.index(stkallBern, "BernFit")$dat
Pi <- inv.logit(Gbern$summary.fitted.values$mean[idb.prd]) # il n'y a pas d'index pour prediction, j'ai setté à NA.
# sd.prd   <- I$summary.fitted.values$sd[id.prd]
mu <- Ggamma$summary.fitted.values$mean[idb.prd]
dat$Fit <- Pi * mu
my.formula <- y ~ x # for smooth in ggplot

inlastation1 <- dat %>%
  ungroup() %>%
  filter(sample_id %in% eggcv$sample_id) %>%
  dplyr::rename(inlastation = Fit) %>%
  dplyr::select(sample_id, inlastation)


pcv <- full_join(egg %>% dplyr::select(sample_id, year, DEP) %>% filter(sample_id %in% inlastation1$sample_id), inlastation1) %>%
  ggplot(aes(x = DEP+1, y = inlastation+1)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_few() +
  stat_poly_eq(
    formula = my.formula,
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
    parse = TRUE, vstep = 0.08
  ) +geom_abline(slope=1, intercept=0, col="red", lty=2)
##attention il faut sauvegarder eggcv pour que cette figure soit reproductible.

pcvEN <- pcv + scale_x_continuous(name = "Observed daily egg production (N/m²/day)", trans = "log10") + scale_y_continuous(name = "Predicted daily egg production (N/m²/day)", trans = "log10")
ggsave(paste0("figures/",year_of_asses,"/INLA/validation", year_to_report, "EN.png"), width = 6, height = 4, units = "in", dpi = 600)
pcvFR <- pcv + scale_x_continuous(name = "Production journali\u00E9re d'oeufs observ\u00E9e (N/m²/jour)",  trans = "log10") + scale_y_continuous(name = "Production journali\u00E9re d'oeufs pr\u00E9dite (N/m²/jour)", trans = "log10")
ggsave(paste0("figures/",year_of_asses,"/INLA/validation", year_to_report, "FR.png"), width = 6, height = 4, units = "in", dpi = 600)
pcvBI <- pcv + scale_x_continuous(name = "Production journali\u00E9re d'oeufs observ\u00E9e (N/m²/jour)\nObserved daily egg production (N/m²/day)",  trans = "log10") + scale_y_continuous(name = "Production journali\u00E9re d'oeufs pr\u00E9dite (N/m²/jour)\nPredicted daily egg production (N/m²/day)", trans = "log10")
ggsave(paste0("figures/",year_of_asses,"/INLA/validation", year_to_report, "BI.png"), width = 6, height = 5, units = "in", dpi = 600)
# end of INLA cross-validation

########INLA get results##########
#run model diagnostics, export random field, get predictions 
source(paste0("R/",year_of_assess,"R/INLA/INLA_results.R")) 
#no sensitivity. this is the one to run every year
INLA_results(res_mesh = 150, res_pred = "station", trans = NULL, myfamily = "ZAG", minw = -10, maxw = 10, smoother = F, Rvar = "DEP", year_to_report = year_to_report, year_of_assess=year_of_assess,trajet = 1)
#sensitivty to outliers
INLA_results(res_mesh = 150, res_pred = "station", trans = NULL, myfamily = "ZAG", minw = -10, maxw = 10, smoother = F, Rvar = "DEP", year_to_report = year_to_report, year_of_assess=year_of_assess,trajet = 1, serie="outliers_rm")
#sensitivity to development time
INLA_results(res_mesh = 150, res_pred = "station", trans = NULL, myfamily = "ZAG", minw = -10, maxw = 10, smoother = F, Rvar = "DEPM", year_to_report = year_to_report, year_of_assess=year_of_assess,trajet = 1)
#sensitivity to timing of survey
INLA_results(res_mesh = 150, res_pred = "station", trans = NULL, myfamily = "ZAG", minw = -10, maxw = 10, smoother = F, Rvar = "DEP", year_to_report = year_to_report, year_of_assess=year_of_assess,trajet = 2)
#combination not done
#INLA_results(res_mesh = 150, res_pred = "station", trans = NULL, myfamily = "ZAG", minw = -10, maxw = 10, smoother = F, Rvar = "DEPM", year_to_report = year_to_report, trajet = 2)
#on density
INLA_results(res_mesh = 150, res_pred = "station", trans = NULL, myfamily = "ZAG", minw = -10, maxw = 10, smoother = F, Rvar = "maq_eggs_stage1_5", year_to_report = year_to_report, year_of_assess=year_of_assess,trajet = 1)
#otuliers in depth estimated
INLA_results(res_mesh = 150, res_pred = "station", trans = NULL, myfamily = "ZAG", minw = -10, maxw = 10, smoother = F, Rvar = "DEP", year_to_report = year_to_report, year_of_assess=year_of_assess,trajet = 1, serie="outliers_depth_rm")



########DEP#########
#INLA predictions substitute missing value to yearly DEP
#function to add a variance column that includes the variance on the INLA estimations.
source(paste0("R/",year_of_assess,"/INLA/getvar.R"))
#####test with density######
predd <- read.delim(paste0("results/",year_of_assess,"INLA/predictions/predictions_maq_eggs_stage1_5_INLA_mesh150km_pred_station_ZAG", year_to_report, "_tajet1.txt"))
DEPda <- full_join(eggt1, left_join(predd, lookup %>% dplyr::select(-c(depth, latitude, longitude)))) %>% mutate(new_egg = coalesce(maq_eggs_stage1_5, Fit),
                                                                                                                 newDEP= new_egg/ I * 24,)

DEPd <- DEPda %>%
                     mutate(
                       trajet = 1) %>%
                     group_by(year, trajet) %>%
                     dplyr::summarize(
                       #DEP.wm = weighted.mean(DEP, stratum_area, na.rm = T),
                       DEP.p = mean(DEP, na.rm = T),
                       )%>%  mutate(test_sens="Density")


####No sensitivty baseline ####
predt1 <- read.delim(paste0("results/",year_of_assess,"/INLA/predictions/predictions_DEP_INLA_mesh150km_pred_station_ZAG", year_to_report, "_tajet1.txt"))
DEPt1a <- full_join(eggt1, left_join(predt1, lookup %>% dplyr::select(-c(depth, latitude, longitude)))) %>% mutate(DEP = coalesce(DEPbackup, Fit),
                                                                                                                   var1 = ifelse(!is.na(DEPbackup), 0, var))

test1a <- DEPt1a %>%  group_by(year) %>%  summarize(var3= mean(var))  

vart1a<- getvar(DEPt1a)
DEPt1 <- full_join(DEPt1a %>%
  mutate(
    trajet = 1) %>%
  group_by(year, trajet) %>%
  dplyr::summarize(
    #DEP.wm = weighted.mean(DEP, stratum_area, na.rm = T),
    DEP.p = mean(DEP, na.rm = T),
    DEP.t = mean(Fit, na.rm = T),
  ),
  vart1a) %>%  mutate(test_sens="No sensitivity")



#check density
check_density<- full_join(DEPd, DEPt1 %>%  dplyr::select(colnames(DEPd)))
ggplot(check_density,aes(x=year, y=DEP.p))+geom_point(aes(col=test_sens))+theme(legend.position = c(0.7, 0.9))

ggsave(paste0("figures/",year_of_assess,"/DEP_test_avec_density_INLA.png"), width=6, height=4, dpi=600, units="in")


#####Sensitivity outliers#########
predto <- read.delim(paste0("results/",year_of_assess,"/INLA/predictions/predictions_DEP_INLA_mesh150km_pred_station_ZAG", year_to_report, "_tajet1outliers_rm.txt"))
DEPtoa <- full_join(eggt1 %>%  mutate(DEPbackup2 =DEPbackup, # for figure below
                                      DEPbackup=ifelse(outliers!="N", NA, DEPbackup)), left_join(predto, lookup %>% dplyr::select(-c(depth, latitude, longitude)))) %>% mutate(DEP = coalesce(DEPbackup, Fit),
                                                                                                                                                                               var1 = ifelse(!is.na(DEPbackup), 0, var))
vartoa<- getvar(DEPtoa)
DEPto <- full_join(DEPtoa %>%
  mutate(
    trajet = 1) %>%
  group_by(year, trajet) %>%
  dplyr::summarize(
    DEP.p = mean(DEP, na.rm = T),
    DEP.t = mean(Fit, na.rm = T),
   ),
  vartoa) %>%  mutate(test_sens="Outliers volume estimated")

my.formula <- y ~ x # for smooth in ggplot


pof<- DEPtoa %>% filter(outliers=="Y") %>%  ggplot(aes(x=DEPbackup2, y=Fit))+ geom_point()+theme_few()+
  geom_smooth(method="lm") +stat_poly_eq(
  formula = my.formula,
  aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
  parse = TRUE, vstep = 0.08
) + geom_abline(slope=1, intercept=0.0000000001, col="red", lty=2)

pofEN <-  pof + scale_x_continuous(name = "Observed daily egg production (N/m²/day)", trans =scales::pseudo_log_trans(base=10), breaks=c(0,1,10,100,1000,2000)) + scale_y_continuous(name = "Predicted daily egg production (N/m²/day)", trans = scales::pseudo_log_trans(base=10), breaks=c(0,1,10,100,1000,2000))
ggsave(paste0("figures/",year_of_assess,"/estimated_outliersEN.png"),  width=6, height=4, dpi=600, units="in")
pofFR <-  pof + scale_x_continuous(name = "Production journali\u00E9re d'oeufs observ\u00E9e (N/m²/jour)", trans = scales::pseudo_log_trans(base=10), breaks=c(0,1,10,100,1000,2000)) + scale_y_continuous(name = "Production journali\u00E9re d'oeufs pr\u00E9dite (N/m²/jour)", trans = scales::pseudo_log_trans(base=10), breaks=c(0,1,10,100,1000,2000))
ggsave(paste0("figures/",year_of_assess,"/estimated_outliersFR.png"),  width=6, height=4, dpi=600, units="in")




#####Sensitivity outliers depth #########
predtod <- read.delim(paste0("results/",year_of_assess,"/INLA/predictions/predictions_DEP_INLA_mesh150km_pred_station_ZAG", year_to_report, "_tajet1outliers_depth_rm.txt"))
DEPtoda <- full_join(eggt1filter %>%  mutate(DEPbackup2 =DEPbackup, # for figure below
                                      DEPbackup=ifelse(outliers_depth!="Depth within CI", NA, DEPbackup)), left_join(predtod, lookup %>% dplyr::select(-c(depth, latitude, longitude)))) %>% mutate(DEP = coalesce(DEPbackup, Fit),
                                                                                                                                                                               var1 = ifelse(!is.na(DEPbackup), 0, var))
vartoda<- getvar(DEPtoda)
DEPtod <- full_join(DEPtoda %>%
                     mutate(
                       trajet = 1) %>%
                     group_by(year, trajet) %>%
                     dplyr::summarize(
                       DEP.p = mean(DEP, na.rm = T),
                       DEP.t = mean(Fit, na.rm = T),
                     ),
                   vartoa) %>%  mutate(test_sens="Outliers depth estimated")

my.formula <- y ~ x # for smooth in ggplot


podf<- DEPtoda %>% filter(outliers_depth!="Depth within CI") %>%  ggplot(aes(x=DEPbackup2, y=Fit, col=outliers_depth))+ geom_point()+theme_few()+
  geom_smooth(method="lm") +stat_poly_eq(
    formula = my.formula,
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
    parse = TRUE, vstep = 0.08
  ) + geom_abline(slope=1, intercept=0.0000000001, col="red", lty=2)

podfEN <-  podf + scale_x_continuous(name = "Observed daily egg production (N/m²/day)", trans =scales::pseudo_log_trans(base=10), breaks=c(0,1,10,100,1000,2000)) + scale_y_continuous(name = "Predicted daily egg production (N/m²/day)", trans = scales::pseudo_log_trans(base=10), breaks=c(0,1,10,100,1000,2000)) + scale_color_manual(name="", values=c("dodgerblue2", "darkred")) + theme(legend.position=c(0.8, 0.15))
ggsave(paste0("figures/",year_of_assess,"/estimated_outliers_depthEN.png"),  width=6, height=4, dpi=600, units="in")
podfFR <-  podf + scale_x_continuous(name = "Production journali\u00E9re d'oeufs observ\u00E9e (N/m²/jour)", trans = scales::pseudo_log_trans(base=10), breaks=c(0,1,10,100,1000,2000)) + scale_y_continuous(name = "Production journali\u00E9re d'oeufs pr\u00E9dite (N/m²/jour)", trans = scales::pseudo_log_trans(base=10), breaks=c(0,1,10,100,1000,2000)) + scale_color_manual(name="", values=c("dodgerblue2", "darkred"), breaks=c("Depth < CI", "Depth > CI"), labels=c("Profondeur < IC", "Profondeur > IC"))+ theme(legend.position=c(0.8, 0.15))
ggsave(paste0("figures/",year_of_assess,"/estimated_outliers_depthFR.png"),  width=6, height=4, dpi=600, units="in")


#########Sensitivity timing#########
predt2 <- read.delim(paste0("results/",year_of_assess,"/INLA/predictions/predictions_DEP_INLA_mesh150km_pred_station_ZAG", year_to_report, "_tajet2.txt"))
DEPt2a <- full_join(eggt2, left_join(left_join(predt2, trajet2filter %>%  dplyr::select(-trajet2_filter)), lookup %>% dplyr::select(-c(depth, latitude, longitude)))) %>% mutate(DEP = coalesce(DEPbackup, Fit),
                                                                                                                                                                                 var1 = ifelse(!is.na(DEPbackup), 0, var))
vart2a<- getvar(DEPt2a)
DEPt2 <- full_join(DEPt2a %>%
  mutate(
    trajet = 2) %>%
  group_by(year, trajet) %>%
  dplyr::summarize(
    DEP.p = mean(DEP, na.rm = T),
    DEP.t = mean(Fit, na.rm = T),
  ),
  vart2a) %>%  mutate(test_sens="Survey timing")


#######Sensitivity developement########
predm <- read.delim(paste0("results/",year_of_assess,"/INLA/predictions/predictions_DEPM_INLA_mesh150km_pred_station_ZAG", year_to_report, "_tajet1.txt"))
DEPMa <- full_join(eggt1, left_join(predm,lookup %>% dplyr::select(-c(depth, latitude, longitude)))) %>% mutate(DEPM = coalesce(DEPM, Fit),
                                                                                                                var1 = ifelse(!is.na(DEPbackup), 0, var))

vartMa<- getvar(DEPMa)
DEPM1 <- full_join(DEPMa %>%
  mutate(
    trajet = 1) %>%
  group_by(year, trajet) %>%
  dplyr::summarize(
    DEP.p = mean(DEPM, na.rm = T),
    DEP.t = mean(Fit, na.rm = T),
  ),
  vartMa)  %>%  mutate(test_sens="Egg development")

#plot of DEP without sensitivity for presentation##

####Effect of incubation time on DEP.
all_years= data.frame(year=seq(1979, year_to_report, 1))
# to do in TEP
right_join(list(DEPt1 ,DEPM1) %>%
  reduce(full_join), expand_grid(all_years, test_sens=c("No sensitivity", "Egg development"), trajet=1)) %>% 
  ungroup() %>%
  mutate(test_sens = recode_factor(test_sens, `No sensitivity` = "Lockwood and Nichols (1977)", `Egg development` = "Mendiola et al. (2006)")) %>%
  mutate(yearjit=ifelse(test_sens=="Lockwood and Nichols (1977)", year, year+0.3)) %>% 
  ggplot(aes(x = yearjit, y = DEP.p, col = test_sens)) +
  geom_point() +
 # geom_errorbar(aes(ymin=DEP.p- 1.96 *sep, ymax=DEP.p+1.96*sep))+
  geom_line() +
  theme_few() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.8, 0.9)
  ) +
  ylab("DEP") +
  xlab("Year") +
  scale_color_manual(values = c("black", "grey55"))
#ggsave("figures/exploration/effect_of_incubation_time_on_DEP.png", width = 7, height = 5, dpi = 600, units = "in")



right_join(DEPt1, expand_grid(all_years)) %>% 
  ungroup() %>%
   ggplot(aes(x = year, y = DEP.p)) +
  geom_point() +
 # geom_errorbar(aes(ymin=DEP.p- 1.96 *sep, ymax=DEP.p+1.96*sep), width=0.8)+
  geom_line() +
  theme_few() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.8, 0.9)
  ) +
  ylab("DEP") +
  xlab("Year") +
  scale_color_manual(values = c("black", "grey55"))+
  ylab("Production journalière d'oeufs (N/m²/jour)\nDaily egg production (N/m²/day)")+
  xlab("Ann\u00E9e\nYear")

ggsave(paste0("figures/",year_of_assess,"/DEP_trajet1_year_BI.png"),width = 8, height = 5, unit = "in", dpi = 600)




DEP <- full_join(lookup %>%  dplyr::rename(lat=latitude, lon=longitude) %>% dplyr::select(station, lat, lon),
                 DEPt1a %>% mutate(estimated = ifelse(is.na(DEPbackup) & !is.na(DEP), "estimated", "observed"))) %>% 
     mutate(latitude=coalesce(latitude, lat),
                          longitude=coalesce(longitude,lon))
basemap2 + 
  geom_point(data = DEP %>%  filter(DEP > 0), aes(x = longitude, y = latitude, size = DEP, col = estimated), shape = 1) + 
  geom_point(data = DEP %>%  filter(DEP == 0), aes(x = longitude, y = latitude), shape = 3, size=0.8) + 
  facet_wrap(~year) + 
  scale_color_manual(values = c("red", "black"), guide = "none") + 
  theme_few()+ theme(legend.position="bottom")+
  scale_size_continuous(name="",breaks=c(seq(200, ceiling(max(DEP$DEP)), 200)),range=c(0.25,5))

ggsave(paste0("figures/",year_of_assess,"/DEP.png"), width=12, height=13, dpi=600, units="in")

basemap2 + 
  geom_point(data = DEP %>%  filter(DEP >0, year<2003), aes(x = longitude, y = latitude, size = DEP, col = estimated), shape = 1) + 
  geom_point(data = DEP%>%  filter(DEP ==0, year<2003), aes(x = longitude, y = latitude), shape = 3, size=0.8) + 
  facet_wrap(~year, ncol=4) + 
  scale_color_manual(values = c("red", "black"), guide = "none") + 
  theme_few()+ theme(legend.position="bottom")+
  scale_size_continuous(name="",breaks=c(seq(200, ceiling(max(DEP$DEP)), 200)),range=c(0.25,5), limits=c(0,1500))

ggsave(paste0("figures/",year_of_assess,"/DEP_1979_2002.png"), width=9, height=12, dpi=600, units="in")


basemap2 + 
  geom_point(data = DEP %>%  filter(DEP >0, year>=2003), aes(x = longitude, y = latitude, size = DEP, col = estimated), shape = 1) + 
  geom_point(data = DEP%>%  filter(DEP ==0, year>=2003), aes(x = longitude, y = latitude), shape = 3, size=0.8) + 
  facet_wrap(~year, ncol=4) + 
  scale_color_manual(values = c("red", "black"), guide = "none") + 
  theme_few()+ theme(legend.position="bottom")+
  scale_size_continuous(name="",breaks=c(seq(200, ceiling(max(DEP$DEP)), 200)),range=c(0.25,5), limits=c(0,1500))

ggsave(paste0("figures/",year_of_assess,"/DEP_2003_",year_to_report,".png"), width=9, height=12, dpi=600, units="in")




dep2017<- DEP %>%  filter(year>year_to_report-6)

basemap2 + 
  geom_point(data = dep2017 %>%  filter(DEP >0), aes(x = longitude, y = latitude, size = DEP, col = estimated), shape = 1) + 
  geom_point(data = dep2017 %>%  filter(DEP ==0), aes(x = longitude, y = latitude), shape = 3, size=0.8) + 
  facet_wrap(~year, ncol=3) + 
  scale_color_manual(values = c("red", "black"), guide = "none") + 
  theme_few()+ theme(legend.position="bottom")+
  scale_size_continuous(name="",breaks=c(seq(50, ceiling(max(dep2017$DEP)), 50)),range=c(0.25,5))

ggsave(paste0("figures/",year_of_assess,"/DEP_last_5_years.png"), width=6, height=6, dpi=600, units="in")

basemap2 + 
  geom_point(data = dep2017 %>%  filter(Fit >0), aes(x = longitude, y = latitude, size = Fit), col = "red", shape = 1) + 
  geom_point(data = dep2017 %>%  filter(Fit ==0), aes(x = longitude, y = latitude), shape = 3, size=0.8) + 
  facet_wrap(~year, ncol=3) + 
  scale_color_manual(values = c("red", "black"), guide = "none") + 
  theme_few()+ theme(legend.position="bottom")+
  scale_size_continuous(name="",breaks=c(seq(50, ceiling(max(dep2017$Fit)), 50)),range=c(0.25,5))

ggsave("figures/DEP_last_5_years_FIT.png", width=6, height=6, dpi=600, units="in")


#exploration for survey timing
eggtrajet<- egg %>%  filter(year %in% c(1988,1989,1990,1992,1994,1996,1998)) 

c1<- basemap2 + geom_point(data = eggtrajet, aes(x = longitude, y = latitude, fill = temperature0_10), shape = 21) +
  facet_grid(year ~ trajet) + 
  scale_fill_viridis_c(option = "turbo", name = "") + 
  theme(legend.position = "bottom")

DEP_checkt <- full_join(DEP %>% dplyr::select(-lat, -lon), full_join(lookup %>%  dplyr::rename(lat=latitude, lon=longitude) %>% dplyr::select(station, lat, lon),
                 DEPt2a %>% filter(trajet==2) %>% mutate(estimated = ifelse(is.na(DEPbackup) & !is.na(DEP), "estimated", "observed")))) %>% 
  mutate(latitude=coalesce(latitude, lat),
         longitude=coalesce(longitude,lon)) %>% filter(year %in% c(1988,1989,1990,1992,1994,1996,1998), trajet %in% c(1,2)) 

c2<- basemap2 + 
  geom_point(data = DEP_checkt %>%  filter(DEP >0), aes(x = longitude, y = latitude, size = DEP, col = estimated), shape = 1) + 
  geom_point(data = DEP_checkt%>%  filter(DEP ==0), aes(x = longitude, y = latitude), shape = 3, size=0.8) + 
  facet_grid(year ~ trajet) + 
  scale_color_manual(values = c("red","black"), guide = "none") + 
   theme(legend.position="bottom")+
  scale_size_continuous(name="",breaks=c(seq(200, ceiling(max(DEP$DEP)), 200)),range=c(0.25,5))
#il manque les valeurs estimées sur cette figure

ggarrange(c1, c2, ncol=2)
#ggsave("figures/exploration/pass_Temp_DEP.png", width=7, height=10, units="in", dpi=600)


#figure for sensitivity
sensM<- full_join(DEPt1a %>%  dplyr::select(sample_id,trajet,year, maq_eggs_stage1_5, I, DEP),
          DEPMa %>%  dplyr::select(sample_id,trajet, year, maq_eggs_stage1_5, IM, DEPM))

summary(sensM %>%  filter(trajet==1, maq_eggs_stage1_5!=0) %>%  dplyr::select(I, IM))
sensM %>%  filter(trajet==1, maq_eggs_stage1_5!=0) %>%  dplyr::select(I, IM) %>%  summarize(sd(I),sd(IM))

summary(sensM %>%  filter(trajet==1) %>%  dplyr::select(DEP, DEPM))
sensM %>%  filter(trajet==1) %>%  dplyr::select(DEP, DEPM) %>%  summarize(sd(DEP),sd(DEPM))

sensM %>%  filter(year %in% c(2021, 2022)) %>%  group_by(year)%>%  summarize(max(DEP),max(DEPM))
sensM %>%  filter(year %in% c(2021, 2022)) %>%  group_by(year)%>%  summarize(mean(DEP),mean(DEPM))


b1<- full_join(expand_grid(year=c(1980:1982,1995,1997,2020), I=NA, IM=NA), sensM %>% filter(trajet==1, maq_eggs_stage1_5!=0) %>%  dplyr::select(sample_id, year,I, IM) %>%  pivot_longer(I:IM) %>%  
                 mutate(name=recode_factor(name, I="Lockwood et al. (1977)", IM="Mendiola et al. (2006)"))) %>%  
  ggplot()+ geom_boxplot(aes(x=as.factor(year), y=value, col=name), outlier.size=0.5) + scale_x_discrete(breaks=seq(1980, year_to_report,5))+theme_few()+ theme(legend.position="bottom", legend.title=element_blank())

b1EN <-  b1 + xlab("Year") +ylab("Incubation time (hours)")

b1FR <-  b1 + xlab("Ann\u00E9e") +ylab("Temps d'incubation (heures)")



b2<- full_join(expand_grid(year=c(1980:1982,1995,1997,2020), DEP=NA, DEPM=NA), sensM %>% filter(trajet==1) %>%  dplyr::select(sample_id, year,DEP, DEPM) %>%  pivot_longer(DEP:DEPM) %>%  
                 mutate(name=recode_factor(name, DEP="Lockwood et al. (1977)", DEPM="Mendiola et al. (2006)"))) %>%  
  ggplot()+ geom_boxplot(aes(x=as.factor(year), y=value, col=name), outlier.size=0.5) + scale_x_discrete(breaks=seq(1980, year_to_report,5))+theme_few()+ theme(legend.position="bottom", legend.title=element_blank())

b2EN <-  b2 + xlab("Year") +ylab("Daily egg production (N/m²/day)")
ggarrange(b1EN, b2EN, nrow=2, common.legend = T, legend="bottom")
ggsave(paste0("figures/",year_of_assess,"/DEP_boxplotEN.png", width=8, height=9, dpi=600, units="in"))  

b2FR <-  b2 + xlab("Ann\u00E9e") +ylab("Production journalière d'oeufs (N/m²/jour)")
ggarrange(b1FR, b2FR, nrow=2, common.legend = T, legend="bottom")
ggsave(paste0("figures/",year_of_assess,"/DEP_boxplotFR.png"), width=8, height=9, dpi=600, units="in")  



######### TEP###############
source(paste0("R/",year_of_assess,"/prop_spawn.R"))
prop_spawn(year_to_report,year_of_assess=year_of_assess, dat = egg, cv = T) # dat is used to calculate median date and in split with trajet.
prop_spawn(year_to_report, year_of_assess=year_of_assess,dat = egg, cv = T,include4X=T) # dat is used to calculate median date and in split with trajet.

#check effect of 4X
#spawn.t<- full_join(read.delim("results/spawning/table_nlmeAR_2022include4X.txt") %>%  mutate(region="4T-4X"),
#                   read.delim("results/spawning/table_nlmeAR_2022.txt") %>%  mutate(region="4T-4W"))
#ggplot(spawn.t, aes(x=year, y=prob_at_med, col=region)) + geom_point()
#extreme values, not kept

########comparaison f.g#############
library(readxl)
spawn.t<- read.delim(paste0("results/",year_of_assess,"/spawning/table_nlmeAR_2022.txt"))

compfg<- full_join(read_excel("data/FG2013_Table8.xlsx") %>%  dplyr::select(-date.med) %>%  mutate(trajet=1),
          spawn.t %>%  filter(trajet==1))
my.formula <- y ~ x # for smooth in ggplot

fg2<- ggplot(compfg , aes(x=prop_spawn_FG*100, y=prob_at_med*100))+geom_point() +theme_few() + geom_smooth(method="lm") +
  geom_abline(col="red", slope=1, intercept=0, lty=2)+
  stat_poly_eq(
  formula = my.formula,
  aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
  parse = TRUE, vstep = 0.08 +geom_abline(col="red", lty=2, slope=1, intercept=0)
) + xlab("Grégoire et al. 2013") 
 fg2FR<- fg2 + ylab("Nouvelle m\u00E9thode") 
 fg2EN<- fg2 + ylab("New method")
 fg2BI<- fg2 + ylab("Nouvelle m\u00E9thode / New method")
 
datfg1<- full_join(compfg %>%  pivot_longer(c("prop_spawn_FG", "prob_at_med")) %>% 
       mutate(name =recode_factor(name, prop_spawn_FG= "Grégoire et al. 2013", prob_at_med="Nouvelle méthode")), 
       expand_grid(all_years,name=c("Grégoire et al. 2013", "Nouvelle méthode"))) %>%  
       mutate(nameEN= recode_factor(name, `Grégoire et al. 2013`="Grégoire et al. 2013", `Nouvelle méthode`="New method"),
              nameBI= recode_factor(name, `Grégoire et al. 2013`="Grégoire et al. 2013", `Nouvelle méthode`="Nouvelle méthode / New method"),
              yearjit=ifelse(name=="Grégoire et al. 2013", year, year+0.2),
              propIClwr = ifelse(name=="Grégoire et al. 2013", NA, propIClwr),
              propICupr = ifelse(name=="Grégoire et al. 2013", NA, propICupr))
        

      fg1FR<-  ggplot(datfg1,aes(x=yearjit, y=value*100, col=name))+geom_point() +geom_errorbar(aes(ymin=propIClwr*100, ymax=propICupr*100), width = 0.8)+
       geom_line() + theme_few()+theme(legend.title=element_blank(),legend.position=c(0.3, 0.85))+
        scale_color_manual(values=c("grey55", "black"))+
        ylab("Pourcentage quotidien d'oeufs pondus") +
       scale_x_continuous(breaks=seq(1980, year_to_report, 5), name="Ann\u00E9e") 
      fg1EN <-  ggplot(datfg1,aes(x=yearjit, y=value*100, col=nameEN))+geom_point()+ geom_errorbar(aes(ymin=propIClwr*100, ymax=propICupr*100), width = 0.8)+
        geom_line() + theme_few()+theme(legend.title=element_blank(),legend.position=c(0.3, 0.85))+
        scale_color_manual(values=c("grey55", "black"))+
        ylab("Percentage of eggs spawned daily") +
        scale_x_continuous(breaks=seq(1980, year_to_report, 5), name="Year") 
      fg1BI<-  ggplot(datfg1,aes(x=yearjit, y=value*100, col=nameBI))+geom_point() + geom_errorbar(aes(ymin=propIClwr*100, ymax=propICupr*100), width = 0.8)+
        geom_line() + theme_few()+theme(legend.title=element_blank(),legend.position=c(0.3, 0.85))+
        scale_color_manual(values=c("grey55", "black"))+
        ylab("Pourcentage quotidien d'oeufs pondus\nPercentage of eggs spawned daily") +
        scale_x_continuous(breaks=seq(1980, year_to_report, 5), name="Année\nYear") 
      
    
ggarrange(fg1FR,fg2FR)
ggsave(paste0("figures/",year_of_assess,"/spawning/ComparaisonFG_FR.png"), width=11, height=4, dpi=600, units="in")
ggarrange(fg1EN,fg2EN)
ggsave(paste0("figures/",year_of_assess,"/spawning/ComparaisonFG_EN.png"), width=11, height=4, dpi=600, units="in")
ggarrange(fg1BI,fg2BI)
ggsave(paste0("figures/",year_of_assess,"/spawning/ComparaisonFG_BI.png"), width=11, height=4, dpi=600, units="in")


area <- sum(distinct(egg %>% dplyr::select(stratum_area, stratum))$stratum_area)

# pour toutes les dates ######
load(paste0("results/",year_of_assess,"/spawning/predicted_propspawn_doy",year_to_report,".RData"))
spawn_st <- full_join(DEPt1a, nlme.fixpredsAR1) %>%
  mutate(DEP.spawn = DEP / prob) %>%
  group_by(year) %>%
  summarize(DEP = mean(DEP.spawn, na.rm = T)
            ) %>%
  mutate(TEP = DEP * area,
         )

#spawn_sttrajet <- full_join(DEPt2a %>%  filter(year %in% trajet2filter[which(trajet2filter$trajet==2),"year"]), nlme.fixpredsAR1) %>%
#  mutate(DEP.spawn = DEP / prob) %>%
#  group_by(year) %>%
#  summarize(DEP.m = mean(DEP.spawn, na.rm = T)) %>%
#  mutate(TEP = DEP.m * area)
#did not do anything more


# pour date med
spaw.med <- read.delim(paste0("results/",year_of_assess,"/spawning/table_nlmeAR_", year_to_report, ".txt")) %>%  mutate(propvar= propsd^2)
spawn.out <- list(DEPt1, DEPt2, DEPM1, DEPto,DEPtod,spaw.med) %>%
  reduce(full_join) %>%
  ungroup() %>%
  mutate(
    TEP = (DEP.p * area) / prob_at_med,
    TEPt = (DEP.t *area)/ prob_at_med) %>%
  filter(!is.na(test_sens))
#covariance
check_covar<- spawn.out %>%  filter(!is.na(prob_at_med)) 
covxy<- abs(cov(check_covar$DEP.p, check_covar$prob_at_med))
covx2y2<- abs(cov(check_covar$DEP.p^2, check_covar$prob_at_med^2))
covxyt<- abs(cov(check_covar$DEP.t, check_covar$prob_at_med))
covx2y2t<- abs(cov(check_covar$DEP.t^2, check_covar$prob_at_med^2))

spawn_error <-  spawn.out #%>%  mutate(#TEPvar= area^2* (covxy + (DEP.p ^2 / prob_at_med^2) -  (DEP.p/ prob_at_med)^2),
                                #      TEPvar=area^2 * (covx2y2 + (varp + DEP.p^2) * (propvar +prob_at_med ^2) -(covxy + DEP.p *prob_at_med)^2),
                                  #    TEPtvar=area^2 * (covx2y2t + (vart + DEP.t^2) * (propvar +prob_at_med ^2) -(covxyt + DEP.t *prob_at_med)^2))
                                      #test_indep_TEPvar= area^2 * (DEP.p^2 * propvar) + (prob_at_med^2 * varp) +( varp*propvar)) wiki page on indepednent.
                                      
#TEPvar= area^2*((DEP.p^2 / propvar) + (prob_at_med^2 / varp) + 2*(varp / propvar)),#from<http://falkenblog.blogspot.com/2008/07/formula-for-varxy.html
                                      #TEPtvar= area^2*((DEP.t^2 / propvar) + (prob_at_med^2 / vart) + (varp / propvar)))#from<http://falkenblog.blogspot.com/2008/07/formula-for-varxy.html
                              #TEPvar = area^2 *((DEP.p / prob_at_med)^2 * ((varp/(DEP.p^2)) + (propvar/(prob_at_med^2)) )),
                              #TEPtvar = area^2 *((DEP.t / prob_at_med)^2 * ((vart/(DEP.t^2)) + (propvar/(prob_at_med^2)))),) # gregoir et al. part page 270
   ## J'ai enlevé la covariance parce que valeurs indépendantes.                                    )


spawn_error<- full_join(spawn_error,
          spawn_error %>%  dplyr::select(year, trajet, test_sens,TEPt) %>%
            filter(test_sens=="No sensitivity") %>%  rename(TEP=TEPt) %>%  mutate(test_sens="Predicted grid"))


spawn.out2 <- full_join(spawn_error, left_join(spaw.med,spawn_st %>% mutate(test_sens = "Spawning at date")))
                        
spawn.out2 <- right_join(spawn.out2, expand_grid(all_years, test_sens=unique(spawn.out2$test_sens)))  %>%
  mutate(
    test_sensEN = recode_factor(test_sens, `Outliers volume estimated`="Volume filtered",`Outliers depth estimated`="Maximum depth sampled",`Egg development` = "Egg development",  `Predicted grid`="Predicted grid",`Survey timing` = "Timing of the survey",`Spawning at date` = "Station Sy,d=x",`Smith et al. 2022`="Smith et al. 2022", `No sensitivity`="No sensitivity"),
    test_sensFR = recode_factor(test_sens, `Outliers volume estimated` = "Volume filtré", `Outliers depth estimated` = "Profondeur maximum échantillonnée",`Egg development` = "Développement des oeufs",  `Predicted grid`="Grille de prédiction", `Survey timing` = "Date du relevé", `Spawning at date` = "Sy,d=x à la station", `Smith et al. 2022`="Smith et al. 2022", `No sensitivity`="Sans sensibilité"),
    yearjit= year+(as.numeric(test_sensEN)*0.05))#,
#  TEPlwr= TEP - 2 * (sqrt(TEPvar)),
 # TEPupr= TEP + 2 * (sqrt(TEPvar)))


##export for assessment
TEP_output<- full_join(all_years, spawn.out2 %>%  filter(test_sens=="No sensitivity", trajet==1) %>%  dplyr::select(year, DEP.p, prob_at_med, TEP))

write.table(TEP_output, file=paste0("results/",year_of_assess,"/TEP_", year_to_report,".txt"), row.names=F, sep="\t", dec=".")



gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


ptEN1 <- ggplot(data = spawn.out2, aes(x = yearjit, y = TEP, col = test_sensEN)) +
  geom_point() +
  #  geom_errorbar(aes(ymin=TEPlwr,ymax=TEPupr))+
  geom_line() +
  theme_few() +
  theme(legend.position = c(0.8, 0.7), legend.title = element_blank()) +
  scale_color_manual(values = c(gg_color_hue(7), "black"))+#, labels = ~par2se(text = .x)) +
  scale_x_continuous(name = "Year", breaks = seq(1980, max(spawn.out2$year), 5)) +
  ylab("Total egg production")
ggsave("figures/",year_of_assess,"/Total_egg_production_depth_outliers_EN.png", width = 8, height = 5, unit = "in", dpi = 600)

ptFR1 <-ggplot(data = spawn.out2, aes(x = year, y = TEP, col = test_sensFR)) +
  geom_point() +
  geom_line() +
  theme_few() +
  theme(legend.position = c(0.8, 0.7), legend.title = element_blank()) +
  scale_color_manual(values = c(gg_color_hue(7), "black")) +
  scale_x_continuous(name = "Ann\u00E9e", breaks = seq(1980, max(spawn.out2$year), 5)) +
  ylab("Production total d'oeufs")
ggsave("figures/",year_of_assess,"/Total_egg_production_depth_outliers_FR.png", width = 8, height = 5, unit = "in", dpi = 600)


#valeurs relatives sensibitlité
#X on axis 
nodatEN<- spawn.out2 %>% group_by(year, test_sens) %>%  mutate(sampled=ifelse(is.na(TEP), "N","Y" )) %>%  filter(sampled=="N", test_sens!="No sensitivity") %>%  
   mutate(name = recode_factor(test_sensEN, `Volume filtered`="Volume filtered",`Maximum depth sampled`="Maximum depth sampled",`Egg development` = "Egg development",  `Predicted grid`="Predicted grid",`Timing of the survey` = "Timing of the survey",`Station Sy,d=x` = "Station Sy,d=x")) %>% 
ungroup() %>% dplyr::select(year, name, sampled)
  
nodatFR<- spawn.out2 %>% group_by(year, test_sens) %>%  mutate(sampled=ifelse(is.na(TEP), "N","Y" )) %>%  filter(sampled=="N", test_sens!="No sensitivity") %>%  
  mutate(name = recode_factor(test_sensFR, `Volume filtré` = "Volume filtré", `Profondeur maximum échantillonnée` = "Profondeur maximum échantillonnée",`Développement des oeufs` = "Développement des oeufs",  `Grille de prédiction`="Grille de prédiction", `Date du relevé` = "Date du relevé", `Sy,d=x à la station` = "Sy,d=x à la station")) %>%
  ungroup() %>% dplyr::select(year, name, sampled)


rel.effect<- spawn.out2 %>%
  dplyr::select(year,trajet, TEP, test_sensEN) %>% filter(!(test_sensEN=="Timing of the survey" & trajet==1), !(test_sensEN!="Timing of the survey" & trajet!=1)) %>%  dplyr::select(-trajet) %>% 
  pivot_wider(names_from = test_sensEN, values_from = TEP) %>%
  pivot_longer(`Timing of the survey`:`Station Sy,d=x`) %>%
  mutate(name = recode_factor(name, `Volume filtered`="Volume filtered",`Maximum depth sampled`="Maximum depth sampled",`Egg development` = "Egg development",  `Predicted grid`="Predicted grid",`Timing of the survey` = "Timing of the survey",`Station Sy,d=x` = "Station Sy,d=x", `Smith et al. 2022`="Smith et al. 2022", `No sensitivity`="No sensitivity")) %>% 
  mutate(diff.rel = (value - `No sensitivity`) / `No sensitivity` * 100)

write.table(rel.effect, "results/",year_of_assess,"/rel_effect_sensitivity_runs.txt", row.names=F, sep="\t", dec=".")


####with depth outliers
dEN <-  ggplot(data=rel.effect , aes(x = year, y = diff.rel)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  theme_few() +
  scale_x_continuous(name = "Year", breaks = seq(1980, max(spawn.out2$year), 5)) +
  ylab("Relative change in index (%)")+
  geom_text(data=nodatEN,aes(x=year, y=0, label="X" ), size=3)+
  facet_wrap(~name, ncol = 2)

ggsave("figures/",year_of_assess,"/Sensitivity_contribution_depthEN.png", width=8, height=8, dpi=600, units="in")


dFR <- spawn.out2 %>%
  dplyr::select(year,trajet, TEP, test_sensFR) %>% filter(!(test_sensFR=="Date du relevé" & trajet==1), !(test_sensFR!="Date du relevé" & trajet!=1)) %>%  dplyr::select(-trajet) %>% 
  pivot_wider(names_from = test_sensFR, values_from = TEP) %>%
  pivot_longer(`Date du relevé`:`Sy,d=x à la station`) %>%
  mutate(name = recode_factor(name, `Volume filtré` = "Volume filtré", `Profondeur maximum échantillonnée` = "Profondeur maximum échantillonnée",`Développement des oeufs` = "Développement des oeufs",  `Grille de prédiction`="Grille de prédiction", `Date du relevé` = "Date du relevé", `Sy,d=x à la station` = "Sy,d=x à la station",`Sans sensibilité`="Sans sensibilité")) %>%
  mutate(diff.rel = (value - `Sans sensibilité`) / `Sans sensibilité` * 100) %>%
  ggplot(aes(x = year, y = diff.rel)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  theme_few() +
  scale_x_continuous(name = "Ann\u00e9e", breaks = seq(1980, max(spawn.out2$year), 5)) +
  ylab("Changement relatif de l'indice (%)")+
  geom_text(data=nodatFR %>%   droplevels(),aes(x=year, y=0, label="X" ), size=3)+
  facet_wrap(~name, ncol = 2)


ggsave("figures/",year_of_assess,"/Sensitivity_contribution_depthFR.png", width = 8, height = 8, dpi = 600, units = "in")



library(CCAM)
year <- 2020
repo <- "https://github.com/iml-mackerel/0.0_model/blob/master/"
smith<- as.data.frame(CCAM::read.ices(url(paste0(repo,'data/',year,'/tep.dat',"?raw=true")))$TEP )%>%  rename(TEP=`-1`) 
smith$year <-  as.numeric(row.names(smith))
smith <-  smith %>%  mutate(trajet=1, test_sens= "Smith et al. 2022", test_sensEN=test_sens, test_sensFR=test_sens,
                            TEP= TEP* (10^12))


   teps<- full_join(smith %>%  mutate(yearjit=year+(7*0.05)), spawn.out2) %>% mutate(test_sensEN = recode_factor(test_sensEN,`Volume filtered`="Volume filtered",`Maximum depth sampled`="Maximum depth sampled",`Egg development` = "Egg development",  `Predicted grid`="Predicted grid",`Timing of the survey` = "Timing of the survey",`Station Sy,d=x` = "Station Sy,d=x", `Smith et al. 2022`="Smith et al. 2022", `No sensitivity` ="Baseline (No sensitivity)"),
                                           test_sensFR = recode_factor(test_sensFR, `Volume filtré` = "Volume filtré", `Profondeur maximum échantillonnée` = "Profondeur maximum échantillonnée",`Développement des oeufs` = "Développement des oeufs",  `Grille de prédiction`="Grille de prédiction", `Date du relevé` = "Date du relevé", `Sy,d=x à la station` = "Sy,d=x à la station", `Smith et al. 2022`="Smith et al. 2022", `Sans sensibilité`="Base de référence (Sans sensibilité)"),
                                           test_sensBI= recode_factor(test_sensFR , `Volume filtré` = "Volume filtré / Volumed filtered", `Profondeur maximum échantillonnée` = "Profondeur maximum échantillonnée / Maximum depth sampled",`Développement des oeufs` = "Développement des oeufs / Egg development",  `Grille de prédiction`="Grille de prédiction / Predicted grid", `Date du relevé` = "Date du relevé / Timing of the survey", `Sy,d=x à la station` = "Station Sy,d=x"))
                                                                  
  ptFR <- ggplot(data=teps ,aes(x = yearjit, y = TEP, col = test_sensFR)) +
  geom_point() +  geom_line() +  theme_few() +
  theme(legend.position = c(0.8, 0.7), legend.title = element_blank()) +
  scale_color_manual(values = c(gg_color_hue(6),"grey55" ,"black")) +
  scale_x_continuous(name = "Ann\u00E9e", breaks = seq(1980, max(spawn.out2$year), 5)) +
  ylab("Production total d'oeufs")
ggsave(paste0("figures/",year_of_assess,"/Total_egg_production_Smith",year_to_report,"FR.png"), width = 8, height = 5, unit = "in", dpi = 600)

ptEN <- ggplot(data=teps ,aes(x = yearjit, y = TEP, col = test_sensEN)) +
  geom_point() +  geom_line() +  theme_few() +
  theme(legend.position = c(0.8, 0.7), legend.title = element_blank()) +
  scale_color_manual(values = c(gg_color_hue(6),"grey55" ,"black")) +
  scale_x_continuous(name = "Year", breaks = seq(1980, max(spawn.out2$year), 5)) +
  ylab("Total egg production")
ggsave(paste0("figures/",year_of_assess,"/Total_egg_production_Smith",year_to_report,"EN.png"), width = 8, height = 5, unit = "in", dpi = 600)

ptBI <- ggplot(data=teps,aes(x = yearjit, y = TEP, col = test_sensBI)) +
  geom_point() +  geom_line() +  theme_few() +
  theme(legend.position = c(0.68, 0.7), legend.title = element_blank()) +
  scale_color_manual(values = c(gg_color_hue(6),"grey55","black")) +
  scale_x_continuous(name = "Ann\u00E9e / Year", breaks = seq(1980, max(spawn.out2$year), 5)) +
  ylab("Production total d'oeufs\nTotal egg production")
ggsave(paste0("figures/",year_of_assess,"/Total_egg_production_Smith",year_to_report,"BI.png"), width = 8, height = 5, unit = "in", dpi = 600)

###no sensitivitity +smith for presentation
ptFR <- ggplot(data=teps %>% filter(test_sens %in% c("Smith et al. 2022", "No sensitivity")) ,aes(x = year, y = TEP, col = test_sensFR)) +
  geom_point() +  geom_line() +  theme_few() +
  theme(legend.position = c(0.8, 0.7), legend.title = element_blank()) +
  scale_color_manual(values = c("grey55", "black")) +
  scale_x_continuous(name = "Ann\u00E9e", breaks = seq(1980, max(spawn.out2$year), 5)) +
  ylab("Production total d'oeufs")
ggsave(paste0("figures/",year_of_assess,"/Total_egg_production_Smith",year_to_report,"2FR.png"), width = 8, height = 5, unit = "in", dpi = 600)

ptEN <- ggplot(data=teps %>% filter(test_sens %in% c("Smith et al. 2022", "No sensitivity")),aes(x = year, y = TEP, col = test_sensEN)) +
  geom_point() +  geom_line() +  theme_few() +
  theme(legend.position = c(0.8, 0.7), legend.title = element_blank()) +
  scale_color_manual(values = c("grey55", "black")) +
  scale_x_continuous(name = "Year", breaks = seq(1980, max(spawn.out2$year), 5)) +
  ylab("Total egg production")
ggsave(paste0("figures/",year_of_assess,"/Total_egg_production_Smith",year_to_report,"2EN.png"), width = 8, height = 5, unit = "in", dpi = 600)

ptBI<- ggplot(data=teps %>% filter(test_sens %in% c("Smith et al. 2022", "No sensitivity")),aes(x = year, y = TEP, col = test_sensBI)) +
  geom_point() +  geom_line() +  theme_few() +
  theme(legend.position = c(0.68, 0.7), legend.title = element_blank()) +
  scale_color_manual(values = c("grey55", "black")) +
  scale_x_continuous(name = "Ann\u00E9e / Year", breaks = seq(1980, max(spawn.out2$year), 5)) +
  ylab("Production total d'oeufs\nTotal egg production")
ggsave(paste0("figures/",year_of_assess,"/Total_egg_production_Smith",year_to_report,"2BI.png"), width = 8, height = 5, unit = "in", dpi = 600)




#####check for spawing at station
#factor levels not changed
check_par<- left_join(spaw.med %>%  filter(trajet==1),spawn.out2 %>%
                        dplyr::select(year,trajet, TEP, test_sensFR) %>% filter(!(test_sensFR=="Date du relevé" & trajet==1), !(test_sensFR!="Date du relevé" & trajet!=1)) %>%  dplyr::select(-trajet) %>% 
                        pivot_wider(names_from = test_sensFR, values_from = TEP) %>%
                        pivot_longer(`Développement des oeufs`:`Sy,d=x à la station`) %>%
                        mutate(name = recode_factor(name, `Volumes aberrants estimés` = "Volumes aberrants estimés", `Développement des oeufs` = "Développement des oeufs", `Sy,d=x à la station` = "Sy,d=x à la station", `Date du relevé` = "Date du relevé")) %>%
                        mutate(diff.rel = (value - `Sans sensibilité`) / `Sans sensibilité` * 100)) %>%  filter(name=="Sy,d=x à la station") %>%  filter(!is.na(value), trajet==1) %>% 
  mutate(doydiff = date.med - xmid)

ggplot(check_par, aes(x=doydiff, y=diff.rel))+geom_point()+geom_label(aes(label=year))
ggplot(check_par, aes(x=duration, y=diff.rel))+geom_point()+geom_label(aes(label=year))

ggplot(check_par, aes(x=prob_at_med, y=diff.rel))+geom_point()+geom_label(aes(label=year))

ggplot(check_par, aes(x=prob_at_med, y=doydiff))+geom_point()+geom_label(aes(label=year, fill=diff.rel))+scale_fill_viridis_c(option="turbo",begin = 0.1, end = 0.9)

change<- ggplot(check_par, aes(x=duration, y=doydiff))+ geom_hline(yintercept=0, lty=2, col="red")+geom_label(aes(label=year, fill=diff.rel))+scale_fill_gradient2(name="",guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) + theme_few()
  
changeEN <- change +  xlab("Sapwning duration (days)") +  ylab("Survey median date - Spawning peak date") 
ggsave("figures/",year_of_assess,"/spawning/explain_changes_with_spawning_at_stationEN.png", dpi=600, units="in", width=6, height=5)

  changeFR <- change +  xlab("Dur\u00E9e du pic de ponte (jours)") +  ylab("Date m\u00E9diane du relev\u00E9 - Date du pic de ponte") 
  ggsave("figures/",year_of_assess,"/spawning/explain_changes_with_spawning_at_stationFR.png", dpi=600, units="in", width=6, height=5)

  changeBI <- change +  xlab("Dur\u00E9e du pic de ponte (jours)\n Sapwning duration (days)") +  ylab("Date m\u00E9diane du relev\u00E9 - Date du pic de ponte\nSurvey median date - Spawning peak date") 
  ggsave("figures/",year_of_assess,"/spawning/explain_changes_with_spawning_at_stationBI.png", dpi=600, units="in", width=6, height=5)
  
  

  #############VALIDATON BY PROPORTION OF STAGE 6
  
  PESd<- read.delim("results/",year_of_assess,"/spawning/table_nlmeAR_2022.txt")  %>%  filter(trajet==1) %>% dplyr::select(year, xmid)
  
  propstages<- read.delim("results/",year_of_assess,"/spawning/proportion_stages_doy.txt")
  
  
  valid_align<- left_join(PESd, propstages)
  
  ggplot(valid_align, aes(x=doy, y= prop_stage6)) +geom_bar(stat="identity") +geom_vline(aes(xintercept=xmid), col="red", lty=2)+facet_wrap(~year)
  
  propstages<- read.delim("results/",year_of_assess,"/spawning/proportion_stages_doy.txt") %>%  dplyr::select(-n, -prop_stage6) %>%  pivot_longer(X5:X8, names_to="stages", values_to="n") %>% 
     mutate(stages=recode_factor(stages, X5 = "5",X6 = "6",X7 = "7",X8 = "8"))
  
  valid_align<- left_join(PESd, propstages)
  
  ggplot(valid_align, aes(x=doy, y= n)) +geom_bar(stat="identity", aes(fill=stages, col=stages)) +geom_vline(aes(xintercept=xmid), col="red", lty=2)+facet_wrap(~year, scale="free_y") +theme_few()+scale_fill_viridis_d(direction=-1) +scale_color_viridis_d(direction=-1)
  
  ggsave("figures/",year_of_assess,"/spawning/validation_stages_prop.png", width=12, height=10, dpi=600, units="in")
      
  
  
  #2. proportion of stages weighted for all n for a year

newprop <- left_join(propstages %>%  dplyr::select(-n),
  propstages %>%  group_by(year) %>%  tally(n)) %>%  
  mutate(prop_stage6 = X6/n) %>% ungroup() %>% dplyr::filter(prop_stage6>0.05)


newprop <- propstages %>% ungroup() %>% dplyr::filter(prop_stage6>0.05)


list(newprop %>% group_by(year) %>%  slice(which.max(prop_stage6))  %>%  dplyr::select(year,doy) %>%  dplyr::rename(maxstage6=doy),
     newprop %>%  group_by(year) %>%  summarize(start6=min(doy, na.rm=T), end6=max(doy, na.rm=T)),
     PESd %>%  group_by(year)) %>% purrr::reduce(left_join) %>% 
    ggplot(aes(x=xmid, y=maxstage6)) +geom_label(aes(label=year))
    
list(newprop %>% group_by(year) %>%  slice(which.max(prop_stage6))  %>%  dplyr::select(year,doy) %>%  dplyr::rename(maxstage6=doy),
     newprop %>%  group_by(year) %>%  summarize(start6=min(doy, na.rm=T), end6=max(doy, na.rm=T)),
     PESd %>%  group_by(year)) %>% purrr::reduce(left_join) %>% 
  ggplot(aes(x=year, y=maxstage6)) + geom_errorbar(aes(ymin=start6, ymax=end6)) +geom_point(aes(y=xmid, x=year))

#weighted average of date, weight = proportion or N?
wa<- list(newprop %>% group_by(year) %>%  summarize(maxstage6 =weighted.mean(doy, w=prop_stage6))  ,
     newprop %>%  group_by(year) %>%  summarize(start6=min(doy, na.rm=T), end6=max(doy, na.rm=T)),
     PESd %>%  group_by(year)) %>% purrr::reduce(left_join) %>% 
  mutate(for_fill= ifelse(year %in% c(1999, 1991), "fit", 
                      ifelse(year %in% c(2006, 2017,2019), "timing", 
                          ifelse(year %in% c(2022), "uncertain", "OK"))))

sd_diff<- wa %>%  filter(year!=2022)%>%  reframe(meandiff_xmid6= sd(xmid-maxstage6, na.rm=T) )
 bias<-  ggplot(data=wa,aes(x=xmid, y=maxstage6)) +geom_abline(intercept=0,slope=1, col="red", lty=2)+geom_abline(intercept=-1.96*sd_diff$meandiff_xmid6,slope=1, col="black", lty=2)+geom_abline(intercept=1.96*sd_diff$meandiff_xmid6,slope=1, col="black", lty=2)+geom_label(aes(label=year, fill=for_fill), alpha=0.5)+
    scale_fill_manual(values=c("#ff6a6a","white", "dodgerblue", "seagreen3"), guide="none") +theme_few()
  
 
 biasEN<- bias + ylab("Day of year of maximum proportion of stage 6") + xlab("Day of year of peak spawning (xmid)")
 ggsave("figures/",year_of_assess,"/spawning/bias_indicator_spawning_EN.png", width=6, height=5, dpi=600, units="in")
 
 biasFR <- bias + ylab("Jour de l'année de la proportion maximale de stade 6") + xlab("Jour de l'année du pic de ponte (xmid)")
 ggsave("figures/",year_of_assess,"/spawning/bias_indicator_spawning_FR.png", width=6, height=5, dpi=600, units="in")
 
  #weighted by gsi
  #wa<- list(newprop %>% group_by(year) %>%  summarize(maxstage6 =weighted.mean(doy, w=mgsi))  ,
  #          newprop %>%  group_by(year) %>%  summarize(start6=min(doy, na.rm=T), end6=max(doy, na.rm=T)),
  #          PESd %>%  group_by(year)) %>% purrr::reduce(left_join) %>% 
  #  mutate(for_fill= ifelse(year %in% c(1999, 1991), "fit", 
  #                          ifelse(year %in% c(2006, 2017,2019), "timing", 
  #                                 ifelse(year %in% c(2022), "uncertain", "OK"))))
#  
#  sd_diff<- wa %>%  filter(year!=2022)%>%  reframe(meandiff_xmid6= sd(xmid-maxstage6, na.rm=T) )
#  ggplot(data=wa,aes(x=xmid, y=maxstage6)) +geom_abline(intercept=0,slope=1, col="red", lty=2)+geom_abline(intercept=-1.96*sd_diff$meandiff_xmid6,slope=1, col="black", lty=2)+geom_abline(intercept=1.96*sd_diff$meandiff_xmid6,slope=1, col="black", lty=2)+geom_label(aes(label=year, fill=for_fill), alpha=0.5)+
#    scale_fill_manual(values=c("#ff6a6a","white", "dodgerblue", "seagreen3"), guide="none") 
  
  
  #############LARVAE##############
  
  load("results/",year_of_assess,"/spawning/predicted_propspawn_doy2022.RData")
  
  tabAR<- read.delim("results/",year_of_assess,"/spawning/table_nlmeAR_2022.txt") %>%  dplyr::select(year, date.med, trajet)
  
  prop<- left_join(tabAR, nlme.fixpredsAR1) %>%  group_by(year, trajet) %>%  filter(doy < date.med) %>%  summarize(prop=sum(prob))
  
  
  df_prop_larvae<- full_join(prop, egg %>%  mutate(prop_larvae =maq_larvae /(maq_eggs_stage1_5 +maq_larvae)) %>%  
                               group_by(year, trajet) %>%  summarize(mean_prop_larvae=mean(prop_larvae, na.rm=T),
                                                                     sd_prop_larvae = sd(prop_larvae, na.rm=T))) %>% 
    mutate(for_fill= ifelse(year %in% c(1999, 1991), "fit", 
                            ifelse(year %in% c(2006, 2017,2019), "timing", 
                                   ifelse(year %in% c(2022), "uncertain", "OK"))))
  
  pl2<- ggplot(data=df_prop_larvae, aes(x=prop, y=mean_prop_larvae))+geom_label(aes(label=year, col=as.factor(trajet), fill=for_fill),alpha=0.5)  + theme_few() +
    scale_color_manual(values=c("black", "purple")) + 
    scale_fill_manual(values=c("#ff6a6a","white", "dodgerblue", "seagreen3"), guide="none")+
    theme(legend.title=element_blank())
  
  pl1 <- ggplot(df_prop_larvae %>%  mutate(year=ifelse(trajet==2, year+0.5, year)) , aes(x=year, y=mean_prop_larvae, col=as.factor(trajet))) +geom_point()+geom_errorbar(aes(ymin=ifelse(mean_prop_larvae-sd_prop_larvae <0, 0,mean_prop_larvae-sd_prop_larvae),
                                                                                                                                                                             ymax=mean_prop_larvae+sd_prop_larvae)) +theme_few() + theme(legend.title=element_blank())+scale_color_manual(values=c("black", "purple"))
  
  pl1EN<- pl1 + scale_x_continuous(name="Year", breaks=seq(1980, year_to_report, 5)) + ylab("Proportion of larvae")
  pl2EN<- pl2 + scale_x_continuous(name="Proportion of spawning before median date") + ylab("Proportion of larvae")
  
  ggarrange(pl1EN, pl2EN, common.legend=T, legend="bottom", ncol=1)
  ggsave("figures/",year_of_assess,"/spawning/Valdiation_spawning_larvesEN.png", width=9, height=8)
  
  
  pl1FR<- pl1 + scale_x_continuous(name="Année", breaks=seq(1980, year_to_report, 5)) + ylab("Proportion de larves")
  pl2FR<- pl2 + scale_x_continuous(name="Proportion de ponte avant la date m\u00E9diane") + ylab("Proportion de larves")
  
  ggarrange(pl1FR, pl2FR, common.legend = T, legend="bottom", ncol=1)
  ggsave("figures/",year_of_assess,"/spawning/Valdiation_spawning_larvesFR.png", width=9, height=8)
  
  pl1BI<- pl1 + scale_x_continuous(name="Année / Year", breaks=seq(1980, year_to_report, 5)) + ylab("Proportion de larves")
  pl2BI<- pl2 + scale_x_continuous(name="Proportion de ponte avant la date m\u00E9diane \n Proportion of spawning before median date") + ylab("Proportion de larves \nProportion of larvae")
  
  ggarrange(pl1BI, pl2BI, common.legend = T, legend="bottom", align="hv", ncol=1)
  ggsave("figures/",year_of_assess,"/spawning/Valdiation_spawning_larvesBI.png", width=9, height=8)
  
  
  
  dur<- read.delim("results/",year_of_assess,"/spawning/table_nlmeAR_2022.txt") %>%  dplyr::filter(!year %in% c(1980:1982,1991, 1999,2020))
  
  my.formula <- y ~ x # for smooth in ggplot
  
  pdu<-   ggplot(dur, aes(x=year, y=duration))+geom_point() + geom_smooth(method="lm")+
    stat_poly_eq(
      formula = my.formula,
      aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
      parse = TRUE, label.x="right"
    ) +theme_few()
  
  lmdur<- lm(duration ~ year, data=dur)
  par(mfcol=c(2,2));plot(lmdur)
  
  pduEN<- pdu +scale_x_continuous(breaks=seq(1980, year_to_report, 5), name="Year")+ ylab("Spawning duration (days)") 
  ggsave("figures/",year_of_assess,"/spawning/durationEN.png", width=7, height=4, dpi=600, units="in")
  pduFR<- pdu +scale_x_continuous(breaks=seq(1980, year_to_report, 5), name="Ann\u00E9e")+ ylab("Durée de la ponte (jours)") 
  ggsave("figures/",year_of_assess,"/spawning/durationFR.png", width=7, height=4, dpi=600, units="in")
  pduBI<- pdu +scale_x_continuous(breaks=seq(1980, year_to_report, 5), name="Ann\u00E9e / Year")+ ylab("Durée de la ponte (jours)\nSpawning duration (days)") 
  ggsave("figures/",year_of_assess,"/spawning/durationBI.png", width=7, height=4, dpi=600, units="in")
  
  
  ####median date
  date.med <- egg %>% filter(year > 1982) %>% 
    ungroup() %>%
    group_by(year, trajet) %>% mutate(date=as.Date(paste(year, month, day, sep="-"))) %>% 
    dplyr::summarize(
      doy.med = round(median(doy, na.rm = T)),
      mindoy = min(doy, na.rm = T),
      maxdoy = max(doy, na.rm = T),
      date.med =paste0("'",substring(median(date, na.rm = T), 6,10)), #' for excel compatibiltiy
      mindate=paste0("'",substring(min(date), 6,10)),
  maxdate=paste0("'",substring(max(date), 6,10))
    )
  
  write.table(date.med, paste0("results/",year_of_assess,"/Survey_dates.txt"), row.names=F, sep="\t", dec=".")
  
  