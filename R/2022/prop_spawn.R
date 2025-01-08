# Load bio data
library(tidyverse)
library(magrittr)
# devtools::install_github("iml-assess/DFOdata")
library(DFOdata)
library(nlme)
library(ggthemes)
library(nlme)
library(nls.multstart)
library(ggpubr)
source(paste0("R/",year_of_assess,"/prop_spawn_CI.R"))


prop_spawn <- function(year_to_report, year_of_assess, dat, cv=F, include4X=F, N=100) {
  if(!include4X)  nameincl4X= ""
  if(include4X)  nameincl4X= "include4X"
  
  #########DATA##########
  
    bio <- get.bio(species = "maquereau", user = my.env$bio.username, password = my.env$bio.password)

    write.table(bio, paste0("data/Maquereau_bio_data_",Sys.Date(),".txt"), row.names=F, sep="\t", dec=".",na="")
    
    #temporaire
    bio[which(bio$year==2022 & bio$nafo=="4W" & bio$matur ==6 & bio$month==11), "day"] <- 11
    bio[which(bio$year==2022 & bio$nafo=="4W" & bio$matur ==6 & bio$month==11), "month"] <- 7
    
if(include4X){
bio <- bio %>% mutate(doy = as.numeric(format(as.Date(paste(year, month, day, sep = "/"), format = "%Y/%m/%d"), format = "%j"))) %>%
   dplyr::filter(
    !is.na(weight), !is.na(wgonad), !is.na(doy),
    weight > 0, weight > wgonad, matur != 0, (nafo %in% c("4T", "4V", "4W") | (nafo=="4X" & matur <7))
  ) 
}
  
    if(!include4X){
      bio <- bio %>% mutate(doy = as.numeric(format(as.Date(paste(year, month, day, sep = "/"), format = "%Y/%m/%d"), format = "%j"))) %>%
        dplyr::filter(
          !is.na(weight), !is.na(wgonad), !is.na(doy),
          weight > 0, weight > wgonad, matur != 0, nafo %in% c("4T", "4V", "4W")
        ) 
    }
      
bio %<>% mutate(gsi = (wgonad / weight) * 100)

#remove outliers
df <- bio %>%
  dplyr::filter(!is.na(gsi), matur > 4, year > 1978, !(gsi > 5 & doy > 250), !(gsi > 10 & doy > 200 & year %in% c(1994, 1999, 2018, 2019, 2020)), !(year == 1999 & matur > 6 & nafo == "4V"), !(nafo!="4T" & matur > 6), !(nafo!="4T" & matur>6 & month >7), length.frozen > 0.2, weight > 0.03, gsi < 40)  %>%  
  group_by(agef) %>% filter(weight>quantile(weight,0.0001), weight<quantile(weight,0.9999))%>%
  arrange(year)


df <- df %>%
  arrange(year) %>%
  mutate(year = as.numeric(as.character(year))) %>%  filter(!is.na(gsi), !is.na(doy), !is.na(year)) %>% mutate(maturcat = ifelse(!is.na(matur) & matur == 6, "6",
                                                                                                                                 ifelse(is.na(matur), "other", "other"))) %>%  ungroup()


ggplot(df %>%  group_by(year, nafo, matur, month) %>% tally() , aes(x=matur, y=year, fill=n))+ geom_tile() + facet_grid(nafo~ month)

date.med <- dat %>%
  ungroup() %>%
  group_by(year, trajet) %>%
  dplyr::summarize(
    date.med = round(median(doy, na.rm = T)),
    mindoy = min(doy, na.rm = T),
    maxdoy = max(doy, na.rm = T)
  )
date.med <- date.med %>% mutate(date.med = ifelse(year == 1979, 166, date.med)) # from Ouellet et al. 1987
date.med <- as.data.frame(date.med)




#doy tiles
df %>%  group_by(year, doy) %>%  tally() %>% 
ggplot(aes(x=as.factor(doy), y=year))+geom_tile(aes(fill=n))



table_prov<- df %>%  group_by(year, nafo) %>%  tally() %>%  pivot_wider(names_from=nafo, values_from=n) %>%  mutate(total=sum(`4T`,`4W`,`4V`, na.rm=T),
                                                                                                                    prop4T=round(`4T`/total,3),
                                                                                                                    prop4W=round(`4W`/total,3),
                                                                                                                    prop4V=round(`4V`/total, 3))

 write.table(table_prov, paste0("results/",year_of_assess,"/spawning/table_provenance_ech",nameincl4X,".txt"), row.names=F, sep="\t", dec=".")

if(!include4X){
ns<- table_prov %>% dplyr::select(-total) %>%  pivot_longer(2:7) %>%  mutate(type=ifelse(grepl(name, pattern="prop"), "Proportion", "n"),
                                                                        region=ifelse(grepl(name, pattern="4T"), "4T",
                                                                                      ifelse(grepl(name, pattern="4V"), "4V",
                                                                                             ifelse(grepl(name, pattern="4W"), "4W", NA)))) %>% dplyr::select(-name) %>% 
  pivot_wider(names_from=type, values_from = value) %>% 
ggplot(aes(x=year, y=n))+geom_bar( stat="identity", aes(fill=region)) + theme_few() + scale_fill_viridis_d(name="", direction=-1)+
  theme(legend.position=c(0.85, 0.85)) #+ geom_text(data=table_prov , aes(x=year, y=1, label=total), angle=90)
}

if(include4X){
  ns<- table_prov %>% dplyr::select(-total) %>%  pivot_longer(2:7) %>%  mutate(type=ifelse(grepl(name, pattern="prop"), "Proportion", "n"),
                                                                               region=ifelse(grepl(name, pattern="4T"), "4T",
                                                                                             ifelse(grepl(name, pattern="4V"), "4V",
                                                                                                    ifelse(grepl(name, pattern="4W"), "4W",
                                                                                                           ifelse(grepl(name, pattern="4X"), "4X",NA))))) %>% dplyr::select(-name) %>% 
    pivot_wider(names_from=type, values_from = value) %>% 
    ggplot(aes(x=year, y=n))+geom_bar( stat="identity", aes(fill=region)) + theme_few() + scale_fill_viridis_d(name="", direction=-1)+
    theme(legend.position=c(0.85, 0.85)) #+ geom_text(data=table_prov , aes(x=year, y=1, label=total), angle=90)
}


nsEN <-  ns + ylab("Number of fish") + scale_x_continuous(name="Year", breaks=seq(1980, max(df$year), 5))
ggsave(paste0("figures/",year_of_assess,"/spawning/Nsamples",nameincl4X,"EN.png"), width=6, height=5, dpi=600, units="in")

nsFR <-  ns + ylab("Nombre de poissons") + scale_x_continuous(name="Ann\u00E9e", breaks=seq(1980, max(df$year), 5))
 ggsave(paste0("figures/",year_of_assess,"/spawning/Nsamples",nameincl4X,"FR.png"), width=6, height=5, dpi=600, units="in")


ns<- df %>%  group_by(year, matur) %>%  tally() %>%  pivot_wider(names_from=matur, values_from=n) %>% 
  pivot_longer(`5` :`8`) %>% 
  ggplot(aes(x=year, y=value))+geom_bar( stat="identity", aes(fill=name)) + theme_few() + scale_fill_viridis_d(name="",direction=-1)+
  theme(legend.position=c(0.85, 0.85),
        legend.background = element_blank()) #+ geom_text(data=table_prov , aes(x=year, y=1, label=total), angle=90)

nsEN <-  ns + ylab("Number of fish") + scale_x_continuous(name="Year", breaks=seq(1980, max(df$year), 5))
ggsave(paste0("figures/",year_of_assess,"/spawning/Nsamples_stages",nameincl4X,"EN.png"), width=6, height=5, dpi=600, units="in")

nsFR <-  ns + ylab("Nombre de poissons") + scale_x_continuous(name="Ann\u00E9e", breaks=seq(1980, max(df$year), 5))
ggsave(paste0("figures/",year_of_assess,"/spawning/Nsamples_stages",nameincl4X,"FR.png"), width=6, height=5, dpi=600, units="in")


#export n stages by doy
propdoy<- left_join(left_join(df %>%  group_by(year,doy, matur) %>%  tally() %>%  pivot_wider(names_from=matur, values_from=n, values_fill=0),
df %>%  group_by(year,doy) %>%  tally()) %>%  mutate(prop_stage6= `6` /n),
df %>%  group_by(year, doy) %>%  filter(matur=="6") %>%  summarize(mgsi=mean(gsi, na.rm=T)))

write.table(propdoy, paste0("results/",year_of_assess,"/spawning/proportion_stages_doy",nameincl4X,".txt"), row.names=F, sep="\t", dec=".")



nsy<- df %>%  group_by(doy, year, matur) %>%  tally() %>% 
  ggplot(aes(x=doy, y=n))+facet_wrap(~year) +geom_bar(linewidth=0.1,stat="identity", aes(fill=as.factor(matur),col=as.factor(matur)))+ theme_few() + scale_fill_viridis_d(name="",direction=-1)+ scale_color_viridis_d(name="",direction=-1)+
  theme(legend.position=c(0.4, 0.05),legend.direction="horizontal",
        legend.background = element_blank()) #+ geom_text(data=table_prov , aes(x=year, y=1, label=total), angle=90)


nsyEN <-  nsy + ylab("Number of fish") + xlab("Day of year")
ggsave(paste0("figures/",year_of_assess,"/spawning/Nsamples_stages",nameincl4X,"_yearEN.png"), width=10, height=6, dpi=600, units="in")

nsyFR <-  nsy + ylab("Nombre de poissons") + xlab("Jour de l'ann\u00E9e")
ggsave(paste0("figures/",year_of_assess,"/spawning/Nsamples_stages",nameincl4X,"_yearFR.png"), width=10, height=6, dpi=600, units="in")


mdoy6<- df %>%  group_by(year) %>% filter(matur==6) %>%  summarize(mdoy_stage6=mean(doy),
                                                                      sddoy_stage6=sd(doy)) 


############MODEL############
#evaluate starting values
start.val <- nls.multstart::nls_multstart(gsi ~ SSlogis(doy, Asym, xmid, scal),
                                          data = df,
                                          start_lower = c(Asym = min(df %>% filter(matur == 5) %>% dplyr::select(gsi)), xmid = min(df$doy), scal = -100),
                                          start_upper = c(Asym = max(df$gsi), xmid = max(df$doy), scal = 100),
                                          iter = 500,
                                          supp_errors = "Y"
)


nlme.fixAR <- nlme(gsi ~ SSlogis(doy, Asym, xmid, scal),
  data = df,
  fixed = Asym + xmid + scal ~ 1,
  start = coefficients(start.val),
  correlation = corAR1(),
  random = Asym + xmid + scal ~ 1 | year
)

nlme_coefAR = as_tibble(coef(nlme.fixAR), rownames = 'year')

#keep starting values in case of needed for update
save(nlme.fixAR,start.val, file = paste0("results/", year_of_assess,"/spawning/nlmeAR", year_to_report, nameincl4X,".RData"))

predtimes <- expand_grid(year = unique(df$year), doy = seq(0.5, 364.5)) # doy aligned to calculate probability at doy

################CROSS VALIDATION##############
if(cv==T){
predcv <-  data.frame()
#for(y in 1984:2019){
  
  for(y in 1984:2019){
  
    mindoy= df %>%  filter(year==2022) %>%  summarize(mindoy=min(doy, na.rm=T))
    
    nlmecv <- tryCatch(
      {nlmecv <- nlme(gsi ~ SSlogis(doy, Asym, xmid, scal),
                    data = df %>%  filter(!(year==y & doy  < mindoy$mindoy)),
                    fixed = Asym + xmid + scal ~ 1,
                    start = coefficients(start.val),
                    correlation = corAR1(),
                    random = Asym + xmid + scal ~ 1 | year
 )}, error = function(e) {return(NULL)})
    

  if(!is.null(nlmecv)){
  nlmepredcv <- left_join(date.med, predtimes %>%
    dplyr::mutate(fitted = predict(nlmecv, newdata = predtimes)) %>%
    group_by(year) %>%
    mutate(
      slope = c(diff(fitted, lag = 1), NA),
      prob = slope / sum(slope, na.rm = T)
    ) %>%
    filter(!is.na(prob)) %>%
    mutate(doy = ceiling(doy))) %>% filter(year==y) # get back the real doy
  
predcv<-  bind_rows(predcv, nlmepredcv)
  }
}

write.table(predcv, paste0("results/",year_of_assess,"/spawning/prediction_CV",year_to_report,nameincl4X,".txt"), row.names=F, sep="\t", dec=".")

tablecv<- predcv %>% filter(doy==date.med) %>% dplyr::rename(probcv=prob) %>%  dplyr::select(year,probcv, trajet)

}

#########MODEL validation##############
valid<- data.frame(resid=residuals(nlme.fixAR , type="pearson"), mu=fitted(nlme.fixAR), obs=df$gsi, year=df$year)

h<- ggplot(valid, aes(x=mu, y=resid))+geom_point(col="grey35", size=0.5) +theme_few() + geom_hline(yintercept=0, col="red", lty=2)
f<- ggplot(valid, aes(x=df$gsi, y=mu))+geom_point(col="grey35", size=0.5)+theme_few() +geom_abline(slope=1, intercept=0, col="red", lty=2)+ geom_smooth(method="lm")

hEN <-  h+ ylab("Pearson's residuals") + xlab("Fitted values")
hFR <-  h+ ylab("R\u00E9sidus de Pearson") + xlab("Valeurs calibr\u00E9es")

fEN <-  f+ ylab("Fitted values") + xlab("Observed values")
fFR <-  f+ ylab("Valeurs calibr\u00E9es") + xlab("Valeurs observ\u00E9es")

ggarrange(hEN, fEN, ncol=2)
ggsave(paste0("figures/",year_of_assess,"/spawning/nlme_valid", year_to_report,nameincl4X,"EN.png"), width=8, height=4, units="in", dpi=600)

ggarrange(hFR, fFR, ncol=2)
ggsave(paste0("figures/",year_of_assess,"/spawning/nlme_valid", year_to_report,nameincl4X,"FR.png"), width=8, height=4, units="in", dpi=600)




h<- ggplot(valid, aes(x=mu, y=resid))+geom_point(col="grey35", size=0.5) +theme_few() + facet_wrap( ~ year) +geom_hline(yintercept=0, col="red", lty=2)
f<- ggplot(valid, aes(x=df$gsi, y=mu))+geom_point(col="grey35", size=0.5)+theme_few() + facet_wrap( ~ year) +geom_abline(slope=1, intercept=0, col="red", lty=2)+ geom_smooth(method="lm")

hEN <-  h+ ylab("Pearson's residuals") + xlab("Fitted values")
hFR <-  h+ ylab("R\u00E9sidus de Pearson") + xlab("Valeurs calibr\u00E9es")

fEN <-  f+ ylab("Fitted values") + xlab("Observed values")
fFR <-  f+ ylab("Valeurs calibr\u00E9es") + xlab("Valeurs observ\u00E9es")

ggarrange(hEN, fEN, ncol=2)
 ggsave(paste0("figures/",year_of_assess,"/spawning/nlme_validyear", year_to_report,nameincl4X,"EN.png"), width=20, height=15, units="in", dpi=600)

ggarrange(hFR, fFR, ncol=2)
ggsave(paste0("figures/",year_of_assess,"/spawning/nlme_validyear", year_to_report,nameincl4X,"FR.png"), width=20, height=15, units="in", dpi=600)


##################results##############

load(paste0("results/",year_of_assess,"/spawning/nlmeAR", year_to_report,nameincl4X, ".RData"))


nlme.fixpredsAR1 <- predtimes %>%
  mutate(fitted =predict(nlme.fixAR, newdata = predtimes,interval = 'confidence')) %>%
  group_by(year) %>%
  mutate(
    slope = c(diff(fitted, lag = 1), NA),
    prob = slope / sum(slope, na.rm = T)
  ) %>%
  filter(!is.na(prob)) %>%
  mutate(
    probcum0025 = cumsum(prob),
    #probcum9725 = rev(cumsum(rev(prob))),
    doy = ceiling(doy)
  ) # get back the real doy



start<- Sys.time()
yvals2 <- replicate(n=N, nlme_boot(fitted=nlme.fixAR, data=df, newframe=nlme.fixpredsAR1)) # temporary fix for problem in data length. save df inside Rdata to avoid.   # ncluster cannot be larger than 2, because mremory is exceeded. 
#started at 14:44 for 100 iterations
end<- Sys.time()
end-start

yout<-as.data.frame(apply(yvals2,MARGIN=c(1, 2) ,FUN=function(x) unlist(x))) 

gsi_pred<- bind_cols(yout[,seq(1,N*2, 2)],nlme.fixpredsAR1) %>%  pivot_longer(1:N) %>%  group_by(year, doy) 

#figure pour la prsésentation
fp<- ggplot(gsi_pred %>%  filter(year==1999)) + geom_line(aes(x=doy, y=value, group=name)) +geom_line(aes(x=doy, y=fitted),lwd=1.5, col="blue")   + geom_vline(xintercept = 10, lty=2, col="red") +ggtitle("1999")

fpEN<- fp + xlab("Day of year") +ylab("Gonadosomatic index")
fpFR<- fp + xlab("Jour de l'ann\u00E9e") +ylab("Indice gonado-somatique")
fpBI<- fp + xlab("Jour de l'ann\u00E9e\n Day of year") +ylab("Indice gonado-somatique\nGonadosomatic index")

#check when sd stabilize
checkn<- gsi_pred %>%  filter(doy ==10)#, year==1999) # remove days where there is no variation
checkn<- checkn %>%  mutate(name=as.numeric(gsub(name, pattern="V", replacement="")), name=ceiling(name/2))

checko  <- list()
for(i in 2:max(unique(checkn$name))){
 checki <-  checkn %>%  filter(name <= i) %>%  dplyr::group_by(year, doy) %>%  dplyr::summarize(sd=sd(value),
                                                                                                IClwr=quantile(value,0.025), ICupr=quantile(value, 0.975)) %>%  dplyr::mutate(it=i)
 checko<- bind_rows(checko, checki)
}

cp<- ggplot(checko , aes(x=it, y=sd, group=it)) + geom_boxplot() +scale_x_continuous(n.breaks = 10)

#ggplot(checko , aes(x=it, y=IClwr, group=it)) + geom_boxplot() +scale_x_continuous(n.breaks = 10)
#ggplot(checko , aes(x=it, y=ICupr, group=it)) + geom_boxplot() +scale_x_continuous(n.breaks = 10)

cpEN <-  cp + xlab("Number of bootstrap samples") + ylab("SD of predictions at day 10")
cpFR <-  cp + xlab("Nombre de r\u00E9\u00E9chantillonage") + ylab("ET des pr\u00E9dictions au jour 10")
cpBI <-  cp + xlab("Nombre de r\u00E9\u00E9chantillonage\nNumber of bootstrap samples") + ylab("ET des pr\u00E9dictions au jour 10\nSD of predictions at day 10")

ggarrange(fpEN, cpEN, align="hv")
ggsave(paste0("figures/",year_of_assess,"/spawning/SD_bootEN.png"), width=9, height=4, dpi=600, unit="in")
ggarrange(fpFR, cpFR, align="hv")
ggsave(paste0("figures/",year_of_assess,"/spawning/SD_bootFR.png"), width=9, height=4, dpi=600, unit="in")
ggarrange(fpBI, cpBI, align="hv")
ggsave(paste0("figures/",year_of_assess,"/spawning/SD_bootBI.png"), width=9, height=4, dpi=600, unit="in")

#is already grouped bu year and doy
gsi_pred<- gsi_pred %>%  summarize(gsiIClwr=quantile(value,0.025), gsiICupr=quantile(value, 0.975)) %>% mutate(doy = ceiling(doy))
prob_pred<- bind_cols(yout[,seq(2,N*2, 2)],nlme.fixpredsAR1) %>%  pivot_longer(1:N) %>%  filter(!is.na(value)) %>%  group_by(year, doy) %>%  summarize(propIClwr=quantile(value,0.025), propICupr=quantile(value, 0.975),
                                                                                                                                                       propsd=sd(value)) %>% mutate(doy = ceiling(doy))

nlme.fixpredsAR1<- list(nlme.fixpredsAR1, gsi_pred, prob_pred) %>%  purrr::reduce(full_join)

save(nlme.fixpredsAR1, file = paste0("results/",year_of_assess,"/spawning/predicted_propspawn_doy", year_to_report, nameincl4X,".RData"))


###########tables for resdoc#############

tabl1 <- left_join(nlme_coefAR %>% mutate(year = as.numeric(year)), full_join(
  nlme.fixpredsAR1 %>% group_by(year) %>% slice(which.max(prob)) %>% rename(doytop = doy) %>% dplyr::select(year, doytop),
  list(
    nlme.fixpredsAR1 %>% group_by(year) %>% filter(probcum0025 > 0.025 & probcum0025 < 0.975) %>% slice(which.min(doy)) %>% rename(doymin = doy) %>% dplyr::select(year, doymin),
    nlme.fixpredsAR1 %>% group_by(year) %>% filter(probcum0025 > 0.025 & probcum0025 <  0.975) %>% slice(which.max(doy)) %>% rename(doymax = doy) %>% dplyr::select(year, doymax),
    nlme.fixpredsAR1 %>% group_by(year) %>% filter(probcum0025 > 0.15 & probcum0025 <  0.85) %>% slice(which.min(doy)) %>% rename(doymin70 = doy) %>% dplyr::select(year, doymin70),
    nlme.fixpredsAR1 %>% group_by(year) %>% filter(probcum0025 > 0.15 & probcum0025 <  0.85) %>% slice(which.max(doy)) %>% rename(doymax70 = doy) %>% dplyr::select(year, doymax70)
    
  ) %>%  purrr::reduce(full_join), 
) %>%
  mutate(duration = doymax - doymin))



tabl2 <- left_join(date.med, nlme.fixpredsAR1) %>%
  filter(date.med == doy) %>%
  mutate(prob_at_med = prob) %>%
  dplyr::select(year,trajet, date.med, prob_at_med, propIClwr,propICupr, propsd)

tabl_out <- left_join(tabl1, tabl2) %>%  mutate(survey_date= ifelse(date.med < doymin70, "Early survey",
                                                                    ifelse(date.med > doymax70, "Late survey", "During 70 % of spawning")))

write.table(tabl_out, paste0("results/",year_of_assess,"/spawning/table_nlmeAR_", year_to_report,nameincl4X, ".txt"), row.names = F, sep = "\t", dec = ".")


nlme.fixpredsAR <- full_join(df, nlme.fixpredsAR1)


if(cv==T){
  table_medcv<- full_join(read.delim(paste0("results/",year_of_assess,"/spawning/table_nlmeAR_", year_to_report,nameincl4X, ".txt")), #%>% 
                            #dplyr::select(-propIClwr,-propICupr)
                          tablecv) %>%  
  pivot_longer(c("prob_at_med","probcv")) %>%  
  mutate(nameEN=recode_factor(name, prob_at_med="all data", probcv="cropped data"),
         nameFR=recode_factor(name, prob_at_med="toutes les données", probcv="données tronquées"),
         nameBI=recode_factor(name, prob_at_med="toutes les données / all data", probcv="données tronquées / cropped data")) %>%  dplyr::select(-name) %>% 
    filter(trajet==1) %>%  mutate(propIClwr=ifelse(nameEN!="all data", NA,propIClwr),
                                  propICupr=ifelse(nameEN!="all data", NA,propICupr))
  
 
  my.formula <- y ~ x # for smooth in ggplot
  all_years= data.frame(year=seq(1979, year_to_report, 1))
  table_medcv<- right_join(table_medcv, expand_grid(all_years,
                                      data.frame(nameEN=unique(table_medcv$nameEN),
                                      nameFR=unique(table_medcv$nameFR),
                                      nameBI=unique(table_medcv$nameBI)
                                      ))) %>%  arrange(year)
  
  
  v1EN <- ggplot(data = table_medcv, aes(x = year, y = value*100, col = nameEN)) +
    geom_point() + geom_line() +theme_few() +geom_errorbar(aes(ymin=propIClwr*100, ymax=propICupr*100), width = 0.8)+
    theme(legend.title = element_blank(), legend.position = c(0.45, 0.88), legend.background = element_blank())+
    ylab("Percentage of eggs spawned at median survey date") +
    scale_color_manual(values=c("black", "red"))+
    scale_x_continuous(name="Year", breaks=seq(min(table_medcv$year), max(table_medcv$year), 5))
  
  v1FR <- ggplot(data = table_medcv, aes(x = year, y = value*100, col = nameFR)) +
    geom_point() +geom_line() +  theme_few() +geom_errorbar(aes(ymin=propIClwr*100, ymax=propICupr*100), width = 0.8)+
    theme(legend.title = element_blank(), legend.position = c(0.45, 0.88), legend.background = element_blank())+
    ylab("Pourcentage d'oeufs pondus à la date médianne du relevé")+
    scale_color_manual(values=c("black", "red"))+
  scale_x_continuous(name="Ann\u00E9e", breaks=seq(min(table_medcv$year), max(table_medcv$year), 5))
  
  
  v1BI <- ggplot(data = table_medcv, aes(x = year, y = value*100, col = nameBI)) +
    geom_point() +  geom_line() +  theme_few() +geom_errorbar(aes(ymin=propIClwr*100, ymax=propICupr*100), width = 0.8)+
    theme(legend.title = element_blank(), legend.position = c(0.45, 0.88), legend.background = element_blank())+
    ylab("Pourcentage d'oeufs pondus à la date médianne du relevé\nPercentage of eggs spawned at median survey date")+
    scale_color_manual(values=c("black", "red"))+
    scale_x_continuous(name="Ann\u00E9e / Year", breaks=seq(min(table_medcv$year), max(table_medcv$year), 5))
  
  
  v2 <- ggplot(data = table_medcv %>% dplyr::select(-nameFR, -nameBI, -propIClwr, -propICupr) %>% pivot_wider(names_from = nameEN, values_from = value), aes(x = `all data`*100, y = `cropped data`*100)) +
    geom_point() +
    geom_smooth(method = "lm") +
    stat_poly_eq(
      formula = my.formula,
      aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")),
      parse = TRUE, vstep = 0.08
    ) +
    geom_abline(slope = 1, intercept = 0, col = "red", lty = 2) +
    theme_few()
  
   v2EN <- v2 + xlab("Percentage of eggs spawned at median survey date\nall data") + ylab("Percentage of eggs spawned at median survey date\ncropped data")
   v2FR <-  v2+xlab("Pourcentage d'oeufs pondus à la date médianne du relevé\ntoutes les données") + ylab("Pourcentage d'oeufs pondus à la date médianne du relevé\ndonnées tronquées")
   v2BI <-  v2+xlab("Pourcentage d'oeufs pondus à la date médianne du relevé (toutes les données)\nPercentage of eggs spawned at median survey date (all data)") + ylab("Pourcentage d'oeufs pondus à la date médianne du relevé (données tronquées)\nPercentage of eggs spawned at median survey date (cropped data)")
   

  ggarrange(v1EN, v2EN, ncol=2)
ggsave(paste0("figures/",year_of_assess,"/spawning/spawning_CV_",nameincl4X,"EN.png"), width=12, height=5, dpi=600, units="in")

ggarrange(v1FR, v2FR, ncol=2)
ggsave(paste0("figures/",year_of_assess,"/spawning/spawning_CV_",nameincl4X,"FR.png"), width=12, height=5, dpi=600, units="in")

ggarrange(v1BI, v2BI, ncol=2)
ggsave(paste0("figures/",year_of_assess,"/spawning/spawning_CV_",nameincl4X,"BI.png"), width=14, height=7, dpi=600, units="in")


}
##########plots of fit with proportion###########
scaling.factor <- 4 # for second axis

sp <- ggplot(nlme.fixpredsAR) +
  theme_few() + theme(legend.position = "none", axis.text.x = element_text(angle = 90)) + 
  coord_cartesian(xlim=c(100, 300))+
  geom_point(aes(x = doy, y = gsi, col = maturcat), size = 0.3) +
  scale_color_manual(values = c("green3", "grey")) +
  geom_line(aes(x = doy, y = fitted), size = 0.5) +
  geom_ribbon(aes(ymin=gsiIClwr, ymax=gsiICupr, x=doy), alpha=0.5, col="black", fill="black")+
  geom_line(aes(x = doy, y = prob * 100 * scaling.factor), col = "red",size=0.5) +
  geom_ribbon(aes(ymin=propIClwr* 100 * scaling.factor, ymax=propICupr* 100 * scaling.factor, x=doy), alpha=0.5, col = "red", fill="red")+
  geom_vline(data=date.med, aes(xintercept=date.med), lty=2) +
  facet_wrap(~year, ncol = 6) 
 
 
spEN <-  sp +
   scale_y_continuous(
    name = "Gonadosomatic index",
    sec.axis = sec_axis(~ . / scaling.factor, name = "Percentage of eggs spawned daily")
  ) +
   xlab("Day of year") 

ggsave(paste0("figures/",year_of_assess,"/spawning/nlme_spawning_proportion_", year_to_report,nameincl4X, "EN.png"), width = 9, height = 12, dpi = 600, units = "in")


spFR <-  sp + scale_y_continuous(
  name = "Indice gonado-somatique",
  sec.axis = sec_axis(~ . / scaling.factor, name = "Pourcentage quotidien d'oeufs pondus")
) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
  xlab("Jour de l'ann\u00E9e")

ggsave(paste0("figures/",year_of_assess,"/spawning/nlme_spawning_proportion_", year_to_report,nameincl4X, "FR.png"), width = 9, height = 12, dpi = 600, units = "in")


sp1 <- ggplot(nlme.fixpredsAR %>% filter(year < 2003)) +
  theme_few() + theme(legend.position = "none", axis.text.x = element_text(angle = 90)) + 
  coord_cartesian(xlim=c(100, 300))+
  geom_point(aes(x = doy, y = gsi, col = maturcat), size = 0.3) +
  scale_color_manual(values = c("green3", "grey")) +
  geom_line(aes(x = doy, y = fitted), size = 0.5) +
  geom_ribbon(aes(ymin=gsiIClwr, ymax=gsiICupr, x=doy), alpha=0.5, col="black", fill="black")+
  geom_line(aes(x = doy, y = prob * 100 * scaling.factor), col = "red", size = 0.5) +
  geom_ribbon(aes(ymin=propIClwr* 100 * scaling.factor, ymax=propICupr* 100 * scaling.factor, x=doy), alpha=0.5, col = "red", fill="red")+
  geom_vline(data=date.med %>% filter(year < 2003), aes(xintercept=date.med), lty=2) +
  facet_wrap(~year, ncol = 6) 

spEN1<- sp1 +
  scale_y_continuous(
    name = "Gonadosomatic index",
    sec.axis = sec_axis(~ . / scaling.factor, name = "Percentage of eggs spawned daily")
  ) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
  xlab("Day of year")
 ggsave(paste0("figures/",year_of_assess,"/spawning/nlme_spawning_proportion_", year_to_report,nameincl4X, "EN1.png"), width = 9, height = 9, dpi = 600, units = "in")

 
 spFR1<- sp1 + scale_y_continuous(
   name = "Indice gonado-somatique",
   sec.axis = sec_axis(~ . / scaling.factor, name = "Pourcentage quotidien d'oeufs pondus")
 ) +
   theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
   xlab("Jour de l'ann\u00E9e")
 
 ggsave(paste0("figures/",year_of_assess,"/spawning/nlme_spawning_proportion_", year_to_report,nameincl4X, "FR1.png"), width = 9, height =9, dpi = 600, units = "in")
 
 
 spBI1<- sp1 + scale_y_continuous(
   name = "Indice gonado-somatique\nGonadosomatic index",
   sec.axis = sec_axis(~ . / scaling.factor, name = "Percentage of eggs spawned daily\nPourcentage quotidien d'oeufs pondus")
 ) +
   theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
   xlab("Jour de l'ann\u00E9e\nDay of year")
 
 ggsave(paste0("figures/",year_of_assess,"/spawning/nlme_spawning_proportion_", year_to_report,nameincl4X, "BI1.png"), width = 9, height = 9, dpi = 600, units = "in")
 
 
 
 sp2 <- ggplot(nlme.fixpredsAR %>% filter(year >= 2003)) +
   theme_few() + theme(legend.position = "none", axis.text.x = element_text(angle = 90)) + 
   coord_cartesian(xlim=c(100, 300))+
   geom_point(aes(x = doy, y = gsi, col = maturcat), size = 0.3) +
   scale_color_manual(values = c("green3", "grey")) +
   geom_line(aes(x = doy, y = fitted), size = 0.5) +
   geom_ribbon(aes(ymin=gsiIClwr, ymax=gsiICupr, x=doy), alpha=0.5, col="black", fill="black")+
   geom_line(aes(x = doy, y = prob * 100 * scaling.factor), col = "red", size = 0.5) +
   geom_ribbon(aes(ymin=propIClwr* 100 * scaling.factor, ymax=propICupr* 100 * scaling.factor, x=doy), alpha=0.5, col = "red", fill="red")+
   geom_vline(data=date.med %>% filter(year >= 2003), aes(xintercept=date.med), lty=2) +
   facet_wrap(~year, ncol = 6) 
 
 
 spEN2 <-  sp2 + scale_y_continuous(
    name = "Gonadosomatic index",
    sec.axis = sec_axis(~ . / scaling.factor, name = "Percentage of eggs spawned daily")
  ) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
  xlab("Day of year")
 ggsave(paste0("figures/",year_of_assess,"/spawning/nlme_spawning_proportion_", year_to_report, nameincl4X,"EN2.png"), width = 9, height = 9, dpi = 600, units = "in")

 spFR2 <-  sp2 + scale_y_continuous(
   name = "Indice gonado-somatique",
   sec.axis = sec_axis(~ . / scaling.factor, name = "Pourcentage quotidien d'oeufs pondus")
 ) +
   theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
   xlab("Jour de l'ann\u00E9e")
 ggsave(paste0("figures/",year_of_assess,"/spawning/nlme_spawning_proportion_", year_to_report, nameincl4X, "FR2.png"), width = 9, height = 9, dpi = 600, units = "in")
 
 spBI2<- sp2 + scale_y_continuous(
   name = "Indice gonado-somatique\nGonadosomatic index",
   sec.axis = sec_axis(~ . / scaling.factor, name = "Percentage of eggs spawned daily\nPourcentage quotidien d'oeufs pondus")
 ) +
   theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
   xlab("Jour de l'ann\u00E9e\nDay of year")
 
 ggsave(paste0("figures/",year_of_assess,"/spawning/nlme_spawning_proportion_", year_to_report,nameincl4X, "BI2.png"), width = 9, height = 9, dpi = 600, units = "in")
 
 
 
 
 
sp3 <- ggplot(nlme.fixpredsAR %>% filter(year >= year_to_report-5)) +
  theme_few() + theme(legend.position = "none", axis.text.x = element_text(angle = 90)) + 
  coord_cartesian(xlim=c(100, 300))+
  geom_point(aes(x = doy, y = gsi, col = maturcat), size = 0.3) +
  scale_color_manual(values = c("green3", "grey")) +
  geom_line(aes(x = doy, y = fitted), size = 0.5) +
  geom_ribbon(aes(ymin=gsiIClwr, ymax=gsiICupr, x=doy), alpha=0.5, col="black", fill="black")+
  geom_line(aes(x = doy, y = prob * 100 * scaling.factor), col = "red", size = 0.5) +
  geom_ribbon(aes(ymin=propIClwr* 100 * scaling.factor, ymax=propICupr* 100 * scaling.factor, x=doy), alpha=0.5, col = "red", fill="red")+
  geom_vline(data=date.med %>%filter(year >= year_to_report-5), aes(xintercept=date.med), lty=2) +
  facet_wrap(~year, ncol = 3) 


spEN3 <-  sp3 +
  scale_y_continuous(
    name = "Gonadosomatic index",
    sec.axis = sec_axis(~ . / scaling.factor, name = "Percentage of eggs spawned daily")
  ) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
  xlab("Day of year")
ggsave(paste0("figures/",year_of_assess,"/spawning/nlme_spawning_proportion_", year_to_report, nameincl4X,"EN3.png"), width = 8, height = 6, dpi = 600, units = "in")



spFR3 <-  sp3 +
  scale_y_continuous(
    name = "Indice gonado-somatique",
    sec.axis = sec_axis(~ . / scaling.factor, name = "Pourcentage quotidien d'oeufs pondus")
  ) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
  xlab("Jour de l'ann\u00E9e")
ggsave(paste0("figures/",year_of_assess,"/spawning/nlme_spawning_proportion_", year_to_report, nameincl4X, "FR3.png"), width = 8, height = 6, dpi = 600, units = "in")


spBI3<- sp3 + scale_y_continuous(
  name = "Indice gonado-somatique\nGonadosomatic index",
  sec.axis = sec_axis(~ . / scaling.factor, name = "Percentage of eggs spawned daily\nPourcentage quotidien d'oeufs pondus")
) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) +
  xlab("Jour de l'ann\u00E9e\nDay of year")

ggsave(paste0("figures/",year_of_assess,"/spawning/nlme_spawning_proportion_", year_to_report,nameincl4X, "BI3.png"), width = 8, height = 6, dpi = 600, units = "in")



date.medfig<- date.med  %>% mutate(mindoy=ifelse(year==1979, date.med, mindoy), maxdoy=ifelse(year==1979, date.med,maxdoy)) %>%  dplyr::rename(doytop=date.med)
tabl_outfig <-  unique( tabl_out %>%  dplyr::select(year, doymin, doymax,doytop, doymin70, doymax70)) %>%  mutate(trajet=1)


pt2 <-    ggplot(data=tabl_outfig, aes(y = doytop, x = year)) +
       geom_line()+
       geom_ribbon(aes(ymin=doymin70, ymax=doymax70), col="goldenrod2", fill="goldenrod2", alpha=0.3) +
      geom_ribbon(aes(ymin=doymin, ymax=doymax), col="gold", fill="gold", alpha=0.3) +
  geom_point(data=date.medfig ,aes(y=doytop,x=year, col=as.factor(trajet))) +
       geom_errorbar(data=date.medfig,aes(ymin = mindoy, ymax = maxdoy, y=doytop,x=year, col=as.factor(trajet))) +
       
       scale_x_continuous(breaks=seq(1975, year_to_report, 5))+
       #scale_y_continuous(limits = c(150, 200)) +
       theme_few()+
       theme(legend.position="bottom",
           legend.key.width = unit(30, "pt"))

pt2EN <-  pt2+ylab("Day of year") + xlab("Year") + scale_color_manual(name="Pass", values=c("black","purple"))
ggsave(paste0("figures/",year_of_assess,"/spawning/spawning_season_", year_to_report,nameincl4X, "EN.png"), width = 7, height = 5, dpi = 600, units = "in")

pt2FR <-  pt2+ylab("Jour de l'ann\u00E9e") + xlab("Ann\u00E9e") + scale_color_manual(name="Trajet", values=c("black","purple"))
ggsave(paste0("figures/",year_of_assess,"/spawning/spawning_season_", year_to_report,nameincl4X, "FR.png"), width = 7, height = 5, dpi = 600, units = "in")

pt2BI <-  pt2+ylab("Jour de l'ann\u00E9e / Day of year") + xlab("Ann\u00E9e / Year")+ scale_color_manual(name="Trajet / Pass", values=c("black","purple"))
ggsave(paste0("figures/",year_of_assess,"/spawning/spawning_season_", year_to_report,nameincl4X, "BI.png"), width = 7, height = 5, dpi = 600, units = "in")


y <- c(120,160,200,240)
Sys.setlocale("LC_ALL", "English")
pt2dayEN <-  pt2 + scale_y_continuous(breaks=y, labels=format(as.Date(y, origin = "2001-01-01"),'%B-%d')) +ylab("Day of year") + xlab("Year")+ scale_color_manual(name="Pass", values=c("black","purple"))
ggsave(paste0("figures/spawning/spawning_season_", year_to_report,nameincl4X, "EN2.png"), width = 7, height = 5, dpi = 600, units = "in")
Sys.setlocale("LC_ALL", "French_Canada.utf8")
pt2dayFR <-  pt2 + scale_y_continuous(breaks=y, labels=format(as.Date(y, origin = "2001-01-01"),'%d-%B')) +ylab("Jour de l'ann\u00E9e") + xlab("Ann\u00E9e") + scale_color_manual(name="Trajet", values=c("black","purple"))
ggsave(paste0("figures/spawning/spawning_season_", year_to_report,nameincl4X, "FR2.png"), width = 7, height = 5, dpi = 600, units = "in")
}


