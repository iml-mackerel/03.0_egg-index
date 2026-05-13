library(dplyr)
library(tidyr)
library(ggthemes)
#windowsFonts(A=windowsFont("Arial"))
#remotes::install_github("ropensci/rnaturalearthhires")
library(rnaturalearth)
library(rnaturalearthhires)
library(rnaturalearthdata)

nec<- ne_countries(country='canada', scale="large")
basetheme <- theme_few()+ theme(axis.title= element_text(face="bold", colour="black", size=13),
                                axis.text = element_text(color="black", size=11),
                                strip.text =element_text(face="bold",color="black", size=13),
                                legend.justification = "left")

basemap <-  ggplot()+
  geom_sf(data = nec, fill="bisque2", col="bisque4") +
  coord_sf(xlim=c(-66.5,-60), ylim=c(45.5,49.5))+ #scale_y_continuous(breaks=c(46,47,48,49))+
  xlab("Longitude")+ ylab("Latitude") + basetheme
 
basemap2 <-  ggplot()+
  geom_sf(data = nec, fill="bisque2", col="bisque4") +
  coord_sf(xlim=c(-66.5,-60), ylim=c(45.5,49.5))+ #scale_y_continuous(breaks=c(48,49,50))+
  xlab("")+ ylab("") + basetheme

