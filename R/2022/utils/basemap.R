library(dplyr)
library(tidyr)
library(ggthemes)
#windowsFonts(A=windowsFont("Arial"))
#remotes::install_github("ropensci/rnaturalearthhires")
library(rnaturalearth)
library(rnaturalearthhires)
library(rnaturalearthdata)
library(ggsn)

nec<- ne_countries(country='canada', scale="large")
basetheme <- theme_few()+ theme(axis.title= element_text(face="bold", colour="black", size=13),
                                axis.text = element_text(color="black", size=11),
                                strip.text =element_text(face="bold",color="black", size=13),
                                legend.justification = "left")

basemap <-  ggplot()+
  geom_polygon(data = fortify(nec), aes(x=long, y=lat, group=group), fill="bisque2", col="bisque4") +
  coord_quickmap(xlim=c(-66.5,-60), ylim=c(45.5,49.5))+ #scale_y_continuous(breaks=c(46,47,48,49))+
  xlab("Longitude")+ ylab("Latitude") + basetheme+
  ggsn::scalebar(dist = 50, dist_unit = "km",
         transform = TRUE, model = "WGS84", st.dist=0.04, st.size=2.5,st.color="grey20",height = 0.04,border.size = 0.3, x.min=-66.3, x.max=-60, y.min=46, y.max=49.5, location="topleft") +
  ggsn::north(symbol=3, scale = 0.1, x.min=-66.3, x.max=-60, y.min=46, y.max=49.5, location="topright")
  
basemap2 <-  ggplot()+
  geom_polygon(data = fortify(nec), aes(x=long, y=lat, group=group), fill="bisque2", col="bisque4") +
  coord_quickmap(xlim=c(-66.5,-60), ylim=c(45.5,49.5))+ #scale_y_continuous(breaks=c(48,49,50))+
  xlab("")+ ylab("") + basetheme

