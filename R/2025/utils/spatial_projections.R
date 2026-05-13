library(raster)

lcc<-CRS("+proj=lcc 
               +lat_1=46
               +lat_2=60
               +lat_0=44
               +lon_0=-68.5
               +x_0=0
               +y_0=0
               +datum=NAD83")
# correspond au point le plus central de la zone de données (pour obtenir les mesures les moins déformées possible)

projdeg<-CRS("+proj=longlat +datum=WGS84")


create_coord<- function (xy, crsin = projdeg, crsout = lcc) 
{
  coord <- as.data.frame(xy)
  colnames(coord) <- c("x", "y")
  sp::coordinates(coord) <- ~x + y
  sp::proj4string(coord) <- crsin
  sp::spTransform(coord, crsout)
}

#library(rgdal)
#batall<-readOGR("data/isobathsJMC/contour10m.shp")
#bat50<-subset(batall, Profond=="50")
#bat60<-subset(batall, Profond=="60")
#bat80<-subset(batall, Profond=="80")
#bat100n<-subset(batall, Profond=="100")
#bat200n<-subset(batall, Profond=="200")
#not availbale anymore use marmap. and geom_controu