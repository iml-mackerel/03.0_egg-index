source(paste0("R/",year_of_assess,"/functions/spatial_projections.R"))


library(ggplot2)
library(raster)
#library(ggmap)
library(rgeos)
library(INLA)
library(mgcv)
library(fields)
library(sp)
library(gstat)
library(maps)
library(maptools)
library(mapdata)
data("worldHiresMapEnv")


make_boundary <- function(dat, year_to_report, trajet){
        
  
  coord_sp <- create_coord(dat[, c("longitude", "latitude")], crsin = projdeg, crsout = lcc)
  dat$X.m <-coord_sp$x
  dat$Y.m <- coord_sp$y
  Loc <- cbind(dat$X.m, dat$Y.m)
  
buf1 <- gBuffer(coord_sp, width=40000, byid=TRUE)
buf2 <- gUnaryUnion(buf1)
buf <-  gBuffer(buf2, width=-20000)
outerRings = Filter(function(f){f@ringDir==1},buf@polygons[[1]]@Polygons)
AreaSP= SpatialPolygons(list(Polygons(outerRings,ID=1)))
proj4string(AreaSP) <- lcc

#copier de l'appendix A
#simplier coastline pour cr?er un polygon avec l'eau seulement
CoastPoly <- maps::map("worldHires", 
                 regions = c("Canada"), 
                 #exact = TRUE,
                 fill = TRUE, 
                 col = "transparent",
                 plot = TRUE,
                 ylim = c(45,53),
                 xlim = c(-70, -60))


# Convert the polygon into metric coordinates
IDs <- sapply(strsplit(CoastPoly$names, ":"), function(x) x[1])
CoastPoly_sp <- map2SpatialPolygons(CoastPoly, 
                                    IDs = IDs,
                                    proj4string = CRS("+proj=longlat +datum=NAD83"))

#buffer around coastline, otherwise point intersect or are outside coastline
buf_coast <- spTransform(gBuffer(spTransform(CoastPoly_sp, lcc), width=-2000, byid=TRUE), projdeg)


CoastSP.smooth <- spTransform(thinnedSpatialPoly(buf_coast, tolerance=0.001, minarea=0.0001,
                                     topologyPreserve = T,
                                     avoidGEOS=F), proj4string(AreaSP))


Water.lcc <- gDifference(AreaSP, CoastSP.smooth)
library(sf)
Watersf<- st_as_sf(Water.lcc)


##waterlcc is the study area ->save the study area for the prediction grid

if(file.exists(paste0("data/",year_of_assess,"/INLA/prediction_grid_area",year_to_report,"_tajet",trajet,".shp"))) st_write(Watersf, dsn=paste0("data/",year_of_assess,"/INLA/prediction_grid_area",year_to_report,"_tajet",trajet,".shp"),delete_dsn=T)
if(!file.exists(paste0("data/",year_of_assess,"/INLA/prediction_grid_area",year_to_report,"_tajet",trajet,".shp"))) st_write(Watersf, dsn=paste0("data/",year_of_assess,"/INLA/prediction_grid_area",year_to_report,"_tajet",trajet,".shp"),delete_dsn=F)

return(Loc)
}



make_mesh<- function(Loc, range=5, year_to_report, trajet){
  library(rgdal)
  Water.lcc<- readOGR(paste0("data/",year_of_assess,"/INLA/prediction_grid_area",year_to_report,".shp"))
  
# 22.6.2 Mesh
#boundary is coastline and liimits of samples
Boundary <- inla.sp2segment(Water.lcc)


Range <- range * 1000
MaxEdge <- Range/ 5

MyCutoff = MaxEdge/2

mesh <- inla.mesh.2d(boundary = Boundary, 
                     max.edge = c(1, 5) * MaxEdge,  #For a paper
                     cutoff = MyCutoff)


png(paste0("figures/",year_of_assess,"/INLA/mesh",range,"km",year_to_report,"_tajet",trajet,".png"), res=600, height=6, width=6, units="in")
par(mfrow = c(1,1), mar=c(1, 1, 1, 1))
plot(mesh, asp = 1)
points(Loc, col = 1, pch = 16, cex = 0.2)
dev.off()

saveRDS(mesh, file=paste0("data/",year_of_assess,"/INLA/mesh",range,"km_",year_to_report,"_tajet",trajet,".rds"), version=2)
}







