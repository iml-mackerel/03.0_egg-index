
make_boundary <- function(dat, year, trajet){
  
  dat<- st_as_sf(dat, coords=c("longitude", "latitude"), crs=4326)
  dat<- st_transform(dat, lcc)
  
  Loc<- st_coordinates(dat)
  colnames(Loc) <-  c("X.m", "Y.m")
  
  
buf1<- st_buffer(dat, dist=40000)  
buf2<- st_union(buf1)  
AreaSP <-  st_buffer(buf2, dist=-20000) %>%  st_simplify(dTolerance=3000)

CoastPoly_sp<- nec %>% dplyr::select(geometry, name_en) %>%  st_transform(lcc)
#buffer around coastline, otherwise point intersect or are outside coastline
buf_coast <- st_buffer(CoastPoly_sp, dist=-5000)
Watersf <- st_difference(AreaSP, CoastPoly_sp)

##waterlcc is the study area ->save the study area for the prediction grid
if(file.exists(paste0("data/",year,"/INLA/prediction_grid_area",year,"_trajet",trajet,".shp"))) st_write(Watersf, dsn=paste0("data/",year,"/INLA/prediction_grid_area",year,"_trajet",trajet,".shp"),delete_dsn=T)
if(!file.exists(paste0("data/",year,"/INLA/prediction_grid_area",year,"_trajet",trajet,".shp"))) st_write(Watersf, dsn=paste0("data/",year,"/INLA/prediction_grid_area",year,"_trajet",trajet,".shp"),delete_dsn=F)

return(Loc)
}



make_mesh<- function(Loc, range=5, year, trajet){
  
  Water.lcc<- st_read(paste0("data/",year,"/INLA/prediction_grid_area",year,"_trajet",trajet,".shp"))
  
# 22.6.2 Mesh
#boundary is coastline and liimits of samples
Boundary <- inla.sp2segment(Water.lcc)

Range <- range * 1000
MaxEdge <- Range/ 5

MyCutoff = MaxEdge/2

mesh <- inla.mesh.2d(boundary = Boundary, 
                     max.edge = c(1, 5) * MaxEdge,  #For a paper
                     cutoff = MyCutoff)

png(paste0("img/",year,"/INLA/mesh",range,"km_",year,"_trajet",trajet,".png"), res=600,units="in", width=6, height=6, bg="transparent")
plot(mesh, asp = 1)
points(Loc, col = 1, pch = 16, cex = 0.2)
dev.off()

saveRDS(mesh, file=paste0("data/",year,"/INLA/mesh",range,"km_",year,"_trajet",trajet,".rds"), version=2)
}







