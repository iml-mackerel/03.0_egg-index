library(INLA)
library(inlabru)

# Get the spatial parameters of an INLA model
MySpatialParams <- function(Model, ThisSpde, wname="w") {
  SpFi <- inla.spde2.result(inla = Model, name = wname, spde = ThisSpde, do.transfer = TRUE) 
  Kappa  <- inla.emarginal(function(x) x, SpFi$marginals.kappa[[1]] )
  sigmau <- inla.emarginal(function(x) sqrt(x),SpFi$marginals.variance.nominal[[1]] )
  Range <- inla.emarginal(function(x) x, SpFi$marginals.range.nominal[[1]] )
  Out <- c(Kappa, sigmau, Range)
  names(Out) <- c("Kappa", "Sigma_u", "Range")
  Out
}


plotSpatialFieldCL <- function(dir, name, minw, maxw, my.year){
library(inlabru)
library(sf)
library(ggpubr)
  load(paste0(dir,name,".Rdata"))
years <- levels(as.factor(dat$year))
 second_panel<-ceiling((length(years)+1)/3)
 
 if(grepl(name, pattern="gaussian")){
  
w     <- mod$summary.random$w$mean

#set.panel(2,second_panel) # 2X2 matrix of plots
p <-  list()
rloc.projdeg <- list()
for (i in 1:length(years)){
  w.pm <- w[wRepl.index$w.repl == i]
  MyTitle <- years[i]
  
  stopifnot(length(w.pm) == mesh$n)
  
# inla.mesh.projector: it creates a lattice using the mesh and specified ranges. 
  proj <- inla.mesh.projector(mesh, 
                              xlim = range(mesh$loc[,1]), 
                              ylim = range(mesh$loc[,2]), 
                              dims = c(400, 400))
  # The function inla.mesh.project can then 
  # be used to project the w's on this grid.
  dat.loc = expand_grid(x=proj$x, y=proj$y)
  dat.loc$field <- inla.mesh.project(mesh,
                                      loc = as.matrix(dat.loc),
                                      field = w.pm)
  coordinates(dat.loc) <- c("x", "y")
  dat.loc<- as(dat.loc, "SpatialPixelsDataFrame"); proj4string(dat.loc) <- lcc
  r.loc<- raster(dat.loc); names(r.loc) <- "layer"
  rloc.projdeg[[i]] <-  projectRaster(r.loc, crs=projdeg)
  
 p[[i]]<-  basemap+gg(mask(rloc.projdeg[[i]], nec, inverse=T))+
           ggtitle(MyTitle)+
           #scale_fill_viridis_c(option="turbo", name = "w",guide=guide_colorbar(frame.colour = "black", ticks.colour="black"), limits=c(minw, maxw))
           scale_fill_gradient2(low="blue", high="red", mid="white", limits=c(minw, maxw), name="", guide=guide_colorbar(frame.colour = "black", ticks.colour="black"))+
          theme(legend.position="top", legend.key.width=unit(1.5,"cm"),
           legend.key.height=unit(0.7, "cm"),legend.text = element_text(size=12))
   
}


#plot correlation
# Here are the spatial parameters
SpatialParams <- MySpatialParams(Model = mod, ThisSpde = spde)
Kappa   <- SpatialParams[1]  #Parameter kappa for the Mattern correlation function
Sigma.u <- SpatialParams[2]  #Sigma of the u
Range   <- SpatialParams[3]  #expressed in metres

# Show correlation structure.
# Obtain the locations of each point of the mesh.
LocMesh <- mesh$loc[,1:2]

# Calculate the distance between each vertex.
D <- as.matrix(dist(LocMesh))

# Using the estimated parameters from the model (see above), we can 
# calculate the imposed Matern correlation values.
d.vec <- seq(0, max(D), length = 500)      
Cor.M <- (Kappa * d.vec) * besselK(Kappa * d.vec, 1) 
Cor.M[1] <- 1

Cor.data <- as.data.frame(cbind(d.vec,Cor.M))

# Which we plot here:
p[[length(years)+1]]<- ggplot(Cor.data,aes(x=d.vec/1000, y=Cor.M))+
  geom_line()+basetheme +theme(plot.background = element_blank()) +xlab("Distance (km)") + ylab("Correlation")+ geom_vline(xintercept = Range/1000, lty=2)+ggtitle("")#+
 # scale_x_continuous(breaks=seq(0, Range*3/1000, 10), limits=c(0, Range*3/1000))

to_return<- ggarrange(plotlist = p, ncol=7, nrow=second_panel, legend="top", common.legend = T)
#ggsave(paste0("figures/",year_of_assess,"/INLA/random_field/",name,".png"), dpi=600, units="in", width=2*7, height=2*second_panel)
}


if(grepl(name, pattern="ZAG")){

  w     <- Gbern$summary.random$w$mean
  wGamma     <- Ggamma$summary.random$wGamma$mean
  
  g <-  list()
  b<-  list()
  
  
  rlocg.projdeg <- list()
  rlocb.projdeg <- list()
  for (i in 1:length(years)){
    w.pm <- w[wReplBern.index$w.repl == i]
    wGamma.pm <- wGamma[wReplGamma.index$wGamma.repl == i]
    MyTitle <- years[i]
    
    stopifnot(length(w.pm) == mesh$n)
    stopifnot(length(wGamma.pm) == mesh$n)
    
    # inla.mesh.projector: it creates a lattice using the mesh and specified ranges. 
    proj <- inla.mesh.projector(mesh, 
                                xlim = range(mesh$loc[,1]), 
                                ylim = range(mesh$loc[,2]), 
                                dims = c(400, 400))
    # The function inla.mesh.project can then 
    # be used to project the w's on this grid.
    dat.loc = expand_grid(x=proj$x, y=proj$y)
    dat.loc$field <- inla.mesh.project(mesh,
                                       loc = as.matrix(dat.loc),
                                       field = w.pm)
    dat.loc$fieldGamma <- inla.mesh.project(mesh,
                                       loc = as.matrix(dat.loc),
                                       field = wGamma.pm)
    
    coordinates(dat.loc) <- c("x", "y")
    dat.loc<- as(dat.loc, "SpatialPixelsDataFrame"); proj4string(dat.loc) <- lcc
    r.locb<- raster(dat.loc, layer="field"); names(r.locb) <- "layer"
    r.locg<- raster(dat.loc, layer="fieldGamma"); names(r.locg) <- "layer"
    rlocb.projdeg[[i]] <-  projectRaster(r.locb, crs=projdeg)
    rlocg.projdeg[[i]] <-  projectRaster(r.locg, crs=projdeg)
    
    b[[i]]<-  basemap2+gg(mask(rlocb.projdeg[[i]], nec, inverse=T))+
      ggtitle(MyTitle)+
      #scale_fill_viridis_c(option="turbo", name = "w",guide=guide_colorbar(frame.colour = "black", ticks.colour="black"), limits=c(minw, maxw))
      scale_fill_gradient2(low="blue", high="red", mid="white", limits=c(minw, maxw), name="", guide=guide_colorbar(frame.colour = "black", ticks.colour="black"))+
      theme(legend.position="top", legend.key.width=unit(1.5,"cm"),
            legend.key.height=unit(0.7, "cm"),legend.text = element_text(size=12))
    g[[i]]<-  basemap2+gg(mask(rlocg.projdeg[[i]], nec, inverse=T))+
      ggtitle(MyTitle)+
      #scale_fill_viridis_c(option="turbo", name = "w",guide=guide_colorbar(frame.colour = "black", ticks.colour="black"), limits=c(minw, maxw))
      scale_fill_gradient2(low="blue", high="red", mid="white", limits=c(minw, maxw), name="", guide=guide_colorbar(frame.colour = "black", ticks.colour="black"))+
      theme(legend.position="top", legend.key.width=unit(1.5,"cm"),
            legend.key.height=unit(0.7, "cm"),legend.text = element_text(size=12))
    
  }
  
  
  #plot correlation
  # Here are the spatial parameters
  SpatialParamsb <- MySpatialParams(Model = Gbern, ThisSpde = spde)
  Kappab   <- SpatialParamsb[1]  #Parameter kappa for the Mattern correlation function
  Sigma.ub <- SpatialParamsb[2]  #Sigma of the u
  Rangeb   <- SpatialParamsb[3]  #expressed in metres
  
  
  SpatialParamsg <- MySpatialParams(Model = Ggamma, ThisSpde = spde, wname="wGamma")
  Kappag   <- SpatialParamsg[1]  #Parameter kappa for the Mattern correlation function
  Sigma.ug <- SpatialParamsg[2]  #Sigma of the u
  Rangeg   <- SpatialParamsg[3]  #expressed in metres
  
  # Show correlation structure.
  # Obtain the locations of each point of the mesh.
  LocMesh <- mesh$loc[,1:2]
  
  # Calculate the distance between each vertex.
  D <- as.matrix(dist(LocMesh))
  
  # Using the estimated parameters from the model (see above), we can 
  # calculate the imposed Matern correlation values.
  d.vec <- seq(0, max(D), length = 500)      
  Cor.Mb <- (Kappab * d.vec) * besselK(Kappab * d.vec, 1) 
  Cor.Mb[1] <- 1
  Cor.datab <- as.data.frame(cbind(d.vec,Cor.Mb))
  
  Cor.Mg <- (Kappag * d.vec) * besselK(Kappag * d.vec, 1) 
  Cor.Mg[1] <- 1
  Cor.datag <- as.data.frame(cbind(d.vec,Cor.Mg))
  
  maxb <- ifelse(Rangeb*3/1000 > 100, Rangeb*2/1000, Rangeb*3/1000)
  maxg <- ifelse(Rangeg*3/1000 > 100, Rangeg*2/1000, Rangeg*3/1000)
  
  
  # Which we plot here:
  b[[length(years)+1]]<- ggplot(Cor.datab,aes(x=d.vec/1000, y=Cor.Mb))+
    geom_line()+basetheme +xlab("Distance (km)") + ylab("Correlation")+ geom_vline(xintercept = Rangeb/1000, lty=2)+ggtitle("")+
    scale_x_continuous(breaks=seq(0, maxb, ifelse(Rangeb > 100, 100,50)), limits=c(0, maxb))+ theme(axis.text.x = element_text(angle=90))
  
  g[[length(years)+1]]<- ggplot(Cor.datag,aes(x=d.vec/1000, y=Cor.Mg))+
    geom_line()+basetheme +xlab("Distance (km)") + ylab("Correlation")+ geom_vline(xintercept = Rangeg/1000, lty=2)+ggtitle("")+
    scale_x_continuous(breaks=seq(0, maxg, ifelse(Rangeg > 100, 100,50)), limits=c(0, maxg)) + theme(axis.text.x = element_text(angle=90))
  
  
 bgg<-  ggarrange(plotlist = b[(length(years)-1):(length(years)+1)], ncol=3, nrow=1, legend="top", common.legend = T)
  ggsave(paste0("img/",my.year,"/INLA/random_field/",name,"_ZAG_bernouilli.png"), dpi=600, units="in", width=5*3, height=5)
  
  
 ggg<- ggarrange(plotlist = g[(length(years)-1):(length(years)+1)], ncol=3, nrow= 1, legend="top", common.legend = T)
  ggsave(paste0("img/",my.year,"/INLA/random_field/",name,"_ZAG_gamma.png"), dpi=600, units="in",width=5*3, height=5)
   to_return<- list(bgg, ggg)
 }

 
 if(grepl(name, pattern="tw")){
   
   w     <- Gtw$summary.random$w$mean
  ptw<-  list()
   
   
   rloct.projdeg <- list()
   for (i in 1:length(years)){
     w.pm <- w[wRepltw.index$w.repl == i]
     MyTitle <- years[i]
     
     stopifnot(length(w.pm) == mesh$n)
    
     # inla.mesh.projector: it creates a lattice using the mesh and specified ranges. 
     proj <- inla.mesh.projector(mesh, 
                                 xlim = range(mesh$loc[,1]), 
                                 ylim = range(mesh$loc[,2]), 
                                 dims = c(400, 400))
     # The function inla.mesh.project can then 
     # be used to project the w's on this grid.
     dat.loc = expand_grid(x=proj$x, y=proj$y)
     dat.loc$field <- inla.mesh.project(mesh,
                                        loc = as.matrix(dat.loc),
                                        field = w.pm)
        
     coordinates(dat.loc) <- c("x", "y")
     dat.loc<- as(dat.loc, "SpatialPixelsDataFrame"); proj4string(dat.loc) <- lcc
     r.locb<- raster(dat.loc, layer="field"); names(r.locb) <- "layer"
    rloct.projdeg[[i]] <-  projectRaster(r.locb, crs=projdeg)

     
     ptw[[i]]<-  basemap2+gg(mask(rloct.projdeg[[i]], nec, inverse=T))+
       ggtitle(MyTitle)+
       #scale_fill_viridis_c(option="turbo", name = "w",guide=guide_colorbar(frame.colour = "black", ticks.colour="black"), limits=c(minw, maxw))
       scale_fill_gradient2(low="blue", high="red", mid="white", limits=c(minw, maxw), name="", guide=guide_colorbar(frame.colour = "black", ticks.colour="black"))+
       theme(legend.position="top", legend.key.width=unit(1.5,"cm"),
             legend.key.height=unit(0.7, "cm"),legend.text = element_text(size=12))
   }
   
   
   #plot correlation
   # Here are the spatial parameters
   SpatialParamsb <- MySpatialParams(Model = Gtw, ThisSpde = spde)
   Kappab   <- SpatialParamsb[1]  #Parameter kappa for the Mattern correlation function
   Sigma.ub <- SpatialParamsb[2]  #Sigma of the u
   Rangeb   <- SpatialParamsb[3]  #expressed in metres
   
  
   # Show correlation structure.
   # Obtain the locations of each point of the mesh.
   LocMesh <- mesh$loc[,1:2]
   
   # Calculate the distance between each vertex.
   D <- as.matrix(dist(LocMesh))
   
   # Using the estimated parameters from the model (see above), we can 
   # calculate the imposed Matern correlation values.
   d.vec <- seq(0, max(D), length = 500)      
   Cor.Mb <- (Kappab * d.vec) * besselK(Kappab * d.vec, 1) 
   Cor.Mb[1] <- 1
   Cor.datab <- as.data.frame(cbind(d.vec,Cor.Mb))
   
   maxt <- ifelse(Rangeb*3/1000 > 100, Rangeb*2/1000, Rangeb*3/1000)
   
   
   # Which we plot here:
   ptw[[length(years)+1]]<- ggplot(Cor.datab,aes(x=d.vec/1000, y=Cor.Mb))+
     geom_line()+basetheme +xlab("Distance (km)") + ylab("Correlation")+ geom_vline(xintercept = Rangeb/1000, lty=2)+ggtitle("")+
     scale_x_continuous(breaks=seq(0, maxt, ifelse(Rangeb > 100, 100,50)), limits=c(0, maxt))+ theme(axis.text.x = element_text(angle=90))
   
 
   #returns 2 years
   bgg<-  ggarrange(plotlist = ptw[(length(years)-1):(length(years)+1)], ncol=3, nrow=1, legend="top", common.legend = T)
   #ggsave(paste0("img/",my.year,"/INLA/random_field/",name,"_tw.png"), dpi=600, units="in", width=2*2, height=2)
   
   to_return<- list(bgg)
 }
 
 return(to_return) 
}

