library(INLA)
library(inlabru)
library(raster)
get_prediction_grid<- function(dir, dir.out,name, trans="1_5"){
  
  if(grepl(name, pattern="gaussian")){
  if(trans=="1_5"){ inverse_trans <-  function(x) x^5}
 load(paste0(dir,name,".Rdata")) 
  
  if(!grepl(mod$call[1], pattern="family = \"gaussian\"")) stop("extraction other than gaussian not implemented yet")

 id.fit <- inla.stack.index(stkall.Repl, 'Fit')$dat
 id.prd <- inla.stack.index(stkall.Repl, 'Pred')$dat
 
if(!exists("df.grid")) stop("df.grid not loaded with model, please load a df.grid into the environment")
  
# Extract mean from fitted values
mean.prd <- mod$summary.fitted.values$mean[id.prd] # not correct for lognormal
#sd.prd   <- I$summary.fitted.values$sd[id.prd]
ic1.prd <- mod$summary.fitted.values$`0.025quant`[id.prd] # not correct for lognormal
ic2.prd <-mod$summary.fitted.values$`0.975quant`[id.prd] # not correct for lognormal

tau<-mod$marginals.hyperpar$`Precision for the Gaussian observations`
MySqrt <- function(x) { 1 / sqrt(x) }
sigma <- inla.emarginal(MySqrt, tau)
sigma


df.grid$Fit <- inverse_trans(mean.prd + sigma^2 /2)
df.grid$IC1 <- inverse_trans(ic1.prd + sigma^2 /2)
df.grid$IC2 <- inverse_trans(ic2.prd + sigma^2 /2)



rast_pred <- stack()
for(y in unique(df.grid$year)){
  cat(y)
  predy<- df.grid %>%  filter(year==y, Y.m >=460631) %>% dplyr::select(X.m, Y.m, Fit)#pour régler un problème pour 2012 plus petit
  
  gridded(predy) = ~X.m+Y.m
  proj4string(predy) <- lcc
  
  rasty<-raster::mask(projectRaster(raster(predy),crs=projdeg), nec, inverse=T); names(rasty) <- y
  
  rast_pred<- stack(rast_pred, rasty)

}
  writeRaster(rast_pred, paste0(dir.out,"predictions_", name, ".grd"), overwrite=T)
}
  
  if(grepl(name, pattern="ZAG")){
    load(paste0(dir,name,".Rdata")) 
    
    idb.prd <- inla.stack.index(stkallBern, 'Pred')$dat
    
    idg.prd <- inla.stack.index(stkallGamma, 'Pred')$dat
    
    
    if(!exists("df.grid")) stop("df.grid not loaded with model, please load a df.grid into the environment")
     library(boot)
       # Extract mean from fitted values
    Pi <- inv.logit(Gbern$summary.fitted.values$mean[idb.prd])
    #sd.prd   <- I$summary.fitted.values$sd[id.prd]
    mu<- Ggamma$summary.fitted.values$mean[idg.prd]
    
   df.grid$Fit <- Pi * mu
   r  <- Ggamma$summary.hyperpar[1,"mean"]
   df.grid$var <- ((Pi * r + pi - pi * r) / r) * mu^2
   
   if(grepl(name, pattern="station")){
     write.table(df.grid, paste0(dir.out,"predictions_", name, ".txt"), row.names=F, sep="\t", dec=".")
   }
     
     
   if(!grepl(name, pattern="station")){
    rast_pred <- stack()
    for(y in unique(df.grid$year)){
       
      predy<- df.grid[which(df.grid$year==y),c("X.m", "Y.m", "Fit")]
      
     # df.grid %>% dplyr::filter(year == y) %>% dplyr::select(X.m, Y.m, Fit)# je ne sais pas pourquoi mais ne fonctionne pas.
       
      gridded(predy) = ~X.m+Y.m
      proj4string(predy) <- lcc
      
      rasty<-raster::mask(projectRaster(raster(predy),crs=projdeg), nec, inverse=T); names(rasty) <- y
      
      rast_pred<- stack(rast_pred, rasty)
    }
    writeRaster(rast_pred, paste0(dir.out,"predictions_", name, ".grd"), overwrite=T)
    }
  }  

}