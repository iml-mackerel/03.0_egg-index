#NOT REALLY WORKING!!!###


library(sdmTMB)
library(sdmTMBextra)

library(sf)
library(tidyverse)
library(INLA)
Water.lcc<- st_read(paste0("data/2024/INLA/prediction_grid_area2024_trajet1.shp"))

# 22.6.2 Mesh
#boundary is coastline and liimits of samples
Boundary <- inla.sp2segment(Water.lcc)

range=150
Range <- range * 1000
MaxEdge <- Range/ 2

MyCutoff = MaxEdge/2


meshin <- inla.mesh.2d(boundary = Boundary, 
                     max.edge = c(1, 5) * MaxEdge,  #For a paper
                     cutoff = MyCutoff)

load(file=paste0("data/2024/eggt1.RData"))
eggt1a<- eggt1 %>%  st_as_sf(coords=c("longitude", "latitude"), crs=4326, remove=F) %>% 
  st_transform(crs=st_crs(Water.lcc))

eggt1b<- bind_cols(eggt1a %>%  st_drop_geometry(), st_coordinates(eggt1a)/1000) %>% 
  mutate(fstation =as.factor(station),
         fyear=as.factor(year))

mesh <- make_mesh(eggt1b, c("X", "Y"), mesh = meshin)
plot(mesh)

plot_pc_matern(range_gt = 150, sigma_lt = 1)

m1<- sdmTMB(DEP ~ (1 | fyear),
        time="year",
       data=eggt1b,
        spatial="on",
         mesh=mesh,  family = delta_gamma(),
       priors = sdmTMBpriors(
         matern_s = pc_matern(range_gt = 150, sigma_lt = 1))
       )

m1<- sdmTMB(DEP ~1,
            time="year",
            data=eggt1b,
            extra_time = c(1980, 1981, 1982, 1995, 1997, 2020),
            spatial="on",
            spatiotemporal = "AR1",
            mesh=mesh,  family = delta_gamma(),
            priors = sdmTMBpriors(
              matern_s = pc_matern(range_gt = 150, sigma_lt = 1)),
)

sanity(m1)
#anisotropie impossible when using priors. 

newgrid = expand_grid(X=seq(floor(min(eggt1b$X)), ceiling(max(eggt1b$X)), 5), 
                      Y=seq(floor(min(eggt1b$Y)), ceiling(max(eggt1b$Y)), 5))
grid_yrs <- replicate_df(newgrid, "year", unique(eggt1b$year))

#library(ggeffects)
predictions <- predict(m1, newdata = grid_yrs) %>%  mutate(est=rlogis(est1)*exp(est2))
#ggpredict(set_delta_model(m1, model=NA), newdata=grid_yrs)


plot_map <- function(dat, column) {
  ggplot(dat, aes(X, Y, fill = {{ column }})) +
    geom_raster() +
    coord_fixed()
}

plot_map(predictions, est) +
  scale_fill_viridis_c(
  ) +
  facet_wrap(~year) +
  ggtitle("Prediction (fixed effects + all random effects)",
          subtitle = paste("maximum estimated biomass density =", round(max(exp(predictions$est))))
  )

