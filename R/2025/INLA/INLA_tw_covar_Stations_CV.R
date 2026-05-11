
INLA_tw_CV <- function(dat, res_mesh, res_pred, Rvar="DEP", trajet, year){
 
  dat$DEP <-  dat[, Rvar]
  
 dat<- dat %>%  mutate(PA = ifelse(DEP > 0,1,0),
                       BiomassPos= ifelse(PA >0, DEP, NA)) #pour éviter na.omit

if(length(unique(dat$PA))==1) stop("Presence/Absence include only one value")

dat <-  dat %>%  mutate(year=as.factor(year), station=as.factor(station)) %>%  
                        dplyr::select(sample_id, DEP,PA,BiomassPos, station, latitude, longitude, year) 


dat<- dat %>%  st_as_sf(coords=c("longitude", "latitude"), crs=4326) %>%  st_transform(crs=lcc)

dat<- bind_cols(dat, st_coordinates(dat) %>% as.data.frame() %>%  rename(X.m=X, Y.m=Y))

Loc <- dat %>%  dplyr::select(X.m, Y.m)

mesh <- readRDS(paste0("data/",year,"/INLA/mesh",res_mesh,"km_",year,"_trajet",trajet,".rds"))


Repl <- as.numeric(as.factor(dat$year))
table(Repl)
NRepl <- length(unique(Repl))
NRepl #Number of replicated random fields we will get

# Define the weight factors aik
#Repl: spatial correlation that changes randomly between years
A.Repl <- inla.spde.make.A(mesh, 
                               loc = Loc, 
                               repl = Repl)

# Define the SPDE 
Range0      <- res_mesh * 1000  
AlphaRange0 <- 0.05  
Sigma0      <- 1      
AlphaSigma0 <- 0.5  

spde <- inla.spde2.pcmatern(mesh,  prior.range = c(Range0, AlphaRange0), 
                                  prior.sigma = c(Sigma0, AlphaSigma0))


# Define the spatial field
wRepltw.index <- inla.spde.make.index('w', 
                                      n.spde = mesh$n,
                                      n.repl = NRepl)

XYear <- model.matrix(~ year, data = dat)




MyStd <- function(x) { (x - mean(x)) / sd(x)}

N <- nrow(dat)
Covariates.Repl <- data.frame(
  Intercept= rep(1,N),
  #years=XYear[,2:ncol(XYear)], # pas besoin si on le mets pas dans le modèle juste dans le spatiotemporal
  station = dat$station
)


StackRepltw <- inla.stack(
  tag  = "twFit",
  data = list(DEP = dat$DEP),  
  A    = list(1, A.Repl),#, A.Repl),                      
  effects = list(                 
    Covariates.Repl,#Intercept
    w   = wRepltw.index))  #spatial correlation
# Shared       = wRepl.index))  #Shared spatial correlation  



#do the formula
fGamtw <- DEP ~ -1 + Intercept + 
  station +
  f(w, model = spde, replicate = w.repl)


stkalltw<- inla.stack(StackRepltw) 


Gtw <- inla(fGamtw,
              family = "tweedie", 
              #control.family = list(link = "log"),
              #lincomb = All.lcs,
              data=inla.stack.data(stkalltw),
              control.compute = list(dic = TRUE, config=T),
              control.predictor = list(compute=T,A = inla.stack.A(stkalltw), link=1),
              control.inla=list(strategy="gaussian"))


exi<- exists(paste0("results/",year,"/INLA/",Rvar,"_INLA_mesh", res_mesh,"km_pred_",res_pred,"_tw",year,"_CV_trajet",trajet,".RData"))
if(!exi | new) save(Gtw, dat, mesh,wRepltw.index, stkalltw, spde,file=paste0("results/",year,"/INLA/",Rvar,"_INLA_mesh", res_mesh,"km_pred_",res_pred,"_tw",year,"_CV_trajet",trajet,".RData"))
if(exi & !new) stop("Preventing overwite. The model was not saved")
}


