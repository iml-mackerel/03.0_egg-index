
INLA_ZAG <- function(dat, res_mesh, res_pred, Rvar="DEP", year_to_report, year_of_assess, trajet, serie=""){
 
  dat$DEP <-  dat[, Rvar]
  
 dat<- dat %>%  mutate(PA = ifelse(DEP > 0,1,0),
                       BiomassPos= ifelse(PA >0, DEP, -999)) #avoid to remove in next na.omit. return to NA after. 

if(length(unique(dat$PA))==1) stop("Presence/Absence include only one value")

dat <-  na.omit(dat %>%  mutate(year=as.factor(year), station=as.factor(station)) %>%  
                        dplyr::select(sample_id,DEP,PA,BiomassPos, station, latitude, longitude, year) ) 

dat <- dat %>%  mutate(BiomassPos=ifelse(BiomassPos ==-999, NA, BiomassPos))

dat<- dat %>%  st_as_sf(coords=c("longitude", "latitude"), crs=4326) %>%  st_transform(crs=lcc)

dat<- bind_cols(dat, st_coordinates(dat) %>% as.data.frame() %>%  rename(X.m=X, Y.m=Y))

Loc <- dat %>%  dplyr::select(X.m, Y.m)

mesh <- readRDS(paste0("data/",year_of_assess,"/INLA/mesh",res_mesh,"km_",year_to_report,"_trajet",trajet,".rds"))


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
wReplBern.index <- inla.spde.make.index('w', 
                                        n.spde = mesh$n,
                                        n.repl = NRepl)
wReplGamma.index <- inla.spde.make.index('wGamma', 
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


StackReplBern <- inla.stack(
  tag  = "BernFit",
  data = list(PA = dat$PA),  
  A    = list(1, A.Repl),#, A.Repl),                      
  effects = list(                 
    Covariates.Repl,#Intercept
    w   = wReplBern.index))  #spatial correlation
# Shared       = wRepl.index))  #Shared spatial correlation  

StackReplGamma <- inla.stack(
  tag  = "GammaFit",
  data = list(BiomassPos = dat$BiomassPos),  
  A    = list(1, A.Repl),#, A.Repl),                      
  effects = list(                 
    Covariates.Repl,#Intercept
    wGamma        = wReplGamma.index))  #spatial correlation
# Shared       = wRepl.index))  #Shared spatial correlation  


#do the formula
  fGambern <- PA ~ -1 + Intercept + 
    station +
    f(w, model = spde, replicate = w.repl)
  
  fGampos <- BiomassPos ~ -1 + Intercept + 
    station +
    f(wGamma, model = spde, replicate = wGamma.repl)


#prediciton data and stack
df.grid <- readRDS(paste0("data/",year_of_assess,"/INLA/prediction_grid_station",year_to_report,".RDS"))

df.grid<-df.grid %>%  filter(year %in% (unique(dat$year))) %>%  mutate(year=as.factor(year))
df.grid <-  na.omit(df.grid %>%  dplyr::select(station, X.m, Y.m, year))


loc.prd <- cbind(df.grid$X.m, df.grid$Y.m)

Repl.prd <- as.numeric(as.factor(df.grid$year))

ARepl.prd<-inla.spde.make.A(mesh, loc.prd, repl=Repl.prd)


XYear.pred <- model.matrix(~ year, data = df.grid)
library(mgcv)


XYear.pred <- model.matrix(~ year, data = df.grid)

N <- nrow(df.grid)
Covariates.pred.Repl <- data.frame(
  Intercept= rep(1,N),
  # years= XYear.pred[,2:ncol(XYear.pred)],
  station=df.grid$station
)

StackReplBern.pred <- inla.stack(
  tag = "Pred",
  data = list(PA = NA),  
  A = list(1, ARepl.prd),                  
  effects = list(   
    Covariates = Covariates.pred.Repl,
    w          = wReplBern.index))

StackReplGamma.pred <- inla.stack(
  tag = "Pred",
  data = list(BiomassPos = NA),  
  A = list(1, ARepl.prd),                  
  effects = list(   
    Covariates = Covariates.pred.Repl,
    wGamma          = wReplGamma.index))


stkallBern<- inla.stack(StackReplBern, StackReplBern.pred) 
stkallGamma<- inla.stack(StackReplGamma, StackReplGamma.pred) 


Gbern <- inla(fGambern,
              family = "binomial", 
              #lincomb = All.lcs,
              data=inla.stack.data(stkallBern),
              control.compute = list(dic = TRUE, config=T),
              control.predictor = list(compute=T,A = inla.stack.A(stkallBern)),
              control.inla=list(strategy="gaussian"))


Ggamma <- inla(fGampos,
               family = "gamma", 
               #lincomb = All.lcs,
               data=inla.stack.data(stkallGamma),
               control.compute = list(dic = TRUE, config=T),
               control.predictor = list(link=1,compute=T,A = inla.stack.A(stkallGamma)),
               control.inla=list(strategy="gaussian"))

exi<- exists(paste0("results/",year_of_assess,"/INLA/",Rvar,"_INLA_mesh", res_mesh,"km_pred_",res_pred,"_ZAG",year_to_report,"_trajet",trajet,serie,".RData"))
overw =NA
if(exi) overw<- menu(c("Yes", "No"), title = "Model already exist. Do you want to overwrite?")
if(!exi | overw ==1) save(Gbern, Ggamma, dat,df.grid, mesh,wReplBern.index ,wReplGamma.index,stkallBern, stkallGamma, spde,file=paste0("results/",year_of_assess,"/INLA/",Rvar,"_INLA_mesh", res_mesh,"km_pred_",res_pred,"_ZAG",year_to_report,"_trajet",trajet,serie,".RData"))
if(exi & overw ==2) stop("Preventing overwite. The model was not saved")

}
