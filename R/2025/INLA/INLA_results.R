#scripts adapted from the krill acoustic project

source(paste0('R/',year_of_assess,'/INLA/plotSpatialFieldCL.R'))
source(paste0('R/',year_of_assess,'/INLA/plotSmoother.R'))
source(paste0('R/',year_of_assess,'/INLA/model_validation.R'))
source(paste0('R/',year_of_assess,'/INLA/get_prediction_grid.R'))
library(dplyr)
library(tidyr)

INLA_results <- function(res_mesh, res_pred, trans, myfamily="gaussian", minw=-10, maxw=10, smoother=T, Rvar="DEP", year_to_report,year_of_assess, trajet, serie=""){
  
   name = paste0(Rvar, "_INLA_mesh", res_mesh,"km_pred_",res_pred,"_",myfamily,year_to_report,"_tajet",trajet, serie)
 
  #plot random field#
  plotSpatialFieldCL(dir=paste0("results/",year_of_assess,"/INLA/"),name=name ,minw=minw, maxw=maxw, year_of_assess=year_of_assess)
  
  
   #plot smoother#
  if(myfamily=="gaussian"){
  plotSmoother(dir=paste0("results/",year_of_assess,"/INLA/"),name=name, smoother=smoother)
  }
  if(myfamily=="ZAG"){
  plotSmoother(dir=paste0("results/",year_of_assess,"/INLA/"),name=name, subfamily="binomial",smoother=smoother)
  plotSmoother(dir=paste0("results/",year_of_assess,"/INLA/"),name=name, subfamily="gamma",smoother=smoother)
  }
  
  #range#
  if(myfamily=="gaussian"){
  myrange(dir=paste0("results/",year_of_assess,"/INLA/"),name=name, subfamily="gaussian")
  }
  if(myfamily=="ZAG"){
  myrange(dir=paste0("results/",year_of_assess,"/INLA/"),name=name, subfamily="binomial")
  myrange(dir=paste0("results/",year_of_assess,"/INLA/"),name=name, subfamily="gamma")
  }
  
  if(myfamily=="gaussian"){
  validation_gaussian(dir=paste0("results/",year_of_assess,"/INLA/"),trans="1/5",name, 
                      rvarpos=2, varpos=5)
  }
  if(myfamily=="ZAG"){
    validation_ZAG(dir=paste0("results/",year_of_assess,"/INLA/"),name=name,
                        rvarpos=2, varpos=5)
  }
   
  #extract predictions#
  get_prediction_grid(dir=paste0("results/",year_of_assess,"/INLA/"), dir.out=paste("results/",year_of_assess,"/INLA/predictions/"),name=name, trans=trans)
  
}