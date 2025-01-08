##### my packages ################################################################################
## CRAN
cran.packages <- c('tidyverse','boot','magrittr','ggpmisc','ggpubr','ggthemes','mgcv',
                   'fields', 'inlabru', 'sf', 'PresenceAbsence', 'verification', 'raster',
                   'scales', 'nlme','nls.multstart', 'stringr', 'ggforce',"readxl", "marmap")



install.this <- cran.packages[!(cran.packages %in% utils::installed.packages()[,"Package"])]
if(length(install.this)>=1) install.packages(install.this)
dummy <- lapply(cran.packages, require, character.only = TRUE)

## github
git.packages <- c('catchR','DFOdata','CCAM', 'INLA')
install.this <- git.packages[!(git.packages %in% utils::installed.packages()[,"Package"])]
if('catchR' %in% install.this)  devtools::install_github("iml-assess/catchR@eli_parallel")
if('DFOdata' %in% install.this)  devtools::install_github("iml-assess/DFOdata")
if('CCAM' %in% install.this)  devtools::install_github("elisvb/CCAM")
if('INLA' %in% install.this)install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

dummy <- lapply(git.packages, require, character.only = TRUE)

##### source R directory  ############################################################################
#invisible(sapply(list.files(pattern="[.]R$", path="R/", full.names=TRUE), source))

##### my ggplot theme ################################################################################
theme_set(theme_mackerel())             # theme_mackerel from catchR
update_geom_defaults("line", list(size = 1))  # no idea why somethimes I get fat lines otherwise

##### passwords databases #############################################################################
source("../../bdOracle.R")

source(paste0("utils/basemap.R")) #if error "plot new has not been call", restart R, package compatibility issues
source(paste0("utils/mackerel_fun_incubation.R")) # Mackerel incubation .
source(paste0("utils/spatial_projections.R"))#
source(paste0("utils/extract_biochem.R"))
source(paste0("utils/extract_T0_10.R"))
source(paste0("INLA/Mesh.R"))
source(paste0("INLA/INLA_ZAG_covar_Stations.R"))
source(paste0("INLA/INLA_ZAG_covar_Stations_CV.R"))
source(paste0("utils/nlme_boot.R"))
source(paste0("INLA/INLA_tw_covar_Stations.R"))
source(paste0("INLA/INLA_tw_covar_Stations_CV.R"))



source(paste0('INLA/plotSpatialFieldCL.R'))
source(paste0('INLA/plotSmoother.R'))
source(paste0('INLA/model_validation.R'))
source(paste0('INLA/get_prediction_grid.R'))
source(paste0('../biochem/PL_Get_SampleID_Batch.R'))
source(paste0('../biochem/PL_Get_Counts_Batch.R'))
source(paste0('../biochem/PL_Read_Filter.R'))
source(paste0('../biochem/PL_Taxonomic_Grouping.R'))


log10p1_trans = function() scales::trans_new("log10p1", transform=function(x) log10(x+1), inverse=function(x) (10^x)-1)#inverse function is necessary for legend

#source(paste0("R/",year_to_report,"/INLA/getvar.R")) # needs to be retought
