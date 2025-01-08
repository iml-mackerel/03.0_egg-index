library(INLA)
library(inlabru)
library(gstat)
library(ggpmisc)

validation_bernouilli <- function(Data, mod, varname){
N=nrow(Data)
Pi <- mod$summary.fitted.values[1:N,"mean"]

Data_valid <- cbind(plotID= seq(1,dim(Data)[1],1),Observed=Data[,varname],Predicted=Pi)

Tresh_opt <- optimal.thresholds(Data_valid,na.rm=TRUE,opt.methods=3)[1,2]
CMX1 <- cmx(Data_valid, na.rm=TRUE, threshold=Tresh_opt)
sens1 <- PresenceAbsence::sensitivity(CMX1, st.dev=FALSE)
spec1 <- PresenceAbsence::specificity(CMX1, st.dev=FALSE)
TSS <- round((sens1+spec1)-1 ,3)
AUC<-round(roc.area(Data_valid[,2], Pi)$A,3)
res<-as.vector(c(paste0("TSS=",TSS),paste0("AUC=", AUC)))

return(res)
}


myrange <- function(dir, name, subfamily="gaussian"){
  load(paste0(dir,name,".Rdata"))
  if(subfamily=="gaussian"){
  SpFi.w <- inla.spde2.result(inla = mod,
                              name = "w",
                              spde = spde,
                              do.transfer = TRUE)
  }
  if(subfamily=="binomial"){
    SpFi.w <- inla.spde2.result(inla = Gbern,
                                name = "w",
                                spde = spde,
                                do.transfer = TRUE)
  }
  if(subfamily=="gamma"){
    SpFi.w <- inla.spde2.result(inla = Ggamma,
                                name = "wGamma",
                                spde = spde,
                                do.transfer = TRUE)
  }  
  
  if(subfamily=="tw"){
    SpFi.w <- inla.spde2.result(inla = Gtw,
                                name = "w",
                                spde = spde,
                                do.transfer = TRUE)
  }  
  Kappa <- inla.emarginal(function(x) x, 
                          SpFi.w$marginals.kappa[[1]] )
  
  sigmau <- inla.emarginal(function(x) sqrt(x), 
                           SpFi.w$marginals.variance.nominal[[1]] )
  
  r <- inla.emarginal(function(x) x, 
                      SpFi.w$marginals.range.nominal[[1]] )
  
  n.r<-as.vector(c(mean=inla.emarginal(function(x) x, SpFi.w$marginals.range[[1]]), 
                   q=inla.hpdmarginal(0.95, SpFi.w$marginals.range[[1]]))[c(1,2,3)])
  
  
  res_range<-as.vector(c(paste0("Kappa=",signif(Kappa,3)), paste0("sigmau=",round(sigmau,3)), paste0("range (km)=",round(r/1000,3)), paste0("erreur around range (km)=",round(n.r[2]/1000,3),"-" ,round(n.r[3]/1000,3))))
  return(res_range)
  }  


validation_gamma <- function(Data, mod, varname, name){
  ####extract residuals Gamma####
  N     <- nrow(Data)
  mu1   <- mod$summary.fitted.values[1:N,"mean"] 
  r1    <- mod$summary.hyperpar[1, "mean"]
  VarY1 <- mu1^2 / r1 
  E1    <- (Data[, varname] - mu1) / sqrt(VarY1)
  
  
  Data$mu1 <- mu1
  Data$E1 <- E1
  Data.Pos <- Data[which(Data[,varname]>0),]
  
  # Plot fitted values versus observed data
  png(paste0("valid_", name,".png"), width=8, height=4, res=600, units="in")
  par(mfrow = c(1,2), mar = c(5,5,2,2), cex.lab = 1.5)
  plot(y = Data.Pos[,varname],
       x = Data.Pos$mu1,
       main = "",
       xlab = "Fitted values",
       ylab = "Observed data",
       xlim = base::range(Data.Pos[,varname]),
       ylim = base::range(Data.Pos[,varname]))
  abline(a=0, b=1, col="red")
  text(mean(base::range(Data.Pos[,varname])),base::range(Data.Pos[,varname])[2]/10,paste0("corr= ",round(cor(Data.Pos[,varname], Data.Pos$mu1),2)),  cex = 1)
  
  plot(y = log1p(Data.Pos[,varname]),
       x = log1p(Data.Pos$mu1),
       main = "",
       xlab = "Fitted values",
       ylab = "Observed data",
       xlim = base::range(log1p(Data.Pos[,varname])),
       ylim = base::range(log1p(Data.Pos[,varname])))
  abline(a=0, b=1, col="red")
  dev.off()
  
  
    # Check the Pearson residuals
  # any spatial dependency using a
  # variogram.
  Data.Pos$Xkm <- Data.Pos$X.utm / 1000
  Data.Pos$Ykm <- Data.Pos$Y.utm / 1000
  
  MyData1 <- data.frame(E1 = Data.Pos$E1, 
                        Xkm = Data.Pos$Xkm, 
                        Ykm = Data.Pos$Ykm)
  coordinates(MyData1) <- c("Xkm", "Ykm")
  V1 <- variogram(E1 ~ 1, 
                  MyData1, 
                  cressie = TRUE)
  plot(V1)
}


validation_ZAG <-  function(dir, name, rvarpos, varpos){
  load(paste0(dir,name,".Rdata"))
   
  dat<- dat %>%  st_drop_geometry()
  
  rvarname= colnames(dat)[rvarpos]
  varnames= colnames(dat)[varpos]
  
  
  #* Section 6.1: ZAG mean and variance ----
  #' Following the expressions for the mean and the variance of the ZAG,
  #' we need to calculate:
  #'   E[Biomass]   = Pi * mu
  #'   var[Biomass] = ( (Pi * r + Pi - Pi *r) / r) * mu^2
  #id.fit <- inla.stack.index(stkall.Repl, 'Fit')$dat

  #' Get the Pi from the Bernoulli model
  N  <- nrow(dat)
  Pi <- Gbern$summary.fitted.values[1:N,"mean"]
  
  #' Get the mu and r from the Gamma model
  mu <- Ggamma$summary.fitted.values[1:N, "mean"]
  r  <- Ggamma$summary.hyperpar[1,"mean"]
  
  #varB <- Pi *(1-Pi)
  #varG <- mu^2/r
  
  #' Calculate the ZAG mean and variance
  ExpY <- Pi * mu
 
  varY <- ( (Pi * r + Pi - Pi *r) / r) * mu^2
  #* Section 6.2: ZAG Pearson residuals ----
  #if(any(varY <0)) Ezag <- (dat[,rvarpos] - ExpY)
   Ezag <- (dat[,rvarpos] - ExpY)  / sqrt(varY)

  #' A Gamma GLM cannot be overdispersed. Therefore, a ZAG cannot be overdispersed.
  my.formula <- y ~ x
 
  
  g1<-ggplot(dat, aes(x=ExpY , y=Ezag))+ ggtitle("")+
    geom_point(col="grey35", size=0.5)+
    #geom_smooth() +
    basetheme + geom_hline(yintercept=0, col="red", lty=2)
  g2<-ggplot(dat, aes(y=ExpY, x=get(rvarname))) + 
    geom_point(col="grey35", size=0.5)+
    geom_smooth(method="lm")+
    geom_abline(slope=1, intercept=0, col="red", lty=2)+
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                 parse = TRUE, vstep=0.08)+
        geom_abline(slope=1, intercept=0, col="red", lty=2)+
    basetheme
  
  g1EN <-  g1 +scale_x_continuous(name="Fitted values") + ylab("Pearson's residuals")
  g1FR <-  g1 +scale_x_continuous(name="Valeurs calibr\u00E9es" , trans="log10p1") + ylab("R\u00E9sidus de Pearson")
  
  g2EN <-  g2 +scale_y_continuous(name="Fitted values" , trans="log10p1", breaks=c(0,10,30,100,300,1000,2000)) + scale_x_continuous(name="Calculated values", trans="log10p1", breaks=c(0,10,30,100,300,1000,2000))
  g2FR <-  g2 +scale_y_continuous(name="Valeurs calibr\u00E9es" , trans="log10p1", breaks=c(0,10,30,100,300,1000,2000)) + scale_x_continuous(name="Valeurs calcul\u00E9es", trans="log10p1", breaks=c(0,10,30,100,300,1000,2000))
  g2BI <-  g2 +scale_y_continuous(name="Valeurs calibr\u00E9es / Fitted values" , trans="log10p1", breaks=c(0,10,30,100,300,1000,2000)) + 
               scale_x_continuous(name="Valeurs calcul\u00E9es / Calculated values", trans="log10p1", breaks=c(0,10,30,100,300,1000,2000))
  
  dir.img <- gsub(dir, pattern="results", replacement="img")
  
  dat  %>% dplyr::mutate(ExpY=ExpY) %>% ungroup() %>%  dplyr::filter(PA ==0) %>% ggplot(aes(x=ExpY))+geom_histogram() + basetheme +
    ylab("N") + xlab("PJO prédites pour PJO calculé de 0\nFitted DEP for calculated DEP of 0"
    )
  ggsave(paste0(dir.img,"/models/histogram_P.png"), dpi=600, units="in", width=5, height=4,)

  #ggarrange(ggarrange(plotlist = all_plots, nrow=1, ncol=3, widths=c(0.2,0.2,0.4)), ggarrange(v), nrow=2, ncol=1) 
 #  ggsave(paste0("figures/INLA/models/model_validation_",name,".png"), dpi=600, units="in", width=10, height=6)
  pEN<- ggarrange(g1EN ,g2EN, nrow=1, ncol=2, align="hv")
  
  ggsave(paste0(dir.img,"/models/model_validation_",name,"EN.png"), dpi=600, units="in", width=9, height=4)
  pFR <- ggarrange(g1FR ,g2FR, nrow=1, ncol=2, align="hv")
  ggsave(paste0(dir.img,"/models/model_validation_",name,"FR.png"), dpi=600, units="in", width=9, height=4)
   
  ggsave(paste0(dir.img,"/models/model_validation_",name,"BI.png"), dpi=600, units="in", width=6, height=4, plot = g2BI)
  
  
 return(pEN) 
}




validation_tw <-  function(dir, name, rvar){
  load(paste0(dir,name,".Rdata"))
  
  dat<- dat %>%  st_drop_geometry()
  
  rvarname= grep(colnames(dat), pattern=rvar, value=T)
   
  idb.prd <- inla.stack.index(stkalltw, "twFit")$data
  # sd.prd   <- I$summary.fitted.values$sd[id.prd]
  ExpY<- Gtw$summary.fitted.values$mean[idb.prd]
  N  <- nrow(dat)
  
  p.power  <- Gtw$summary.hyperpar[1,"mean"]
  sigma  <- Gtw$summary.hyperpar[2,"mean"]
  
  tweedie_var <- sigma^2*(ExpY)*p.power

  Ezag <- (dat[,rvar] - ExpY)  / sqrt(tweedie_var)
 
  my.formula <- y ~ x
  
  
  
  g1<-ggplot(dat, aes(x=ExpY , y=Ezag))+ ggtitle("")+
    geom_point(col="grey35", size=0.5)+
    #geom_smooth() +
    basetheme + geom_hline(yintercept=0, col="red", lty=2)
  g2<-ggplot(dat, aes(y=ExpY, x=get(rvarname))) + 
    geom_point(col="grey35", size=0.5)+
    geom_smooth(method="lm")+
    geom_abline(slope=1, intercept=0, col="red", lty=2)+
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                 parse = TRUE, vstep=0.08)+
    geom_abline(slope=1, intercept=0, col="red", lty=2)+
    basetheme
  
  g1EN <-  g1 +scale_x_continuous(name="Fitted values") + ylab("Pearson's residuals")
  g1FR <-  g1 +scale_x_continuous(name="Valeurs calibr\u00E9es") + ylab("R\u00E9sidus de Pearson")
  g1BI <-  g1 +scale_x_continuous(name="Valeurs calibr\u00E9es / Fitted values") + ylab("R\u00E9sidus de Pearson / Pearson's residuals")
  
  g2EN <-  g2 +scale_y_continuous(name="Fitted values" , trans="log10p1", breaks=c(0,10,30,100,300,1000,2000)) + scale_x_continuous(name="Observed values", trans="log10p1", breaks=c(0,10,30,100,300,1000,2000))
  g2FR <-  g2 +scale_y_continuous(name="Valeurs calibr\u00E9es" , trans="log10p1", breaks=c(0,10,30,100,300,1000,2000)) + scale_x_continuous(name="Valeurs observ\u00E9es", trans="log10p1", breaks=c(0,10,30,100,300,1000,2000))
  g2BI <-  g2 +scale_y_continuous(name="Valeurs calibr\u00E9es / Fitted values" , trans="log10p1", breaks=c(0,10,30,100,300,1000,2000)) + scale_x_continuous(name="Valeurs observ\u00E9es / Observed values", trans="log10p1", breaks=c(0,10,30,100,300,1000,2000))
  
  #ggarrange(ggarrange(plotlist = all_plots, nrow=1, ncol=3, widths=c(0.2,0.2,0.4)), ggarrange(v), nrow=2, ncol=1) 
  #  ggsave(paste0("figures/INLA/models/model_validation_",name,".png"), dpi=600, units="in", width=10, height=6)
  pEN<- ggarrange(g1EN ,g2EN, nrow=1, ncol=2, align="hv")
  dir.img <- gsub(dir, pattern="results", replacement="img")
    ggsave(paste0(dir.img,"/models/model_validation_",name,"EN.png"), dpi=600, units="in", width=9, height=4)
  
  pFR <- ggarrange(g1FR ,g2FR, nrow=1, ncol=2, align="hv")
  ggsave(paste0(dir.img,"/models/model_validation_",name,"FR.png"), dpi=600, units="in", width=9, height=4)
  
  pBI <- ggarrange(g1BI ,g2BI, nrow=1, ncol=2, align="hv")
  ggsave(paste0(dir.img,"/models/model_validation_",name,"BI.png"), dpi=600, units="in", width=9, height=4)
  
  
  return(pEN) 
}



validation_gaussian <- function(dir,name, trans="log", rvarpos, varpos, res_agg){
  library(gstat)
  load(paste0(dir,name,".Rdata"))
  id.fit <- inla.stack.index(stkall.Repl, 'Fit')$dat
  id.prd <- inla.stack.index(stkall.Repl, 'Pred')$dat
  
  rvarname= colnames(dat)[rvarpos]
  varnames= colnames(dat)[varpos]
  
  if(trans=="1/5") inverse_trans <-  function(x) x^5
  if(trans=="1/5") transfct <-  function(x) x^(1/5)
  
  
  ####extract residuals####
  N     <- nrow(dat)
  tau <- mod$marginals.hyperpar$`Precision for the Gaussian observations`
  MySqrt <- function(x) { 1 / sqrt(x) }
  sigma <- inla.emarginal(MySqrt, tau)
  fit   <- mod$summary.fitted.values[id.fit,"mean"] 
  
  
 
  VarY1 <- sigma^2
 
  E1    <- (transfct(dat[, rvarname]) - fit) / sqrt(VarY1)
  
  dat$mu <- inverse_trans(fit)
  dat$E1 <- E1
  
 
  g1<-ggplot(dat, aes(x=transfct(mu), y=E1))+ ggtitle("")+
    geom_point()+
    geom_smooth() +basetheme
  g2<-ggplot(dat, aes(x=mu, y=get(rvarname))) + 
      geom_point()+
      geom_smooth(method="lm")+
      xlab("Fitted values")+
      ylab("Observed data")+
      geom_abline(slope=1, intercept=0, col="red", lty=2)+
      ggtitle(paste0("corr= ",round(cor(dat[,rvarname], dat$mu),2)))+basetheme
  
  
  datvar<- dat %>%  dplyr::select(E1, varnames)
  datvar<- datvar %>%  pivot_longer(2:ncol(datvar))
  
    v<- ggplot(datvar, aes(x=value, y=E1))+ xlab("")+
      geom_point()+
      geom_smooth() +
      facet_wrap(~ name, nrow=1,strip.position="bottom", scales="free_x") + basetheme +theme(strip.placement = "outside")
  
  
 test<-as.data.frame(cbind(dat[,rvarname],dat$mu)); colnames(test)<- c("Observations", "Fitted Values")
 testg<-gather(test, key="Type", value="ResponseVar")
 
  h1<-ggplot(testg, aes(x=ResponseVar, fill=Type))+
    geom_histogram(bins=200, position="identity", alpha=0.5)+
   basetheme + theme(legend.position="top")
 
  # Check the Pearson residuals
  # any spatial dependency using a
  # variogram.
  
  dat$Xkm <- dat$X.m / 1000
  dat$Ykm <- dat$Y.m / 1000
  
  Mydat1 <- data.frame(E1 = dat$E1, 
                        Xkm = dat$Xkm, 
                        Ykm = dat$Ykm)
  coordinates(Mydat1) <- c("Xkm", "Ykm")
  V1 <- variogram(E1 ~ 1, 
                  Mydat1, 
                  cressie = TRUE)
  
  pv1<-plot(V1)

  all_plots<- list(g1,g2,h1)
  
  ggarrange(ggarrange(plotlist = all_plots, nrow=1, ncol=3, widths=c(0.2,0.2,0.4)), ggarrange(v), nrow=2, ncol=1) 
  ggsave(paste0("figures/INLA/models/model_validation_",name,".png"), dpi=600, units="in", width=10, height=6)
  
  
  }



