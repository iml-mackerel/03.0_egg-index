library(INLA)
library(inlabru)
library(gstat)
library(ggpmisc)

validation_bernouilli <- function(Data, mod, varname){
library(PresenceAbsence)
library(verification)
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
  
  Kappa <- inla.emarginal(function(x) x, 
                          SpFi.w$marginals.kappa[[1]] )
  
  sigmau <- inla.emarginal(function(x) sqrt(x), 
                           SpFi.w$marginals.variance.nominal[[1]] )
  
  r <- inla.emarginal(function(x) x, 
                      SpFi.w$marginals.range.nominal[[1]] )
  
  n.r<-as.vector(c(mean=inla.emarginal(function(x) x, SpFi.w$marginals.range[[1]]), 
                   q=inla.hpdmarginal(0.95, SpFi.w$marginals.range[[1]]))[c(1,2,3)])
  
  
  res_range<-as.vector(c(paste0("Kappa=",signif(Kappa,3)), paste0("sigmau=",round(sigmau,3)), paste0("range (km)=",round(r/1000,3)), paste0("erreur around range (km)=",round(n.r[2]/1000,3),"-" ,round(n.r[3]/1000,3))))
  write.table(as.data.frame(res_range), file=paste0(dir,"range",name,subfamily, ".txt"), row.names=T, sep="\t", dec=".")
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
   
  rvarname= colnames(dat)[rvarpos]
  varnames= colnames(dat)[varpos]
  
  
  #* Section 6.1: ZAG mean and variance ----
  #' Following the expressions for the mean and the variance of the ZAG,
  #' we need to calculate:
  #'   E[Biomass]   = Pi * mu
  #'   var[Biomass] = ( (Pi * r + Pi - Pi *r) / r) * mu^2
  #id.fit <- inla.stack.index(stkall.Repl, 'Fit')$dat
  library(boot)
  #' Get the Pi from the Bernoulli model
  N  <- nrow(dat)
  Pi <- inv.logit(Gbern$summary.fitted.values[1:N,"mean"]) 
  
  #' Get the mu and r from the Gamma model
  N  <- nrow(dat)
  mu <- Ggamma$summary.fitted.values[1:N, "mean"]
  r  <- Ggamma$summary.hyperpar[1,"mean"]
  
  
  #' Calculate the ZAG mean and variance
  ExpY <- Pi * mu
  varY <- ((Pi * r + pi - pi * r) / r) * mu^2
  
  #varY <- ((Pi * r) + (pi - (pi * r)) / r) * mu^2
  
  
  #* Section 6.2: ZAG Pearson residuals ----
  if(any(varY <0)) Ezag <- (dat[,rvarpos] - ExpY)
  if(!any(varY <0)) Ezag <- (dat[,rvarpos] - ExpY)  / sqrt(varY)
  
  #' A Gamma GLM cannot be overdispersed. Therefore, a ZAG cannot be overdispersed.
  
  my.formula <- y ~ x
  
  g1<-ggplot(dat, aes(x=mu , y=Ezag))+ ggtitle("")+
    geom_point(col="grey35", size=0.5)+
    #geom_smooth() +
    basetheme + geom_hline(yintercept=0, col="red", lty=2)
  g2<-ggplot(dat, aes(y=mu +1, x=get(rvarname) +1)) + 
    geom_point(col="grey35", size=0.5)+
    geom_smooth(method="lm")+
    geom_abline(slope=1, intercept=0, col="red", lty=2)+
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
                 parse = TRUE, vstep=0.08)+
        geom_abline(slope=1, intercept=0, col="red", lty=2)+
    basetheme
  
  g1EN <-  g1 +scale_x_continuous(name="Fitted values" ) + ylab("Pearson's residuals")
  g1FR <-  g1 +scale_x_continuous(name="Valeurs calibr\u00E9es" ) + ylab("R\u00E9sidus de Pearson")
  
  g2EN <-  g2 +scale_y_continuous(name="Fitted values" , trans="log10") + scale_x_continuous(name="Observed values", trans="log10")
  g2FR <-  g2 +scale_y_continuous(name="Valeurs calibr\u00E9es" , trans="log10") + scale_x_continuous(name="Valeurs observ\u00E9es", trans="log10")
  
  
  
  
  datvar<- dat %>% mutate(Ezag=Ezag) %>%  dplyr::select(Ezag, varnames)
  datvar<- datvar %>%  pivot_longer(2:ncol(datvar))
  
  v<- ggplot(datvar, aes(x=value, y=Ezag))+ xlab("")+
    geom_point()+
    geom_smooth() +
    facet_wrap(~ name, nrow=1,strip.position="bottom", scales="free_x") + basetheme +theme(strip.placement = "outside")
  
  
  test<-as.data.frame(cbind(dat[,rvarname],mu)); colnames(test)<- c("Observations", "Fitted Values")
  testg<-gather(test, key="Type", value="ResponseVar")
  
  h1<-ggplot(testg, aes(x=ResponseVar, fill=Type))+
    geom_histogram(bins=200, position="identity", alpha=0.5)+
    basetheme + theme(legend.position="top")
  
  # Check the Pearson residuals
  # any spatial dependency using a
  # variogram.
  
  dat$Xkm <- dat$X.m / 1000
  dat$Ykm <- dat$Y.m / 1000
  
  Mydat1 <- na.omit(data.frame(Ezag = Ezag, 
                       Xkm = dat$Xkm, 
                       Ykm = dat$Ykm))
  coordinates(Mydat1) <- c("Xkm", "Ykm")
  V1 <- variogram(Ezag ~ 1, 
                  Mydat1, 
                  cressie = TRUE)
  
  pv1<-plot(V1)
  
  all_plots<- list(g1,g2,h1)
  
  #ggarrange(ggarrange(plotlist = all_plots, nrow=1, ncol=3, widths=c(0.2,0.2,0.4)), ggarrange(v), nrow=2, ncol=1) 
 #  ggsave(paste0("figures/INLA/models/model_validation_",name,".png"), dpi=600, units="in", width=10, height=6)
  ggarrange(g1EN ,g2EN, nrow=1, ncol=2, align="hv")
  ggsave(paste0("figures/INLA/models/model_validation_",name,"EN.png"), dpi=600, units="in", width=9, height=4)
  ggarrange(g1FR ,g2FR, nrow=1, ncol=2, align="hv")
  ggsave(paste0("figures/INLA/models/model_validation_",name,"FR.png"), dpi=600, units="in", width=9, height=4)
  
  
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



