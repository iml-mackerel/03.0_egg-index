plotSmoother <- function(dir,name,subfamily="gaussian", smoother){

  load(paste0(dir,name,".Rdata"))
  dir.fig <-  gsub(dir, pattern="results", replacement="img")
  
   if(subfamily=="gaussian") mod=mod
   if(subfamily=="binomial") mod=Gbern
   if(subfamily=="gamma") mod=Ggamma
  if(subfamily=="tw") mod=Gtw
  
   
   if(smoother){
   
  #get automatically the number of smoother
  smoothname <- gsub(rownames(mod$summary.lincomb.derived), pattern="lc", replacement="")
  smoothname <- unique(gsub(smoothname, pattern='[0-9]+', replacement=""))
  Nsmoother=length(smoothname)
  
 if(Nsmoother >3) stop("cette fonction ne permet pas plus de 3 smoothers pour l'instant")
  
Ns <- nrow(dat)
if(Nsmoother>=1){
f.smooth1    <- mod$summary.lincomb.derived[1:Ns + 0 * Ns, "mean"] 
SeLo.smooth1 <- mod$summary.lincomb.derived[1:Ns + 0 * Ns,"0.025quant"] 
SeUp.smooth1 <- mod$summary.lincomb.derived[1:Ns + 0 * Ns,"0.975quant"]
Ismooth1 <- order(dat[,smoothname[1]])

}

if(Nsmoother >=2){
f.smooth2    <- mod$summary.lincomb.derived[1:Ns + 1 * Ns, "mean"] 
SeLo.smooth2 <- mod$summary.lincomb.derived[1:Ns + 1 * Ns,"0.025quant"] 
SeUp.smooth2 <- mod$summary.lincomb.derived[1:Ns + 1 * Ns,"0.975quant"]
Ismooth2 <- order(dat[,smoothname[2]])
}

if(Nsmoother >=3){
  f.smooth3    <- mod$summary.lincomb.derived[1:Ns + 2 * Ns, "mean"] 
  SeLo.smooth3 <- mod$summary.lincomb.derived[1:Ns + 2 * Ns,"0.025quant"] 
  SeUp.smooth3 <- mod$summary.lincomb.derived[1:Ns + 2 * Ns,"0.975quant"]
  Ismooth3 <- order(dat[,smoothname[3]])
}



if(Nsmoother==1){
  Mydat <- data.frame(
    mu   = c(f.smooth1[Ismooth1]), 
    SeUp = c(SeUp.smooth1[Ismooth1]), 
    SeLo = c(SeLo.smooth1[Ismooth1]), 
    Xaxis = c(sort(dat[,smoothname[1]])),
    ID    = factor(rep(c(smoothname[1]), each = nrow(dat))))
  my.ggp.yrange <- c(0)
  }

if(Nsmoother==2){
Mydat <- data.frame(
  mu   = c(f.smooth1[Ismooth1],    f.smooth2[Ismooth2]), 
  SeUp = c(SeUp.smooth1[Ismooth1], SeUp.smooth2[Ismooth2]), 
  SeLo = c(SeLo.smooth1[Ismooth1], SeLo.smooth2[Ismooth2]), 
  Xaxis = c(sort(dat[,smoothname[1]]), sort(dat[,smoothname[2]])),
  ID    = factor(rep(c(smoothname[1],smoothname[2]), each = nrow(dat))))
my.ggp.yrange <- c(0, 0)
}

if(Nsmoother==3){
  Mydat <- data.frame(
    mu   = c(f.smooth1[Ismooth1],    f.smooth2[Ismooth2],  f.smooth3[Ismooth3]), 
    SeUp = c(SeUp.smooth1[Ismooth1], SeUp.smooth2[Ismooth2] , SeUp.smooth3[Ismooth3]), 
    SeLo = c(SeLo.smooth1[Ismooth1], SeLo.smooth2[Ismooth2], SeLo.smooth3[Ismooth3]), 
    Xaxis = c(sort(dat[,smoothname[1]]), sort(dat[,smoothname[2]]), sort(dat[,smoothname[3]])),
    ID    = factor(rep(c(smoothname[1],smoothname[2],smoothname[3]), each = nrow(dat))))
  my.ggp.yrange <- c(0, 0,0)
  }
Mydat$Y <- 0



# Figure 21.11 and for Figure 21.12:  add 3x exp 
p <- ggplot(data = Mydat) + xlab("") + ylab("Partial effect")+ basetheme+
       geom_line(aes(x = Xaxis, y = mu))+
       geom_ribbon(data = Mydat, aes(x = Xaxis, ymax = SeUp, ymin = SeLo),alpha = 0.6)+
       facet_wrap(~ID, scales = "free", ncol = 3,strip.position="bottom") +
       geom_text(aes(y = Y,x = Xaxis, label = "|"), size = 1) +theme(strip.placement = "outside")


png(paste0(dir.fig,"/smoother_",name,subfamily,".png"), width=3*Nsmoother, height=4, res=600, units="in")
print(p)
dev.off()
}

fixed_effect <- mod$summary.fixed[,c("mean","sd", "0.025quant", "0.975quant")]
#plot is not exported.
fixed_effect %>%  mutate(cof= row.names(fixed_effect)) %>%  ggplot(aes(x=cof, y=mean))+geom_point()+geom_errorbar(aes(ymin=`0.025quant`, ymax=`0.975quant`))

write.table(fixed_effect, file=paste0(dir,"/fixed_effects_",name,subfamily, ".txt"), row.names=T, sep="\t", dec=".")

#save(p, Mydat, file=paste0(dir.fig,"smoother_",name,subfamily,".Rdata"))
#takes too much place.
 
}