
#https://stats.stackexchange.com/questions/231074/confidence-intervals-on-predictions-for-a-non-linear-mixed-model-nlme
#see cited paper at the end of the page

nlme_boot <-function(fitted, data, newframe){   
  
  pp <- predict(fitted,levels=1)
  rr <- residuals(fitted, type="pearson")
  dd <- data.frame(data,pred=pp,res=rr)
  
   
      bsamp2<- dd %>% #dplyr::group_by(year) %>%
        dplyr::mutate( gsi= pred + sample(res, replace=T)) %>%  dplyr::arrange(year) #changed! my version
    #bsamp2$gsi[bsamp2$gsi<0] <- 0
    fm=tryCatch(
      {update(fitted,data=bsamp2)},
      error = function(e){return(NULL)}
    )
    
    if(!is.null(fm)){
      out<- predict(fm,newdata=newframe,level=1) #changed to 1
      frameout<- bind_cols(doy=newframe$doy, year=newframe$year, out=out)
      
      frameout<- as_tibble(frameout %>% group_by(year) %>%  mutate(slope = c(diff(out, lag = 1), NA),
                                                         prob = slope / sum(slope, na.rm = T)))
    }
    
    if(is.null(fm)) {
      out <- rep(NA, nrow(newframe)) 
      slope <- rep(NA, nrow(newframe)) 
      prob <- rep(NA, nrow(newframe)) 
      frameout<- as_tibble(bind_cols(doy=newframe$doy, year=newframe$year, out=out,slope=slope, prob=prob))
      
    }
    
    return(list(frameout$out, frameout$prob))
 
}      
  
 






