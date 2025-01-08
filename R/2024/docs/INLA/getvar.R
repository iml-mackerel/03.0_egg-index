#basé sur : 
#x <- rnorm(100)
#g <- gl(10,10)
#mns <- tapply(x, g, mean)
#vs <- tapply(x, g, var)
#9/99*(sum(vs) + 10*var(mns))


getvar<- function(pred.data){
#total if using all predicted values
vs <- tapply(pred.data$Fit, pred.data$year, var)
vs<- as.data.frame(vs)  
vs<- bind_cols(vs,year=as.numeric(row.names(vs)))

nl_nt<- (length(unique(pred.data$station))-1)/(nrow(pred.data)-1)

pred.out1<- full_join(pred.data %>%  group_by(year)  %>%  summarize(partial_var=sum(var)), vs)%>%  
       group_by(year) %>% summarize(vart = (nl_nt* (partial_var +  length(unique(pred.data$station))*vs))) %>%  
       mutate(sdt=sqrt(vart),
              set=sdt/sqrt(66))

#partial if using only values that were estimated
vs <- tapply(pred.data$Fit, pred.data$year, var)
vs<- as.data.frame(vs)  
vs<- bind_cols(vs,year=as.numeric(row.names(vs)))

nl_nt<- (length(unique(pred.data$station))-1)/(nrow(pred.data)-1)

pred.out2<- full_join(pred.data %>%  group_by(year)  %>%  summarize(partial_var=sum(var1)), vs)%>%  #variance is 0 for not estimated
  group_by(year) %>% summarize(varp = (nl_nt* (partial_var +  length(unique(pred.data$station))*vs))) %>%  
  mutate(sdp=sqrt(varp),
         sep=sdp/sqrt(66))
pred.out<- full_join(pred.out1, pred.out2)

return(pred.out)
}
