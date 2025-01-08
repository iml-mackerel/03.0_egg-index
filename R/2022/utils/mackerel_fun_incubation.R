# Mackerel egg incubation functions to assist in the stock assessment of Atlantic mackerel (northern contingent)

#### Lockwood, S. J., J. H. Nichols, and S. H. Coombs. 1977. The development rates of mackerel (Scomber scombrusL.) eggs over a range of temperatures. ICES CM 1977/J:13, 13 p 
#' Mackerel egg incubation - Lockwood 1977 - NEA mackerel
#' Bay of Biscay
#' @param t t = temperature in degrees Celsisus. 
#' @return incubation time in hours
#' @export
#' @examples I = I_Lockwood(10)
I_lockwood <- function(t){
  I = exp(-1.614 * log(t) + 7.759)
  I
}

### Mendiola, D., P. Alvarez, U. Cotano, E. Etxebeste, and A. M. de Murguia. 2006. Effects of temperature on development and mortality of Atlantic mackerel fish eggs. Fish. Res. 80:158-168.
#' Mackerel egg incubation - Mendiola et al., 2006 - NEA mackerel
#' for egg stage 1b which is roughly equivalent to our stage 1 (see ICES documents and Girard 2000)
#' @param t temperature in degrees Celsius
#' @return incubation time in hours
#' @export
#' @examples 
#' I_mendiola(10)
I_mendiola <- function(t){
  I = exp(-1.313 * log(t) + 6.902)
  I
}

library(magrittr);library(tidyverse); library(ggthemes)

Temp <- c(seq(5,25,0.5))

Temp <- as.data.frame(Temp)
p1<- Temp %>% mutate(`Lockwood et al. (1977)` = I_lockwood(Temp), # notice that this is the same as I_la3
                `Mendiola et al. (2006)` = I_mendiola(Temp)) %>%
  pivot_longer(cols = 2:3, names_to = "model", values_to = "Incubation") %>%
  ggplot(aes(Temp, Incubation))+
  theme_few()+theme(legend.title=element_blank(), 
                    legend.position=c(0.75, 0.8),
                    legend.background = element_blank())+
  geom_rect(aes(xmin=6.7, xmax=14.9, ymin=Inf, ymax=-Inf), col="black", fill="grey85")+
    geom_line(size = 4, alpha = 0.7, aes(colour = model))+geom_point(aes(colour = model))
 
p1EN <-  p1 + xlab("Temperature (\u00B0C)") +ylab("Incubation time (hours)")
ggsave(paste0("figures/",year_of_assess,"/Incubation_timeEN.png"), width=6, height=4, dpi=600, units="in")  
  
p1FR <-  p1 + xlab("Temp\u00E9rature (\u00B0C)") +ylab("Temps d'incubation (heures)")
ggsave(paste0("figures/",year_of_assess,"/Incubation_timeFR.png"), width=6, height=4, dpi=600, units="in")  

p1BI <-  p1 + xlab("Temp\u00E9rature (\u00B0C)") +ylab("Temps d'incubation (heures)\n Incubation time (hours)")
ggsave(paste0("figures/",year_of_assess,"/Incubation_timeBI.png"), width=6, height=4, dpi=600, units="in")  

#quantile(eggt1$temperature0_10, probs=c(0.025, 0.975))
  