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


#quantile(eggt1$temperature0_10, probs=c(0.025, 0.975))
  