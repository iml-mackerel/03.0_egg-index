PL_Abs2SqMeter <- function(value, start_depth, end_depth, split_fraction, volume) {
  
  # Function for converting absolute quantities to per unit are basis based on
  # sample parameters.
  #
  # Input: value : variable to convert - numeric array
  #        start_depth : tow start depth - numeric array
  #        end_depth : tow end depth - numeric array
  #        split_fraction : sample split fraction - numeric array
  #        volume : tow filtered volume - numeric array
  #
  # Output: var_out : converted variable (e.g in #/m2) - numeric array
  #
  # Last update: 20151123
  # Benoit.Casault@dfo-mpo.gc.ca
  
  return(value*abs(start_depth-end_depth)/split_fraction/volume)
}