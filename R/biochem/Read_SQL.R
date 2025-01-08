Read_SQL <- function(sql_file) {
  
  # Function to read a sql query from a text file.
  #
  # Input: sql_file : sql filename - string
  #
  # Output: sql : sql text string - string
  #
  # Last update: 20141015
  # Benoit.Casault@dfo-mpo.gc.ca
  
  # read sql query
  sql <- scan(file=sql_file, what=ls(), sep="\n")
  sql <- paste(sql, collapse=" ")

  return(sql)
}