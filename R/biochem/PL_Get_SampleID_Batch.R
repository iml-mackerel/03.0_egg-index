PL_Get_SampleID_Batch <- function(sql_file, output_file) {
  
  # Function for the batch extraction of PL sampleID from Biochem.
  #
  # Input: none
  #
  # Output: None
  #
  # Last update: 20141015
  # Benoit.Casault@dfo-mpo.gc.ca
  
  source("R/biochem/Write_Database_Table.R")
  source("R/biochem/PL_BioChem_Data.R")
  
  # get the data
  PL_Biochem_Data(sql_file, output_file, my.env$host, my.env$port, my.env$sid, my.env$username, my.env$password)
  
  
}


