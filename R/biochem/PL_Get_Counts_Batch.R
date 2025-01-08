PL_Get_Counts_Batch <- function(ifile, output_file, sql_file) {
  
  # Function for the batch extraction of PL counts from Biochem.
  #
  # Input: none
  #
  # Output: None
  #
  # Last update: 20141015
  # Benoit Casault
  # Update: 20160914 by Lehoux, Caroline. Using a more general formulation for ifile and deal with multiple files
  #       + remove quote="" when opening df. filter, changed flag by meters-sqd_flag 
  #       + REG
  #Lehoux C 20160915 removed the filter. To ensure that all data needs to be there, open sampleID file and adjust sql querry accordingly.
  
  source("R/biochem/Write_Database_Table.R")
  source("R/biochem/PL_BioChem_Data.R")

  df <- readRDS(ifile)
  Write_Database_Table(df, "TMP", my.env$host, my.env$port, my.env$sid, my.env$username, my.env$password)
  # sql file
  # get the data
  PL_Biochem_Data(sql_file, output_file, my.env$host, my.env$port, my.env$sid, my.env$username, my.env$password)
    
  return()
}