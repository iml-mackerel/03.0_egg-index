Run_Database_Query <- function(sql_file, host, port, sid, username, password) {
  
  # Function to extract data from a database (e.g. Biochem).  The function reads a sql query, opens a
  # connection to a database, extracts the data and returns the extracted data and the associated
  # sql code.
  #
  # Input: sql_file : sql filename - string
  #        data_source : data source name for database - string
  #        username : username - string
  #        password: user password - string
  #
  # Output: data : list of two elements
  #         data[[1]] : sql code for the query - string
  #         data[[2]] : data extracted from the database - dataframe
  #
  # Last update: 20141015
  # Benoit.Casault@dfo-mpo.gc.ca
  
  # load odbc library
  library(ROracle)
  
  # source custom functions
  source("R/biochem/Read_SQL.R")
  
  # declare empty list to store outputs
  data <- list()
  
  # read sql file
  data[[1]] <- Read_SQL(sql_file)
  
  # create an Oracle Database instance and create connection
  drv <- dbDriver("Oracle")
  
  # connect string specifications
  connect.string <- paste(
    "(DESCRIPTION=",
    "(ADDRESS=(PROTOCOL=tcp)(HOST=", host, ")(PORT=", port, "))",
    "(CONNECT_DATA=(SID=", sid, ")))", sep = "")

  # use username/password authentication
  conn <- ROracle::dbConnect(drv, username=username, password=password, dbname = connect.string, encoding="UTF-8")
  #"windows-1252"
  # run SQL statement by creating first a resultSet object
  rs <- dbSendQuery(conn, data[[1]])
  
  # fetch records from the resultSet into a dataframe
  data[[2]] <- fetch(rs)
  
  # close database connection
  conn <- dbDisconnect(conn)
  
  return(data)
}
