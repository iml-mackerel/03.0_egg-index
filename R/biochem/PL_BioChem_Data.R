PL_Biochem_Data <- function(sql_file, output_file, host, port, sid, username, password) {
    
  # Function to extract data from database, perform basic quality control and save the data in .Rdata and
  # tsv (tab seperated value) formats.
  #
  # Input: sql_file : sql filename - string
  #        data_source : data source name for database - string
  #        username : username - string
  #        password: user password - string
  #        output_file : output filename - string
  #
  # Output: None
  #
  # Last update: 20141015
  # Benoit.Casault
  
  # source custom functions
  source("R/biochem/Run_Database_Query.R")
  
  # get the data
  data <- Run_Database_Query(sql_file, host, port, sid, username, password)

  # rename variables
  sql_str <- data[[1]]
  df <- data[[2]]
  rm(data)

  # rename df variables to lower case
  col_names <- tolower(names(df))
  names(df) <- col_names
  
  # find data type (i.e. whether sample_id's, counts or weight data)
  samples <- FALSE
  weights <- FALSE
  counts <- FALSE
  if ("wet_weight" %in% col_names | "dry_weight" %in% col_names) {
    weights <- TRUE
  } else if ("counts" %in% col_names) {
    counts <- TRUE    
  } else {
    samples <- TRUE
  }
  
  # quality control checks
  # checks on sample_id's
  if (samples) {
    # check depth values
    index <- df$start_depth==df$end_depth
    if (any(index)) {
      index <- which(index)
      sample_str <- paste(df$custom_sample_id[index], collapse=',')
      warning("Line 38: cases where start_depth=end_depth were found\n",
              "  Sample_ID: ", sample_str,"\n",
              "  These samples were removed from further analysis")
      # remove samples from dataframe
      df <- df[-index,]
    }
    # check volume values
    index <- is.na(df$volume)
    if (any(index)) {
         warning("cases where volume=NA were found\n"
              )
    }
    # checks on data
  } else if (weights | counts) {
    # check split fraction
    index <- df$split_fraction==0 | is.na(df$split_fraction)
    if (any(index)) {
      index <- which(index)
      sample_str <- paste(df$custom_sample_id[index], collapse=',')
      warning("Line 61: cases where split_fraction=0 or NA were found\n",
              "  Sample_ID: ", sample_str,"\n",
              "  These samples were removed from further analysis")
      # remove samples from dataframe
      df <- df[-index,]
    }
    # check sieve sizes  Need some work because of NA is sieve entries
    #       index <- df$min_sieve==df$max_sieve
    #       if (any(index)) {
    #         index <- which(index)
    #         sample_str <- paste(df$custom_sample_id[index], collapse=',')
    #         warning("Line 61: cases where min_sieve==max_sieve were found\n",
    #                 "  Sample_ID: ", sample_str,"\n",
    #                 "  These samples were removed from further analysis")
    #         # remove samples from dataframe
    #         df <- df[-index,]
    #       }
    # There should be another check to detect multiple values of split fraction for given parameter
  }
  
    # print to tsv (tab seperated values) file
 
  write.table(df, file=gsub(output_file, pattern="RData", replacement="tsv"), dec=".", sep="\t")
  write_rds(df, file=output_file)
  
  return()
}