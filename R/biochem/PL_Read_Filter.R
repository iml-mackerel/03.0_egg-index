PL_Read_Filter <- function(filter_file) {
  
  # Function to read a taxonomic filter text file.
  #
  # Input: filter_file : filter filename - string
  #
  # Output: filter : df.filter - dataframe
  #
  # Last update: 20151123
  # Benoit.Casault@dfo-mpo.gc.ca
  
  # required package
  library(dplyr)

  ## read filter file
  df.filter <- read.table(filter_file, header=TRUE, sep="\t", stringsAsFactors=FALSE, na.strings="",
                          comment.char="#", check.names=FALSE)
  df.filter<-df.filter[, -6]  #enle`ve la colonne code 1
  
  ## reshape filter to long format
  group_levels <- base::setdiff(names(df.filter), c("taxonomic_name", "stage", "molt_number", "sex"))
  df.filter <- df.filter %>%
    tidyr::gather(., group_level, group_name, one_of(group_levels)) %>%      # melt the data
    dplyr::select(., -group_level) %>%                                       # subset variables
    dplyr::filter(., !is.na(group_name)) %>%                                 # remove NA group
    dplyr::distinct(.)                                                       # eliminate duplicate rows
  
  # output
  return(df.filter)
}