PL_Taxonomic_Grouping <- function(df.data, df.filter) {

  ## required packages
  library(dplyr)
  library(tidyr)

  ## subset input data to taxonomic_name appearing in filter
  target_names <- base::unique(df.filter$taxonomic_name)
  df.data <- df.data %>%
    dplyr::filter(., taxonomic_name %in% target_names)

  ## unique groups
  group_names <- base::sort(base::unique(df.filter$group_name))

  ## indexing - match data vs filter
  for (i_group in group_names) {

    ## initialize index for given group
    index <- rep(FALSE, nrow(df.data))

    ## subset entries for given group
    tmp <- df.filter %>%
      dplyr::filter(., group_name==i_group) %>%                      # filter for given group
      dplyr::select(., taxonomic_name, stage, molt_number, sex) %>%  # select columns for matching with data
      dplyr::mutate(., index=TRUE)                                   # add index column (for join operation below)

    for (i_row in seq(1, nrow(tmp))) {

      ## select columns to match for given entry
      column_names <- base::setdiff(names(tmp[i_row,!is.na(tmp[i_row,])]), "index")

      ## update index for given group
#       if (length(column_names==0)) {
#         index <- index | rep(TRUE, nrow(df.data))  # default case: matching everything (i.e. name,stage,molt,sex are empty)
#       } else {
        index <- index |
          as.vector(!(dplyr::left_join(dplyr::select(df.data, column_names),                  # select only columns to match in data
                                       dplyr::select(tmp[i_row,], c(column_names, "index")),  # select only columns to match in filter
                                       by=column_names)  %>%                                         # join based on all selected columns
                        dplyr::select(., index) %>%                                                  # select index variable (desired result)
                        is.na))                                                                      # assign 1's to TRUE and NA's to FALSE
#       }
    }
    ## add group index to data as new variable
    df.data <- cbind(df.data, index) %>%                         # append index as new column
      dplyr::rename(setNames("index", i_group))  # rename index variable to given group name
  }

  ## calculate totals for each group
  df.grouped <- df.data %>%
    tidyr::gather(., group_name, index, one_of(group_names)) %>%    # melt the data
    dplyr::filter(., index==TRUE) %>%                               # keep only entries belonging to given groups
    dplyr::select(., sample_id, n_m2, group_name) %>%                 # subset variables
    dplyr::group_by(., sample_id, group_name) %>%                   # group by sample and group_name
    dplyr::summarise(., total=sum(n_m2)) %>%                          # calculate total
    tidyr::spread(., group_name, total, fill=0)                     # cast the data

  ## output
  return(df.grouped)
}
