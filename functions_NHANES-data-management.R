# Funcitons for NHANES data management
# Created: 2018-10-14
# Last edited: 2018-10-14
# Required packages: RNHANES, dplyr
# Previously called: NHANES_data-management-functions.R
library(RNHANES)
library(dplyr)

################################################################################
# Function: nhanes_yr_combiner()
# Description: Function pulls and combines data from a specified dataset
# across multiple cycles. Takes string vector dat and vector yrs where dat is the
# name of the dataset or a vector (or a vector of said dataset including cycle codes) 
# and yrs are the cycles.
# Requires: RNHANES
nhanes_yr_combiner <- function(dat, yrs){
  sets <- list()
  for(i in 1:length(yrs)){
    #pull data
    # If dat is a single string, rnhanes will automatically append cycle code to the end
    if(length(dat) == 1){
      sets[[i]] <- nhanes_load_data(dat, yrs[i], demographics=F) 
      #note: double brackets required for certain nested data structures in lists
      #combine data
    # If dat is a vector of multiple strings, rnhanes won't try to append cycle code
    # this helps for some edge cases where rnahnes append doesn't work.
    }else{
      sets[[i]] <- nhanes_load_data(dat[i], yrs[i], demographics=F) 
    }
    if(i == 1){
      d <- as.data.frame(sets[[i]])
      names(d)
    }else{
      d <- bind_rows(d, as.data.frame(sets[[i]]))
    }
  }#closefor
  return(d)
}#closefunc
###

################################################################################
# Function: nhanes_merge()
# Description: merges a list of NHANES datasets with one another. This is 
# intended to be used with different dataset types (ex. bpq and demo), NOT
# the same type across different cycles.
# The variable dat is a list of data frames.
# Requires: dplyr
nhanes_merge <- function(dat){
  byvars <- c("SEQN", "begin_year", "end_year", "cycle")
  
  # Initialize merged dataset and strip file_name (confuses merge)
  mergedat <- subset(dat[[1]], select=-file_name)
  
  # Merge remaining data frames
  for(i in 2:length(dat)){
    m <- subset(dat[[i]], select=-file_name)
    mergedat <- merge(mergedat, m, by=byvars, all = T)
  } #closefor
  return(mergedat)
} #closefunc
###
