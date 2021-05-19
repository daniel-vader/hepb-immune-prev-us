### Load data from NHANES and merge
### Author: Daniel Vader

#library(devtools)
#devtools::install_github("silentspringinstitute/RNHANES") # Install latest v of RNHANES
library(tidyverse)
library(RNHANES)

source("../functions_NHANES-data-management.R")


### Initialize local funcitons #################################################
# Function to factorize standard yes/no NHANES variables
ynfactor <- function(var){
  factor(var, levels=c(1,2,7,9,NA),
         labels=c("Yes","No","Refused","Don't Know","NA"),
         exclude=NULL)
}


### Pull data from NHANES ######################################################

# Pull 2015-2016 data
# Sets: DEMO (demographic questionnaire), IMQ (immunization questionnaire), 
#       HEPBD (surface antigen, core antibody), HEPB_S (surface antibody)

#create vector of years
cycles <- c("2015-2016") # can add cycles to the vector as necessary. I = 2015-2016

#create list of nhanes datasets to merge. Only include sets that exist in all 
#   three cycles.
data.list.all <- c("DEMO", "IMQ", "HEPBD", "DUQ", "SXQ", "HUQ")

# Pull data and merge across cycles
sets.merger <- list()
for(i in 1:length(data.list.all)){
  print(data.list.all[i])
  sets.merger[[i]] <- nhanes_yr_combiner(data.list.all[i], cycles)
}

# RNHANES doesn't like appending the cycle code to the end of HEPB_S
hbs_codes <- c("HEPB_S_I")
sets.merger[[length(sets.merger) + 1]] <- nhanes_yr_combiner(hbs_codes, cycles)

# merge across NHANES data category (aka, elements of data.list.all) and save
nhanes.vac <- nhanes_merge(sets.merger)


### Create indicator and weight variable(s) #####################################

# HepB Vaccination/Immunity indicator variables (do we want to do this in this file or in the analysis?)
# Subject is considered immune if HBsAg negative and anti-HBs positive.

# nhanes.vac$hepbv <- ifelse(nhanes.vac$LBXHBS > 1, 0, # IF not pos for anti-HBs, not immune
#                      ifelse(nhanes.vac$LBXHBC > 1, 1, # IF anti-HBc negative, immune
#                             ifelse(nhanes.vac$LBDHBG > 1, 1, 0) # IF hbsag negative, immune
#                             )
#                      )

# same result when looking only at anti-hbs
nhanes.vac$hepbv <- ifelse(nhanes.vac$LBXHBS == 1, 1, 0) # If anti-HBs postitive, immune. Else, not.

# Self report. If refused or don't know (7,9) then combine cat (4). Otherwise keep
#   original cats (1 = 3+ doses, 2 = less than 3 doses, 3 = no doses)
nhanes.vac$hepbs <- ifelse(nhanes.vac$IMQ020 >= 7, 4, nhanes.vac$IMQ020)


# Create weights. If subjects do not fall into our target population 
# (age>=18), then set to a near zero weight. If they do, divide assigned weight
# by number of cycles. Ref: https://wwwn.cdc.gov/nchs/data/nhanes/2011-2012/analyticguidelines/analytic_guidelines_11_16.pdf
num.cyc <- length(cycles)

nhanes.vac$wt_mec_n <- ifelse(nhanes.vac$RIDAGEYR < 18, 1*10^(-6),
                              nhanes.vac$WTMEC2YR/num.cyc) 

### Save data ##################################################################

# Drop subjects who did not complete the MEC exam
# nhanes.vac.mec <- nhanes.vac %>% filter(RIDSTATR == 2) 

# saveRDS(nhanes.vac.mec, "data/nhanes_vac1.rds")
saveRDS(nhanes.vac, "data/nhanes_vac1.rds")


### Create factors #############################################################
nhanes.vac <- readRDS("data/nhanes_vac1.rds")

# gender
nhanes.vac$female <- nhanes.vac$RIAGENDR - 1 # 0 if male, 1 if female


# country of birth (note, US born category in NHANES only includes the 50 states
#   and not US territories). Classify subjects who refused to respond or responded
#   that they "don't know" if they were born in the US as NA.
nhanes.vac$fborn <- ifelse(nhanes.vac$DMDBORN4 > 3, NA, nhanes.vac$DMDBORN4 - 1) # shift to 0 (US Born) 1 (foreign born)

# education
nhanes.vac$educ <- ifelse(nhanes.vac$RIDAGEYR < 20, # Subjects under 20 use a different ed classification var
                          ifelse(nhanes.vac$DMDEDUC3 > 70, NA, #refused or don't know set to NA
                                 ifelse(nhanes.vac$DMDEDUC3 > 50, 1, #Subjects with less than 9th grade ed
                                    ifelse(nhanes.vac$DMDEDUC3 >= 13, 3, #Subjects with GED
                                           ifelse(nhanes.vac$DMDEDUC3 >= 9, 2, 1) #Subjects with 9th-11th grade ed, anyone else less than 9th
                                           )
                                 )
                          ), ifelse(nhanes.vac$DMDEDUC2 > 5, NA, nhanes.vac$DMDEDUC2))
# nhanes.vac$educ <- factor(nhanes.vac$educ, levels = c(1,2,3,4,5), 
#                           labels = c("< 9th grade", "9th-11th grade", "High school/GED",
#                                      "Some college/AA", "College graduate"))

# income (use income-poverty ratio)

# Drug use. Dataset: DUQ
    # Have you ever, even once, used a needle to inject a drug not 
        # prescribed by a doctor? (DUQ370)
    # 1 = yes, 2 = no, 7 = refused, 9 = don't know. 
nhanes.vac$idu <- ifelse(nhanes.vac$DUQ370 > 3, NA, nhanes.vac$DUQ370) %>%
  factor(levels = c(1,2), labels = c("Yes", "No"))

# sexual activity. Dataset: SXQ 
    # Num F sex partners / lifetime (M) - SXD171
    # Num M anal/oral sex partners / lifetime (M) - SXQ410
    # Had sex with new partner past 12mo - SXQ648
    # Num M sex partners lifetime (F) - SXD101
    # Num F sex partners lifetime (F) - SXQ130
    # Num times had vaginal or anal sex/year - SXQ610 ** 
# Notes: is this important enough to use? This is a risk factor
    # for infection but not necessarily vaccination.
    # SXQ610 has a lot of missing due to age (<19 and >59). May also have missing
    # due to prereq variables.
nhanes.vac$sexp <- ifelse(nhanes.vac$SXQ610 > 70, NA, 
                          ifelse(nhanes.vac$SXQ610 > 4, 4, 
                                 nhanes.vac$SXQ610)) %>%
  factor(levels=c(0,1,2,3,4), labels=c("Never", "Once", "2-11 times", 
                                       "12-51 times", ">51 times"))


# Doctor's visit. Dataset: HUQ
    # How long since last healthcare visit? (HUQ061)
    # 1 = <= 6mo, 2 = >6mo and <=1yr, 3 = >1yr and <= 2yrs, 4 = >2yrs and <= 5yrs
        # 5 = > 5yrs, 6 = Never?, 77 = refused, 99 = don't know
    # also requires HUQ051: # times received health care in past year
nhanes.vac$hvisit <- ifelse(nhanes.vac$HUQ051 > 70, NA, 
                            ifelse(nhanes.vac$HUQ051 > 0, 1,
                                   ifelse(nhanes.vac$HUQ061 <= 2, 1,
                                          ifelse(nhanes.vac$HUQ061 == 3, 2,
                                                 ifelse(nhanes.vac$HUQ061 == 4, 3,
                                                        ifelse(nhanes.vac$HUQ061 <= 6, 4, NA
                                                               )
                                                        )
                                                 )
                                          )
                                   )

                            )
nhanes.vac$hvisit <- factor(nhanes.vac$hvisit, levels=c(1,2,3,4),
                            labels=c("<=1 year", "1-2 years", "2-5 years", ">5 years"))

# Save dataset
saveRDS(nhanes.vac, "data/nhanes_vac2.rds")
