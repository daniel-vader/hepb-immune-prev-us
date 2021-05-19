### Descriptive statistics ###
# Author: Daniel Vader

library(meta)
library(tidyr)
library(ggplot2)
library(ggpubr)

### Import data and setup variables ###
hbsdat <- read.csv("data/antihbs_spec_dat.csv") 
names(hbsdat)[1] <- c("authors") # CSV file saved this column name badly

# Combine author and year into one var for labeling
hbsdat$authyr <- paste(hbsdat$authors, " (", hbsdat$year, ")", sep="")

hbsdat$prepost <- ifelse(hbsdat$postboost == 1, "Post-booster", "Pre-booster")

### Data format as entered is too split up, collapse to one line per "study"
hbsdat_pre <- hbsdat %>% dplyr::filter(postboost == 0)
hbsdat_post <- hbsdat %>% dplyr::filter(postboost > 0)

hbsdat_rc <- hbsdat_pre
names(hbsdat_rc)[names(hbsdat_rc) == "n_pos"] <- "n_pos_pre"
names(hbsdat_rc)[names(hbsdat_rc) == "n_total"] <- "n_total_pre"
hbsdat_rc$n_pos_post <- hbsdat_post$n_pos
hbsdat_rc$n_total_post <- hbsdat_post$n_total
hbsdat_rc$postboost <- hbsdat_post$postboost


################################################################################
### code adjusted postboost estimates for ltfu
# postboost = 1, full pre and post booster data.
# postboost = 2, postboost data only for subjects who were - on initial test.
################################################################################
# Pre-booster total = a
# Pre-booster(+) = b
# Pre-booster(-) = a - b = c
# Post-booster total = d
# Post-booster(+) = e
# Post-booster missing | only-negative subjects boosted = c-d = m
# Missing adjustment factor = d/c = f
# Adjusted post-booster proporiton = (e + f*b) / (d + f*b)
hbsdat_rc$adjf <- ifelse(hbsdat_rc$postboost == 2, hbsdat_rc$n_total_post / (hbsdat_rc$n_total_pre - hbsdat_rc$n_pos_pre), 1)
hbsdat_rc$n_pos_post_adj <- 
  ifelse(hbsdat_rc$postboost == 2, 
         round(hbsdat_rc$n_pos_pre * hbsdat_rc$adjf + hbsdat_rc$n_pos_post, digits = 0),
         hbsdat_rc$n_pos_post)
hbsdat_rc$n_total_post_adj <- 
  ifelse(hbsdat_rc$postboost == 2,
         round(hbsdat_rc$n_pos_pre * hbsdat_rc$adjf + hbsdat_rc$n_total_post, digits = 0),
         hbsdat_rc$n_total_post
  )
hbsdat_rc$n_pos_pre_adj <- ifelse(hbsdat_rc$postboost == 2, 
                                  round(hbsdat_rc$n_pos_pre * hbsdat_rc$adjf, digits = 0),
                                  hbsdat_rc$n_pos_pre) 
saveRDS(hbsdat_rc, "data/hbslit.rds")


################################################################################
### Create a table to display basic info on each study
################################################################################
studydet <- hbsdat_rc %>%
  dplyr::filter(vac_type != "unvaccinated", pre_neg_only != 1) %>% 
  dplyr::select(authors, year, age, yr_grp)

studydet$ages <- ifelse(studydet$age == "Infants", 1,
                        ifelse(studydet$age == "Children", 2,
                               ifelse(studydet$age == "Infants and Children", 1.5,
                               ifelse(studydet$age == "Adolescents", 3,
                                      ifelse(studydet$age == "Adults", 4, 5)))))
studydet$yr_grps <- ifelse(studydet$yr_grp == "5 to 9", 1,
                           ifelse(studydet$yr_grp == "10 to 14", 2,
                                  ifelse(studydet$yr_grp == "15 to 19", 3, 4)))
studydet <- studydet[order(studydet$yr_grps, studydet$ages),] # sort data
write.csv(studydet, "study_details.csv")



################################################################################
### Summarize data in forest plots for paper 
# Use the clopper pearson "exact" binomal interval to generate 95% CIs
# Random effects estimates calculated using the 
   # Hartung-Knapp method for random effects meta-analysis
   # https://cran.r-project.org/web/packages/meta/meta.pdf (pg 4)
################################################################################
lit <- readRDS("data/hbslit.rds") %>% 
  dplyr::filter(vac_type != "unvaccinated", 
                pre_neg_only != 1) #%>% 
  #dplyr::filter(age == "Adults" | age == "Adolescents" | 
  #                yr_grp == "15 to 19" | yr_grp == "20 to 24")
lit <- lit[order(lit$year),]

lit$ageclass <- ifelse(lit$age %in% c("Adults", "Adolescents"), "Vaccinated as adults or adolescents (age 13+)",
                       ifelse(lit$age %in% c("Infants", "Infants and Children", "Children"), "Vaccinated as infants or children (age <13)",
                              ifelse(lit$age == "All ages", "Vaccinated at any age", NA)))


# Vector for relabeling left column of forest plots
llab <- c("Study", "Pre-boost +", "Total +")
rlab <- c("Sensitivity", "95% CI")

### generate forest plots ###
# Use the clopper pearson "exact" binomial interval to generate 95% CIs

# Function to draw forest plots with options set.
pforest <- function(d){
  forest(d, xlim = c(0,1), 
         lty.fixed = NULL, lty.random = NULL, 
         leftlabs = llab,
         comb.fixed = F, comb.random = T, prediction=F, 
         pooled.events = T, pooled.totals = T, hetstat = F, 
         text.random.w = "Subgroup summary",
         text.random = "Summary",
         col.by = "#666666",
         weight.study = "random")
}

# Figure directory
fdir <- "C:/Users/Daniel/Google Drive/PhDEpi/Research/Thesis/Publication prep/Aim 2/figures/"


# 5 to 14 years post vaccination ###############################################
lit2 <- lit %>% dplyr::filter(yr_grp == "5 to 9" | yr_grp == "10 to 14", 
                              lit$ageclass != "Vaccinated as adults or adolescents (age 13+)")

# Pre-booster
mp1 <- metaprop(n_pos_pre, n_total_pre, studlab=authyr, 
                data=lit2, 
                method.ci = "cp") 
pforest(mp1)

# Post-booster
mp2 <- metaprop(n_pos_post_adj, n_total_post_adj, studlab=authyr, 
                data=lit2, 
                method.ci = "cp") 
pforest(mp2)

# Pre-booster sensitivity
mp.5_14 <- metaprop(n_pos_pre_adj, n_pos_post_adj, studlab=authyr, 
                     data=lit2, print.byvar = F,
                     method.ci = "cp") 
tiff(paste0(fdir, "sens_0514.tiff"), units="in", width=10, height=5.5, res=300)
pforest(mp.5_14)
dev.off()

# 15 to 19 years post vaccination ##############################################
lit3 <- lit %>% 
  dplyr::filter(yr_grp == "15 to 19", 
                lit$ageclass != "Vaccinated as adults or adolescents (age 13+)")

# Pre-booster
mp3 <- metaprop(n_pos_pre, n_total_pre, studlab=authyr, 
                data=lit3, 
                method.ci = "cp") 
pforest(mp3)

# Post-booster
mp4 <- metaprop(n_pos_post_adj, n_total_post_adj, studlab=authyr, 
                data=lit3, 
                method.ci = "cp") 
pforest(mp4)

# Pre-booster sensitivity
mp.15_19 <- metaprop(n_pos_pre_adj, n_pos_post_adj, studlab=authyr, 
                    data=lit3, print.byvar = F,
                    method.ci = "cp") 
tiff(paste0(fdir, "sens_a1519.tiff"), units="in", width=10, height=5.5, res=300)
pforest(mp.15_19)
dev.off()


# 20 to 24 years or Adults #####################################################
# (These data used to inform priors.)

# Pre-booster
lit4 <- lit %>% dplyr::filter(yr_grp == "20 to 24" | lit$ageclass == "Vaccinated as adults or adolescents (age 13+)")
mp5 <- metaprop(n_pos_pre, n_total_pre, studlab=authyr, 
                data=lit4, 
                method.ci = "cp") 
pforest(mp5)

# Post-booster
mp6 <- metaprop(n_pos_post_adj, n_total_post_adj, studlab=authyr, 
                data=lit4, 
                method.ci = "cp") 
pforest(mp6)

# Pre-booster sensitivity (export figure)
tiff(paste0(fdir, "sens_a20plus.tiff"), units="in", width=10, height=5.5, res=300)
  mp.20_24 <- metaprop(n_pos_pre_adj, n_pos_post_adj, studlab=authyr, 
                      data=lit4, byvar = ageclass, print.byvar = F,
                      method.ci = "cp") 
  pforest(mp.20_24)
dev.off()

################################################################################
### Create frequency table ###
################################################################################
library(tidyr)

# Setup variables (Create factors in desired order)
nhanes.vac <- readRDS("data/nhanes_vac2.rds") %>% dplyr::filter(RIDAGEYR >= 19, !is.na(fborn), !is.na(hepbv))

nhanes.vac$reth <- ifelse(nhanes.vac$RIDRETH3 <= 2, 4, # Hispanic
                   ifelse(nhanes.vac$RIDRETH3 == 3, 1, # White
                   ifelse(nhanes.vac$RIDRETH3 == 4, 2, # Black
                   ifelse(nhanes.vac$RIDRETH3 == 6, 3, # Asian American
                   ifelse(nhanes.vac$RIDRETH3 == 7, 5, # Other
                          NA))))) %>%
  factor(levels=c(1,2,3,4,5), labels=c("White", "Black","Asian American",
                                       "Hispanic","Other"))

nhanes.vac$srvac <- ifelse(nhanes.vac$hepbs == 1, 1,
                    ifelse(nhanes.vac$hepbs <=3, 2, 4)) %>%
  factor(levels=c(1,2,4), 
         labels=c(">= 3 doses", "< 3 doses", "Don't know/refused"))

nhanes.vac$agegrp <- ifelse(nhanes.vac$RIDAGEYR < 30, 1,
                     ifelse(nhanes.vac$RIDAGEYR <50, 2, 3)) %>%
  factor(levels=c(1,2,3), labels=c("< 30", "30 to 49", ">= 50"))

nhanes.vac$incpov <- ifelse(is.na(nhanes.vac$INDFMPIR), 4,
          ifelse(nhanes.vac$INDFMPIR <= 1.3, 1, 2)) %>%
  factor(levels=c(1,2,4), labels=c("<= 130%", "> 130%", "Missing"))

nhanes.vac$sex <- ifelse(nhanes.vac$RIAGENDR == 1, 1, 2) %>%
  factor(levels=c(1,2), labels=c("Male", "Female"))

nhanes.vac$fborn2 <- ifelse(nhanes.vac$fborn == 0, 1, 2) %>%
  factor(levels=c(1,2), labels=c("Born in US", "Not born in US"))

nhanes.vac$hepbv <- ifelse(nhanes.vac$hepbv == 1, 1, 2) %>%
  factor(levels=c(1,2), labels=c("antiHBs_pos", "antiHBs_neg"))

# Build table
tabvars <- c("reth", "agegrp", "fborn2", "srvac", "incpov", "sex")
for(i in 1:length(tabvars)){
  t <- table(nhanes.vac[,tabvars[i]], nhanes.vac[,"hepbv"]) # Get freqs
  p <- chisq.test(t) # Estimate chi-squared p-value
  tp <- t %>% prop.table(2) %>% as.data.frame.matrix() # Get percentages
  t <- t %>% as.data.frame.matrix() # Convert freq table to dataframe
  t$total <- t[,1] + t[,2]
  t$total_per <- t$total/sum(t$total) %>% as.vector()
  t$p <- p$p.value # Add p-value to freq dataframe
  
  # Bind freqs, percentages, and p-value together in order.
  tpp <- cbind(t[1], tp[1], t[2], tp[2], t[3], t[4], t[5])
  
  # Add this chunk to the main frequency table dataframe
  if(i == 1){
    freqtab <- tpp
  }else{
    freqtab <- rbind(freqtab, tpp)
  }
}
names(freqtab) <- c( "antiHBs_pos_n", "antiHBs_pos_per", "antiHBs_neg_n",
                    "antiHBs_neg_per", "total_n", "total_per", "p")
# Write table to file
write.csv(freqtab, "data/freqtab.csv")
