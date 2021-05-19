# Set up correction for misclassificaiton of Hep B immune status 
#   when measured with anti-HBs using NHANES data.
# Author: Daniel Vader
# Dependencies: 
    # m1.txt - file containing JAGS models built in 3-3_misclass_models.R
    # 3-1_misclass_functions.R
    # 3-2_misclass_model_init_functions.R

# Load packages
library(tidyverse)
library(coda)

################################################################################
### Step 1: Intialize indicators 
################################################################################

# Load data prepared in data step
nhanes.vac <- readRDS("data/nhanes_vac2.rds") %>% dplyr::filter(RIDAGEYR >= 19)

### Setup indicator variables ##################################################
nhanes.vac$srvac_gt3 <- ifelse(nhanes.vac$hepbs == 1, 1, 0) # 3+ doses
nhanes.vac$srvac_lt3 <- ifelse(nhanes.vac$hepbs == 2 | nhanes.vac$hepbs == 3, 1, 0) # <3 doses
nhanes.vac$srvac_dnr <- ifelse(nhanes.vac$hepbs == 4, 1, 0) # don't know/refused

# race/eth, ref = white
nhanes.vac$reth_hisp <- ifelse(nhanes.vac$RIDRETH3 <= 2, 1, 0)
nhanes.vac$reth_black <- ifelse(nhanes.vac$RIDRETH3 == 4, 1, 0)
nhanes.vac$reth_asian <- ifelse(nhanes.vac$RIDRETH3 == 6, 1, 0)
nhanes.vac$reth_other <- ifelse(nhanes.vac$RIDRETH3 == 7, 1, 0)

# Age groups, ref = 50+
nhanes.vac$age_lt30 <- ifelse(nhanes.vac$RIDAGEYR < 30, 1, 0)
nhanes.vac$age_3049 <- ifelse(nhanes.vac$RIDAGEYR < 50 & nhanes.vac$RIDAGEYR >= 30, 1, 0)
nhanes.vac$age_gte50 <- ifelse(nhanes.vac$RIDAGEYR >= 50, 1, 0)

# Sensitivity splitter
nhanes.vac$age_lte25 <- ifelse(nhanes.vac$RIDAGEYR <= 25, 1, 0)

# Income-poverty ratio <=130%
nhanes.vac$incpov_lt130 <- ifelse(is.na(nhanes.vac$INDFMPIR), 0,
                                  ifelse(nhanes.vac$INDFMPIR <= 1.3, 1, 0))
nhanes.vac$incpov_na <- ifelse(is.na(nhanes.vac$INDFMPIR), 1, 0)

# Sex
nhanes.vac$male <- ifelse(nhanes.vac$RIAGENDR == 2, 0, 1)

# Switch strata to easy indicators
nhanes.vac$strata <- factor(nhanes.vac$SDMVSTRA)
nhanes.vac$strata <- droplevels(nhanes.vac$strata)
nhanes.vac$strata <- as.integer(nhanes.vac$strata)


# Save abbreviated data set with indicator variables
t <- nhanes.vac %>% dplyr::select(srvac_gt3, srvac_lt3, srvac_dnr, reth_hisp, reth_black,
                           reth_asian, reth_other, RIDAGEYR, male, age_lt30, age_3049, age_lte25,
                           fborn, WTMEC2YR, hepbv, SDMVPSU, SDMVSTRA, incpov_lt130, 
                           incpov_na, strata, RIDSTATR) %>%
  dplyr::filter(RIDAGEYR >= 19, !is.na(fborn), !is.na(hepbv)) # 252 no MEC interview, 451 missing hepbv, 1 missing fborn
saveRDS(t, "data/nhanes_vac3.rds") 


################################################################################
### Step 1.5: Start here to run models ###
################################################################################

### Model 1: Distribution centered at 65% ######################################
source("3-2_misclass_model-init-functions.R")

# Set priors for sensitivity and specificity
beta.sens.a <- 105+1 # set beta distribution for sn0 in JAGS model
beta.sens.b <- 109+1
beta.spec.a <- 265+1
beta.spec.b <- 18+1
beta.sens1.a <- 210+1
beta.sens1.b <- 34+1

# initialize jags data with new beta distribution
jags.dat <- init.jags.dat() 

# Run JAGS model with parallel processing. 
  # 3 chains. 20000 iterations. 2000 burn in. Thin=1
jf1 <- jags.parallel(data=jags.dat, #inits=jags.inits, 
                     parameters.to.save = jags.params, model.file = "m1.txt",
                     n.chains=3, n.iter=15000, n.burnin = 2000, n.thin=1)

# Save data 
saveRDS(jf1, "D:/researchdat/aim2/jf1.rds") 


# Run a traceplot of key parameters to check chain mixing
R2jags::traceplot(jf1, mfrow = c(2, 2), varname = c("m_pvac", "m_pvacstar", "sn0", "sp0", "sn1", "sp1"))


### Model 2: Sens distribution centered at 87%; Same spec. ########################
  #(Trivello et al. 1995); Orig: 816/942 -> (816+1, 126+1)

# .rs.restartR() # Restart R session to completely clean out memory

# Initialize resources needed to run this model
source("3-2_misclass_model-init-functions.R")

# Set initial values (new values for new sens/spec)

beta.sens.a <- 81+1 # set beta distribution for sn0 in JAGS model
beta.sens.b <- 215+1
beta.spec.a <- 265+1
beta.spec.b <- 18+1
beta.sens1.a <- 210+1
beta.sens1.b <- 34+1

# initialize jags data with new beta distribution
jags.dat <- init.jags.dat() 

# Run JAGS model with parallel processing. 
# 3 chains. 20000 iterations. 2000 burn in. Thin = 2
jf2 <- jags.parallel(data=jags.dat, 
                     parameters.to.save = jags.params, model.file = "m1.txt",
                     n.chains=3, n.iter=15000, n.burnin = 2000, n.thin=1)

# Save data 
saveRDS(jf2, "D:/researchdat/aim2/jf2.rds") 

# Run a traceplot of key parameters to check chain mixing
R2jags::traceplot(jf2, mfrow = c(2, 2), varname = c("m_pvac", "m_pvacstar", "sn0", "sp0", "sn1", "sp1"))


### Model 3: Distribution centered at 38% ######################################
   # (Poovorawan et al. 2010); Orig: 81/215 -> (81+1,134+1)
   # TP/P: Hammit (2/23), Bassal (17/45), Lin (215/580), Su (160/398) -> (2+17+215+160+81+1, 21+28+365+238+134+1) -> 
     # FP/N: 17/236 FP (Ismail), 1/47 FP (Huzley) -> (219+46, 1+17)
     # Laboratory Evaluation of a Fully Automated Chemiluminescence Immunoassay for Rapid Detection of HBsAg, Antibodies to HBsAg, and Antibodies to Hepatitis C Virus  
 .rs.restartR() # Restart R session to completely clean out memory

# Initialize resources needed to run this model
source("3-2_misclass_model-init-functions.R")

beta.sens.a <- 105+1 # set beta distribution for sn0 in JAGS model
beta.sens.b <- 109+1
beta.spec.a <- 265+1
beta.spec.b <- 18+1
beta.sens1.a <- 95+1
beta.sens1.b <- 6+1

# initialize jags data with new beta distribution
jags.dat <- init.jags.dat() 

# Run JAGS model with parallel processing. 
# 3 chains. 20000 iterations. 2000 burn in. 
jf3 <- jags.parallel(data=jags.dat, 
                     parameters.to.save = jags.params, model.file = "m1.txt",
                     n.chains=3, n.iter=15000, n.burnin = 2000, n.thin=1)

# Save data 
saveRDS(jf3, "D:/researchdat/aim2/jf3.rds") # save RJags data


# Run a traceplot of key parameters to check chain mixing
R2jags::traceplot(jf3, mfrow = c(2, 2), varname = c("m_pvac", "m_pvacstar", "sn0", "sp0"))


### Further model / parameter evaluation #######################################
library(tidyr)
library(bayesplot)
library(R2jags)
library(ggplot2)

color_scheme_set("viridis")
ps <- c("m_pvacstar", "m_pvac", "m_psrvac", "sn0", "sp0", "sn1")

# Model 1
jf1.mcmc <- readRDS("D:/researchdat/aim2/jf1.rds") %>% as.mcmc()
gc()

tiff("jf1-density.tiff", units="in", width=7, height=4, res=300)
jf1.density <- mcmc_dens_overlay(jf1.mcmc, pars=ps) + theme_dark()
jf1.density
dev.off()

tiff("jf1-trace.tiff", units="in", width=7, height=3, res=300)
jf1.trace <- mcmc_trace(jf1.mcmc, pars=ps) + theme_dark()
jf1.trace
dev.off()

# Model 2
jf2.mcmc <- readRDS("D:/researchdat/aim2/jf2.rds") %>% as.mcmc()
gc()

tiff("jf2-density.tiff", units="in", width=7, height=4, res=300)
jf2.density <- mcmc_dens_overlay(jf2.mcmc, pars=ps) + theme_dark()
jf2.density
dev.off()

tiff("jf2-trace.tiff", units="in", width=7, height=3, res=300)
jf2.trace <- mcmc_trace(jf2.mcmc, pars=ps) + theme_dark()
jf2.trace
dev.off()

# Model 3
jf3.mcmc <- readRDS("D:/researchdat/aim2/jf3.rds") %>% as.mcmc()
gc()

tiff("jf3-density.tiff", units="in", width=7, height=4, res=300)
jf3.density <- mcmc_dens_overlay(jf3.mcmc, pars=ps) + theme_dark()
jf3.density
dev.off()

tiff("jf3-trace.tiff", units="in", width=7, height=3, res=300)
jf3.trace <- mcmc_trace(jf3.mcmc, pars=ps) + theme_dark()
jf3.trace
dev.off()


################################################################################
### Step 2: Build poststratification table from ACS 2018 5-year public 
###   microdata surveys
################################################################################
library(tidyr)

# Read in data and setup variables so that names match those in the JAGS model
ps.demo <- read.csv("poststrat_table_cleaned_fb.csv") 
ps.demo$male<- ifelse(ps.demo$sex == "Male", 1, 0)
ps.demo$fborn <- ifelse(ps.demo$fb == "no", 0, 1)
ps.demo$age1 <- ifelse(ps.demo$age == "19 to 29", 1, 0)
ps.demo$age2 <- ifelse(ps.demo$age == "30 to 49", 1, 0)
ps.demo$incpov1 <- ifelse(is.na(ps.demo$pov_ratio), 0,
                         ifelse(ps.demo$pov_ratio == "<=130%", 1, 0))
ps.demo$incpov2 <- ifelse(is.na(ps.demo$pov_ratio), 1, 0)
ps.demo$reth_hisp <- ifelse(ps.demo$hispanic == "yes", 1, 0)
ps.demo$reth_black <- ifelse(ps.demo$hispanic == "no" & ps.demo$race == "black", 1, 0)
ps.demo$reth_asian <- ifelse(ps.demo$hispanic == "no" & ps.demo$race == "asian", 1, 0)
ps.demo$reth_other <- ifelse(ps.demo$hispanic == "no" & ps.demo$race == "other", 1, 0)
names(ps.demo)[1] <- "n" # remove name artifact from csv

# Remove na values from income/poverty
#ps.demo2 <- ps.demo %>% dplyr::filter(!is.na(incpov))

# Sum across duplicate hispanic / race rows created in recode
ps.demo2 <- ps.demo %>% 
  dplyr::group_by(male, fborn, age1, age2, incpov1, incpov2, reth_hisp, reth_black, 
           reth_asian, reth_other) %>%
  dplyr::summarise(n=sum(n)) %>% as.data.frame()
 
# Reorder variables and remove extraneous            
ps.demo2 <- ps.demo2 %>% 
  dplyr::select(age1, age2, male, incpov1, incpov2, reth_hisp, reth_black, reth_asian, 
                reth_other, fborn, n) 

# Calculate population weight
ps.demo2$pwt <- ps.demo2$n / sum(ps.demo2$n)

# Sort data
ps.demo2 <- ps.demo2[order(ps.demo2$age1, ps.demo2$age2, ps.demo2$male, ps.demo2$incpov1, ps.demo2$incpov2,
               ps.demo2$reth_hisp, ps.demo2$reth_black, ps.demo2$reth_asian,
               ps.demo2$reth_other),]

# Add strata identifiers to data
ps.demo2$strata <- 1:nrow(ps.demo2)

saveRDS(ps.demo2, "data/psdemo.rds")



################################################################################
### Step 3: Add poststratification weights to posterior samples and calculate 
###    prevalence. You can start here if Step 1 and Step 2 have already been run
################################################################################

library(tidyr)
source("3-1_misclass_functions.R")

jf1.list <- jf1$BUGSoutput$sims.list
jf1.poststrat <- calc.poststrat.prev(jf1.list)
saveRDS(jf1.poststrat, "D:/researchdat/aim2/jf1_poststrat-prev.rds")

jf2.list <- jf2$BUGSoutput$sims.list
jf2.poststrat <- calc.poststrat.prev(jf2.list)
saveRDS(jf2.poststrat, "D:/researchdat/aim2/jf2_poststrat-prev.rds")

jf3.list <- jf3$BUGSoutput$sims.list
jf3.poststrat <- calc.poststrat.prev(jf3.list)
saveRDS(jf3.poststrat, "D:/researchdat/aim2/jf3_poststrat-prev.rds")


# Make table of results ########################################################
library(tidyr)

prev.tab <- function(dat){
  cri <- c(0.5, 0.025, 0.975)
  prev.tab <- data.frame(
    p_vac=numeric(),
    p_vac_lci=numeric(),
    p_vac_uci=numeric(),
    ps_vac=numeric(),
    ps_vac_lci=numeric(),
    ps_vac_uci=numeric(),
    p_vac_star=numeric(),
    p_vac_star_lci=numeric(),
    p_vac_star_uci=numeric(),
    ps_vac_star=numeric(),
    ps_vac_star_lci=numeric(),
    ps_vac_star_uci=numeric(),
    p_srvac=numeric(),
    p_srvac_lci=numeric(),
    p_srvac_uci=numeric(),
    ps_srvac=numeric(),
    ps_srvac_lci=numeric(),
    ps_srvac_uci=numeric()
  )
  g1 <- NULL
  for(i in 1:length(dat)){
    q <- quantile(dat[,i], cri, na.rm = T)
    g1 <- c(g1, q)
    if(i%%6 == 0){
      index <- i/6
      prev.tab[index,] <- g1
      g1 <- NULL
    }
  }
  
  prev.tab <- prev.tab[,c(1:3,7:9,13:15,4:6,10:12,16:18)]
  
  rownames(prev.tab) <- c("marginal", "white", "black", "asian", "other", "hispanic",
                          "age1", "age2", "age3", "fb","fbn")
  prev.tab <- prev.tab %>% tibble::rownames_to_column()
  return(prev.tab)
}

# Calculate absolute bias of anti-HBs+ and self-report against 
# adjusted estimates of immunity
bias.tab <- function(dat){
  cri <- c(0.5, 0.025, 0.975)
  bias.tab <- data.frame(
    b_vacstar=numeric(),
    b_vacstar_lci=numeric(),
    b_vacstar_uci=numeric(),
    b_srvac=numeric(),
    b_srvac_lci=numeric(),
    b_srvac_uci=numeric()
  )
  for(i in 1:length(dat)){
    q1 <- quantile(dat[[i]][4]-dat[[i]][[2]], cri, na.rm = T)
    q2 <- quantile(dat[[i]][6]-dat[[i]][[2]], cri, na.rm = T)
    r <- c(q1, q2)
    bias.tab[i,] <- r
  }
  rownames(bias.tab) <- c("marginal", "white", "black", "asian", "other", "hispanic",
                          "age1", "age2", "age3", "fb", "fbn")
  bias.tab <- bias.tab %>% tibble::rownames_to_column()
  return(bias.tab)
}

ps1 <- readRDS("D:/researchdat/aim2/jf1_poststrat-prev.rds") %>% as.data.frame() %>% dplyr::select(-contains("tweight"))
ps2 <- readRDS("D:/researchdat/aim2/jf2_poststrat-prev.rds") %>% as.data.frame() %>% dplyr::select(-contains("tweight"))
ps3 <- readRDS("D:/researchdat/aim2/jf3_poststrat-prev.rds") %>% as.data.frame() %>% dplyr::select(-contains("tweight"))

prev.tab.ps1 <- prev.tab(ps1)
prev.tab.ps2 <- prev.tab(ps2)
prev.tab.ps3 <- prev.tab(ps3)

write.csv(prev.tab.ps1, "data/jf1_results_table2.csv")
write.csv(prev.tab.ps2, "data/jf2_results_table2.csv")
write.csv(prev.tab.ps3, "data/jf3_results_table2.csv")

ps1.list <- readRDS("D:/researchdat/aim2/jf1_poststrat-prev.rds") 
ps2.list <- readRDS("D:/researchdat/aim2/jf2_poststrat-prev.rds") 
ps3.list <- readRDS("D:/researchdat/aim2/jf3_poststrat-prev.rds") 

bias.tab.ps1 <- bias.tab(ps1.list)
bias.tab.ps2 <- bias.tab(ps2.list)
bias.tab.ps3 <- bias.tab(ps3.list)

write.csv(bias.tab.ps1, "data/jf1_biastab.csv")
write.csv(bias.tab.ps2, "data/jf2_biastab.csv")
write.csv(bias.tab.ps3, "data/jf3_biastab.csv")

################################################################################
### Plot bias
################################################################################
library(ggplot2)

# Put prevalence table data in long format by estimator 
# to make it easier to work with ggplot
longdat <- function(dat){
  t1.1 <- dat %>% dplyr::select(rowname, ps_vac, ps_vac_lci, ps_vac_uci)
  t1.1$stat <- "Immunity"
  names(t1.1) <- c("rowname", "est", "lci", "uci", "stat")
  
  t1.2 <- dat %>% dplyr::select(rowname, ps_vac_star, ps_vac_star_lci, ps_vac_star_uci)
  t1.2$stat <- "Anti-HBs +"
  names(t1.2) <- c("rowname", "est", "lci", "uci", "stat")
  
  t1.3 <- dat %>% dplyr::select(rowname, ps_srvac, ps_srvac_lci, ps_srvac_uci)
  t1.3$stat <- "Self-report"
  names(t1.3) <- c("rowname", "est", "lci", "uci", "stat")
  
  tc <- rbind(t1.1, t1.2, t1.3)
  
  # Turn rowname into a factor and give nicer labels
  tc$rowname <- factor(tc$rowname, levels=c("marginal", "white", "black", "hispanic", "asian", "other", 
                                            "age1", "age2", "age3", "fb","fbn"),
                       labels=c("Marginal", "White", "Black", "Hispanic", "Asian", "Mixed race / other",
                                 "Age 19 to 29", "Age 30 to 49", "Age 50+", 
                                "Not born in US", "Born in US"))
  return(tc)
}

t1.p <- longdat(prev.tab.ps1)
t2.p <- longdat(prev.tab.ps2)
t3.p <- longdat(prev.tab.ps3)

t1.p$model <- "Model 1"
t2.p$model <- "Model 2"
t3.p$model <- "Model 3"
tt <- rbind(t1.p, t2.p, t3.p)W

# Plot prevalence by model
p1 <- ggplot(data=tt, aes(y=est, x=stat, fill=stat)) + 
  geom_bar(stat="identity") + 
  geom_point() +
  geom_linerange(aes(ymin=lci, ymax=uci)) + 
  coord_flip() + 
  scale_fill_viridis_d() + 
  ylab("Prevalence") +
  facet_grid(rowname~model) + 
  xlab("") + 
  theme(axis.text.y = element_blank(), 
        strip.text.y = element_text(angle=0),
        axis.ticks.y = element_blank(),
        legend.position = "bottom") + 
  labs(fill = "Estimator")

tiff("m1-prev.tiff", units="in", width=7, height=6, res=300)
p1
dev.off()

# Plot bias by model and measure ###############################################

# Swap data to long format. Function same as longdat but for absolute bias tables.
longdat.b <- function(dat){
  
  t1.2 <- dat %>% dplyr::select(rowname, b_vacstar, b_vacstar_lci, b_vacstar_uci, model)
  t1.2$stat <- "Anti-HBs +"
  names(t1.2) <- c("rowname", "est", "lci", "uci", "model", "stat")
  
  t1.3 <- dat %>% dplyr::select(rowname, b_srvac, b_srvac_lci, b_srvac_uci, model)
  t1.3$stat <- "Self-report"
  names(t1.3) <- c("rowname", "est", "lci", "uci", "model", "stat")
  
  tc <- rbind(t1.2, t1.3)
  tc$rowname <- factor(tc$rowname)
  return(tc)
}

bbl <- longdat.b(bb)

p7 <- ggplot(data=bbl, aes(y=est, x=stat, fill=stat)) + 
  geom_bar(stat="identity") + 
  geom_point() +
  geom_linerange(aes(ymin=lci, ymax=uci)) + 
  coord_flip() + 
  scale_fill_viridis_d(option="D") + 
  ylab("Absolute bias") +
  facet_grid(rowname~model) + 
  xlab("") + 
  theme(axis.text.y = element_blank(), 
        strip.text.y = element_text(angle=0),
        axis.ticks.y = element_blank(),
        legend.position = "bottom") + 
  labs(fill = "Estimator")

tiff("ab-allmod.tiff", units="in", width=6, height=5, res=300)
p7
dev.off()

