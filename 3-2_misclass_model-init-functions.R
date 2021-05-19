# Model init functions
# These need to run before running R2jags

library(R2jags)

# parameters to monitor
jags.params <- c("m_pvac", 
                 "m_pvacstar", 
                 "m_psrvac",
                 "sn0", 
                 "sp0",
                 "sn1",
                 "sp1",
                 #"p_vac", 
                 #"p_incpov",
                 "p_vac_star",
                 "a0", "b0", "a[1:9]","b[1:9]"
)

# Function for initializing data to be read by JAGS model
init.jags.dat <- function(){
  jags.dat <- 
    list(
      "vac_star"=nhanes.vac$hepbv,
      "srvac1"=nhanes.vac$srvac_gt3,
      "srvac2"=nhanes.vac$srvac_dnr,
      "male"=nhanes.vac$male,
      "reth_hisp"=nhanes.vac$reth_hisp,
      "reth_black"=nhanes.vac$reth_black,
      "reth_asian"=nhanes.vac$reth_asian,
      "reth_other"=nhanes.vac$reth_other,
      "age1"=nhanes.vac$age_lt30,
      "age2"=nhanes.vac$age_3049,
      "age25"=nhanes.vac$age_lte25,
      "fborn"=nhanes.vac$fborn,
      "incpov1"=nhanes.vac$incpov_lt130,
      "incpov2"=nhanes.vac$incpov_na,
      "n"=nrow(nhanes.vac),
      "sens_a"=beta.sens.a,
      "sens_b"=beta.sens.b,
      "spec_a"=beta.spec.a,
      "spec_b"=beta.spec.b,
      "sens1_a"=beta.sens1.a,
      "sens1_b"=beta.sens1.b)
  return(jags.dat)
}

# Load data (2015-2016 NHANES cycle, ages 19+)
nhanes.vac <- readRDS("data/nhanes_vac3.rds")
