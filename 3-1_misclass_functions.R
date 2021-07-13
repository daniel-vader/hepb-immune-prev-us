### Init functions for poststratification ###

# Setup data
nhanes.vac.poststrat <- readRDS("data/nhanes_vac3.rds")
nhanes.vac.poststrat$incpov1 <- ifelse(is.na(nhanes.vac.poststrat$incpov_lt130), 0, 
                                       nhanes.vac.poststrat$incpov_lt130)
nhanes.vac.poststrat$incpov2 <- ifelse(is.na(nhanes.vac.poststrat$incpov_lt130), 1, 0)
nhanes.vac.poststrat$srvac1 <- nhanes.vac.poststrat$srvac_gt3
nhanes.vac.poststrat$srvac2 <- nhanes.vac.poststrat$srvac_dnr
nhanes.vac.poststrat$age1 <- nhanes.vac.poststrat$age_lt30
nhanes.vac.poststrat$age2 <- nhanes.vac.poststrat$age_3049
nhanes.vac.poststrat$age25 <- nhanes.vac.poststrat$age_lte25
nhanes.vac.poststrat <- nhanes.vac.poststrat %>% 
  dplyr::select(srvac1, srvac2, male, reth_hisp, reth_black, reth_asian, reth_other, 
                age1, age2, age25, incpov1, incpov2, fborn)
#ps.demo2 <- readRDS("data/psdemo.rds")
ps.demo2 <- readRDS("data/psdemo2.rds") %>%
  mutate(age1 = ifelse(age18to25 + age26to29 > 0, 1, 0),
         age2 = ifelse(age30to49 == 1, 1, 0),
         age25 = ifelse(age18to25 == 1, 1, 0)) %>%
  select(-age18to25, -age30to49, -age26to29)

# Assign strata identifier and strata probability weights to data
nvp <- dplyr::left_join(nhanes.vac.poststrat, ps.demo2) %>% 
  dplyr::select(-n)

# Set up filtering functions to automate subgrouping
init.filter.func <- function(){
  filter.white <- function(dat){
    dat[dat$reth_black==0 & dat$reth_asian==0 & dat$reth_other==0 & dat$reth_hisp==0,]
  }
  filter.black <- function(dat){
    dat[dat$reth_black==1,]
  }
  filter.asian <- function(dat){
    dat[dat$reth_asian==1,]
  }
  filter.other <- function(dat){
    dat[dat$reth_other==1,]
  }
  filter.hisp <- function(dat){
    dat[dat$reth_hisp==1,]
  }
  filter.age1 <- function(dat){
    dat[dat$age1==1,]
  }
  filter.age2 <- function(dat){
    dat[dat$age2==1,]
  }
  filter.age3 <- function(dat){
    dat[dat$age1==0 & dat$age2==0,]
  }
  filter.fb <- function(dat){
    dat[dat$fborn==1,]
  }
  filter.fbn <- function(dat){
    dat[dat$fborn==0,]
  }
  
  # put filtering functions in a list
  filter.list <- list(function(x) filter.white(x), function(x) filter.black(x), 
                      function(x) filter.asian(x), 
                      function(x) filter.other(x), function(x) filter.hisp(x),
                      function(x) filter.age1(x), function(x) filter.age2(x),
                      function(x) filter.age3(x),
                      function(x) filter.fb(x), function(x) filter.fbn(x))
  return(filter.list)
}

# set vector of values to be used to rescale subgroup weights
init.scales <- function(filter.list, dat){
  scales <- 1:length(filter.list)
  for(j in 1:length(filter.list)){
    scales[j] <- filter.list[[j]](dat) %>% dplyr::select(pwt) %>% sum()
  }
  return(scales)
}

# Set up dataframes to store poststratification prev estimates
init.post.frames <- function() {
  vac.marginal <- data.frame(p_vac=numeric(), ps_vac=numeric(), p_vac_star=numeric(),
                             ps_vac_star=numeric(), p_srvac=numeric(),
                             ps_srvac=numeric(), tweight=numeric())
  vac.white <- vac.marginal
  vac.black <- vac.marginal
  vac.asian <- vac.marginal
  vac.other <- vac.marginal
  vac.hisp <- vac.marginal
  vac.age1 <- vac.marginal
  vac.age2 <- vac.marginal
  vac.age3 <- vac.marginal
  vac.fb <- vac.marginal
  vac.fbn <- vac.marginal
  
  # Put dataframes in list
  ps.dat.list <- list(vac.marginal, vac.white, vac.black, vac.asian, vac.other, 
                      vac.hisp, vac.age1, vac.age2, vac.age3, vac.fb, vac.fbn)
  names(ps.dat.list) <- c("marginal", "white", "black", "asian", "other", "hisp", 
                          "age1", "age2", "age3", "fb", "fbn")
  return(ps.dat.list)
}

# Calculate poststratified estimates of prevalence in each simulation, both marginal and conditional.
library(svMisc) # Load library that lets is have a progress bar
calc.poststrat.prev <- function(dat){
  filter.list <- init.filter.func()
  ps.dat.list <- init.post.frames() # initialized list of poststratified data frames
  scales <- init.scales(filter.list, ps.demo2)
  nsims <- length(dat$a[,1])
  
  for(i in 1:nsims){
    progress(i, nsims)
    tmp <- nvp
    
    # self-reported vaccination
    srvac.odds <- exp(
      dat$a0[i] + 
        tmp$reth_hisp*dat$a[i,1] + 
        tmp$reth_black*dat$a[i,2] +
        tmp$reth_asian*dat$a[i,3] +
        tmp$reth_other*dat$a[i,4] +
        tmp$age1*dat$a[i,5] +
        tmp$age2*dat$a[i,6] +
        tmp$fborn*dat$a[i,7] +
        tmp$incpov1*dat$a[i,8] +
        tmp$incpov2*dat$a[i,9]
    )
    tmp$psrvac <- srvac.odds / (1 + srvac.odds)
    
    # immunity
    vac.odds <- exp(
      dat$b0[i] + 
        tmp$reth_hisp*dat$b[i,1] + 
        tmp$reth_black*dat$b[i,2] +
        tmp$reth_asian*dat$b[i,3] +
        tmp$reth_other*dat$b[i,4] +
        tmp$age1*dat$b[i,5] +
        tmp$age2*dat$b[i,6] +
        tmp$fborn*dat$b[i,7] +
        tmp$incpov1*dat$b[i,8] +
        tmp$incpov2*dat$b[i,9] +
        tmp$psrvac*dat$b[i,10]
    )
    tmp$pvac <- vac.odds / (1 + vac.odds)
    
    # anti-HBs
    tmp$pvacstar <- tmp$age25*(dat$sn0[i]*tmp$pvac + (1-tmp$pvac)*(1-dat$sp0[i])) +
      (1-tmp$age25)*(dat$sn1[i]*tmp$pvac + (1-tmp$pvac)*(1-dat$sp1[i]))
    
    
    # Calculate poststratified p_vac, p_vac_star, and psrvac
    tmp2 <- tmp %>% dplyr::group_by(strata) %>% 
      dplyr::summarize(pwt = mean(pwt), #pwt not a mean, just a way to reduce dimension
                       ps_vac = mean(pvac)*mean(pwt),
                       ps_vac_star = mean(pvacstar)*mean(pwt),
                       ps_srvac = mean(psrvac)*mean(pwt),
                       .groups="drop") %>%
      ungroup()
    ps.dat.list[[1]] <- 
      dplyr::add_row(ps.dat.list[[1]], 
                     p_vac = mean(tmp$pvac),
                     ps_vac = sum(tmp2$ps_vac),
                     p_vac_star = mean(tmp$pvacstar), 
                     ps_vac_star = sum(tmp2$ps_vac_star),
                     p_srvac = mean(tmp$psrvac),
                     ps_srvac = sum(tmp2$ps_srvac),
                     tweight = sum(tmp2$pwt))
    
    # Calculate poststratified p_vac, p_vac_star, and srvac by demographic subgroup
    for(j in 1:length(filter.list)){
      tmp2 <- filter.list[[j]](tmp) 
      tmp3 <- tmp2 %>% dplyr::group_by(strata) %>% 
        dplyr::summarize(pwt = mean(pwt/scales[j]), #pwt not a mean, just a way to reduce dimension
                         ps_vac = mean(pvac)*mean(pwt),
                         ps_vac_star = mean(pvacstar)*mean(pwt),
                         ps_srvac = mean(psrvac)*mean(pwt),
                         .groups="drop"
        ) %>% 
        ungroup()
      
      
      ps.dat.list[[j+1]] <- 
        dplyr::add_row(ps.dat.list[[j+1]], 
                       p_vac = mean(tmp2$pvac),
                       ps_vac = sum(tmp3$ps_vac), 
                       p_vac_star = mean(tmp2$pvacstar), 
                       ps_vac_star = sum(tmp3$ps_vac_star), 
                       p_srvac = mean(tmp2$psrvac),
                       ps_srvac = sum(tmp3$ps_srvac),
                       tweight = sum(tmp2$pwt))
    }
    
  }
  return(ps.dat.list)
}

