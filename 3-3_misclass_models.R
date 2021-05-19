# JAGS Models for 3_misclass.R
# sn0 and sp0 are fed into the model in the data.

### Primary Model ##############################################################
bugsmod <-
  "
model {
  for(i in 1:n){
    # outcome model, log odds of hep b vaccination given predictors
    vac[i] ~ dbern(p_vac[i])
    logit(p_vac[i])<- b0+b[1]*reth_hisp[i]+
      b[2]*reth_black[i]+b[3]*reth_asian[i]+b[4]*reth_other[i]+b[5]*age1[i]+
      b[6]*age2[i]+b[7]*fborn[i]+b[8]*incpov1[i]+b[9]*incpov2[i]
    
    # exposure models
    srvac1[i] ~ dbern(p_srvac[i])
    logit(p_srvac[i]) <- a0+a[1]*reth_hisp[i]+
      a[2]*reth_black[i]+a[3]*reth_asian[i]+a[4]*reth_other[i]+a[5]*age1[i]+
      a[6]*age2[i]+a[7]*fborn[i]+a[8]*incpov1[i]+a[9]*incpov2[i]
    
    # measurement model, hepatitis b vaccination (vac)  is imperfectly measured 
      # by vac_star
    vac_star[i] ~ dbern(p_vac_star[i])
    p_vac_star[i] <- age25[i]*(sn0*vac[i]+(1-vac[i])*(1-sp0)) + 
                      (1-age25[i])*(sn1*vac[i]+(1-vac[i])*(1-sp1))
    
  }

  # Outcome model priors
  b0 ~ dnorm(0,.0001)
  for(j in 1:9){
    b[j] ~ dnorm(0,.0001)
  }
  
  a0 ~ dnorm(0,.0001)
  # Srvac model priors
  for(k in 1:9){
    a[k] ~ dnorm(0,.0001)
  }
  
  # Measurement model priors
  sn0~dbeta(sens_a,sens_b)
  sp0~dbeta(spec_a,spec_b)
  
  sn1~dbeta(sens1_a,sens1_b)
  sp1~dbeta(spec_a,spec_b) # Same prior dist for this spec
  
  m_pvac <- mean(p_vac)
  m_pvacstar <- mean(p_vac_star)
  m_psrvac <- mean(p_srvac)
}
"
# write model to file
writeLines(bugsmod, "m1.txt") 
