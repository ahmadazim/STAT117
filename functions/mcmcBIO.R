mcmcBIO = function(probeset,jagsmodel,small=FALSE){
  BIO = makeBIO(probeset,small)
  BIOmodelRjags = jags.model(textConnection(jagsmodel),
                             data = BIO,
                             n.chains = 4,
                             n.adapt = 100)
  print(jagsmodel)
  if (jagsmodel == BIOmodelRotatedJoint){
    BIOmcmc = coda.samples(BIOmodelRjags,c("mu","phi","pi","alpha","beta","mean.rat","mean.diff","precision.rat","precision.diff"),n.iter = 9000,thin=50)
    return(BIOmcmc)
  }
  # if (jagsmodel == BIOmodelRotatedJointTake1){
    # BIOmcmc = coda.samples(BIOmodelRjags,c("mu","phi","pi","alpha","beta","mean.center","mean.diff","precision.center","precision.diff"),n.iter = 9000,thin=50)
  # }
  if (jagsmodel == BIOmodelDirect){
    BIOmcmc = coda.samples(BIOmodelRjags,c("mu","phi","mean.mu","precision.mu"),n.iter = 9000,thin=100)
    return(BIOmcmc)
  }
  if (jagsmodel == BIOmodelRotated){
    BIOmcmc = coda.samples(BIOmodelRjags,c("mu","phi","mean.center","mean.diff","precision.center","precision.diff"),n.iter = 9000,thin=50)
    return(BIOmcmc)
  }
}
