BayesSlope <- function(x,y,nSubj,ploting=c(TRUE,FALSE),SaveName, Xlab){
  
  #########################
  setwd("~/Documentos/Omar/Literature/Statistics/kruschke/")
  
  source("plotPost.R")
  source("openGraphSaveGraph.R")
  require(rjags)
  #########################
  
  modelstring = "
  model {
  for( i in 1 : Ndata ) {
     y[i] ~ dnorm( mu[i] , tau )
     mu[i] <- beta0 + beta1 * x[i]
    }
  beta0 ~ dnorm( 0 , 1.0E-12 )
  beta1 ~ dnorm( 0 , 1.0E-12 )
  tau ~ dgamma( 0.001 , 0.001 )
  }
  " 

  writeLines(modelstring,con="model.txt")
  
  #########################
  
  xM = mean( x ) ; xSD = sd( x )
  yM = mean( y ) ; ySD = sd( y )
  zx = ( x - xM ) / xSD
  zy = ( y - yM ) / ySD
  
  dataList = list(
    x = zx ,
    y = zy ,
    Ndata = nSubj
  )
  
  #########################
  r = cor(x,y)
  initsList = list(
    beta0 = 0 ,    # because data are standardized
    beta1 = r ,        # because data are standardized
    tau = 1 / (1-r^2)  # because data are standardized
  )
  #########################
  
  # RUN THE CHAINS
  
  parameters = c("beta0" , "beta1" , "tau")  # The parameter(s) to be monitored.
  adaptSteps = 500              # Number of steps to "tune" the samplers.
  burnInSteps = 500            # Number of steps to "burn-in" the samplers.
  nChains = 3                   # Number of chains to run.
  numSavedSteps=50000           # Total number of steps in chains to save.
  thinSteps=1                   # Number of steps to "thin" (1=keep every step).
  nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
  
  # Create, initialize, and adapt the model:
  jagsModel = jags.model( "model.txt" , data=dataList , inits=initsList ,
                          n.chains=nChains , n.adapt=adaptSteps )
  
  
  
  # Burn-in:
  cat( "Burning in the MCMC chain...\n" )
  update( jagsModel , n.iter=burnInSteps )
  # The saved MCMC chain:
  cat( "Sampling final MCMC chain...\n" )
  codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                              n.iter=nPerChain , thin=thinSteps )
  
  str(codaSamples)
  
  head(codaSamples[[1]])
  
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  
  #------------------------------------------------------------------------------
  # EXAMINE THE RESULTS
  
  checkConvergence = FALSE
  if ( checkConvergence ) {
    openGraph(width=7,height=7)
    autocorr.plot( codaSamples[[1]] , ask=FALSE ) 
    show( gelman.diag( codaSamples ) )
    effectiveChainLength = effectiveSize( codaSamples ) 
    show( effectiveChainLength )
  }
  
  # Convert coda-object codaSamples to matrix object for easier handling.
  # But note that this concatenates the different chains into one long chain.
  # Result is mcmcChain[ stepIdx , paramIdx ]
  mcmcChain = as.matrix( codaSamples )
  
  head(mcmcChain)
  
  # Extract chain values:
  z0 = mcmcChain[, "beta0" ]
  z1 = mcmcChain[, "beta1" ]
  zTau = mcmcChain[, "tau" ]
  zSigma = 1 / sqrt( zTau ) # Convert precision to SD
  
  # Convert to original scale:
  b1 = z1 * ySD / xSD
  b0 = ( z0 * ySD + yM - z1 * ySD * xM / xSD )
  sigma = zSigma * ySD
  
  # Posterior prediction:
  # Specify x values for which predicted y's are needed:
  xPostPred = seq(55,80,1)
  # Define matrix for recording posterior predicted y values at each x value.
  # One row per x value, with each row holding random predicted y values.
  postSampSize = length(b1)
  yPostPred = matrix( 0 , nrow=length(xPostPred) , ncol=postSampSize )
  # Define matrix for recording HDI limits of posterior predicted y values:
  yHDIlim = matrix( 0 , nrow=length(xPostPred) , ncol=2 )
  # Generate posterior predicted y values.
  # This gets only one y value, at each x, for each step in the chain.
  for ( chainIdx in 1:postSampSize ) {
    yPostPred[,chainIdx] = rnorm( length(xPostPred) ,
                                  mean = b0[chainIdx] + b1[chainIdx] * xPostPred ,
                                  sd = rep( sigma[chainIdx] , length(xPostPred) ) )
  }
  
  
  source("HDIofMCMC.R")
  for ( xIdx in 1:length(xPostPred) ) {
    yHDIlim[xIdx,] = HDIofMCMC( yPostPred[xIdx,] )
  }
  
  
  #####################################################################
  
  histInfo = plotPost( z1 , xlab = Xlab)
  
  if(ploting==T){
    
    setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")
    
    saveGraph(file=SaveName,type="eps")
  
  }
  
  return(histInfo[1])
  
}

