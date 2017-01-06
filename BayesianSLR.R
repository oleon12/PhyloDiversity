setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

AllInd <- read.csv("AllIndex_Table.csv")


head(AllInd, 5L)

setwd("~/Documentos/Omar/Literature/Statistics/kruschke/")


fileNameRoot="SimpleLinearRegressionJags" # for constructing output filenames
source("openGraphSaveGraph.R")
source("plotPost.R")
require(rjags)         # Kruschke, J. K. (2011). Doing Bayesian Data Analysis:
# A Tutorial with R and BUGS. Academic Press / Elsevier.
#------------------------------------------------------------------------------
# THE MODEL.
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
# close quote for modelstring
writeLines(modelstring,con="model.txt")

#------------------------------------------------------------------------------
# THE DATA.

nSubj = NROW(AllInd)
x = AllInd$f.Rich
y = AllInd$DT
# Re-center data at mean, to reduce autocorrelation in MCMC sampling.
# Standardize (divide by SD) to make initialization easier.
xM = mean( x ) ; xSD = sd( x )
yM = mean( y ) ; ySD = sd( y )
zx = ( x - xM ) / xSD
zy = ( y - yM ) / ySD

# Specify data, as a list.
dataList = list(
  x = zx ,
  y = zy ,
  Ndata = nSubj
)

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.

r = cor(x,y)
initsList = list(
  beta0 = 0 ,    # because data are standardized
  beta1 = r ,        # because data are standardized
  tau = 1 / (1-r^2)  # because data are standardized
)

#------------------------------------------------------------------------------
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


# Display the posterior of the b1:
openGraph(10,4)
par( mar=c(4,4,1,1)+0.1 , mgp=c(2.5,0.8,0) )
layout( matrix(1:2,nrow=1) )
histInfo = plotPost( z1 , xlab="\nStandardized slope PR ~ TD"  )
histInfo = plotPost( b1 , xlab="Slope (pounds per inch)"  )
saveGraph(file=paste(fileNameRoot,"PostSlope",sep=""),type="eps")

# Display data with believable regression lines and posterior predictions.
openGraph()
par( mar=c(3,3,2,1)+0.5 , mgp=c(2.1,0.8,0) )
# Plot data values:
xRang = max(x)-min(x)
yRang = max(y)-min(y)
limMult = 0.25
xLim= c( min(x)-limMult*xRang , max(x)+limMult*xRang )
yLim= c( min(y)-limMult*yRang , max(y)+limMult*yRang )
plot( x , y , cex=1.5 , lwd=2 , col="black" , xlim=xLim , ylim=yLim ,
      xlab="Phylogenetic Richness" , ylab="TD" , cex.lab=1.5 ,
      cex.main=1.33  )
# Superimpose a smattering of believable regression lines:
for ( i in seq(from=1,to=length(b0),length=50) ) {
  abline( b0[i] , b1[i] , col="skyblue" )
}
saveGraph(file=paste(fileNameRoot,"DataLines",sep=""),type="eps")

#------------------------------------------------------------------------------

library(ggplot2)

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

ggplot(data = AllInd)+
  geom_point(aes(x=log(t.Rich),y=log(DT), colour="TD"))+
  geom_point(aes(x=log(t.Rich),y=log(PD), colour="PD"))+
  geom_point(aes(x=log(t.Rich),y=log(AvDT), colour="AvTD"))+
  geom_abline(slope = 0.973, intercept = .2, colour= "blue")+
  geom_abline(slope = 0.67, intercept = -.5, colour="green")+
  geom_abline(slope = 0.53, intercept =  -2, colour="red")



TDcor <- ggplot(data = AllInd)+
  geom_point(aes(x=log(t.Rich),y=log(DT)))+
  geom_abline(aes(slope = 1, intercept = 0, colour="Expected"), size=1)+
  geom_abline(aes(slope = .97, intercept= 0, colour="TR"), size=1)+
  geom_abline(aes(slope=.99, intercept = 0, colour= "PR"), size=1)+
  xlab("log Richness")+ylab("log TD")+
  geom_text(aes(x=5,y=.4),label="Slope PR = 0.99", size=10)+
  geom_text(aes(x=5,y=0),label="Slope TR = 0.97", size=10)+
  scale_color_discrete(name="Slopes")+
  theme(legend.position=c(.8,.3),
        legend.title=element_text(size = 30,face = "bold"),
        legend.text = element_text(size = 30),
        legend.key.size = unit(1, "cm"),
        legend.key = element_rect(fill= "white",colour="black"),
        axis.text.y=element_text(size=15),
        axis.text.x=element_blank(),
        axis.title.y=element_text(size=30),
        axis.title.x=element_blank())+  
  guides(colour = guide_legend(override.aes = list(size=2)))


png("TD_Corr.png", width =2000 , height =1200 ,res = 200)
TDcor
dev.off()

PDcor <- ggplot(data = AllInd)+
  geom_point(aes(x=log(t.Rich),y=log(PD)))+
  geom_abline(aes(slope = 1, intercept = 0, colour="Expected"), size=1)+
  geom_abline(aes(slope = .67, intercept=0, colour="TR"), size=1)+
  geom_abline(aes(slope=.65, intercept = 0, colour= "PR"), size=1)+
  xlab("log Richness")+ylab("log PD")+
  geom_text(aes(x=5,y=-4.3),label="Slope PR = 0.65", size=10)+
  geom_text(aes(x=5,y=-5),label="Slope TR = 0.67", size=10)+
  scale_color_discrete(name="Slopes")+
  theme(legend.position=c(.8,.3),
        legend.title=element_text(size = 30,face = "bold"),
        legend.text = element_text(size = 30),
        legend.key.size = unit(1, "cm"),
        legend.key = element_rect(fill= "white",colour="black"),
        axis.text.y=element_text(size=15),
        axis.text.x=element_blank(),
        axis.title.y=element_text(size=30),
        axis.title.x=element_blank())+  
  guides(colour = guide_legend(override.aes = list(size=2)))


png("PD_Corr.png", width =2000 , height =1200 ,res = 200)
PDcor
dev.off()

AvTDcor <- ggplot(data = AllInd)+
  geom_point(aes(x=log(t.Rich),y=log(AvDT)))+
  geom_abline(aes(slope = 1, intercept = 0, colour="Expected"), size=1)+
  geom_abline(aes(slope = .53, intercept=0, colour="TR"), size=1)+
  geom_abline(aes(slope=.55, intercept = 0, colour= "PR"), size=1)+
  xlab("log Richness")+ylab("log AvTD")+
  geom_text(aes(x=5,y=-5.3),label="Slope PR = 0.55", size=10)+
  geom_text(aes(x=5,y=-6),label="Slope TR = 0.53", size=10)+
  scale_color_discrete(name="Slopes")+
  theme(legend.position=c(.8,.3),
        legend.title=element_text(size = 30,face = "bold"),
        legend.text = element_text(size = 30),
        legend.key.size = unit(1, "cm"),
        legend.key = element_rect(fill= "white",colour="black"),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15),
        axis.title.y=element_text(size=30),
        axis.title.x=element_text(size=30))+  
  guides(colour = guide_legend(override.aes = list(size=2)))


png("AvTD_Corr.png", width =2000 , height =1300 ,res = 200)
AvTDcor
dev.off()
