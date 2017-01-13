library(ggplot2)
library(vegan)
library(shapefiles)
library(maptools)
library(Imap)

## Go to the work directory where are the shape poly 
setwd("/home/omar/Documentos/Omar/Tesis/Scripts/Distribution/shp/Grid/")
# Read the shape file
grid <- readShapePoly("grid_25g.shp")
#
# Now, go to the directory where are the index results
setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")
# Read the results depending on the cell size
x <- read.csv("AvTD.grid25",header=T)
# Extract the index values
area <- read.csv("DT.grid25")$area 
x <- x$Dplus
# Check the values
x
# Due to there are many cell without index values (this is beacuse the species distribution)
# remove the 0 values, or the quantile classification will be biased
x2 <- x[-grep(0,x)]

# No, do the quantile clasification, in this case, given the values index, the intervals will be generated
brks <- classIntervals(x2,n=5,style = "quantile")
brks <- brks$brks
# Here, the intervarls given the index value
brks

# Now, each index values will be classified in a category given the quantile intervals
# Here, the all index values (include those cells with index values equal 0)
# Just, beacuase the final shape file neead a value for each cell
class <- findInterval(x,brks,all.inside = T)
# Now, beacuase the focus is only in the quantile 5 and 4, (those with highest index values)
# The cells out of the Q5 and Q4 will be replace with NO
class[-c(grep(5,class),grep(4,class))] <-  "NO"
# Now, those cell int in the two last quantiles, will be filled with the quotes "Q5" and "Q4"
class[grep(5,class)] <- "Q5"
class[grep(4,class)] <- "Q4"
# See the result
class
#
# Now, the hard part...
# 

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/")

rawGrid <- read.csv("grid_25g.dist.matrix",header = T)

sp <- rawGrid$especie

rawGrid <- rawGrid[,-1]
rownames(rawGrid) <- sp

rawGrid <- t(rawGrid)

out.cells <- c()

sp.len <- length(rawGrid[1,])

for(i in 1:length(rownames(rawGrid))){
  
  len.c <- length(grep(0,rawGrid[i,]))
  
  if(len.c==sp.len){
    
    out.cells <- c(out.cells,rownames(rawGrid)[i])
    
  }
  
}


rawGrid2 <- rawGrid[-which(rownames(rawGrid)%in%out.cells),]

dist2 <- as.matrix(vegdist(rawGrid2,method = 'jaccard',binary = T,na.rm = T))

dist2 <- dist2[which(rownames(dist2)%in%area[grep("Q5",class)]), ]

brks2 <- c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1)

for(i in 1:length(rownames(dist2))){
  
  class.v <- findInterval(dist2[i,], brks2, all.inside = T)
  
  dist2[i, ] <- class.v
}

dist2

CI.cells <- c()

l2 <- length(dist2[,1])

for (i in 1:length(colnames(dist2))){
  
  l1 <- length(grep(10,dist2[ ,i]))
  
  if(l1==l2){
    CI.cells <- c(CI.cells,colnames(dist2)[i])
  }
  
}

dist.final <- dist2[,which(colnames(dist2)%in%CI.cells)]

class[which(area%in%colnames(dist.final))]

class[which(area%in%colnames(dist.final))] <- "CI.1"

class

#########################################################################

rawGrid3 <- rawGrid

rownames(rawGrid3) <- class

Q5CC.Grid <- rawGrid3[c(grep("Q5",class),grep("CI.1",class)),]

rownames(Q5CC.Grid)

#########################################################################

modulus <- function(Val1, Val2){
  
  t1<-floor(Val1/Val2)
  return(Val1-t1*Val2)
  
  
}


XYto1 <- function(X,Y, maxY){
  
  One <- X * maxY + Y
  
  return(One)
}

CoordGrid <- matrix(0, ncol = 2, nrow = length(grid$ID))
colnames(CoordGrid) <- c("X","Y")

rownames(CoordGrid) <- class

for(i in 1:length(CoordGrid[,1])){
  
  CoordGrid[i,1] <- as.numeric(grid2$shp[[i]]$points[1, 1])
  
  CoordGrid[i,2] <- as.numeric(grid2$shp[[i]]$points[1, 2])
  
}

head(CoordGrid)
tail(CoordGrid)

####################################################

Xm <- mean(CoordGrid[,1])
Xsd <- sd(CoordGrid[,1])

Ym <- mean(CoordGrid[,2])
Ysd <- sd(CoordGrid[,2])

Zx <- (CoordGrid[,1]-Xm)/Xsd

Zy <- (CoordGrid[,2]-Ym)/Ysd

maxY <- max(Zy)

UniqVal <- c()

for(i in 1:length(Zx)){
  
  UniqTemp <- XYto1(Zx[i],Zy[i],maxY)
  
  UniqVal <- c(UniqVal, UniqTemp)
  
}

CoordGrid <- as.data.frame(cbind(CoordGrid,UniqVal))

head(CoordGrid)

UniqValData <- as.matrix(UniqVal)

Q5CC.Uniq <- as.matrix(UniqValData[c(grep("Q5",class),grep("CI.1",class)),])

length(Q5CC.Uniq[,1])

length(Q5CC.Grid[,1])

CoordDist <- as.matrix(vegdist(Q5CC.Uniq, method = "manhattan"))

IndDist <- as.matrix(vegdist(Q5CC.Grid, method = "jaccard", binary = T, na.rm = T))


length(CoordDist[,1])

length(IndDist[,1])

mantel(CoordDist,IndDist)

mantel(IndDist, CoordDist)


CoorDist2 <- CoordDist[grep("CI.1",rownames(CoordDist)),grep("Q5",colnames(CoordDist))]

IndDist2 <- IndDist[grep("CI.1",rownames(IndDist)),grep("Q5",colnames(IndDist))]

length(CoorDist2[,1])
length(IndDist2[,1])


PostBayesSlope <- c()

for(i in 134:length(colnames(CoorDist2))){
  
  print(paste(paste("Iteration",i,sep = " "), paste("of", length(colnames(CoorDist2)), sep=" "),sep=" "))
  
  PBS <- BayesSlope(x = CoorDist2[,i], y= IndDist2[,i], nSubj = length(IndDist2[,1]))
  
  PostBayesSlope <- c(PostBayesSlope, PBS)

}

PostBayesSlope

hist(PostBayesSlope)

ggplot()+geom_histogram(aes(PostBayesSlope))


###########################################################################
###########################################################################

rawGrid3 <- rawGrid

rownames(rawGrid3) <- class

Q5CC.Grid <- rawGrid3[c(grep("Q5",class),grep("CI.1",class)),]

rownames(Q5CC.Grid)

#########################################################################

modulus <- function(Val1, Val2){
  
  t1<-floor(Val1/Val2)
  return(Val1-t1*Val2)
  
  
}


CoordGrid <- matrix(0, ncol = 2, nrow = length(grid$ID))
colnames(CoordGrid) <- c("X","Y")

rownames(CoordGrid) <- class

for(i in 1:length(CoordGrid[,1])){
  
  CoordGrid[i,1] <- as.numeric(grid2$shp[[i]]$points[1, 1])
  
  CoordGrid[i,2] <- as.numeric(grid2$shp[[i]]$points[1, 2])
  
}

head(CoordGrid)
tail(CoordGrid)

####################################################


Q5CC.Uniq <- as.matrix(CoordGrid[c(grep("Q5",class),grep("CI.1",class)),])

length(Q5CC.Uniq[,1])

length(Q5CC.Grid[,1])


CoorDist <- matrix(0, ncol = length(Q5CC.Uniq[,1]),nrow=length(Q5CC.Uniq[,1]))


for(i in 1:length(Q5CC.Uniq[,1])){
  
  for(j in 1:length(Q5CC.Uniq[,1])){
    
    distVal <- gdist(lon.1 =Q5CC.Uniq[i,1] ,
                     lat.1 =Q5CC.Uniq[i,2] ,
                     lon.2 = Q5CC.Uniq[j,1],
                     lat.2 = Q5CC.Uniq[j,2])
    
    CoorDist[i,j] <- distVal
    
  }
  
}

colnames(CoorDist) <- rownames(Q5CC.Uniq)
rownames(CoorDist) <- rownames(Q5CC.Uniq)

CoorDist

IndDist <- as.matrix(vegdist(Q5CC.Grid, method = "jaccard", binary = T, na.rm = T))

length(CoordDist[,1])

length(IndDist[,1])

mantel(CoordDist,IndDist)

mantel(IndDist, CoordDist)

CoorDist2 <- CoordDist[grep("CI.1",rownames(CoordDist)),grep("Q5",colnames(CoordDist))]

IndDist2 <- IndDist[grep("CI.1",rownames(IndDist)),grep("Q5",colnames(IndDist))]

length(CoorDist2[,1])
length(IndDist2[,1])

PostBayesSlope2 <- c()

for(i in 134:length(colnames(CoorDist2))){
  
  print(paste(paste("Iteration",i,sep = " "), paste("of", length(colnames(CoorDist2)), sep=" "),sep=" "))
  
  PBS <- BayesSlope(x = CoorDist2[,i], y= IndDist2[,i], nSubj = length(IndDist2[,1]))
  
  PostBayesSlope2 <- c(PostBayesSlope2, PBS)
  
}

PostBayesSlope2

mean(PostBayesSlope2)

hist(PostBayesSlope2)

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

png("CIvsDist.png", width = 1500, height = 800)
ggplot()+geom_histogram(aes(PostBayesSlope2))+
         ylab(NULL)+xlab("Posterior Slope Means")+
         geom_text(aes(.25,15), label=paste("Mean",round(mean(PostBayesSlope2),digits=3),sep = ": "), size=15)+
         geom_text(aes(.25,14), label=paste("Max",round(max(PostBayesSlope2),digits=3),sep = ": "), size=15)+
         geom_text(aes(.25,13), label=paste("Min",round(min(PostBayesSlope2),digits=3),sep = ": "), size=15)+
  theme(axis.text.x=element_text(size=30),
        axis.text.y=element_text(size=30),
        axis.ticks=element_blank(),
        axis.title.y=element_text(size=35),
        axis.title.x=element_text(size=35))
dev.off()
