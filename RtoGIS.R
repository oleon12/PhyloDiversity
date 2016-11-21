# Load libraries
library(vegan)
library(maptools)
library(RColorBrewer)
library(classInt)
library(maps)
library(maptools)
library(dismo)

## Go to the work directory where are the shape poly 
setwd("/home/omar/Documentos/Omar/Tesis/Scripts/Distribution/shp/Grid/")
# Read the shape file
grid <- readShapePoly("grid_1g.shp")
#
# Now, go to the directory where are the index results
setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")
# Read the results depending on the cell size
x <- read.csv("AvTD.grid1",header=T)
# Extract the index values
area <- read.csv("DT.grid1")$area 
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

rawGrid <- read.csv("grid_1g.dist.matrix",header = T)

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


rawGrid <- rawGrid[-which(rownames(rawGrid)%in%out.cells),]

dist2 <- as.matrix(vegdist(rawGrid,method = 'jaccard',binary = T,na.rm = T))

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

grid$Index <- class

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

writePolyShape(grid,fn = "Grid1_AvTD")
