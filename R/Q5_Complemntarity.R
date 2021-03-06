library(maptools)
library(classInt)
library(vegan)
library(shapefiles)


setwd("/home/omar/Documentos/Omar/Tesis/Scripts/Distribution/shp/Grid/")

grid <- readShapePoly("grid_25g.shp")

grid2 <- read.shp("grid_25g.shp")

setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

x <- read.csv("AvTD.grid25",header=T)

area <- read.csv("DT.grid25")$area 
x <- x$Dplus

x

x2 <- x[which(x>0.000000)]

brks <- classIntervals(x2,n=5,style = "quantile")
brks <- brks$brks

brks

class <- findInterval(x,brks,all.inside = T)

class[grep(1,class)] <- "Q1"
class[grep(2,class)] <- "Q2"
class[grep(3,class)] <- "Q3"
class[grep(4,class)] <- "Q4"
class[grep(5,class)] <- "Q5"

class

Q5pos <- grep("Q5",class)

InitClass <- class

#######################################################

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/")

rawGrid <- read.csv("grid_25g.dist.matrix",header = T)

sp <- rawGrid$especie

rawGrid <- rawGrid[,-1]
rownames(rawGrid) <- sp

rawGrid <- t(rawGrid)

#######

out.cells <- c()

sp.len <- length(rawGrid[1,])

for(i in 1:length(rownames(rawGrid))){
  
  len.c <- length(grep(0,rawGrid[i,]))
  
  if(len.c==sp.len){
    
    out.cells <- c(out.cells,rownames(rawGrid)[i])
    
  }
  
}

class2 <- class[-which(rownames(rawGrid)%in%out.cells)]

classOut <- which(rownames(rawGrid)%in%out.cells)

Q5names <- rownames(rawGrid)[Q5pos]

rawGrid <- rawGrid[-which(rownames(rawGrid)%in%out.cells),]


colnames(rawGrid)
#rownames(rawGrid) <- class2
rownames(rawGrid)

########################

dist2 <- as.matrix(vegdist(rawGrid,method = 'jaccard',binary = T,na.rm = T))

dist3 <- dist2[which(rownames(dist2)%in%Q5names),-which(rownames(dist2)%in%Q5names)]

SumQ5 <- apply(dist3, 2, sum)

SumQ5 <- SumQ5/length(dist3[,1])

SumQ5

brks2 <- classIntervals(SumQ5,n=4,style = "quantile")
brks2 <- brks2$brks

brks2

Allclass <- findInterval(SumQ5,brks2,all.inside = T)

names(Allclass) <- names(SumQ5)

Allclass

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/")

rawGrid <- read.csv("grid_25g.dist.matrix",header = T)

rawGrid <- rawGrid[,-1]

ComVal <- rep(NA, length(colnames(rawGrid)))

names(ComVal) <- colnames(rawGrid)


ComVal[which(names(ComVal)%in%names(Allclass))] <- Allclass

ComVal[which(is.na(ComVal)==TRUE)] <- rep("NO", length(which(is.na(ComVal)==TRUE)))

ComVal

###############################

grid$CompVal <- ComVal

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

writePolyShape(grid,fn = "Q5vsAll-CC.shp")

