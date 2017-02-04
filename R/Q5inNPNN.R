# Load libraries
library(vegan)
library(maptools)
library(RColorBrewer)
library(classInt)
library(maps)
library(maptools)
library(dismo)
library(shapefiles)
library(SDMTools)
library(doParallel)

setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")
# Read the shape file
grid <- readShapePoly("Grid25_AvTD.shp")
pos.Q5 <- which(grid$Index=="Q5")
pos.Q5 <- pos.Q5 + 1
#
setwd("~/Documentos/Omar/Tesis/Scripts/Distribution/shp/Grid/")
grid2 <- read.shp("grid_25g.shp")

setwd("~/Documentos/Omar/Tesis/Scripts/Distribution/shp/NAB/")
PNN.names <- readShapePoly("PNN2_NAB.shp")
PNN.names <- PNN.names$NAME
PNN <- read.shp("PNN2_NAB.shp")

CellPerc <- matrix(0,ncol = length(PNN.names), nrow = 1)

for(i in 1:length(pos.Q5)){
  
  pos <- pos.Q5[i]
  
  print(paste(i, paste("of",length(pos.Q5), sep=" "), sep = " "))
  
  cell <- grid2$shp[[pos]]$points
  
  for (j in 1:length(PNN.names)) {
    
    PNN.coord <- PNN$shp[[j]]$points
    plot(PNN.coord)
    
    pnt <- pnt.in.poly(pnts = cell, poly.pnts = PNN.coord)
    
    n <- grep(1,pnt$pip)
    
    if(length(n)==5){
      
      CellPerc[1,j] <- (CellPerc[1,j]+1)
      plot(cell)
    }
    
    if(length(n)==4){
      
      CellPerc[1,j] <- (CellPerc[1,j]+1)
      plot(cell)
    }
    
    if(length(n)==3){
      
      CellPerc[1,j] <- (CellPerc[1,j]+0.75)
      
    }
    
    if(length(n)==2){
      
      CellPerc[1,j] <- (CellPerc[1,j]+0.50)
      
    }
    
    if(length(n)==1){
      
      CellPerc[1,j] <- (CellPerc[1,j]+0.25)
      
    }
    
  }
}

sum(CellPerc)/length(PNN.names)

###############################################################################

setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")
# Read the shape file
grid <- readShapePoly("Grid25_PD.shp")
pos.Q5 <- which(grid$Index=="Q5")
pos.Q5 <- pos.Q5 + 1
#
setwd("~/Documentos/Omar/Tesis/Scripts/Distribution/shp/Grid/")
grid2 <- read.shp("grid_25g.shp")

setwd("~/Documentos/Omar/Tesis/Scripts/Distribution/shp/NAB/")
PNN.names <- readShapePoly("PNN2_NAB.shp")
PNN.names <- PNN.names$NAME
PNN <- read.shp("PNN2_NAB.shp")

CellPerc <- matrix(0,ncol = length(PNN.names), nrow = 1)

cl <- makeCluster(8)
registerDoParallel(cl)

for(i in 1:length(pos.Q5)){
  
  pos <- pos.Q5[i]
  
  print(paste(i, paste("of",length(pos.Q5), sep=" "), sep = " "))
  
  cell <- grid2$shp[[pos]]$points
  
  for (j in 1:length(PNN.names)) {
    
    PNN.coord <- PNN$shp[[j]]$points
    plot(PNN.coord)
    
    pnt <- pnt.in.poly(pnts = cell, poly.pnts = PNN.coord)
    
    n <- grep(1,pnt$pip)
    
    if(length(n)==5){
      
      CellPerc[1,j] <- (CellPerc[1,j]+1)
      plot(cell)
    }
    
    if(length(n)==4){
      
      CellPerc[1,j] <- (CellPerc[1,j]+1)
      plot(cell)
    }
    
    if(length(n)==3){
      
      CellPerc[1,j] <- (CellPerc[1,j]+0.75)
      
    }
    
    if(length(n)==2){
      
      CellPerc[1,j] <- (CellPerc[1,j]+0.50)
      
    }
    
    if(length(n)==1){
      
      CellPerc[1,j] <- (CellPerc[1,j]+0.25)
      
    }
    
  }
}

stopCluster(cl)

sum(CellPerc)/length(PNN.names)

###############################################################################

setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")
# Read the shape file
grid <- readShapePoly("Grid25_TD.shp")
pos.Q5 <- which(grid$Index=="Q5")
pos.Q5 <- pos.Q5 + 1
#
setwd("~/Documentos/Omar/Tesis/Scripts/Distribution/shp/Grid/")
grid2 <- read.shp("grid_25g.shp")

setwd("~/Documentos/Omar/Tesis/Scripts/Distribution/shp/NAB/")
PNN.names <- readShapePoly("PNN2_NAB.shp")
PNN.names <- PNN.names$NAME
PNN <- read.shp("PNN2_NAB.shp")

CellPerc <- matrix(0,ncol = length(PNN.names), nrow = 1)

cl <- makeCluster(8)
registerDoParallel(cl)

for(i in 1:length(pos.Q5)){
  
  pos <- pos.Q5[i]
  
  print(paste(i, paste("of",length(pos.Q5), sep=" "), sep = " "))
  
  cell <- grid2$shp[[pos]]$points
  
  for (j in 1:length(PNN.names)) {
    
    PNN.coord <- PNN$shp[[j]]$points
    plot(PNN.coord)
    
    pnt <- pnt.in.poly(pnts = cell, poly.pnts = PNN.coord)
    
    n <- grep(1,pnt$pip)
    
    if(length(n)==5){
      
      CellPerc[1,j] <- (CellPerc[1,j]+1)
      plot(cell)
    }
    
    if(length(n)==4){
      
      CellPerc[1,j] <- (CellPerc[1,j]+1)
      plot(cell)
    }
    
    if(length(n)==3){
      
      CellPerc[1,j] <- (CellPerc[1,j]+0.75)
      
    }
    
    if(length(n)==2){
      
      CellPerc[1,j] <- (CellPerc[1,j]+0.50)
      
    }
    
    if(length(n)==1){
      
      CellPerc[1,j] <- (CellPerc[1,j]+0.25)
      
    }
    
  }
}

stopCluster(cl)

sum(CellPerc)/length(PNN.names)

