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

## Go to the work directory where are the shape poly 
setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Results/Julio1/RawIndex/")
# Read the shape file
grid <- readShapePoly("Grid25_PD.shp")
pos.Q5 <- grep("Q5", grid$Index)
#
setwd("~/Documentos/Omar/Tesis/Scripts/Distribution/shp/Grid/")
grid2 <- read.shp("grid_25g.shp")

setwd("~/Documentos/Omar/Tesis/Scripts/Distribution/shp/PNN/PNN_NAB/")
PNN.names <- readShapePoly("PNN_NAB.shp")
PNN.names <- PNN.names$NAME
PNN <- read.shp("PNN_NAB.shp")

PNNinCell <- c()
PNNinCell100 <- c()
PNNinCell50 <- c()
CellinPNN <- c()


for(i in 1:length(pos.Q5)){
  
  pos <- pos.Q5[i]
  
  cell <- grid2$shp[[pos]]$points
  
  for (j in 1:length(PNN.names)) {
  
    PNN.coord <- PNN$shp[[j]]$points
    
    pnt <- pnt.in.poly(pnts = cell, poly.pnts = PNN.coord)
    
    if(1%in%pnt$pip){
      
      CellinPNN <- c(CellinPNN, pos)
      PNNinCell <- c(PNNinCell, as.character(PNN.names[j]))
      
    }
    
  }
}

for(i in 1:length(PNN.names)){
  
  PNN.coord <- PNN$shp[[i]]$points

  for (j in 1:length(pos.Q5)) {
    
    pos <- pos.Q5[j]

    cell <- grid2$shp[[pos]]$points  
    
    pnt2 <- pnt.in.poly(pnts = PNN.coord,poly.pnts = cell)
    
    perc50 <- length(PNN.coord[,1])/2
    perc100 <- length(PNN.coord[,1])
    
    if(length(grep(1, pnt2$pip))>=perc50){ PNNinCell50 <- c(PNNinCell50, as.character(PNN.names[i]))}
    if(length(grep(1, pnt2$pip))==perc100){ PNNinCell100 <- c(PNNinCell100, as.character(PNN.names[i]))}
    
  }
}


CellinPNN
length(CellinPNN)/length(pos.Q5)

PNNinCell
unique(PNNinCell50)
unique(PNNinCell100)

length(PNNinCell)/length(PNN.names)


