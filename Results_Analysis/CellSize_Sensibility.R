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
setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")
# Read the shape file
grid <- readShapePoly("Grid1_AvTD.shp")
pos1.Q5 <- which(grid$Index=="Q5")
pos1.Q5 <- pos1.Q5 + 1
grid1 <- read.shp("Grid1_AvTD.shp")
#
grid <- readShapePoly("Grid50_AvTD.shp")
pos50.Q5 <- which(grid$Index=="Q5")
pos50.Q5 <- pos50.Q5 + 1
grid50 <- read.shp("Grid50_AvTD.shp")
#
grid <- readShapePoly("Grid25_AvTD.shp")
pos25.Q5 <- which(grid$Index=="Q5")
pos25.Q5 <- pos25.Q5 + 1
grid25 <- read.shp("Grid25_AvTD.shp")
#

g25.50 <- c()
g25.1 <- c()
g50.1 <- c()

for(i in 1:length(pos25.Q5)){
  
  pos <- pos25.Q5[i]
  cellmin <- grid25$shp[[pos]]$points
  
  for(j in 1:length(pos50.Q5)){
    
    pos2 <- pos50.Q5[j]
    cellmax <- grid50$shp[[pos2]]$points
    
    pnt <- pnt.in.poly(pnts = cellmin, poly.pnts = cellmax)
    
    if (length(grep(1,pnt$pip))) {g25.50 <- c(g25.50, pos25.Q5[i])}
    
  }
}


one <- length(unique(g25.50))/length(pos25.Q5)

#################################################

for(i in 1:length(pos25.Q5)){
  
  pos <- pos25.Q5[i]
  cellmin <- grid25$shp[[pos]]$points
  
  for(j in 1:length(pos1.Q5)){
    
    pos2 <- pos1.Q5[j]
    cellmax <- grid1$shp[[pos2]]$points
    
    pnt <- pnt.in.poly(pnts = cellmin, poly.pnts = cellmax)
    
    if (length(grep(1,pnt$pip))) {g25.1 <- c(g25.1, pos25.Q5[i])}
    
  }
}


two <- length(unique(g25.1))/length(pos25.Q5)

################################################################

for(i in 1:length(pos50.Q5)){
  
  pos <- pos50.Q5[i]
  cellmin <- grid50$shp[[pos]]$points
  
  for(j in 1:length(pos1.Q5)){
    
    pos2 <- pos1.Q5[j]
    cellmax <- grid1$shp[[pos2]]$points
    
    pnt <- pnt.in.poly(pnts = cellmin, poly.pnts = cellmax)
    
    if (length(grep(1,pnt$pip))) {g50.1 <- c(g50.1, pos50.Q5[i])}
    
  }
}


three <- length(unique(g50.1))/length(pos50.Q5)

comparison <- cbind(one*100,two*100,three*100)
colnames(comparison) <- c("0.25° in 0.50°", "0.25° in 1°", "0.50° in 1°")
rownames(comparison) <- "Percentage"

comparison <- as.data.frame(comparison)

comparison

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/")

write.csv(comparison, "AvTD_SensibilityAnalysis", quote = F, row.names =T ,col.names = T)
