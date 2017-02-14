## Autor
# Leon-Alvarado, Omar Daniel.
# leon.alvarado12@gmail.com

## License
# The follow script was created under the GNU/GPLv2. license.
# http://www.gnu.org/licenses/old-licenses/gpl-2.0-standalone.html

## Title
# Cells in PAs

## Description 
# This script compares the results obtained from the grid cells prioritization and the PAs
# Two approaches are implemented:
# 1. Compare how many Q5 cells match with a PA
#  Here, the count of cells are made using the next rules
#
#    -If the 4 points of the cell are inside, then the cell get a value of 1
#    -If 3 points of the cell are inside, then the cell get a value of 0.75
#    -If 2 points of the cell are inside, then the cell get a value of 0.50
#    -If 1 point of the cell are inside, the the cell get a value of 0.25
#
#2. Compare the percentage of PA that are inside a cell
#  Here, only were evaluate if the 50% or the 100% of a PA are inside a cell

############################################################################

# Load libraries
library(maptools)
library(shapefiles)
library(SDMTools)

## Go to the work directory where are the shape poly 
setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")
# Read the shape file and find only the position of the Q5 cells
grid <- readShapePoly("Grid25_AvTD.shp")
pos.Q5 <- which(grid$Index=="Q5")
pos.Q5 <- pos.Q5+1
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

#######################################################################
#                     Approach 1: Q5 cells in PAs                     #
#######################################################################

for(i in 1:length(pos.Q5)){
  
  pos <- pos.Q5[i]
  
  cell <- grid2$shp[[pos]]$points
  
  for(j in 1:length(PNN.names)){
    
    PNN.coord <- PNN$shp[[j]]$points
    
    pnt <- pnt.in.poly(pnt= cell, poly.pnts= PNN.coord)
    
    n <- grep(1, pnt$pip)
    
    if(length(n)==5){CellinPNN <- c(CellinPNN,1)}
    
    if(length(n)==4){CellinPNN <- c(CellinPNN,1)}
    
    if(length(n)==3){CellinPNN <- c(CellinPNN,0.75)}
    
    if(length(n)==2){CellinPNN <- c(CellinPNN,0.50)}
    
    if(length(n)==1){CellinPNN <- c(CellinPNN,0.25)}
    
  }
}

sum(CellinPNN) / length(pos.Q5)

#######################################################################
#                     Approach 2: PAs in Q5 cell                      #
#######################################################################

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


PNNinCell
unique(PNNinCell50)
unique(PNNinCell100)

length(PNNinCell)/length(
