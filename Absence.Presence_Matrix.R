## Autor
# Leon-Alvarado, Omar Daniel.
# leon.alvarado12@gmail.com

## License
# The follow script was created under the GNU/GPLv2. license.
# http://www.gnu.org/licenses/old-licenses/gpl-2.0-standalone.html

## Title
# Absence/Presence Matrix.

## Description
# This R scripts create a absebce/presence matrix for the calculation of Taxonomic Distincness (DT) (Vane-Wrigth et al. 1991), Phylogenetic Diversity (PD) (Faith 1992) and Average Taxonomic Distincness (AvTD) (Clarke & Warwick 1998) indeces
# The script take all occurrences from each species and determines if have at least one point inside an area
# The script requieres a file with the occurrences. This file must have three colums:
# Species_Colum | Latitude_Colum | Longitude_Colum
# The files used for the areas are shapes build in QGIS.
# The outcome is a matrix object with species as rows and areas as colums, filled with 0/1

## Load libraries

library(rgeos)
library(maptools)
library(rgbif)
library(dismo)
library(shapefiles)
library(ape)
library(phangorn)
library(phytools)
library(picante)
library(jrich)
library(SDMTools)
library(doParallel)

## Set thw working directory.

setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Aves/Results/Ene8/")

## Read the occurences data.
sp.dist <- read.csv("Dist_ene8.occ",header = T)
# Create a new data frame pasting the two first original columns.
#sp.dist <- data.frame(Sp=paste(sp.dist[,1],sp.dist[,2],sep="_"), 
#                      Lon=sp.dist$Lon,
#                      Lat=sp.dist$Lat)
# Write the new table.
#write.csv(data,"mammals2.occ",row.names=F,col.names=T,quote=F) 

head(sp.dist,3L)

## Create an empty list object.
area <- list() 
# Read the directory where the endemis area's shape files are.
area.dir <- dir("/home/omar/Documentos/Omar/Tesis/Scripts/Distribution/shp/Morrone/")
# Exctrac only the .shp files from the entire directory.
area.dir <- area.dir[grep(".shp",area.dir)]

# Read each shp file and put it in the list object.
for(i in 1:length(area.dir)){
  setwd("/home/omar/Documentos/Omar/Tesis/Scripts/Distribution/shp/Morrone/")
  area[[i]] <- read.shp(shp.name = area.dir[i])
}

# Create a matrix filled with 0, the rows will be the species and the columns will be the areas.
# This will be our absece/presence matrix.
data.abpre <- matrix(0,nrow = length(levels(sp.dist$especie)),ncol=length(area))
colnames(data.abpre) <- area.dir
rownames(data.abpre) <- levels(sp.dist$especie)

# Now we need to fill the matrix with the presences.
for(i in 1:length(area)){
  for(j in 1:length(data.abpre[,1])){
    # Extrac only the Lat and Long for a j specie given the data.abpre matrix.
    sp.data <- sp.dist[which(sp.dist$especie==rownames(data.abpre)[j]),2:3]
    # Compare if at least one georeference point are inside the area polygon.
    out <- pnt.in.poly(sp.data,area[[i]]$shp[[1]]$points)
    # If a point are inside the area polygon, will fill his respective space with 1.
    if (1%in%out$pip){ data.abpre[j,i]<- 1; print(area.dir[i])} 
  }
}

# Save the matrix in his respective directory.

setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Aves/Results/Ene8/")
write.table(data.abpre,file="Area.dist.matrix",row.names=T, col.names=T,quote=F,sep=",")


## Now do the same, but replacing the endemism areas for a grid.

# Read the grid from his respective directory.

setwd("/home/omar/Documentos/Omar/Tesis/Scripts/Distribution/shp/Grid/")
grid.g <- read.shp("grid_50g.shp")

# Create a matrix filled with 0, the rows will be the species and the columns will be the areas.
# This will be our absece/presence matrix.

data.grid <- matrix(0, nrow = length(levels(sp.dist$especie)),ncol = length(grid.g$shp))
colnames(data.grid) <- 1:length(grid.g$shp)
rownames(data.grid) <- levels(sp.dist$especie)

cl <- makeCluster(8)
registerDoParallel(cl)

for(i in 1:length(data.grid[,1])){
  print(paste("Total species:",length(data.grid[,1])," "))
  print(paste("Specie:",i,sep="")) # Checkpoint
  for (j in 1:length(grid.g$shp)){
    # Extrac only the Lat and Long for a j specie given the data.grid matrix.
    sp.data <- sp.dist[which(sp.dist$especie==rownames(data.grid)[i]),2:3]
    # Compare if at least one georeference point are inside the area polygon.
    out <- pnt.in.poly(sp.data,grid.g$shp[[j]]$points)
    # If a point are inside the area polygon, will fill his respective space with 1.
    if(1%in%out$pip){data.grid[i,j]<-1; print(paste("Cell with presence",j,sep=""))}
  }
}

# Save the created data.
setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Aves/Results/Ene8/")
write.table(data.grid,file="grid_50g.dist.matrix",row.names=T, col.names=T,quote=F,sep=",")


# Read the PNN shp file from his respective directory.

setwd("/home/omar/Documentos/Omar/Tesis/Scripts/Distribution/shp/PNN/")
dir <- dir()[grep(".shp",dir())]
PNN <- list()

for(i in 1:length(dir)){
  PNN[[i]] <- read.shp(shp.name = dir[i])
}

# Create a matrix filled with 0, the rows will be the species and the columns will be the areas.
# This will be our absece/presence matrix.

data.PNN <- matrix(0, nrow = length(levels(sp.dist$especie)),ncol = length(PNN))
colnames(data.PNN) <- dir
rownames(data.PNN) <- levels(sp.dist$especie)

cl <- makeCluster(8)
registerDoParallel(cl)

for(i in 1:length(PNN)){
  print(paste("Total PNN:",length(PNN),sep=""))
  print(paste("PNN #:",i,sep=""))
  for(j in 1:length(data.PNN[,1])){
    # Extrac only the Lat and Long for a j specie given the data.abpre matrix.
    sp.data <- sp.dist[which(sp.dist$especie==rownames(data.PNN)[j]),2:3]
    # Compare if at least one georeference point are inside the area polygon.
    out <- pnt.in.poly(sp.data,PNN[[i]]$shp[[1]]$points)
    # If a point are inside the area polygon, will fill his respective space with 1.
    if (1%in%out$pip){ data.PNN[j,i]<- 1; print(i)} 
  }
}

# Save the created data.
setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Aves/Results/Ene8/")
write.table(data.PNN,file="PNN.dist.matrix",row.names=T, col.names=T,quote=F,sep=",")
