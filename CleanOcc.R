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

setwd("~/Documentos/Omar/Tesis/Taxa/Richness/")

## Read the occurences data.
sp.dist <- read.csv("Richness.occ",header = T)
# Create a new data frame pasting the two first original columns.
#sp.dist <- data.frame(Sp=paste(sp.dist[,1],sp.dist[,2],sep="_"), 
#                      Lon=sp.dist$Lon,
#                      Lat=sp.dist$Lat)
# Write the new table.
#write.csv(data,"mammals2.occ",row.names=F,col.names=T,quote=F) 

head(sp.dist,3L)

setwd("~/Documentos/Omar/Tesis/Scripts/Distribution/shp/NAB/")

NAB.poly <- read.shp("NAB.shp")

out.occ <- c()

for(i in 1:length(sp.dist$Sp)){
  
  print(paste("Time:",length(sp.dist$Sp)-i,sep = " "))
  
  occ.d <- sp.dist[i, c(2,3)]
  
  out <- pnt.in.poly(occ.d,NAB.poly$shp[[1]]$points)

  if(out$pip==0){out.occ <- c(out.occ,i)}
  
}

length(out.occ)

sp.dist <- sp.dist[-out.occ,]

setwd("~/Documentos/Omar/Tesis/Taxa/Richness/")

write.table(sp.dist, "Richness.occ", quote = F, row.names = F, col.names = T, sep = ",")
