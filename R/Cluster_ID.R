library(ape)
library(phytools)
library(picante)
library(phangorn)
library(jrich)
library(fpc)
library(cluster)
library(maptools)
library(rgdal)
library(ggplot2)
library(gridExtra)
library(factoextra)

source("...phylo.id.table.R")
source("...pd2 .b.R")
source("...toNum.R")
source("...match.phy.dist.R")
source("...cstats.table.R")
source("...cluster.stats.2.R")
source("...phylo.beta.R")
source("...get.cluster.n.R")

########################################################################
#                              Areas List                              #
########################################################################

area <- read.csv("../Data/Distributions/grid_25g.dist.matrix")

str(area)

########################################################################
#                                Q5 Cells                              #
########################################################################

#### PD Data
pd <- readOGR("../Data/shp/Grid25_PD.shp")@data

#ID Cells
Cells <- as.character(1:length(pd$SP_ID))

pd2 <- data.frame(Id=Cells[which(pd$Index=="Q5")],
                 Class=pd$Index[which(pd$Index=="Q5")],
                 PD=pd$IndexVals[which(pd$Index=="Q5")])
head(pd2)

########################################################################

sp <- area$Species

area <- area[,-1]

areaQ5 <- area[,which(pd$Index=="Q5")]

########################################################################

distQ5 <- as.data.frame(vegdist(t(areaQ5), method = "jaccard"))

distQ5

########################################################################

clustN <- fviz_nbclust(x = t(areaQ5),kmeans, method = "wss")
clustN

clustN2 <- fviz_nbclust(x = t(areaQ5),kmeans, method = "silhouette")
clustN2

clustN3 <- fviz_nbclust(x = t(areaQ5), kmeans, nstart = 25,  method = "gap_stat", nboot = 50)
clustN3

clustN3$data

########################################################################

stat.clust <- get.cluster.n(dist = distQ5, k = 9)

plot(stat.clust$N.groups, stat.clust$tot.withinss)
lines(stat.clust$N.groups, stat.clust$tot.withinss)

dend <- hclust(distQ5,method = "ward.D2")

#Beta total 3 clusters
#Beta rep 3 clusters
#Beta rich 4 clusters
plot(dend)
rect.hclust(dend, 6)

class.cluster <- cutree(dend, k=6)

##########################################################################

grid.25 <- readOGR("../Data/shp/grid_25g.shp")
area <- read.csv("../Data/Distributions/grid_25g.dist.matrix")

cells.n <- rep("NO", length(area[1,])-1)
names(cells.n) <- colnames(area)[-1]


cells.n[which(names(cells.n)%in%names(class.cluster))] <- class.cluster


grid.25$Cluster <- cells.n

writeOGR(obj = grid.25,driver = "ESRI Shapefile",dsn = "../Data/shp/", layer = "Grid25_Cluster.Jacc6")
