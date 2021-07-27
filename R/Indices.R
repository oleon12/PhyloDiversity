library(ape)
library(phytools)
library(vegan)
library(maptools)
library(RColorBrewer)
library(classInt)
library(maps)
library(maptools)
library(dismo)
library(picante)
library(rgdal)
library(ggplot2)
library(gridExtra)

########################################################################

setwd("../Data/Trees/")

########################################################################

source("....find.root.R")
source("....match.phy.dist.R")
source("....pd2.R")
source("....toNum.R")

########################################################################
#                               Tree List                              #
########################################################################

dir.tree <- as.list(dir()[grep("_tree",dir())])

# Creat a empty list where we will put the trees.
multi.phylo <- lapply(dir.tree, function(x){read.tree(x)})

names(multi.phylo) <- dir.tree

########################################################################
#                              Areas List                              #
########################################################################

setwd("../Data/Richness/")

dir.data <- as.list(dir()[grep(".matrix",dir())])

multi.data <- lapply(dir.data, function(x){read.csv(x)})

# Set names for each data.
names(multi.data) <- dir.data
str(multi.data)

########################################################################
#                   Phylogenetic Diversity Indices                     #
########################################################################

##PD
pd.out <- pd2(Phylogeny = multi.phylo, Distribution = multi.data, pdRoot = T)

names(pd.out) <- c("Area.end","Endemic.sp","Grid.25")

save(pd.out, file = "../Data/Indices/pd.out.rda")

########################################################################
#                       Quantil prioritization                         #
########################################################################

grid <- readShapePoly("../Data/shp/grid_25g.shp")

ind.grid <- data.frame(PD = pd.out$Grid.25$PD)

ind.grid2 <- data.frame(Ricness = pd.out$Grid.25$SR,
                        PD = pd.out$Grid.25$PD)

write.csv(ind.grid2, "..Data/Indices/Indices_Table.csv", quote = F, row.names = F)

PD <- ind.grid$PD
PD2 <- PD[which(PD>0.0000000)]
head(PD2)

brks <- classIntervals(PD2,n=5,style = "quantile")
brks <- brks$brks
brks

class <- findInterval(PD,brks,all.inside = T)

Pos5 <- grep(5,class)
Pos4 <- grep(4,class)
Pos3 <- grep(3,class)
Pos2 <- grep(2,class)
Pos1 <- grep(1,class)
NonVal <- which(PD==0)

#class[-Pos5] <-  "Non-Q5"
#class[Pos5] <- "Q5"
#ClassPD1 <- class # Distiguish Q5 and Non-Q5

class[Pos5] <- "Q5"
class[Pos4] <- "Q4"
class[Pos3] <- "Q3"
class[Pos2] <- "Q2"
class[Pos1] <- "Q1"
class[NonVal] <- "NO"

ClassPD2 <- class #All categories

grid$Index <- ClassPD2
grid$IndexVals <- PD

writePolyShape(grid,fn = "../Data/shp/Grid25_PD")
PD.poly <- grid

save(PD.poly, file = "../Data/Indices/PD.poly.rda")
