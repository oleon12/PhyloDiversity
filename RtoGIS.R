library(maptools)
library(RColorBrewer)
library(classInt)
library(maps)
library(maptools)
library(dismo)


setwd("/home/omar/Documentos/Omar/Tesis/Scripts/Distribution/shp/Grid/")

grid <- readShapePoly("grid_50g.shp")



setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Results/May18/Raw_IndexR/")

x <- read.csv("PD.grid50",header=T)

x <- x$PD
x

x2 <- x[-grep(0,x)]

brks <- classIntervals(x2,n=5,style = "quantile")
brks <- brks$brks
brks

class <- findInterval(x,brks,all.inside = T)

class[-c(grep(5,class),grep(4,class))] <-  "NO"

class[grep(5,class)] <- "Q5"
class[grep(4,class)] <- "Q4"

class

grid$Index <- class

setwd("~/Documentos/Omar/Tesis/Taxa/Results/May18/Raw_IndexR/")

writePolyShape(grid,fn = "Grid50_PD")
