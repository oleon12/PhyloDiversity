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

multi.data <- list(read.csv("../Data/Distributions/PNN.inn_25g.matrix"),read.csv("../Data/Distributions/PNN.out_25g.matrix"))
names(multi.data) <- c("PNN.in","PNN.out")

########################################################################
#                   Phylogenetic Diversity Indices                     #
########################################################################

##PD
pd.out <- pd2(Phylogeny = multi.phylo, Distribution = multi.data, pdRoot = T)

names(pd.out) <- names(multi.data)

save(pd.out, file = "../Data/Indices/pd.out_PNN.rda")

###########################################################

PNN.in <- data.frame(Ricness = pd.out$PNN.in$SR,
                     PD = pd.out$PNN.in$PD,
                     AvTD = avtd.out$PNN.in$Dplus,
                     EDGE = apply(t(edge.out$PNN.in),1,sum))

write.csv(PNN.in, "../Data/Indices/Indices_PNN.in_Table.csv", quote = F, row.names = F)

PNN.out <- data.frame(Ricness = pd.out$PNN.out$SR,
                      PD = pd.out$PNN.out$PD,
                      AvTD = avtd.out$PNN.out$Dplus,
                      EDGE = apply(t(edge.out$PNN.out),1,sum))

write.csv(PNN.out, "../Data/Indices/Indices_PNN.out_Table.csv", quote = F, row.names = F)

###########################################################

inPNN <- read.csv("../Data/Indices/Indices_PNN.in_Table.csv")
outPNN <- read.csv("../Data/Indices/Indices_PNN.out_Table.csv")
Total <- read.csv("../Data/Indices/Indices_Table.csv")

outTotal <- apply(outPNN, 2, sum)
inTotal <- apply(inPNN, 2, sum)
Total <- apply(Total, 2, sum)

Total
outTotal
inTotal

Total; (outTotal+inTotal)

OutTable <- rbind(round(outTotal/Total, digits = 2),round(inTotal/Total, digits = 2))
rownames(OutTable) <- c("Outside_PNN","Within_PNN")
OutTable <- as.data.frame(OutTable)

OutTable

write.csv(OutTable, "../Data/Indices/PNN_Results.csv", quote = F, row.names = T)
