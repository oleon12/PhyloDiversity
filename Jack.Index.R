# This script perfomance a Jacknife support for three phylogenetic diversityy index.
# To use this script only nee to change the working directory
# DT function requiere a very specific data format, check it at the package example

## Load libraries.

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

## Load the phylogenies and distributions


# Set working directory.
setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Aves/Phyml/Ene8/")

# Read the directory.
dir.tree <- dir()[grep("tree.txt",dir())]

# Creat a empty list where we will put the trees.
multi.phylo <- list()

for (i in 1:length(dir.tree)){
  tree <- read.tree(dir.tree[i]) # Read each tree and...
  plot.phylo(tree)
  multi.phylo[[i]] <- tree # Put inside the list, at the enda we create a multiphylo object.
}

# Create a multidata with the ocurrences files created before.

setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Aves/Results/Ene8/")
dir.data <-(dir()[grep(".matrix",dir())])
dir.data <- dir.data[-grep(".matrix~",dir.data)]

multi.data <- list()

for(i in 1:length(dir.data)){
  multi.data[[i]] <- read.csv(dir.data[i]) 
}

# Set names for each data.
names(multi.data) <- dir.data
head(multi.data$Area.dist.matrix,5L)

########################################################################

## Extract the tips from all phylogenies
## Those species will be our pool data species, for the permutation

## Vector where the species will put

dead.pool <- c()

## Extract all species terminals

for (i in 1:length(multi.phylo)){

  tax <- multi.phylo$tip.label
  dead.pool <- c(dead.pool,tax)

}

