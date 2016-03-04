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
setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Squamata/Results/Feb12/")

# Read the directory.
dir.tree <- dir()[grep("_tree",dir())]
dir.tree <- dir.tree[-grep("~",dir.tree)]

# Creat a empty list where we will put the trees.
multi.phylo <- list()

for (i in 1:length(dir.tree)){
  tree <- read.tree(dir.tree[i]) # Read each tree and...
  #plot.phylo(tree)
  multi.phylo[[i]] <- tree # Put inside the list, at the enda we create a multiphylo object.
}

# Create a multidata with the ocurrences files created before.

setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Squamata/Results/Feb12/")
dir.data <-(dir()[grep(".matrix",dir())])
dir.data <- dir.data[-grep(".matrix~",dir.data)]

multi.data <- list()

for(i in 1:length(dir.data)){
  multi.data[[i]] <- read.csv(dir.data[i]) 
}

# Set names for each data.
names(multi.data) <- dir.data
head(multi.data$Area.dist.matrix,5L)

multi.data[[5]]$Yurubí <- rep(0, length(multi.data[[5]]$Yurubí))

########################################################################

## Extract the tips from all phylogenies
## Those species will be our pool data species, for the permutation

## Vector where the species will put

dead.pool <- c()

## Extract all species terminals

for (i in 1:length(multi.phylo)){
  
  tax <- multi.phylo[[i]]$tip.label
  dead.pool <- c(dead.pool,tax)
}

# All species from phylogenies
dead.pool 

# Find the phylogeny species that match with the distribution species
match.sp <- dead.pool[which(dead.pool%in%multi.data[[1]]$especie)]
match.sp
## Create three empty list where the success will put

dt.match <- as.list(rep(NA,length(multi.data)))
pd.match <- as.list(rep(NA,length(multi.data)))
avtd.match <- as.list(rep(NA,length(multi.data)))

## Create empty lists for the original data results
dt.origin <- list()
pd.origin <- list() 
avtd.origin <- list()

for (i in 1:length(multi.phylo)){
  # And for each distribution data
  for(j in 1:length(multi.data)){
    # Both data have the same species ? Extrac from the data distribution only the species shared with the phylogeny
    dist<- multi.data[[j]][which((multi.data[[j]]$especie%in%multi.phylo[[i]]$tip.label)==T),]
    # Sometimes after the below process are more species in the phylogeny than the distribution data
    if(length(dist$especie)<length(multi.phylo[[i]]$tip.label)){
      # Extrac from the phylogeny the specie that miss in the data distribution
      miss.sp <- multi.phylo[[i]]$tip.label[which((multi.phylo[[i]]$tip.label%in%dist$especie)==F)]
      # Create a vector with distribution species and the miss specie, it will be necessary next.
      sp <- c(as.character(dist$especie),miss.sp)
      # Now attach the miss specie(s) to the distribution matrix.
      # Because this specie not present occurences in the area, the row will fill with 0
      for (x in 1:length(miss.sp)){
        dist.tmp <- rep(NA,length(colnames(dist)))
        dist <- rbind(dist,dist.tmp)
        dist[length(dist$especie),2:length(colnames(dist))] <- 0
      }
      # Now replace the species column with the sp vector created before
      dist$especie <- sp
    }
  
    # First, extract colnames and species in different vectors.
    sp.dist0 <- dist$especie
    areas0 <- colnames(dist)
    
    # Remove the species colum, and the transpose the matrix.
    # The transpose is necessary for PD function.
    dist2<- dist[,-1]
    dist2 <- t(dist2)
    # Put the species' names as column names
    colnames(dist2) <- sp 
    
    # DT Calculation
    Ind <- Calculate.Index(tree=multi.phylo[[i]],dist=dist)
    
    name.data <- Ind$area
    if(i == 1){
      dt.origin[[j]] <- Ind
    }else{
      dt.origin[[j]]$area <- NA
      Ind$area <- NA
      dt.origin[[j]] <-dt.origin[[j]] + Ind
      dt.origin[[j]]$area<- name.data}
    
    #PD Calculation
    pd <- pd(samp = dist2, tree = multi.phylo[[i]] , include.root = T)
    
    if (i == 1){ 
      pd.origin[[j]] <- pd
    }else{
      pd.origin[[j]] <- pd.origin[[j]] + pd}
    
    #AvDT calculation
    tree.dist <- cophenetic.phylo(multi.phylo[[i]])
    
    avdt <- taxondive(comm = dist2, dis = tree.dist)
    
    avdt2 <- data.frame(Species=avdt$Species,
                        D=avdt$D,
                        Dstar=avdt$Dstar,
                        Lambda=avdt$Lambda,
                        Dplus=avdt$Dplus,
                        sd.Dplus=avdt$sd.Dplus,
                        SDplus=avdt$SDplus,
                        ED=avdt$ED,
                        EDstar=avdt$EDstar,
                        EDplus=avdt$EDplus)
    
    if(i == 1){
      avtd.origin[[j]] <- avdt2
    }else{
      
      avdt2[is.na(avdt2)] <- 0
      avtd.origin[[j]][is.na(avtd.origin[[j]])] <- 0
      
      avtd.origin[[j]] <- avtd.origin[[j]] + avdt2}
  
  }
}


## Create a list were the results of the resampling will put

dt.jack <- list()
pd.jack <- list()
avdt.jack <- list()


## Jacknife calculation

for (ii in 1:100){
  ## Resampling from the species distribution
  ## 25% of species will resample
  rem.samp <- match.sp[sample(x = 1:length(match.sp), size = 0.25*length(match.sp))]
  
  
  for (i in 1:length(multi.phylo)){
    # And for each distribution data
    for(j in 1:length(multi.data)){
      # Both data have the same species ? Extrac from the data distribution only the species shared with the phylogeny
      dist<- multi.data[[j]][which((multi.data[[j]]$especie%in%multi.phylo[[i]]$tip.label)==T),]
      # Sometimes after the below process are more species in the phylogeny than the distribution data
      if(length(dist$especie)<length(multi.phylo[[i]]$tip.label)){
        # Extrac from the phylogeny the specie that miss in the data distribution
        miss.sp <- multi.phylo[[i]]$tip.label[which((multi.phylo[[i]]$tip.label%in%dist$especie)==F)]
        # Create a vector with distribution species and the miss specie, it will be necessary next.
        sp <- c(as.character(dist$especie),miss.sp)
        # Now attach the miss specie(s) to the distribution matrix.
        # Because this specie not present occurences in the area, the row will fill with 0
        for (x in 1:length(miss.sp)){
          dist.tmp <- rep(NA,length(colnames(dist)))
          dist <- rbind(dist,dist.tmp)
          dist[length(dist$especie),2:length(colnames(dist))] <- 0
        }
        # Now replace the species column with the sp vector created before
        dist$especie <- sp
      }
      
      # First, extract colnames and species in different vectors.
      sp.dist0 <- dist$especie
      areas0 <- colnames(dist)
      
      # Remove the species colum, and the transpose the matrix.
      # The transpose is necessary for PD function.
      dist2<- dist[,-1]
      dist2 <- t(dist2)
      # Put the species' names as column names
      colnames(dist2) <- sp 
      
      # DT Calculation
      Ind <- Calculate.Index(tree=multi.phylo[[i]],dist=dist)
      
      name.data <- Ind$area
      if(i == 1){
        dt.jack[[j]] <- Ind
      }else{
        dt.jack[[j]]$area <- NA
        Ind$area <- NA
        dt.jack[[j]] <-dt.jack[[j]] + Ind
        dt.jack[[j]]$area<- name.data}
      
      #PD Calculation
      pd <- pd(samp = dist2, tree = multi.phylo[[i]] , include.root = T)
      
      if (i == 1){ 
        pd.jack[[j]] <- pd
      }else{
        pd.jack[[j]] <- pd.jack[[j]] + pd}
      
      #AvDT calculation
      tree.dist <- cophenetic.phylo(multi.phylo[[i]])
      
      avdt <- taxondive(comm = dist2, dis = tree.dist)
      
      avdt2 <- data.frame(Species=avdt$Species,
                          D=avdt$D,
                          Dstar=avdt$Dstar,
                          Lambda=avdt$Lambda,
                          Dplus=avdt$Dplus,
                          sd.Dplus=avdt$sd.Dplus,
                          SDplus=avdt$SDplus,
                          ED=avdt$ED,
                          EDstar=avdt$EDstar,
                          EDplus=avdt$EDplus)
      
      if(i == 1){
        avdt.jack[[j]] <- avdt2
      }else{
        
        avdt2[is.na(avdt2)] <- 0
        avdt.jack[[j]][is.na(avdt.jack[[j]])] <- 0
        
        avdt.jack[[j]] <- avdt.jack[[j]] + avdt2}
      
    }
  }
  
  for(xx in 1:length(dt.origin)){
    if(dt.jack[[xx]]$I==dt.origin[[xx]]$I){dt.match[[xx]][ii] <- 1}
    if(pd.jack[[xx]]$PD==pd.origin[[xx]]$PD){pd.match[[xx]][ii] <- 1}
    if(avdt.jack[[xx]]$Dplus==avtd.origin[[xx]]$Dplus){avtd.match[[xx]][ii] <- 1}
  }
  
}
