# This script perfomance the calculation of DT, PD and AvDT index, from differents tree files.
# To use this script online need to change the working directories.
# DT function requiere a very specific data format, check it at the package examples

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

###########################################################################################
###########################################################################################
##                                                                                       ##
##                                        DT CALCULATION                                 ##
##                                                                                       ##
###########################################################################################
###########################################################################################

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

## Now the DT calculation.
#i <- 1; j <- 1
# For each i in multiphylo
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
    # Set the directory where the results will put
    setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Aves/Results/Ene8/")
    # Dt index calculation
    Ind <- Calculate.Index(tree=multi.phylo[[i]],dist=dist)
    # Create the output file name
    file <- paste(paste(dir.tree[i],dir.data[j],sep = "_"),".res")
    # Save the results
    write.csv(Ind,file = file,row.names = F,quote = F)
  }
}

# Sum all DT results.
# Setworking Dircetory


dir <- dir()[grep(".res",dir())]
# Finn each kind of result
enar <- dir[grep("Area",dir)]
g1 <- dir[grep("1g",dir)]
g25 <- dir[grep("25g",dir)]
g50 <- dir[grep("50g",dir)]
pnn <- dir[grep("PNN",dir)]

# Create a list with all type data
# and a empty list where the results will put
sumas <- list(enar,g1,g25,g50,pnn)
total <- list()

for (i in 1:length(sumas)){
  # Extract the names of the areas
  name.data <- read.csv(sumas[[i]][1])$area
  # Sum the two first data
  suma <- read.csv(sumas[[i]][1]) + read.csv(sumas[[i]][2])
  
  for(j in 3:length(sumas[[i]])){
    # Sum the rest of the data
    suma <- suma + read.csv(sumas[[i]][j])
  }
  # Because characters aren't summables, NA are introduced during summatory
  # Add original labels
  suma$area <- name.data
  # Finally each total summatory for each data kind are put in the total list
  total[[i]] <- suma
}

# Add respective names to the list
names(total) <- c("DTenar.total","DTg1.total","DTg25.total","DTg50.total","DTPNN.total")
# Save all results.
for ( i in 1:length(total)){
  write.csv(x=total[[i]],file=names(total)[i],row.names = F,quote = F)
}

system("rm *.res")

############################################################################################################
############################################################################################################
##                                                                                                        ##
##                                           PD CALCULATION                                               ##
##                                                                                                        ##
############################################################################################################
############################################################################################################

## The same directory and the same multi phylo and data will be used
## The script will be too similar to the DT scripts with few variations
## Similar variables to DT CALCULATION

dir.data ## Directory of the distribution data
dir.tree ## Directory of the Phylogenies 
multi.data ## All distribution data in a list
multi.phylo ## All phylogenies in a list

#####################################################################
# Now the calculation of PD, this calculation don't inlcude de root.

## Now the PD calculation.
#i <- 1; j <- 1
# For each i in multiphylo
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
    #At this moment begin all variations from DT calculation.
    
    # First, extract colnames and species in different vectors.
    sp.dist0 <- dist$especie
    areas0 <- colnames(dist)
    
    # Remove the species colum, and the transpose the matrix.
    # The transpose is necessary for PD function.
    dist<- dist[,-1]
    dist <- t(dist)
    # Put the species' names as column names
    colnames(dist) <- sp
    
    # For PD function a rooted tree is required. Root the tree.
    # Always the outgroup is the tip label 1, so this will be the rooting terminal
    # If there are a basal politomy, resolve.root fix it.
    
    #tree.root <- root(multi.phylo[[5]], outgroup = "Corvus_corax",resolve.root = T)
    
    # Create the output file name, set working directory and save th results.
    pd <- pd(samp = dist, tree = multi.phylo[[i]] , include.root = T)
    setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Aves/Results/Ene8/")
    file <- paste(paste(dir.tree[i],dir.data[j],sep = "_"),".res")
    write.csv(pd,file = file, row.names = T, quote = F)
  }
}

# Sum all PD results.
# Setworking Dircetory
#setwd("")

dir <- dir()[grep(".res",dir())]
# Finn each kind of result
enar <- dir[grep("Area",dir)]
g1 <- dir[grep("1g",dir)]
g25 <- dir[grep("25g",dir)]
g50 <- dir[grep("50g",dir)]
pnn <- dir[grep("PNN",dir)]

# Create a list with all type data
# and a empty list where the results will put

sumas <- list(enar,g1,g25,g50,pnn)
total <- list()

for (i in 1:length(sumas)){
  # Extract the names of the areas
  name.data <- read.csv(sumas[[i]][1])$X
  # Sum the two first data
  suma <- read.csv(sumas[[i]][1]) + read.csv(sumas[[i]][2])
  
  for(j in 3:length(sumas[[i]])){
    # Sum the rest of the data
    suma <- suma + read.csv(sumas[[i]][j])
  }
  # Because characters aren't summables, NA are introduced during summatory
  # Add original labels
  suma$X <- name.data
  # Finally each total summatory for each data kind are put in the total list
  total[[i]] <- suma
}

# Add respective names to the list
names(total) <- c("PDenar.total","PDg1.total","PDg25.total","PDg50.total","PDPNN.total")
# Save all results.
for ( i in 1:length(total)){
  write.csv(x=total[[i]],file=names(total)[i],row.names = F,quote = F)
}

system("rm *.res")

###############################################################################################
###############################################################################################
##                                                                                           ##
##                                      AvDT CALCULATION                                     ##
##                                                                                           ##
###############################################################################################
###############################################################################################

## The script are the same to PD CALCULATION
## The only differences is the creation of pairwise distance matrix from the phylogeny
## The four same important variables


dir.data ## Directory of the distribution data
dir.tree ## Directory of the Phylogenies 
multi.data ## All distribution data in a list
multi.phylo ## All phylogenies in a list

## Now the calculation of the AvDT.

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
    dist<- dist[,-1]
    dist <- t(dist)
    # Put the species' names as column names
    colnames(dist) <- sp
    
    #At this moment begin all variations from DT calculation.
    
    # AvDT requiere a pairwise distance matrix from the phylogeny.
    # So a cophenectic distance is used.
    
    tree.dist <- cophenetic.phylo(multi.phylo[[i]])
    
    avdt <- taxondive(comm = dist, dis = tree.dist)
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
    # Create output file name, set workig directory and save results.
    setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Aves/Results/Ene8/")
    file <- paste(paste(dir.tree[i],dir.data[j],sep = "_"),".res")
    write.table(avdt2,file=file,row.names=T,col.names = T,quote=F,sep=",")
    
  }
}


# Sum all PD results.
# Setworking Dircetory
setwd("")

dir <- dir()[[grep(".res",dir())]
# Finn each kind of result
enar <- dir[grep("Area",dir)]
g1 <- dir[grep("1g",dir)]
g25 <- dir[grep("25g",dir)]
g50 <- dir[grep("50g",dir)]
pnn <- dir[grep("PNN",dir)]

# Create a list with all type data
# and a empty list where the results will put
sumas <- list(enar,g1,g25,g50,pnn)
total <- list()

for(i in 1:length(sumas)){
  aa <- read.csv(sumas[[i]][1])
  aa[is.na(aa)] <- 0
  
  bb <- read.csv(sumas[[i]][2])
  bb[is.na(bb)] <- 0
  
  suma <- aa + bb
  
  for(j in 3:length(sumas[[i]])){
    cc <- read.csv(sumas[[i]][j])
    cc[is.na(cc)] <- 0
    
    suma <- suma+cc
  }
  
  total[[i]] <- suma
}

# Add respective names to the list
names(total) <- c("AVDTenar.total","AVDTg1.total","AVDTg25.total","AVDTg50.total","AVDTPNN.total")
# Save all results.
for ( i in 1:length(total)){
  write.csv(x=total[[i]],file=names(total)[i],row.names = F,quote = F)
}

system("rm *.res")
