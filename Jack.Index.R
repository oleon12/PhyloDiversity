## Autor
# Leon-Alvarado, Omar Daniel.
# leon.alvarado12@gmail.com

## License
# The follow script was created under the GNU/GPLv2. license.
# http://www.gnu.org/licenses/old-licenses/gpl-2.0-standalone.html

## Title
# Jackknife implementation for Multiple Phylogenetic Diversity Index 

## Description
# This R scripts perfomance the a Jackknife calculation for Taxonomic Distincness (DT) (Vane-Wrigth et al. 1991), Phylogenetic Diversity (PD) (Faith 1992) and Average Taxonomic Distincness (AvTD) (Clarke & Warwick 1998) index for multiple distribution matrices and phylogenies
# This script requires two list objects: 
# 1. A list object with distribution matrices.
# 2. A list object with phylogenies.
# Species' names must be the same in the distribution matrices and phylogenies.
# The distribution matrices format to use are very specific, the same implemented in the package Jrich (Dmirandae/Jrich), see packages examples
# To use this script must change the working directories
# The numbers of iterations could be change modifying the lenght of the ii variable in the first loop for.
# Also, the remotion percetange could be change modifying the default value (0.25) in the variable rem.samp at the beginning of the first loop for.
# The script returns a matrix with Jackknife support value of the index for each distibution matrix. 

####################################################################################3
toNum <-function(x){
  out <- matrix(NA, nrow = length(rownames(x)), ncol = length(colnames(x)))
  for (i in 1:length(colnames(x))){
    out[,i] <- as.numeric(x[,i])
  }
  colnames(out) <- colnames(x)
  rownames(out) <- rownames(x)
  
  return(out)
}

####################################################################################
avtd.root <- function(Phylo, Dist){
  
  rootName <- find.root(Phylo)
  
  for(i in 1:length(Dist[,1])){
    
    l <- grep(1,Dist[i, ])
    
    if(length(l)==1){
      
      Dist[i,grep(rootName,colnames(Dist))] <- 1
      
    }
    
  }
  
  return(Dist)
}
######################################################################################

find.root <- function(Phylo){
  
  library(ape, verbose = F,quietly = T)
  library(phangorn,verbose = F, quietly = T)
  
  ##########################################
  # Filter
  if(class(Phylo)!="phylo"){
    stop("Your input must be a class phylo")
  }
  if(is.rooted(Phylo)==F){
    stop("Your phylogeny must be rooted")
  }
  ###########################################
  
  RootNode <- length(Phylo$tip.label)+1
  
  Option <- c(1,length(Phylo$tip.label))
  
  RootTip <- c()
  
  for(i in 1:2){
    
    tip <- Option[i]
    
    Anc <- Ancestors(Phylo,tip)
    
    if(length(Anc)==1){
      
      RootTip <- Phylo$tip.label[Option[i]]
    }
  }
  
  if(is.null(RootTip)){
    
    warning("Your root terminal are not a unique specie, thus, just one will choose")
    
    for(i in 1:2){
      
      tip <- Option[i]
      
      Anc <- Ancestors(Phylo,tip)
      
      if(length(Anc)==2){
        
        RootTip <- Phylo$tip.label[Option[i]]
      }
    }
  }
  return(RootTip)
}

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
  setwd("~/Documentos/Omar/Result_TesisOmar/Trees/")

# Read the directory.
dir.tree <- dir()[grep("_tree",dir())]
#dir.tree <- dir.tree[-grep("~",dir.tree)]

# Creat a empty list where we will put the trees.
multi.phylo <- list()

for (i in 1:length(dir.tree)){
  tree <- read.tree(dir.tree[i]) # Read each tree and...
  #plot.phylo(tree)
  multi.phylo[[i]] <- tree # Put inside the list, at the enda we create a multiphylo object.
}

# Create a multidata with the ocurrences files created before.

setwd("~/Documentos/Omar/Result_TesisOmar/Final2/")
dir.data <-(dir()[grep(".matrix",dir())])
#dir.data <- dir.data[-grep(".matrix~",dir.data)]

multi.data <- list()

for(i in 1:length(dir.data)){
  multi.data[[i]] <- read.csv(dir.data[i]) 
}

# Set names for each data.
names(multi.data) <- dir.data
head(multi.data$Area.dist.matrix,5L)

#multi.data[[5]]$Yurubí <- rep(0, length(multi.data[[5]]$Yurubí))

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


## Create empty lists for the original data results
dt.origin <- list()
pd.origin <- list() 
avtd.origin <- list()

for (i in 1:length(multi.phylo)){
  
  print(paste(paste("Phylogeny", i, sep=" "), paste("from", length(multi.phylo), sep=" "), sep=" "))
  
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
    
    if (is.numeric(dist2[1,1])==F){
      
      dist2 <- toNum(dist2)
      
    }
    
    # Put the species' names as column names
    colnames(dist2) <- sp 
    
    # DT Calculation
    Ind <- Calculate.Index(tree=multi.phylo[[i]],dist=dist, verbose = F)
    
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
    distRoot <- avtd.root(Phylo = multi.phylo[[i]],Dist = dist2)
    
    tree.dist <- cophenetic.phylo(multi.phylo[[i]])
    
    avdt <- taxondive(comm = distRoot, dis = tree.dist)
    
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


area.jack <- matrix(0, nrow = length(dt.origin[[1]]$area),ncol = 100)
grid1.jack <- matrix(0, nrow = length(dt.origin[[2]]$area),ncol = 100)
grid25.jack <- matrix(0, nrow = length(dt.origin[[3]]$area),ncol = 100)
grid50.jack <- matrix(0, nrow = length(dt.origin[[4]]$area),ncol = 100)
pnn.jack <- matrix(0, nrow = length(dt.origin[[5]]$area),ncol = 100)

jackTD.check <- list(area.jack, grid1.jack, grid25.jack, grid50.jack, pnn.jack)
jackPD.check <- list(area.jack, grid1.jack, grid25.jack, grid50.jack, pnn.jack)
jackAVTD.check <- list(area.jack, grid1.jack, grid25.jack, grid50.jack, pnn.jack)

## Jacknife calculation
for (ii in 1:100){
  ## Resampling from the species distribution
  ## 25% of species will resample
  rem.samp <- match.sp[sample(x = 1:length(match.sp), size = round(0.25*length(match.sp),digits = 0))]
  print(paste("Iteration",ii,sep = ":"))
  
  for (i in 1:length(multi.phylo)){
    
    print(paste(paste("Phylogeny No",i,sep = ": "), paste("From", length(multi.phylo), sep=": "), sep = " "  ))
    
    # And for each distribution data
    for(j in 1:length(multi.data)){
      
      dist1 <- multi.data[[j]]
      
      dist1[which(dist1$especie%in%rem.samp),2:length(colnames(dist1))] <- 0
      
      # Both data have the same species ? Extrac from the data distribution only the species shared with the phylogeny
      dist<- dist1[which((dist1$especie%in%multi.phylo[[i]]$tip.label)==T),]
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
      
      
      ## Remove species
      
      #dist[which(dist$especie%in%rem.samp),2:length(colnames(dist))] <- 0
      
      # First, extract colnames and species in different vectors.
      sp.dist0 <- dist$especie
      areas0 <- colnames(dist)
      
      # Remove the species colum, and the transpose the matrix.
      # The transpose is necessary for PD function.
      dist2<- dist[,-1]
      dist2 <- t(dist2)
      
      if (is.numeric(dist2[1,1])==F){
        
        dist2 <- toNum(dist2)
        
      }
      
      
      # Put the species' names as column names
      colnames(dist2) <- sp 
      
      # DT Calculation
      Ind <- Calculate.Index(tree=multi.phylo[[i]],dist=dist,verbose = F)
      
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
      distRoot <- avtd.root(Phylo = multi.phylo[[i]],Dist = dist2)
      
      tree.dist <- cophenetic.phylo(multi.phylo[[i]])
      
      avdt <- taxondive(comm = distRoot, dis = tree.dist)
      
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
    
    jackTD.check[[xx]][ ,ii] <- dt.jack[[xx]]$area[order(dt.jack[[xx]]$W, decreasing = T)]
    
    jackPD.check[[xx]][ ,ii] <- rownames(pd.jack[[xx]])[order(pd.jack[[xx]]$PD, decreasing = T)]
    
    jackAVTD.check[[xx]][ ,ii] <- rownames(avdt.jack[[xx]])[order(avdt.jack[[xx]]$Dplus, decreasing = T)]
    
    
  }
  
}


jackTD.check[[1]]

jackPD.check[[1]]

jackAVTD.check[[1]]


ii <- 1
xx <- 1

## Create three empty list where the success will put

dt.match <- as.list(rep(0,length(multi.data)))
pd.match <- as.list(rep(0,length(multi.data)))
avtd.match <- as.list(rep(0,length(multi.data)))


for(ii in 1:length(dt.origin)){
  
  for(xx in 1:length(jackAVTD.check[[1]][1,])){
    
    RankOriginTD <- dt.origin[[ii]]$area[order(dt.origin[[ii]]$W, decreasing = T)] 
    RankOriginPD <- rownames(pd.origin[[ii]])[order(pd.origin[[ii]]$PD, decreasing = T)]
    RankOriginAVTD <- rownames(avtd.origin[[ii]])[order(avtd.origin[[ii]]$Dplus, decreasing = T)]
    
    JackRankTD <- jackTD.check[[ii]][,xx]
    JackRankPD <- jackPD.check[[ii]][,xx]
    JackRankAVTD <- jackAVTD.check[[ii]][,xx]
    
    if(ii==1){
      
      c1 <- length(which((RankOriginTD[1:3]==JackRankTD[1:3])==T))
      if(c1==length(RankOriginTD[1:3])){dt.match[[ii]][xx] <- 1}
      
      c2 <- length(which((RankOriginPD[1:3]==JackRankPD[1:3])==T))
      if(c2==length(RankOriginPD[1:3])){pd.match[[ii]][xx] <- 1}
      
      c3 <- length(which((RankOriginAVTD[1:3]==JackRankAVTD[1:3])==T))
      if(c3==length(RankOriginAVTD[1:3])){avtd.match[[ii]][xx] <- 1}
    
    }
    
    if(ii==5){
      
      c1 <- length(which((RankOriginTD[1:3]==JackRankTD[1:3])==T))
      if(c1==length(RankOriginTD[1:3])){dt.match[[ii]][xx] <- 1}
      
      c2 <- length(which((RankOriginPD[1:3]==JackRankPD[1:3])==T))
      if(c2==length(RankOriginPD[1:3])){pd.match[[ii]][xx] <- 1}
      
      c3 <- length(which((RankOriginAVTD[1:3]==JackRankAVTD[1:3])==T))
      if(c3==length(RankOriginAVTD[1:3])){avtd.match[[ii]][xx] <- 1}
      
    }else{
      
      c1 <- length(which((RankOriginTD[1]==JackRankTD[1])==T))
      if(c1==length(RankOriginTD[1])){dt.match[[ii]][xx] <- 1}
      
      c2 <- length(which((RankOriginPD[1]==JackRankPD[1])==T))
      if(c2==length(RankOriginPD[1])){pd.match[[ii]][xx] <- 1}
      
      c3 <- length(which((RankOriginAVTD[1]==JackRankAVTD[1])==T))
      if(c3==length(RankOriginAVTD[1])){avtd.match[[ii]][xx] <- 1}
      
      
    }
    
    
  }
  
}


r.support <- matrix(0, nrow = length(dt.match), ncol = 3)

colnames(r.support) <- c("DT","PD","AvDT")
rownames(r.support) <- names(multi.data)

for(k in 1:length(rownames(r.support))){
  
  r.support[k,] <- c(length(which((dt.match[[k]]==1)==T))/100,
                     length(which((pd.match[[k]]==1)==T))/100,
                     length(which((avtd.match[[k]]==1)==T))/100) 
  
}

r.support

write.csv(r.support,"Soporte.Jack", quote = F, row.names = T)
