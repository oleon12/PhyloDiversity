## Autor
# Leon-Alvarado, Omar Daniel.
# leon.alvarado12@gmail.com

## License
# The follow script was created under the GNU/GPLv2. license.
# http://www.gnu.org/licenses/old-licenses/gpl-2.0-standalone.html

## Title
# Endemic_MultindexCalculation

## Description
# This R script perfomace a kind of Jackknife focused only in Endemic species. It remove those species from the indeces calculation aa% in ii times.
# an Endemic species it is one that have a distribution restricted only inside and Area of endemism.
# The file requiere for this script is a absence/presence matrix. 
# The outcome file is 3 list corresponding to the three indeces. Each list have (NumberOfMatrices X NumberOfRemove%) slots and each one have ii Index values for each area

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

#################################################################################################


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
library(doParallel)


set.seed(100)

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

dir.data <- dir.data[c(1,3)]

multi.data <- list()

for(i in 1:length(dir.data)){
  multi.data[[i]] <- read.csv(dir.data[i]) 
}

# Set names for each data.
names(multi.data) <- dir.data


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

# Find all endemic species from the distriution data

end.sp <- c()

for (i in 1:length(multi.data[[1]]$especie)){
  
  if(length(which(multi.data[[1]][i,]==1))==1){
    #print(paste(multi.data[[1]]$especie[i],"es endémica"))
    end.sp <- c(end.sp,as.character(multi.data[[1]]$especie[i]))
  }
  
}

end.sp

# Use ony those endemic species thar are in both, distribution and phylogenies

end.sp <- end.sp[which(end.sp%in%match.sp)]
end.sp

# Define the removal percentage

perc.rem <- c(0,.25,.50,.75,1)

# Make a number sequences for the organization of the final outcome

sec <- seq(from=0,to=((length(perc.rem))*length(multi.data))-1,by=length(multi.data))

## Create three list where the results will put it

dt.sum <- list()
pd.sum <- list()
avdt.sum <- list()

## Create three list where the index summatories will put it

dt.origin <- list()
pd.origin <- list()
avdt.origin <- list()


no_cores <- detectCores() 
cl <- makeCluster(no_cores)
registerDoParallel(cl)

## For each removal percentage do...
for(aa in 1:length(perc.rem)){
  ## For each interation do...
  for (ii in 1:30){
    ## Resampling from the species distribution
    # remove the aa percentage of endemic species 
    rem.samp <- end.sp[sample(x = 1:length(end.sp), size = perc.rem[aa]*length(end.sp))]
    print(paste("Iteration",ii,sep = ":"))
    
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
        
        ## Remove species
        # The fastest way is put 0 in all the distribution of each endemic species 
        dist[which(dist$especie%in%rem.samp),2:length(colnames(dist))] <- 0
        
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
          avdt.origin[[j]] <- avdt2
        }else{
          
          avdt2[is.na(avdt2)] <- 0
          avdt.origin[[j]][is.na(avdt.origin[[j]])] <- 0
          
          avdt.origin[[j]] <- avdt.origin[[j]] + avdt2}
        
      }
    }
    
    for(xx in 1:length(dt.origin)){
      
      if (ii==1){
        
        dt.sum[[(xx+sec[aa])]] <- list(rep(NA,1))
        pd.sum[[(xx+sec[aa])]] <- list(rep(NA,1))
        avdt.sum[[(xx+sec[aa])]] <- list(rep(NA,1))
        
        dt.sum[[(xx+sec[aa])]][[1]] <- as.matrix(dt.origin[[xx]]$W)
        pd.sum[[(xx+sec[aa])]][[1]] <- as.matrix(pd.origin[[xx]]$PD)
        avdt.sum[[(xx+sec[aa])]][[1]] <- as.matrix(avdt.origin[[xx]]$Dplus)
      }
      
      dt.sum[[(xx+sec[aa])]][[1]] <- rbind(dt.sum[[(xx+sec[aa])]][[1]],as.matrix(dt.origin[[xx]]$W))
      pd.sum[[(xx+sec[aa])]][[1]] <- rbind(pd.sum[[(xx+sec[aa])]][[1]],as.matrix(pd.origin[[xx]]$PD))
      avdt.sum[[(xx+sec[aa])]][[1]] <- rbind(avdt.sum[[(xx+sec[aa])]][[1]],as.matrix(avdt.origin[[xx]]$Dplus))
      
    }
    
  }
  
  
}

################################################################################

Areas <- rep(dt.origin[[1]]$area,(ii+1))
Areas <- rep(Areas,length(perc.rem))
Areas

PP <- rep(perc.rem, each=length(dt.sum[[1]][[1]]))

sec2 <- seq(from=1,to=length(dt.sum),by=2)
sec2

td.val <- 0
pd.val <- 0
avtd.val <- 0

for (i in 1:length(sec2)){
  
  td.val <- rbind(td.val,dt.sum[[sec2[i]]][[1]])
  pd.val <- rbind(pd.val,pd.sum[[sec2[i]]][[1]])
  avtd.val <- rbind(avtd.val,avdt.sum[[sec2[i]]][[1]])
  
}

end1 <- cbind(td.val,pd.val,avtd.val)
end1 <- end1[-1,]

head(end1,5L)

end1 <- cbind(Areas,end1,PP)

colnames(end1) <- c("Area","TD","PD","AvTD","Percentage")

head(end1,5L)

###########################################################################

Areas <- rep(dt.origin[[2]]$area,(ii+1))
Areas <- rep(Areas,length(perc.rem))
Areas

PP <- rep(perc.rem, each=length(dt.sum[[2]][[1]]))

sec2 <- seq(from=0,to=length(dt.sum),by=2)
sec2 <- sec2[-1]
sec2

td.val <- 0
pd.val <- 0
avtd.val <- 0

for (i in 1:length(sec2)){
  
  td.val <- rbind(td.val,dt.sum[[sec2[i]]][[1]])
  pd.val <- rbind(pd.val,pd.sum[[sec2[i]]][[1]])
  avtd.val <- rbind(avtd.val,avdt.sum[[sec2[i]]][[1]])
  
}

end2 <- cbind(td.val,pd.val,avtd.val)
end2 <- end2[-1,]

head(end2,5L)

end2 <- cbind(Areas,end2,PP)

colnames(end2) <- c("Area","TD","PD","AvTD","Percentage")

head(end2,5L)

############################################################################

setwd("~/Documentos/Omar/Result_TesisOmar/Final2/Endemic/")

write.table(end1,"Area.end", quote = F, row.names = F, col.names = T, sep = ",")
write.table(end2,"Grid.end", quote = F, row.names = F, col.names = T, sep = ",")
