setwd("~/Documentos/Omar/Tesis/Taxa/Results/Julio1/")
# Read absence/presence matrix for the areas of endemism
area.dist <- read.csv("Area.dist.matrix",header = T)
# Total species used

n <- length(levels(area.dist2$especie))

# Check it
head(area.dist)
# empty vector for the summatory
total.sp <- c()
# Make the sumatory fo each area
for(i in 2:length(colnames(area.dist))){
  
  sp <- sum(area.dist[,i])
  
  total.sp <- c(total.sp,sp)
  
}
#Assign the names
names(total.sp) <- colnames(area.dist)[2:length(colnames(area.dist))]
# Check it
total.sp

########################################################################
########################################################################

## Load libraries.

library(ape)
library(phangorn)
library(phytools)

## Load the phylogenies and distributions


# Set working directory.
setwd("~/Documentos/Omar/Tesis/Taxa/Trees/")

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

########################################################################

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Julio1/")

g.info <- read.csv("General.info", header = T)


# Find the phylogeny species that match with the distribution species
match.sp <- area.dist$especie[which(area.dist$especie%in%g.info$Sp)]
match.sp
# Extract from the absence/presence matrix only the species that match with the species phylogeny
area.dist2 <- area.dist[which(match.sp%in%area.dist$especie),]
# empty vector for the summatory
phylo.sp <- c()
# Make the summatory
for(i in 2:length(colnames(area.dist2))){
  
  sp <- sum(area.dist2[,i])
  
  phylo.sp <- c(phylo.sp,sp)
  
}
# Assing the names
names(phylo.sp) <- colnames(area.dist2)[2:length(colnames(area.dist2))]
# Check it
phylo.sp

#################################################################
# Endemic SP
#################################################################

##################################################
# First, find all endemic

end.sp <- c()

for (i in 1:length(area.dist$especie)){
  
  if(length(which(area.dist[i,]==1))==1){
    #print(paste(multi.data[[1]]$especie[i],"es endémica"))
    end.sp <- c(end.sp,as.character(area.dist$especie[i]))
  }
  
}

end.sp

area.dist3 <- area.dist[which(area.dist$especie%in%end.sp),]

str(area.dist3)

end.d <- c()

for(i in 2:length(colnames(area.dist3))){
  
  sp <- sum(area.dist3[,i])
  
  end.d <- c(end.d,sp)
  
}
# Assing the names
names(end.d) <- colnames(area.dist2)[2:length(colnames(area.dist2))]
# Check it
end.d

# Use ony those endemic species thar are in both, distribution and phylogenies

end.match <- end.sp[which(end.sp%in%match.sp)]
end.match

area.dist4 <- area.dist3[which(area.dist3$especie%in%end.match),]

str(area.dist4)

end.dMatch <- c()

for(i in 2:length(colnames(area.dist4))){
  
  sp <- sum(area.dist4[,i])
  
  end.dMatch <- c(end.dMatch,sp)
  
}
# Assing the names
names(end.dMatch) <- colnames(area.dist2)[2:length(colnames(area.dist2))]
# Check it
end.dMatch


################################################################
# Branch Length per Area
################################################################


area.bl <- area.dist[which(area.dist$especie%in%g.info$Sp),]

for(i in 2:length(rownames(area.bl))){
  
  bl <- g.info$BL[grep(area.bl$especie[i],g.info$Sp)]
  
  if(any(1%in%area.bl[i,])){
  
  area.bl[i,which(area.bl[i, ]==1)] <- bl

  }
  
}

area.bl

total.bl <- c()

for(i in 2:length(colnames(area.bl))){
  
  sp <- sum(area.bl[,i])
  
  total.bl <- c(total.bl,sp)
  
}

names(total.bl) <- colnames(area.bl)[2:length(colnames(area.bl))]

total.bl

#########################################

area.bl2 <- area.bl[which(area.bl$especie%in%end.match), ]

end.bl <- c()

for(i in 2:length(colnames(area.bl2))){
  
  sp <- sum(area.bl2[,i])
  
  end.bl <- c(end.bl,sp)
  
}

names(end.bl) <- colnames(area.bl2)[2:length(colnames(area.bl2))]

end.bl

########################################################3

# Calculate the percentage of the species of the phylogeny given the total species
percentage1 <- (phylo.sp/total.sp)*100
percentage1

percentage2 <- (end.dMatch/end.d)*100

percentage3 <- (end.d/total.sp)*100

percentage4 <- (end.dMatch/phylo.sp)*100

percentage5 <- (total.bl/phylo.sp)

# Join the tree vectors in a unique matrix
final.cal <- rbind(total.sp, 
                   phylo.sp, 
                   percentage1,
                   end.d,
                   end.dMatch,
                   percentage2,
                   percentage3,
                   percentage4,
                   total.bl,
                   end.bl,
                   percentage5)
# Assing row names
rownames(final.cal) <- c("TotalSp", 
                         "PhyloSp",
                         "PhyloSp/TotalSp %",
                         "TotalEnd",
                         "PhyloEnd",
                         "PhyloEnd/TotalEnd %",
                         "TotalEnd/TotalSp %",
                         "PhyloEnd/PhyloSp %",
                         "TotalBL",
                         "EndBL",
                         "TotalBL/PhyloSp")
# Transpose the matrix for a easy read
final.cal <- as.data.frame(t(final.cal))

final.cal <- final.cal[order(final.cal$TotalSp, decreasing = T),]

final.cal
# Save the results
setwd("~/Documentos/Omar/Tesis/Taxa/Results/Julio1/")

write.table(final.cal, file = "Areas_sp.info",sep = ",", row.names = T, col.names = T, quote = F)
