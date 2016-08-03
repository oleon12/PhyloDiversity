setwd("~/Documentos/Omar/Tesis/Taxa/Results/Julio1/")
# Read absence/presence matrix for the areas of endemism
area.dist <- read.csv("Area.dist.matrix",header = T)
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
match.sp <- dead.pool[which(dead.pool%in%area.dist$especie)]
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

# Calculate the percentage of the species of the phylogeny given the total species
percentage <- phylo.sp/total.sp
percentage
# Join the tree vectors in a unique matrix
final.cal <- rbind(total.sp, phylo.sp, percentage)
# Assing row names
rownames(final.cal) <- c("Total.sp", "Phylo.sp","Phylo/Total")
# Transpose the matrix for a easy read
final.cal <- t(final.cal)

# Save the results
setwd("~/Documentos/Omar/Tesis/Taxa/Results/Julio1/")

write.table(final.cal, file = "Areas_sp.info",sep = ",", row.names = T, col.names = T, quote = F)
