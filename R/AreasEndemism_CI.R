## Autor
# Leon-Alvarado, Omar Daniel.
# leon.alvarado12@gmail.com

## License
# The follow script was created under the GNU/GPLv2. license.
# http://www.gnu.org/licenses/old-licenses/gpl-2.0-standalone.html

## Title
# Complementarity index with areas of endemism

## Description
# This script calculate the complementarity values for all areas of endemism, and then extract the complementarity
# areas given the most important area of endemism, which in this case is Magdalena

###################################################################################

# Load the libraries

library(vegan)

#Set the working directory and read the absence/presence matrix for the areas of endemism

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final/")

areas <- read.csv("Area.dist.matrix")

# Extract the species names and let it as rownames, then transponse the matrix
# Finally, the colnames will be the species and the rownames will be the areas of endemism

name <- areas$especie 

areas <- areas[,-1] 

rownames(areas) <- name 

areas <- t(areas)

areas

# Calculated the jaccard index using the vegdist function from the "vegan" package

d <- as.matrix(vegdist(areas, method = "jaccard"))

# Extract all colums only from the row which correspond with the Magdalena area

magdalenaCI <- d[grep("Magdalena",rownames(areas)),]

# Transform the vector to a matrix and change the colum names

magdalenaCI <- as.matrix(magdalenaCI[order(magdalenaCI,decreasing = T)])

colnames(magdalenaCI) <- "CI"

magdalenaCI

write.csv(magdalenaCI,"~/Documentos/Omar/Tesis/Taxa/Results/Final/Raw_IndexR/AreaEndemis.CI",
          quote = F, row.names = T)
