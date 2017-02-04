#Load libraries
library(vegan)
library(maptools)
library(RColorBrewer)
library(classInt)
library(maps)
library(maptools)
library(dismo)



setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/")
# Read absence/presence matrix for the areas of endemism
area.dist <- read.csv("grid_25g.dist.matrix",header = T)

# Total species used
n <- length(levels(area.dist$especie))

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

################################################################################
################################################################################

## Go to the work directory where are the shape poly 
setwd("/home/omar/Documentos/Omar/Tesis/Scripts/Distribution/shp/Grid/")
# Read the shape file
grid <- readShapePoly("grid_25g.shp")


#area <- x$area 

x <- total.sp
# Check the values
x
# Due to there are many cell without index values (this is beacuse the species distribution)
# remove the 0 values, or the quantile classification will be biased
x2 <- x[which(x>0)]

# No, do the quantile clasification, in this case, given the values index, the intervals will be generated
brks <- classIntervals(x2,n=4,style = "jenks")
brks <- brks$brks
# Here, the intervarls given the index value
brks

# Now, each index values will be classified in a category given the quantile intervals
# Here, the all index values (include those cells with index values equal 0)
# Just, beacuase the final shape file neead a value for each cell
class <- findInterval(x,brks,all.inside = T)

names(class) <- names(total.sp)

# Now, beacuase the focus is only in the quantile 5 and 4, (those with highest index values)
# The cells out of the Q5 and Q4 will be replace with NO

#class[-c(grep(5,class),grep(4,class))] <-  "NO"

# Now, those cell int in the two last quantiles, will be filled with the quotes "Q5" and "Q4"

class[-which(x>0)] <- "NO"

# See the result
class
#

grid$Richness <- class

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R//")

writePolyShape(grid,fn = "Richness_Grid25_TD")
