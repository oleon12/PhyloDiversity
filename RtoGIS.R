# Load libraries

library(maptools)
library(RColorBrewer)
library(classInt)
library(maps)
library(maptools)
library(dismo)

## Go to the work directory where are the shape poly 
setwd("/home/omar/Documentos/Omar/Tesis/Scripts/Distribution/shp/Grid/")
# Read the shape file
grid <- readShapePoly("grid_1g.shp")
#
# Now, go to the directory where are the index results
setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Results/Julio1/RawIndex/")
# Read the results depending on the cell size
x <- read.csv("DT.grid1",header=T)
# Extract the index values
x <- x$W
# Check the values
x
# Due to there are many cell without index values (this is beacuse the species distribution)
# remove the 0 values, or the quantile classification will be biased
x2 <- x[-grep(0,x)]

# No, do the quantile clasification, in this case, given the values index, the intervals will be generated
brks <- classIntervals(x2,n=5,style = "quantile")
brks <- brks$brks
# Here, the intervarls given the index value
brks

# Now, each index values will be classified in a category given the quantile intervals
# Here, the all index values (include those cells with index values equal 0)
# Just, beacuase the final shape file neead a value for each cell
class <- findInterval(x,brks,all.inside = T)
# Now, beacuase the focus is only in the quantile 5 and 4, (those with highest index values)
# The cells out of the Q5 and Q4 will be replace with NO
class[-c(grep(5,class),grep(4,class))] <-  "NO"
# Now, those cell int in the two last quantiles, will be filled with the quotes "Q5" and "Q4"
class[grep(5,class)] <- "Q5"
class[grep(4,class)] <- "Q4"
# See the result
class
#
# Now, the hard part...
# 

 vegdist()
grid$Index <- class

setwd("~/Documentos/Omar/Tesis/Taxa/Results/May18/Raw_IndexR/")

writePolyShape(grid,fn = "Grid50_PD")
