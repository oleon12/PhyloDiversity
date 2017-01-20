# Load libraries
library(vegan)
library(maptools)
library(RColorBrewer)
library(classInt)
library(maps)
library(maptools)
library(dismo)

## Go to the work directory where are the shape poly 
setwd("/home/omar/Documentos/Omar/Tesis/Scripts/Distribution/shp/Grid/")
# Read the shape file
grid <- readShapePoly("grid_25g.shp")
#
# Now, go to the directory where are the index results
setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")
# Read the results depending on the cell size
x <- read.csv("AvTD.grid25",header=T)
# Extract the index values
area <- read.csv("DT.grid25")$area 
x <- x$Dplus
# Check the values
x
# Due to there are many cell without index values (this is beacuse the species distribution)
# remove the 0 values, or the quantile classification will be biased
x2 <- x[which(x>0.0000000)]

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
Pos5 <- grep(5,class)
Pos4 <- grep(4,class)
Pos3 <- grep(3,class)
Pos2 <- grep(2,class)
NonVal <- which(x==0)

class[-Pos5] <-  "Non-Q5"
# Now, those cell int in the two last quantiles, will be filled with the quotes "Q5" and "Q4"
class[Pos5] <- "Q5"

ClassAvTD1 <- class

class[Pos5] <- "Q5"
class[Pos4] <- "Q4"
class[Pos3] <- "Q3"
class[Pos2] <- "Q2"
class[NonVal] <- "NO"

ClassAvTD2 <- class

grid$Index <- ClassAvTD2

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

writePolyShape(grid,fn = "Grid25_AvTD")

########################################################
########################################################

## Go to the work directory where are the shape poly 
setwd("/home/omar/Documentos/Omar/Tesis/Scripts/Distribution/shp/Grid/")
# Read the shape file
grid <- readShapePoly("grid_25g.shp")
#
# Now, go to the directory where are the index results
setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")
# Read the results depending on the cell size
x <- read.csv("PD.grid25",header=T)
# Extract the index values
area <- read.csv("DT.grid25")$area 
x <- x$PD
# Check the values
x
# Due to there are many cell without index values (this is beacuse the species distribution)
# remove the 0 values, or the quantile classification will be biased
x2 <- x[which(x>0.0000000)]

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
Pos5 <- grep(5,class)
Pos4 <- grep(4,class)
Pos3 <- grep(3,class)
Pos2 <- grep(2,class)
NonVal <- which(x==0)

class[-Pos5] <-  "Non-Q5"
# Now, those cell int in the two last quantiles, will be filled with the quotes "Q5" and "Q4"
class[Pos5] <- "Q5"

ClassPD1 <- class

class[Pos5] <- "Q5"
class[Pos4] <- "Q4"
class[Pos3] <- "Q3"
class[Pos2] <- "Q2"
class[NonVal] <- "NO"

ClassPD2 <- class

grid$Index <- ClassPD2

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

writePolyShape(grid,fn = "Grid25_PD")


#####################################################
#####################################################

## Go to the work directory where are the shape poly 
setwd("/home/omar/Documentos/Omar/Tesis/Scripts/Distribution/shp/Grid/")
# Read the shape file
grid <- readShapePoly("grid_25g.shp")
#
# Now, go to the directory where are the index results
setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")
# Read the results depending on the cell size
x <- read.csv("DT.grid25",header=T)
# Extract the index values
area <- read.csv("DT.grid25")$area 
x <- x$W
# Check the values
x
# Due to there are many cell without index values (this is beacuse the species distribution)
# remove the 0 values, or the quantile classification will be biased
x2 <- x[which(x>0.0000000)]

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
Pos5 <- grep(5,class)
Pos4 <- grep(4,class)
Pos3 <- grep(3,class)
Pos2 <- grep(2,class)
NonVal <- which(x==0)

class[-Pos5] <-  "Non-Q5"
# Now, those cell int in the two last quantiles, will be filled with the quotes "Q5" and "Q4"
class[Pos5] <- "Q5"

ClassTD1 <- class

class[Pos5] <- "Q5"
class[Pos4] <- "Q4"
class[Pos3] <- "Q3"
class[Pos2] <- "Q2"
class[NonVal] <- "NO"

ClassTD2 <- class

grid$Index <- ClassTD2

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

writePolyShape(grid,fn = "Grid25_AvTD")


######################################################################
######################################################################

AllClass <- cbind(as.matrix(ClassAvTD1),
                  as.matrix(ClassPD1),
                  as.matrix(ClassTD1))

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

write.csv(AllClass, "AllIndices_Quantiles.csv", quote = F, row.names = F, col.names = T)
