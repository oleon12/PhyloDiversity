library(maptools)


## Set working directory
setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Results/Final2/")

## Read the absence/presence tables for the total richness

grid.rich <- read.csv("grid_25g.dist.matrix")

## Create two vectors where the total richness will are

grid.t.rich <- c()

## Make the summatory of the richness (This could be made by a apply function, check it later)

grid.zero <- c()

for(i in 2:length(colnames(grid.rich))){
  total <- sum(grid.rich[,i])
  grid.t.rich[i] <- total
  
  if(total==0){
    grid.zero <- c(grid.zero, colnames(grid.rich)[i])
  }
}

grid.t.rich <- grid.t.rich[-1]


# Check the results

grid.t.rich


###################################################################################

setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

index.grid <- data.frame(Cell.Grid=read.csv("DT.grid25")$area,
                         TD=read.csv("DT.grid25")$W,
                         TDs=read.csv("DT.grid25")$Ws,
                         TDe=read.csv("DT.grid25")$We,
                         TDse=read.csv("DT.grid25")$Wse,
                         PD=read.csv("PD.grid25")$PD,
                         AvTD=read.csv("AvTD.grid25")$Dplus,
                         f.Rich=read.csv("DT.grid25")$rich,
                         t.Rich=grid.t.rich)

head(index.grid)

CellClass <- readShapePoly("Grid25_AvTD2.shp")$Index

CellClass

index.grid <- cbind(index.grid,CellClass)

str(index.grid)


#################################################################################
#################################################################################

Xrich <- index.grid$t.Rich[grep("Q5",index.grid$CellClass)]

Yind <- index.grid$PD[grep("Q5",index.grid$CellClass)]

nSubj <- length(Xrich)

BayesSlope(x = Xrich, y = Yind, nSubj = nSubj, ploting = T, SaveName = "Q5PDvsRich")

#################################################################################
################################################################################

Xrich <- index.grid$t.Rich[grep("Q5",index.grid$CellClass)]

Yind <- index.grid$AvTD[grep("Q5",index.grid$CellClass)]

nSubj <- length(Xrich)

BayesSlope(x = Xrich, y = Yind, nSubj = nSubj, ploting = T, SaveName = "Q5AvTDvsRich")


############################################################################################3


#################################################################################
#################################################################################

Xrich <- index.grid$t.Rich[grep("Q2",index.grid$CellClass)]

Yind <- index.grid$PD[grep("Q2",index.grid$CellClass)]

nSubj <- length(Xrich)

BayesSlope(x = Xrich, y = Yind, nSubj = nSubj, ploting = T, SaveName = "Q2PDvsRich2")

#################################################################################
################################################################################

Xrich <- index.grid$t.Rich[grep("Q2",index.grid$CellClass)]

Yind <- index.grid$AvTD[grep("Q2",index.grid$CellClass)]

nSubj <- length(Xrich)

BayesSlope(x = Xrich, y = Yind, nSubj = nSubj, ploting = T, SaveName = "Q2AvTDvsRich2")


######################################################################################

#################################################################################
#################################################################################

nonQ5 <- c(grep("Q4",index.grid$CellClass),grep("Q3",index.grid$CellClass),grep("Q2",index.grid$CellClass),grep("Q1",index.grid$CellClass))

Xrich <- index.grid$t.Rich[nonQ5]

Yind <- index.grid$PD[nonQ5]

nSubj <- length(Xrich)

BayesSlope(x = Xrich, y = Yind, nSubj = nSubj, ploting = T, SaveName = "nonQ5PDvsRich2")

#################################################################################
################################################################################

Xrich <- index.grid$t.Rich[nonQ5]

Yind <- index.grid$AvTD[nonQ5]

nSubj <- length(Xrich)

BayesSlope(x = Xrich, y = Yind, nSubj = nSubj, ploting = T, SaveName = "nonQ5AvTDvsRich2")


################################################################################


################################################################################

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

Non0Grid <- read.csv("AllIndex_Table.csv")$Cell.Grid

index.grid2 <- index.grid[which(index.grid$Cell.Grid%in%Non0Grid),]

Xrich <- index.grid2$t.Rich

Yind <- index.grid2$TDs

nSubj <- length(Xrich)

BayesSlope(x = Xrich, y = Yind, nSubj = nSubj, ploting = T, SaveName = "TDsvsRich2")


Yind <- index.grid2$TDe

BayesSlope(x = Xrich, y = Yind, nSubj = nSubj, ploting = T, SaveName = "TDevsRich2")

