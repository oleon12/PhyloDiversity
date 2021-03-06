## Autor
# Leon-Alvarado, Omar Daniel.
# leon.alvarado12@gmail.com

## License
# The follow script was created under the GNU/GPLv2. license.
# http://www.gnu.org/licenses/old-licenses/gpl-2.0-standalone.html

## Title
# Cells of Grid Comparisons. 

## Description
# The R cuantify the percentage of shared cells of grid given a quantile classification.
# The construction of this script was specific for only three phylogenetic index (DT, PD and AvTDT),
# so, if you want to make comparisons with more or less indexes, you must need modify the script.
#
# First, a classification of the values of the three index are made using quantiles.
# then, looking for those cells who was classified in the two last quantiles (4 and 5) for each index
# and makes a pairwise comparison index, cuantifying the percentage of shared cells.
#
# Pairwise comparisons:
# 1.) DT + PD
# 2.) DT + AvTD
# 3.) PD + AvTD
#
# The results is a matrix with three shared percentage: A%, B% and A+B% for each comparison
# A% = Total shared cells / Total Cells of the Left index
# B% = Total shared cells / Total Cells of the Right index
# A+B% = Total shared cells / (Total Cells of Left index + Total Cells of the Right index)
#
# How to know what is the right or left index ?
# e.g. DT + PD == Left index + Rigth index


## Cell index values comparison 
library(classInt)

## Set working directory
setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

## Create a data frame with index values and cell ID
index.grid <- data.frame(Cells=read.csv("DT.grid25")$area,
                         DT=read.csv("DT.grid25")$W,
                         PD=read.csv("PD.grid25")$PD,
                         AVDT=read.csv("AvTD.grid25")$Dplus)

str(index.grid)

## Matrix outcome where the interval classification will put for each cell index value
r1 <- as.data.frame(matrix(NA, nrow = length(index.grid$Cells),ncol = length(colnames(index.grid))))
colnames(r1) <- colnames(index.grid)
r1$Cells <- index.grid$Cells ##


## Classification
for(i in 2:length(colnames(index.grid))){
  ## Made the quatiles range
  
  ind2 <- index.grid[which(index.grid[,i]>0.000000),i]
  
  brks <- classIntervals(ind2,n=5,style = "quantile")
  brks <- brks$brks
  brks
  ## Classification given quantiles intervals.
  r1[,i] <- findInterval(index.grid[,i],brks)
}

r1


## Identify the important cells (cell with inde values inside the 5 and 4 quantiles)
dt.imp <- c(grep("5",r1$DT),grep("4",r1$DT))
pd.imp <- c(grep("5",r1$PD),grep("4",r1$PD))
avtd.imp <- c(grep("5",r1$AVDT),grep("4",r1$AVDT))


## Comparisons and results

## Create outcome matrix

r2 <- matrix(NA,nrow = 3, ncol = 3)
rownames(r2) <- c("TD+PD","TD+AvDT","PD+AvDT")
colnames(r2) <- c("A%","B%","A+B%")

## DT v PD
r2[1,1] <- round((length(which(dt.imp%in%pd.imp))/ length(dt.imp))*100,digits = 2)# NoShared/TotalDT 
r2[1,2] <- round((length(which(dt.imp%in%pd.imp))/ length(pd.imp))*100,digits = 2)# NoShared/TotalPD
r2[1,3] <- round((length(which(dt.imp%in%pd.imp))/ (length(dt.imp)+length(pd.imp)))*100,digits = 2)# NoShared/Total(DT+PD)
# DT v AvTD
r2[2,1] <- round((length(which(dt.imp%in%avtd.imp))/ length(dt.imp))*100,digits = 2)# NoShared/TotalDT 
r2[2,2] <- round((length(which(dt.imp%in%avtd.imp))/ length(avtd.imp))*100,digits = 2)# NoShared/TotalAvTD
r2[2,3] <- round((length(which(dt.imp%in%avtd.imp))/ (length(dt.imp)+length(avtd.imp)))*100,digits = 2)
# PD v AvTD
r2[3,1] <- round((length(which(pd.imp%in%avtd.imp))/ length(pd.imp))*100,digits = 2)# NoShared/TotalPD
r2[3,2] <- round((length(which(pd.imp%in%avtd.imp))/ length(avtd.imp))*100,digits = 2)# NoShared/TotalAvTD
r2[3,3] <- round((length(which(pd.imp%in%avtd.imp))/ (length(pd.imp)+length(avtd.imp)))*100,digits = 2)# NoShared/Total(PD+AvTD)

r2

write.table(r2,file = "~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/Index.ShareCells",
            quote = F, row.names = T, col.names = T)

##########################################################################################3

r3 <- matrix(NA,nrow = 3, ncol = 3)
rownames(r3) <- c("PD+TD","AvDT+TD","AvDT+PD")
colnames(r3) <- c("A%","B%","A+B%")

## PD v TD
r3[1,1] <- round((length(which(pd.imp%in%dt.imp))/ length(pd.imp))*100,digits = 2)# NoShared/TotalDT 
r3[1,2] <- round((length(which(pd.imp%in%dt.imp))/ length(dt.imp))*100,digits = 2)# NoShared/TotalPD
r3[1,3] <- round((length(which(pd.imp%in%dt.imp))/ (length(dt.imp)+length(pd.imp)))*100,digits = 2)# NoShared/Total(DT+PD)
# AvTD v TD
r3[2,1] <- round((length(which(avtd.imp%in%dt.imp))/ length(avtd.imp))*100,digits = 2)# NoShared/TotalDT 
r3[2,2] <- round((length(which(avtd.imp%in%dt.imp))/ length(dt.imp))*100,digits = 2)# NoShared/TotalAvTD
r3[2,3] <- round((length(which(avtd.imp%in%dt.imp))/ (length(dt.imp)+length(avtd.imp)))*100,digits = 2)
# AvTD v PD
r3[3,1] <- round((length(which(avtd.imp%in%pd.imp))/ length(avtd.imp))*100,digits = 2)# NoShared/TotalPD
r3[3,2] <- round((length(which(avtd.imp%in%pd.imp))/ length(pd.imp))*100,digits = 2)# NoShared/TotalAvTD
r3[3,3] <- round((length(which(avtd.imp%in%pd.imp))/ (length(pd.imp)+length(avtd.imp)))*100,digits = 2)# NoShared/Total(PD+AvTD)

r3

write.table(r3,file = "~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/Index2.ShareCells",
            quote = F, row.names = T, col.names = T)


r2
r3
