library(maptools)
library(ggplot2)

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

CellClass <- read.csv("AllIndices_Quantiles.csv", header = T)

head(CellClass)

index.grid <- cbind(index.grid,CellClass)

str(index.grid)

############################################################################

setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Results/Final2/")

GI <- read.csv("General.info")

str(GI)

GI <- GI[ ,1:2]

str(GI)

GI <- GI[which(GI$Sp%in%grid.rich$especie),]

str(GI)

grid.bl <- grid.rich

pb <- txtProgressBar(min = 0, max = length(GI$Sp), style = 3)
for(i in 1:length(GI$Sp)){
  
  gridRow <- grep(GI$Sp[i],grid.bl$especie)
  
  grid.bl[grep(1,grid.bl[gridRow,])] <- GI$BL[i]
  
  setTxtProgressBar(pb, i)
}

close(pb)

bl <- apply(grid.bl, 2, sum)

bl <- bl[which(names(bl)%in%index.grid$Cell.Grid)]

index.all <- cbind(index.grid,bl)

########################################################################3

ind <- index.all[which(index.all$f.Rich>0.00000), ]

ggplot()+
  geom_point(aes(x=log(ind$f.Rich), y=log(ind$AvTD), shape=ind$AvTDclass, colour=ind$bl))+
  scale_color_gradient(low="blue", high = "red")


ggplot()+
  geom_point(aes(x=log(ind$bl), y=log(ind$AvTD), colour=ind$f.Rich))+
  scale_color_gradient(low="green", high = "blue")



ggplot()+
  geom_point(aes(x=log(ind$bl), y=log(ind$PD), colour=ind$f.Rich))+
  scale_color_gradient(low="green", high = "blue")
###############################################################################

indQ5 <- ind[which(ind$AvTDclass=="Q5"), ]

PRvsInd <-  ggplot()+
  geom_point(aes(x=log(indQ5$f.Rich), y=log(indQ5$AvTD), colour=indQ5$bl), size=4)+
  scale_color_gradient(low="blue", high = "red", name="BL")+
  xlab("Phylogenetic Richness")+ylab("AvTD")+
  theme(legend.position = c(.2,.7),
        legend.title=element_text(size = 30,face = "bold"),
        legend.text = element_text(size = 30),
        legend.key.size = unit(1, "cm"),
        legend.key = element_rect(fill= NULL,colour="black"),
        legend.background = element_rect(fill=NA))+
  theme( axis.text.y=element_text(size=15),
         axis.text.x=element_text(size=15),
         axis.title.y=element_text(size=30),
         axis.title.x=element_text(size=30))
PRvsInd

png("PRvsInd_Q5.png",width = 800, height = 800)
PRvsInd
dev.off()

BLvsInd <- ggplot()+
  geom_point(aes(x=log(indQ5$bl), y=log(indQ5$AvTD), colour=indQ5$f.Rich), size=4)+
  scale_color_gradient(low="green", high = "blue", name="PR")+
  xlab("Accumulated Branch Length")+ylab("AvTD")+
  theme(legend.position = c(.4,.6),
        legend.title=element_text(size = 30,face = "bold"),
        legend.text = element_text(size = 30),
        legend.key.size = unit(1, "cm"),
        legend.key = element_rect(fill= NULL,colour="black"),
        legend.background = element_rect(fill=NA))+
  theme( axis.text.y=element_text(size=15),
         axis.text.x=element_text(size=15),
         axis.title.y=element_text(size=30),
         axis.title.x=element_text(size=30))

BLvsInd

png("BLvsInd_Q5.png",width = 800, height = 800)
BLvsInd
dev.off()

