library(rgdal)
library(rgeos)
library(ggplot2)
library(SDMTools)
library(maptools)
library(maps)

#################################################

nab <- readOGR("../shp/NAB.shp")
grid <- readOGR("../shp/grid_25g.shp")

abpre <- read.csv("../Distributions/grid_25g.dist.matrix")
Sp <- rownames(abpre)
SpCoords <- read.csv("../Distributions/Richness_end.2019.occ", row.names = 1)

nab2 <- fortify(nab)
grid2 <- fortify(grid)

intPol <- over(grid,nab)
posInNab <- which(intPol$id==1)

abpre <- abpre[,posInNab]
totalGrid <- dim(abpre)[2]

##############################################################

col <- rep("blue",length(intPol$id))
col[posInNab] <- "red"
col <- rep(col, each=5)

ggplot()+
  geom_polygon(aes(x=grid2$long,y=grid2$lat,group=grid2$group, fill=col))+
  coord_fixed()

##############################################################

DistGrid <- c()

SpCells <- as.matrix(apply(abpre, 1, function(x){DistGrid <- c(DistGrid,sum(x))}))
Percent <- SpCells/totalGrid

SpCells <- as.data.frame(cbind(SpCells, Percent))
colnames(SpCells) <- c("NoCells","Percentage")

SpCells <- SpCells[which(SpCells$Percentage<0.25),]

centroGrid <- gCentroid(grid, byid = T)
Dist <- pointDistance(centroGrid@coords[1:5,], centroGrid@coords[1:5,], lonlat = T, allpairs = T)
minDist25g <- Dist[2,1]

minDist25g


###############################################################

SpCoords <- SpCoords[which(SpCoords$Sp%in%rownames(SpCells)), ]

End <- c()

for(i in 1:length(SpCells$NoCells)){
  
  Sp1 <- SpCoords[which(SpCoords$Sp%in%rownames(SpCells)[i]),]
  
  if(length(Sp1[,1])==1){
    End <- c(End, as.character(unique(Sp1$Sp)))
  }
  if(length(Sp1[,1])>1){
    
    Spc <- Sp1[,2:3]
    DistP <- pointDistance(Spc,Spc, lonlat = T, allpairs = T)
  
    mean1 <- apply(DistP, 1, mean)
    mean2 <- mean(mean1)
  
    if(mean2<minDist25g*2){
      End <- c(End, as.character(unique(Sp1$Sp)))
    }
  }
}

End

############################################################

for(i in 1:length(End)){

  fileN <- paste(End[i],".png", sep = "")
  fileN <- paste("../Data/Endemica/", fileN, sep = "")
  Coo <- SpCoords[which(SpCoords$Sp%in%End[i]),2:3]
  
  ggplot()+
    geom_polygon(aes(x=nab2$long,y=nab2$lat,group=nab2$group))+
    geom_point(aes(x=Coo[,1],y=Coo[,2]),colour="red")+
    xlab("Lon")+ylab("Lat")+
    coord_fixed()
    
  ggsave(filename = fileN, device = "png")
}

############################################################

abpre <- read.csv("..Data/Distributions/grid_25g.dist.matrix")

abpre <- abpre[which(rownames(abpre)%in%End), ]
 
write.csv(abpre, "..Data/Distributions/Endemic_25g.matrix")
