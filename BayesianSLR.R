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

CellClass <- read.csv("AllIndices_Quantiles.csv", header = T)

head(CellClass)

index.grid <- cbind(index.grid,CellClass)

str(index.grid)


BayesSlope(x = index.grid$t.Rich, y=index.grid$TD, nSubj = length(index.grid[,1]), 
           ploting = T, SaveName = "tRichVsTD", Xlab = "Richness ~ TD")

Q5pos <- which(index.grid$TDclass=="Q5")

BayesSlope(x = index.grid$t.Rich[Q5pos], y=index.grid$TD[Q5pos], 
           nSubj = length(Q5pos), 
           ploting = T, SaveName = "Q5tRichVsTD", Xlab = "Q5 Richness ~ TD")

BayesSlope(x = index.grid$t.Rich[-Q5pos], y=index.grid$TD[-Q5pos], 
           nSubj = length(which(index.grid$TDclass=="Non-Q5")), 
           ploting = T, SaveName = "NonQ5tRichVsTD", Xlab = "Non-Q5 Richness ~ TD")

library(ggplot2)

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

TDcor <- ggplot(data = index.grid)+
  geom_point(aes(x=log(t.Rich),y=log(TD), colour=TDclass))+
  geom_abline(aes(slope = .98, intercept = 0, colour="All"), size=1)+
  geom_abline(aes(slope = .91, intercept= 0, colour="Q5"), size=1)+
  geom_abline(aes(slope=.97, intercept = 0, colour= "Non-Q5"), size=1)+
  xlab("log Total Richness")+ylab("log TD")+
  geom_text(aes(x=5,y=.8),label="Slope All = 0.99", size=7)+
  geom_text(aes(x=5,y=0),label="Slope Non-Q5 = 0.97", size=7)+
  geom_text(aes(x=5,y=.4),label="Slope Q5 = 0.91", size=7)+
  scale_color_discrete(name="Slopes")+
  theme(legend.position=c(.8,.35),
        legend.title=element_text(size = 30,face = "bold"),
        legend.text = element_text(size = 30),
        legend.key.size = unit(1, "cm"),
        legend.key = element_rect(fill= "white",colour="black"),
        axis.text.y=element_text(size=15),
        axis.text.x=element_blank(),
        axis.title.y=element_text(size=30),
        axis.title.x=element_blank())+  
  guides(colour = guide_legend(override.aes = list(size=2)))


png("TDvsTR_Corr.png", width =2000 , height =1200 ,res = 200)
TDcor
dev.off()

###
###

BayesSlope(x = index.grid$f.Rich, y=index.grid$TD, nSubj = length(index.grid[,1]), 
           ploting = T, SaveName = "fRichVsTD", Xlab = "Richness ~ TD")

Q5pos <- which(index.grid$TDclass=="Q5")

BayesSlope(x = index.grid$f.Rich[Q5pos], y=index.grid$TD[Q5pos], 
           nSubj = length(Q5pos), 
           ploting = T, SaveName = "Q5fRichVsTD", Xlab = "Q5 Richness ~ TD")

BayesSlope(x = index.grid$f.Rich[-Q5pos], y=index.grid$TD[-Q5pos], 
           nSubj = length(which(index.grid$TDclass=="Non-Q5")), 
           ploting = T, SaveName = "NonQ5fRichVsTD", Xlab = "Non-Q5 Richness ~ TD")

library(ggplot2)

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

TDcor <- ggplot(data = index.grid)+
  geom_point(aes(x=log(f.Rich),y=log(TD), colour=TDclass))+
  geom_abline(aes(slope = .99, intercept = 0, colour="All"), size=1)+
  geom_abline(aes(slope = .97, intercept= 0, colour="Q5"), size=1)+
  geom_abline(aes(slope=.98, intercept = 0, colour= "Non-Q5"), size=1)+
  xlab("log Phylogenetic Richness")+ylab("log TD")+
  geom_text(aes(x=5,y=.8),label="Slope All = 0.99", size=7)+
  geom_text(aes(x=5,y=0),label="Slope Non-Q5 = 0.98", size=7)+
  geom_text(aes(x=5,y=.4),label="Slope Q5 = 0.97", size=7)+
  scale_color_discrete(name="Slopes")+
  theme(legend.position=c(.8,.35),
        legend.title=element_text(size = 30,face = "bold"),
        legend.text = element_text(size = 30),
        legend.key.size = unit(1, "cm"),
        legend.key = element_rect(fill= "white",colour="black"),
        axis.text.y=element_text(size=15),
        axis.text.x=element_blank(),
        axis.title.y=element_text(size=30),
        axis.title.x=element_blank())+  
  guides(colour = guide_legend(override.aes = list(size=2)))


png("TDvsPR_Corr.png", width =2000 , height =1200 ,res = 200)
TDcor
dev.off()


#####################################################################################
#####################################################################################
#####################################################################################

BayesSlope(x = index.grid$t.Rich, y=index.grid$PD, nSubj = length(index.grid[,1]), 
           ploting = T, SaveName = "tRichVsPD", Xlab = "Richness ~ PD")

Q5pos <- which(index.grid$PDclass=="Q5")

BayesSlope(x = index.grid$t.Rich[Q5pos], y=index.grid$PD[Q5pos], 
           nSubj = length(Q5pos), 
           ploting = T, SaveName = "Q5tRichVsPD", Xlab = "Q5 Richness ~ PD")

BayesSlope(x = index.grid$t.Rich[-Q5pos], y=index.grid$PD[-Q5pos], 
           nSubj = length(which(index.grid$PDclass=="Non-Q5")), 
           ploting = T, SaveName = "NonQ5tRichvsPD", Xlab = "Non-Q5 Richness ~ PD")

library(ggplot2)

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

PDcor <- ggplot(data = index.grid)+
  geom_point(aes(x=log(t.Rich),y=log(PD), colour=PDclass))+
  geom_abline(aes(slope = .77, intercept = 0, colour="All"), size=1)+
  geom_abline(aes(slope = .43, intercept= 0, colour="Q5"), size=1)+
  geom_abline(aes(slope=.93, intercept = 0, colour= "Non-Q5"), size=1)+
  xlab("log Total Richness")+ylab("log PD")+
  geom_text(aes(x=5,y=-3.5),label="Slope All = 0.77", size=7)+
  geom_text(aes(x=5,y=-4.3),label="Slope Non-Q5 = 0.93", size=7)+
  geom_text(aes(x=5,y=-5),label="Slope Q5 = 0.43", size=7)+
  scale_color_discrete(name="Slopes")+
  theme(legend.position=c(.8,.37),
        legend.title=element_text(size = 30,face = "bold"),
        legend.text = element_text(size = 30),
        legend.key.size = unit(1, "cm"),
        legend.key = element_rect(fill= "white",colour="black"),
        axis.text.y=element_text(size=15),
        axis.text.x=element_blank(),
        axis.title.y=element_text(size=30),
        axis.title.x=element_blank())+  
  guides(colour = guide_legend(override.aes = list(size=2)))


png("PDvsTR_Corr.png", width =2000 , height =1200 ,res = 200)
PDcor
dev.off()

###
###

BayesSlope(x = index.grid$f.Rich, y=index.grid$PD, nSubj = length(index.grid[,1]), 
           ploting = T, SaveName = "fRichVsPD", Xlab = "Richness ~ PD")

Q5pos <- which(index.grid$PDclass=="Q5")

BayesSlope(x = index.grid$f.Rich[Q5pos], y=index.grid$PD[Q5pos], 
           nSubj = length(Q5pos), 
           ploting = T, SaveName = "Q5fRichVsPD", Xlab = "Q5 Richness ~ PD")

BayesSlope(x = index.grid$f.Rich[-Q5pos], y=index.grid$PD[-Q5pos], 
           nSubj = length(which(index.grid$PDclass=="Non-Q5")), 
           ploting = T, SaveName = "NonQ5fRichVsPD", Xlab = "Non-Q5 Richness ~ PD")

library(ggplot2)

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

PDcor <- ggplot(data = index.grid)+
  geom_point(aes(x=log(f.Rich),y=log(PD), colour=PDclass))+
  geom_abline(aes(slope = .75, intercept = 0, colour="All"), size=1)+
  geom_abline(aes(slope = .45, intercept= 0, colour="Q5"), size=1)+
  geom_abline(aes(slope=.96, intercept = 0, colour= "Non-Q5"), size=1)+
  xlab("log Total Richness")+ylab("log PD")+
  geom_text(aes(x=4.5,y=-3.5),label="Slope All = 0.75", size=7)+
  geom_text(aes(x=4.5,y=-4.3),label="Slope Non-Q5 = 0.96", size=7)+
  geom_text(aes(x=4.5,y=-5),label="Slope Q5 = 0.45", size=7)+
  scale_color_discrete(name="Slopes")+
  theme(legend.position=c(.8,.37),
        legend.title=element_text(size = 30,face = "bold"),
        legend.text = element_text(size = 30),
        legend.key.size = unit(1, "cm"),
        legend.key = element_rect(fill= "white",colour="black"),
        axis.text.y=element_text(size=15),
        axis.text.x=element_blank(),
        axis.title.y=element_text(size=30),
        axis.title.x=element_blank())+  
  guides(colour = guide_legend(override.aes = list(size=2)))

png("PDvsPR_Corr.png", width =2000 , height =1200 ,res = 200)
PDcor
dev.off()

#####################################################################################
#####################################################################################
#####################################################################################

BayesSlope(x = index.grid$t.Rich, y=index.grid$AvTD, nSubj = length(index.grid[,1]), 
           ploting = T, SaveName = "tRichVsAvTD", Xlab = "Richness ~ AvTD")

Q5pos <- which(index.grid$AvTDclass=="Q5")

BayesSlope(x = index.grid$t.Rich[Q5pos], y=index.grid$AvTD[Q5pos], 
           nSubj = length(Q5pos), 
           ploting = T, SaveName = "Q5tRichVsAvTD", Xlab = "Q5 Richness ~ AvTD")

BayesSlope(x = index.grid$t.Rich[-Q5pos], y=index.grid$AvTD[-Q5pos], 
           nSubj = length(which(index.grid$AvTDclass=="Non-Q5")), 
           ploting = T, SaveName = "NonQ5tRichvsAvTD", Xlab = "Non-Q5 Richness ~ AvTD")

library(ggplot2)

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

AvTDcor <- ggplot(data = index.grid)+
  geom_point(aes(x=log(t.Rich),y=log(AvTD), colour=AvTDclass))+
  geom_abline(aes(slope = .41, intercept = 0, colour="All"), size=1)+
  geom_abline(aes(slope = .34, intercept= 0, colour="Q5"), size=1)+
  geom_abline(aes(slope=.84, intercept = 0, colour= "Non-Q5"), size=1)+
  xlab("log Total Richness")+ylab("log AvTD")+
  geom_text(aes(x=5,y=-3.2),label="Slope All = 0.41", size=7)+
  geom_text(aes(x=5,y=-4.2),label="Slope Non-Q5 = 0.84", size=7)+
  geom_text(aes(x=5,y=-5),label="Slope Q5 = 0.34", size=7)+
  scale_color_discrete(name="Slopes")+
  theme(legend.position=c(.8,.36),
        legend.title=element_text(size = 30,face = "bold"),
        legend.text = element_text(size = 30),
        legend.key.size = unit(1, "cm"),
        legend.key = element_rect(fill= "white",colour="black"),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15),
        axis.title.y=element_text(size=30),
        axis.title.x=element_text(size=30))+  
  guides(colour = guide_legend(override.aes = list(size=2)))


png("AvTDvsTR_Corr.png", width =2000 , height =1200 ,res = 200)
AvTDcor
dev.off()

###
###

BayesSlope(x = index.grid$f.Rich, y=index.grid$AvTD, nSubj = length(index.grid[,1]), 
           ploting = T, SaveName = "fRichVsAvTD", Xlab = "Richness ~ AvTD")

Q5pos <- which(index.grid$AvTDclass=="Q5")

BayesSlope(x = index.grid$f.Rich[Q5pos], y=index.grid$AvTD[Q5pos], 
           nSubj = length(Q5pos), 
           ploting = T, SaveName = "Q5fRichVsAvTD", Xlab = "Q5 Richness ~ AvTD")

BayesSlope(x = index.grid$f.Rich[-Q5pos], y=index.grid$AvTD[-Q5pos], 
           nSubj = length(which(index.grid$AvTDclass=="Non-Q5")), 
           ploting = T, SaveName = "NonQ5fRichVsAvTD", Xlab = "Non-Q5 Richness ~ AvTD")

library(ggplot2)

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

AvTDcor <- ggplot(data = index.grid)+
  geom_point(aes(x=log(f.Rich),y=log(AvTD), colour=AvTDclass))+
  geom_abline(aes(slope = .38, intercept = 0, colour="All"), size=1)+
  geom_abline(aes(slope = .27, intercept= 0, colour="Q5"), size=1)+
  geom_abline(aes(slope=.88, intercept = 0, colour= "Non-Q5"), size=1)+
  xlab("log Phylogenetic Richness")+ylab("log AvTD")+
  geom_text(aes(x=5,y=-3.2),label="Slope All = 0.38", size=7)+
  geom_text(aes(x=5,y=-4.2),label="Slope Non-Q5 = 0.88", size=7)+
  geom_text(aes(x=5,y=-5),label="Slope Q5 = 0.27", size=7)+
  scale_color_discrete(name="Slopes")+
  theme(legend.position=c(.8,.36),
        legend.title=element_text(size = 30,face = "bold"),
        legend.text = element_text(size = 30),
        legend.key.size = unit(1, "cm"),
        legend.key = element_rect(fill= "white",colour="black"),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15),
        axis.title.y=element_text(size=30),
        axis.title.x=element_text(size=30))+  
  guides(colour = guide_legend(override.aes = list(size=2)))

png("AvTDvsPR_Corr.png", width =2000 , height =1200 ,res = 200)
AvTDcor
dev.off()

######################################################################################
######################################################################################

BayesSlope(x = index.grid$TD, y=index.grid$PD, nSubj = length(index.grid[,1]), 
           ploting = T, SaveName = "TDvsPD", Xlab = "TD ~ PD")

BayesSlope(x = index.grid$TD, y=index.grid$AvTD, nSubj = length(index.grid[,1]), 
           ploting = T, SaveName = "TDvsAvTD", Xlab = "TD ~ AvTD")

BayesSlope(x = index.grid$PD, y=index.grid$AvTD, nSubj = length(index.grid[,1]), 
           ploting = T, SaveName = "PDvsAvTD", Xlab = "PD ~ AvTD")

######################################################################################
######################################################################################
