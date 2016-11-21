library(ggplot2)
library(maptools)
library(gridExtra)

setwd("~/Documentos/Omar/Tesis/Scripts/Distribution/shp/NAB/")

PAshp <- as.data.frame(readShapePoly("PNN2_NAB.shp"))
PAshp2 <- PAshp[order(PAshp$NAME), ]

AEshp <- as.data.frame(readShapePoly("A_endemism.NAB.shp"))
AEshp2 <- AEshp[order(AEshp$Name), ]

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

PA <- data.frame(ID=read.csv("TD.PNN")$area,
                 TD=read.csv("TD.PNN")$W,
                 PD=read.csv("PD.PNN")$PD,
                 AvTD=read.csv("AvTD.PNN")$Dplus)

AE <-data.frame(ID=read.csv("DT.Area")$area,
                TD=read.csv("DT.Area")$W,
                PD=read.csv("PD.Area")$PD,
                AvTD=read.csv("AvTD.Area")$Dplus)


PA2 <- PA[order(PA$ID), ]
AE2 <- AE[order(AE$ID), ]

PNNin <- grep("PNNin",PA2$ID)
PNNout <- grep("PNNout",PA2$ID)
NABin <- grep("NABin",PA2$ID)

PAend <- cbind(PA2[-c(PNNin,PNNout,NABin),], Area=PAshp2$Area)
AEend <- cbind(AE2[-c(PNNin,PNNout,NABin),], Area=AEshp2$Area)

cor1 <- cor.test(x=AEend$Area,y=log(AEend$TD), method = "spearman")

p1 <- ggplot(data=AEend)+ geom_point(aes(x=Area,y=log(TD)))+xlab("Size of the Areas of Endemism")+
  ylab("log TD Index")+
  geom_smooth(aes(x=Area,y=log(TD)), method = "lm")+
  geom_label(aes(x=80000000000,y=5.1), label=(paste("R",cor1$estimate,sep = " : ")))+
  theme(plot.margin = unit(c(0,5,1,1),units="points"))
p1 #bottom


cor2 <- cor.test(x=AEend$Area,y=log(AEend$PD), method = "spearman")

p2 <- ggplot(data = AEend)+geom_point(aes(x=Area,y=log(PD)))+xlab(NULL)+ylab("log PD index")+
  geom_smooth(aes(x=Area, y=log(PD)), method = "lm")+
  geom_label(aes(x=80000000000,y=3.1), label=(paste("R", round(cor2$estimate,digits = 1),sep = " : ")))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0,5,1,1),units="points"))
p2 #middle


cor3 <- cor.test(x=AEend$Area,y=log(AEend$AvTD), method = "spearman")

p3 <- ggplot(data = AEend)+geom_point(aes(x=Area,y=log(AvTD)))+xlab(NULL)+ylab("AvTD index")+
  geom_smooth(aes(x=Area, y=log(AvTD)), method = "lm")+
  geom_label(aes(x=80000000000,y=1.5), label=(paste("R", round(cor3$estimate,digits = 1),sep = " : ")))+  
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0,5,1,1),units="points"))
p3 #middle

png("SizeAE_Ind.png")
grid.arrange(p3,p2,p1,heights = c(1/3, 1/3, 1/3))
dev.off()

####################################################################################################

cor1 <- cor.test(x=PAend$Area,y=log(PAend$TD), method = "spearman")

p1 <- ggplot()+ geom_point(aes(x=PAend$Area[-which(PAend$TD==0)],y=log(PAend$TD[-which(PAend$TD==0)])))+
  xlab("Size of the Protected Areas")+
  ylab("log TD Index")+
  geom_smooth(aes(x=PAend$Area[-which(PAend$TD==0)],y=log(PAend$TD[-which(PAend$TD==0)])), method = "lm")+
  geom_label(aes(x=5000000000,y=0.5), label=(paste("R",round(cor1$estimate,digits=1),sep = " : ")))+
  theme(plot.margin = unit(c(0,5,1,1),units="points"))
p1 #bottom


cor2 <- cor.test(x=PAend$Area,y=log(PAend$PD), method = "spearman")

p2 <- ggplot()+geom_point(aes(x=PAend$Area[-which(PAend$PD==0)],y=log(PAend$PD[-which(PAend$PD==0)])))+
  xlab(NULL)+ylab("log PD index")+
  geom_smooth(aes(x=PAend$Area[-which(PAend$PD==0)], y=log(PAend$PD[-which(PAend$PD==0)])), method = "lm")+
  geom_label(aes(x=5000000000,y=-3.5), label=(paste("R", round(cor2$estimate,digits = 1),sep = " : ")))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0,5,1,1),units="points"))
p2 #middle


cor3 <- cor.test(x=PAend$Area,y=log(PAend$AvTD), method = "spearman")

p3 <- ggplot()+geom_point(aes(x=PAend$Area[-which(PAend$AvTD==0)],y=log(PAend$AvTD[-which(PAend$AvTD==0)])))+
  xlab(NULL)+ylab("AvTD index")+
  geom_smooth(aes(x=PAend$Area[-which(PAend$AvTD==0)], y=log(PAend$AvTD[-which(PAend$AvTD==0)])), method = "lm")+
  geom_label(aes(x=5000000000,y=-5.1), label=(paste("R", round(cor3$estimate,digits = 1),sep = " : ")))+  
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0,5,1,1),units="points"))
p3 #middle

png("SizePA_Ind.png")
grid.arrange(p3,p2,p1,heights = c(1/3, 1/3, 1/3))
dev.off()
