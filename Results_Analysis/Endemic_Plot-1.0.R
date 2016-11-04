library(ggplot2)
library(classInt)

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final/Endemic/")

end <- read.csv("Area.end",header = T)

head(end,5L)

col <- c(Cauca="#990000",Choco.Darien="#CCCC33",Guajira="#00CCCC",Imeri="#CCCCFF",Magdalena="#000066",Napo="#CC3366",Paramo="#336600",
         Sabana="#FF9999", Venezuela="#CCFF66",WEcuador="#FF9933")

ord <- unique(end$Area[order(end$TD[which(end$Percentage==0)], decreasing = T)])

png(filename = "DT_area_end.png",width = 1500,height = 1000)

p <- ggplot(data = end, aes(x=factor(Percentage),y=TD))
p + geom_violin() +
  geom_jitter(aes(colour=Area),width=0.7, size=2)+
  #geom_point(aes(x=factor(Percentage),y=mean(Values)),colour="red",shape=17,size=4)+
  ylab("TD Index")+ xlab("Remove %")+
  ggtitle("Endemic species for areas of endemism (TD)")+
  scale_color_manual(breaks=ord,values = col)+
  theme(panel.background = element_rect(fill = "gray97"))+
  theme(axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 30,),
        axis.title.x = element_text(size = 30,),
        plot.title=element_text(size=50),
        legend.title = element_text(size=35),
        legend.text = element_text(size = 35),
        legend.key.size = unit(2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=10)))


dev.off()

#################################################################################

#end <- read.csv("PD_area.endemic",header = T)

head(end,5L)

ord <- unique(end$Area[order(end$PD[which(end$Percentage==0)], decreasing = T)])

png(filename = "PD_area_end.png",width = 1500,height = 1000,)

p <- ggplot(data = end, aes(x=factor(Percentage),y=PD))
p + geom_violin() +
  geom_jitter(aes(colour=Area),width=0.7, size=2)+
  #geom_point(aes(x=factor(Percentage),y=mean(PD)),colour="red",shape=17,size=4)+
  ylab("PD")+ xlab("Remove %")+
  ggtitle("Endemic species for areas of endemism (PD)")+
  scale_color_manual(breaks=ord,values = col)+
  theme(panel.background = element_rect(fill = "gray97"))+
  theme(axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 30,),
        axis.title.x = element_text(size = 30,),
        plot.title=element_text(size=50),
        legend.title = element_text(size=35),
        legend.text = element_text(size = 35),
        legend.key.size = unit(2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=10)))


dev.off()

##################################################################################

#end <- read.csv("AvDT_area.endemic",header = T)

head(end,5L)

ord <- unique(end$Area[order(end$AvTD[which(end$Percentage==0)], decreasing = T)])

png(filename = "AvDT_area_end.png",width = 1500,height = 1000)

p <- ggplot(data = end, aes(x=factor(Percentage),y=AvTD))
p + geom_violin() +
  geom_jitter(aes(colour=Area), width = .7, size=2)+
  #geom_point(aes(x=factor(Percentage),y=mean(Values)),colour="red",shape=17,size=4)+
  ylab("AvTD")+ xlab("Remove %")+
  ggtitle("Endemic species for areas of endemism (AvTD)")+
  scale_color_manual(breaks=ord,values = col)+
  theme(panel.background = element_rect(fill = "gray97"))+
  theme(axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 30,),
        axis.title.x = element_text(size = 30,),
        plot.title=element_text(size=50),
        legend.title = element_text(size=35),
        legend.text = element_text(size = 35),
        legend.key.size = unit(2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=10)))


dev.off()

ord2 <- ord[c(grep("Paramo", ord),grep("Cauca", ord))]
col2 <- col[c(grep("Paramo", names(col)),grep("Cauca", names(col)))]

endAVTD <- end[c(grep("Cauca", end$Area),grep("Paramo", end$Area)), ]
endAVTD

png(filename = "AvDT2_area_end.png",width = 1500,height = 1000)
p <- ggplot(data = endAVTD, aes(x=factor(Percentage),y=AvTD))
p + geom_violin() +
  geom_jitter(aes(colour=Area),width = .7, size=2)+
  #geom_point(aes(x=factor(Percentage),y=mean(Values)),colour="red",shape=17,size=4)+
  ylab("AvTD")+ xlab("Remove %")+
  ggtitle("Endemic species for Cauca and Paramo (AvTD)")+
  scale_color_manual(breaks=ord2,values = col2)+
  theme(panel.background = element_rect(fill = "gray97"))+
  theme(axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 30,),
        axis.title.x = element_text(size = 30,),
        plot.title=element_text(size=50),
        legend.title = element_text(size=35),
        legend.text = element_text(size = 35),
        legend.key.size = unit(2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=10)))

dev.off()

#################################################################################
library(classInt)
library(RColorBrewer)
library(maptools)

end <- read.csv("Grid.end",header = T)

class <- readShapePoly("~/Documentos/Omar/Tesis/Taxa/Results/Final/Raw_IndexR/Grid25_TD.shp")

q5 <- paste("X",class$SP_ID[grep("Q5",class$Index)],sep="")
q5 <- which(end$Area%in%q5)

q4 <- paste("X",class$SP_ID[grep("Q4",class$Index)],sep="")
q4 <- which(end$Area%in%q4)

#brewer.pal.info

col <- rep(brewer.pal(9,"Reds"), length(unique(end$Area[q5])))
col <- col[length(col):1]

png(filename = "DT_grid25Q5_end.png",width = 1500,height = 1000)

p <- ggplot()
p +geom_jitter(aes(x=factor(end$Percentage[q5]),
                   y=end$TD[q5],
                   colour=end$Area[q5]), size=2)+
  theme(legend.position = "none")+ 
  geom_violin(aes(x=factor(end$Percentage[q5]),y=end$TD[q5]))+
  ylab("TD")+ xlab("Remove %")+
  ggtitle("Endemic species for 25° grid  (TD-Q5)")+
  theme(panel.background = element_rect(fill = "gray97"))+
  scale_color_manual(values = col)+
  theme(axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 30,),
        axis.title.x = element_text(size = 30,),
        plot.title=element_text(size=50))

dev.off()


col <- rep(brewer.pal(9,"Blues"), length(unique(end$Area[q4])))
col <- col[length(col):1]

png(filename = "DT_grid25Q4_end.png",width = 1500,height = 1000)

p <- ggplot()
p +geom_jitter(aes(x=factor(end$Percentage[q4]),
                   y=end$TD[q4],
                   colour=end$Area[q4]), size=2)+
  theme(legend.position = "none")+ 
  geom_violin(aes(x=factor(end$Percentage[q4]),y=end$TD[q4]))+
  ylab("TD")+ xlab("Remove %")+
  ggtitle("Endemic species for 25° grid  (TD-Q4)")+
  theme(panel.background = element_rect(fill = "gray97"))+
  scale_color_manual(values = col)+
  theme(axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 30,),
        axis.title.x = element_text(size = 30,),
        plot.title=element_text(size=50))

dev.off()


novoQ <- rep(NA, length(end$Area))

novoQ[q5] <- rep("Q5", length(novoQ[q5]))
novoQ[q4] <- rep("Q4", length(novoQ[q4]))

end$Area <- as.factor(novoQ)

png(filename = "DT_grid25Q4Q5_end.png",width = 1500,height = 1000)

p <- ggplot()
p +geom_jitter(aes(x=factor(end$Percentage[c(q5,q4)]),
                   y=end$TD[c(q5,q4)],
                   colour=end$Area[c(q5,q4)]), size=2)+
  theme(legend.position = "bottom")+ 
  geom_violin(aes(x=factor(end$Percentage[c(q5,q4)]),y=end$TD[c(q5,q4)]))+
  ylab("TD")+ xlab("Remove %")+
  scale_color_manual(name="Cells",breaks=c("Q5","Q4"), values = c(Q5="red",Q4="blue"))+
  ggtitle("Endemic species for 25° grid  (TD-Q5+Q4)")+
  theme(panel.background = element_rect(fill = "gray97"))+
  theme(axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 30,),
        axis.title.x = element_text(size = 30,),
        plot.title=element_text(size=50),
        legend.title = element_text(size=35),
        legend.text = element_text(size = 35),
        legend.key.size = unit(2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=10)))

dev.off()
##################################################################################
##################################################################################

end <- read.csv("Grid.end",header = T)

class <- readShapePoly("~/Documentos/Omar/Tesis/Taxa/Results/Final/Raw_IndexR/Grid25_PD.shp")

q5 <- paste("X",class$SP_ID[grep("Q5",class$Index)],sep="")
q5 <- which(end$Area%in%q5)

q4 <- paste("X",class$SP_ID[grep("Q4",class$Index)],sep="")
q4 <- which(end$Area%in%q4)

#brewer.pal.info

col <- rep(brewer.pal(9,"Reds"), length(unique(end$Area[q5])))
#col <- col[length(col):1]

png(filename = "PD_grid25Q5_end.png",width = 1500,height = 1000)

p <- ggplot()
p +geom_jitter(aes(x=factor(end$Percentage[q5]),
                   y=end$PD[q5],
                   colour=end$Area[q5]), size=2)+
  theme(legend.position = "none")+ 
  geom_violin(aes(x=factor(end$Percentage[q5]),y=end$PD[q5]))+
  ylab("PD")+ xlab("Remove %")+
  ggtitle("Endemic species for 25° grid  (PD-Q5)")+
  theme(panel.background = element_rect(fill = "gray97"))+
  scale_color_manual(values = col)+
  theme(axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 30,),
        axis.title.x = element_text(size = 30,),
        plot.title=element_text(size=50))

dev.off()


col <- rep(brewer.pal(9,"Blues"), length(unique(end$Area[q4])))
col <- col[length(col):1]

png(filename = "PD_grid25Q4_end.png",width = 1500,height = 1000)

p <- ggplot()
p +geom_jitter(aes(x=factor(end$Percentage[q4]),
                   y=end$PD[q4],
                   colour=end$Area[q4]), size=2)+
  theme(legend.position = "none")+ 
  geom_violin(aes(x=factor(end$Percentage[q4]),y=end$PD[q4]))+
  ylab("PD")+ xlab("Remove %")+
  ggtitle("Endemic species for 25° grid  (PD-Q4)")+
  theme(panel.background = element_rect(fill = "gray97"))+
  scale_color_manual(values = col)+
  theme(axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 30,),
        axis.title.x = element_text(size = 30,),
        plot.title=element_text(size=50))

dev.off()


novoQ <- rep(NA, length(end$Area))

novoQ[q5] <- rep("Q5", length(novoQ[q5]))
novoQ[q4] <- rep("Q4", length(novoQ[q4]))

end$Area <- as.factor(novoQ)

png(filename = "PD_grid25Q4Q5_end.png",width = 1500,height = 1000)

p <- ggplot()
p +geom_jitter(aes(x=factor(end$Percentage[c(q5,q4)]),
                   y=end$PD[c(q5,q4)],
                   colour=end$Area[c(q5,q4)]), size=2)+
  theme(legend.position = "bottom")+ 
  geom_violin(aes(x=factor(end$Percentage[c(q5,q4)]),y=end$PD[c(q5,q4)]))+
  ylab("PD")+ xlab("Remove %")+
  scale_color_manual(name="Cells",breaks=c("Q5","Q4"), values = c(Q5="red",Q4="blue"))+
  ggtitle("Endemic species for 25° grid  (PD-Q5+Q4)")+
  theme(panel.background = element_rect(fill = "gray97"))+
  theme(axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 30,),
        axis.title.x = element_text(size = 30,),
        plot.title=element_text(size=50),
        legend.title = element_text(size=35),
        legend.text = element_text(size = 35),
        legend.key.size = unit(2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=10)))

dev.off()
##################################################################################
##################################################################################
##################################################################################

end <- read.csv("Grid.end",header = T)

class <- readShapePoly("~/Documentos/Omar/Tesis/Taxa/Results/Final/Raw_IndexR/Grid25_AvTD.shp")

q5 <- paste("X",class$SP_ID[grep("Q5",class$Index)],sep="")
q5 <- which(end$Area%in%q5)

q4 <- paste("X",class$SP_ID[grep("Q4",class$Index)],sep="")
q4 <- which(end$Area%in%q4)

#brewer.pal.info

col <- rep(brewer.pal(9,"Reds"), length(unique(end$Area[q5])))
col <- col[length(col):1]

end$AvTD <- log2(end$AvTD)

png(filename = "AvDT_grid25Q5_end.png",width = 1500,height = 1000)

p <- ggplot()
p +geom_jitter(aes(x=factor(end$Percentage[q5]),
                   y=end$AvTD[q5],
                   colour=end$Area[q5]), size=2)+
  theme(legend.position = "none")+ 
  geom_violin(aes(x=factor(end$Percentage[q5]),y=end$AvTD[q5]))+
  ylab("AvTD")+ xlab("Remove %")+
  ggtitle("Endemic species for 25° grid  (AvTD-Q5)")+
  theme(panel.background = element_rect(fill = "gray97"))+
  scale_color_manual(values = col)+
  theme(axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 30,),
        axis.title.x = element_text(size = 30,),
        plot.title=element_text(size=50))

dev.off()


col <- rep(brewer.pal(9,"Blues"), length(unique(end$Area[q4])))
col <- col[length(col):1]

png(filename = "AvDT_grid25Q4_end.png",width = 1500,height = 1000)

p <- ggplot()
p +geom_jitter(aes(x=factor(end$Percentage[q4]),
                   y=end$AvTD[q4],
                   colour=end$Area[q4]), size=2)+
  theme(legend.position = "none")+ 
  geom_violin(aes(x=factor(end$Percentage[q4]),y=end$AvTD[q4]))+
  ylab("AvTD")+ xlab("Remove %")+
  ggtitle("Endemic species for 25° grid  (AvTD-Q4)")+
  theme(panel.background = element_rect(fill = "gray97"))+
  scale_color_manual(values = col)+
  theme(axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 30,),
        axis.title.x = element_text(size = 30,),
        plot.title=element_text(size=50))

dev.off()


novoQ <- rep(NA, length(end$Area))

novoQ[q5] <- rep("Q5", length(novoQ[q5]))
novoQ[q4] <- rep("Q4", length(novoQ[q4]))

end$Area <- as.factor(novoQ)

png(filename = "AvDT_grid25Q4Q5_end.png",width = 1500,height = 1000)

p <- ggplot()
p +geom_jitter(aes(x=factor(end$Percentage[c(q5,q4)]),
                   y=end$AvTD[c(q5,q4)],
                   colour=end$Area[c(q5,q4)]), size=2)+
  theme(legend.position = "bottom")+ 
  geom_violin(aes(x=factor(end$Percentage[c(q5,q4)]),y=end$AvTD[c(q5,q4)]))+
  ylab("AvTD")+ xlab("Remove %")+
  scale_color_manual(name="Cells",breaks=c("Q5","Q4"), values = c(Q5="red",Q4="blue"))+
  ggtitle("Endemic species for 25° grid  (AvTD-Q5+Q4)")+
  theme(panel.background = element_rect(fill = "gray97"))+
  theme(axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 30,),
        axis.title.x = element_text(size = 30,),
        plot.title=element_text(size=50),
        legend.title = element_text(size=35),
        legend.text = element_text(size = 35),
        legend.key.size = unit(2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=10)))

dev.off()
##################################################################################
##################################################################################
