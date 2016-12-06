library(ggplot2)
library(cowplot)
library(gridExtra)

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/")

general.info <- read.csv("General.info",header = T)

match.sp <- as.factor(read.csv("Match.sp")$x)

general.info <- general.info[which(general.info$Sp%in%match.sp),]

general.info$Sp[which(general.info$Ende.WD=="End")]

dt1 <- ggplot(data = general.info) + 
  geom_violin(aes(x=factor("1"),y=DT.area),fill="gray97")+
  geom_jitter(aes(x=factor("1"),y=DT.area,colour=Ende.WD),size=2)+
  #geom_point(aes(x=1,y=mean(DT.area)),size=5)+
  ylab("TD")+ xlab("")+
  ggtitle("Areas of endemism")+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=32),
        axis.title.y = element_text(size = 30,face="bold"))


dt2 <- ggplot(data = general.info) + 
  geom_violin(aes(x=factor("1"),y=DT.grid),fill="gray97")+
  geom_jitter(aes(x=factor("1"),y=DT.grid,colour=Ende.WD),size=2)+
  #geom_point(aes(x=1,y=mean(DT.area)),size=5)+
  ylab("")+ xlab("")+
  ggtitle("Grid of 0.25° cells")+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=32),
        axis.title.y = element_text(size = 30,face="bold"))


pd1 <- ggplot(data = general.info)+ 
  geom_violin(aes(x=factor("1"),y=log(PD.area)),fill="gray97")+
  geom_jitter(aes(x=factor("1"),y=log(PD.area),colour=Ende.WD),size=2)+
  #geom_point(aes(x=1,y=mean(DT.area)),size=5)+
  ylab("log PD")+ xlab("")+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 30,face="bold"))

pd2 <- ggplot(data = general.info) + 
  geom_violin(aes(x=factor("1"),y=log(PD.grid)),fill="gray97")+
  geom_jitter(aes(x=factor("1"),y=log(PD.grid),colour=Ende.WD),size=2)+
  #geom_point(aes(x=1,y=mean(DT.area)),size=5)+
  ylab("")+ xlab("")+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 30,face="bold"))

avdt1 <- ggplot(data = general.info)+
  geom_violin(aes(x=factor("1"),y=log(AvDT.area)),fill="gray97")+
  geom_jitter(aes(x=factor("1"),y=log(AvDT.area),colour=Ende.WD),size=2)+
  #geom_point(aes(x=1,y=mean(DT.area)),size=5)+
  scale_color_discrete(name="Species", labels=c("Endemic","Widespread"))+
  ylab("log AvTD")+ xlab("")+
  theme(legend.title=element_text(face="bold", size = 25),
        legend.position="bottom",
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 30,face="bold"),
        legend.text = element_text(size = 20),
        legend.key.size = unit(2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=5)))

avdt2 <- ggplot(data = general.info) + 
  geom_violin(aes(x=factor("1"),y=log(AvDT.grid)),fill="gray97")+
  geom_jitter(aes(x=factor("1"),y=log(AvDT.grid),colour=Ende.WD),size=2)+
  #geom_point(aes(x=1,y=mean(DT.area)),size=5)+
  scale_color_discrete(name="Species", labels=c("Endemic","Widespread"))+
  ylab("")+ xlab("")+
  theme(legend.title=element_text(face="bold", size = 25),
        legend.position="bottom",
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 30,face="bold"),
        legend.text = element_text(size = 20),
        legend.key.size = unit(2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=5)))


grid.arrange(dt1, dt2, pd1, pd2,avdt1,avdt2, ncol=2, nrow =3)

png(filename = "~/Documentos/Omar/Tesis/Taxa/Results/Final2/Endemic/grid.png",
    width = 1001,
    height=1000)
grid.arrange(dt1, dt2, pd1, pd2,avdt1,avdt2, ncol=2, nrow =3)
dev.off()


End1 <- general.info[grep("End", general.info$Ende.WD), ]
End1 <- End1[order(End1$BL, decreasing = T), ]
WD1 <- general.info[grep("WD", general.info$Ende.WD), ]
WD1 <- WD1[order(WD1$BL, decreasing = T), ]

GI.end <- as.data.frame(rbind(End1, WD1))
299
head(GI.end)

bl.plot <- ggplot()+geom_point(aes(x=(1:length(WD1$BL))+1,y=log(WD1$BL),colour=WD1$Ende.WD),size=6)+
  geom_point(aes(x=1:length(End1$BL),y=log(End1$BL),colour=End1$Ende.WD),size=6)+
  ylab("log Branch Lengths")+xlab("Species")+ 
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=30),
        axis.ticks=element_blank(),
        axis.title.y=element_text(size=40),
        axis.title.x=element_text(size=40),
        panel.background = element_rect(fill = "gray97"),
        panel.grid.major = element_line(colour = "white"))+
  scale_colour_discrete(name="Species",
                        labels=c("Endemic","Widespread"))+
  theme(legend.position=c(.8,.9),
        legend.direction="horizontal",
        legend.title=element_text(size = 35,face = "bold"),
        legend.text = element_text(size = 35),
        legend.key.size = unit(1, "cm"),
        legend.key = element_rect(fill= "white",colour="black"),
        legend.background = element_rect(fill = "white"))+  
  guides(colour = guide_legend(override.aes = list(size=10)))+
  geom_text(aes(x=740,y=1.25),label="Endemic = 299 spp",size=15)+
  geom_text(aes(x=720,y=0),label="Total = 1255 spp",size=15)

png(filename = "~/Documentos/Omar/Tesis/Taxa/Results/Final/Endemic/BL.png",
    width = 1500,
    height=1000)
bl.plot
dev.off()
