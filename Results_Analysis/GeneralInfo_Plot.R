library(ggplot2)
library(cowplot)
library(gridExtra)

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Abril1/")

general.info <- read.csv("General.info",header = T)

general.info <- general.info[which(general.info$Sp%in%match.sp),]

general.info$Sp[general.info$Ende.WD=="End"]%in%end.sp

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
        plot.title = element_text(size=22))
  

dt2 <- ggplot(data = general.info) + 
  geom_violin(aes(x=factor("1"),y=DT.grid),fill="gray97")+
  geom_jitter(aes(x=factor("1"),y=DT.grid,colour=Ende.WD),size=2)+
  #geom_point(aes(x=1,y=mean(DT.area)),size=5)+
  ylab("")+ xlab("")+
  ggtitle("Grid of 25° cells")+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size=22))
 

pd1 <- ggplot(data = general.info)+ 
  geom_violin(aes(x=factor("1"),y=PD.area),fill="gray97")+
    geom_jitter(aes(x=factor("1"),y=PD.area,colour=Ende.WD),size=2)+
    #geom_point(aes(x=1,y=mean(DT.area)),size=5)+
    ylab("PD")+ xlab("")+
    theme(legend.title=element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks = element_blank())

pd2 <- ggplot(data = general.info) + 
  geom_violin(aes(x=factor("1"),y=PD.grid),fill="gray97")+
    geom_jitter(aes(x=factor("1"),y=PD.grid,colour=Ende.WD),size=2)+
    #geom_point(aes(x=1,y=mean(DT.area)),size=5)+
    ylab("")+ xlab("")+
    theme(legend.title=element_blank(),
        legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks = element_blank())

avdt1 <- ggplot(data = general.info)+
  geom_violin(aes(x=factor("1"),y=AvDT.area),fill="gray97")+
  geom_jitter(aes(x=factor("1"),y=AvDT.area,colour=Ende.WD),size=2)+
  #geom_point(aes(x=1,y=mean(DT.area)),size=5)+
  ylab("AvTD")+ xlab("")+
  theme(legend.title=element_blank(),
        legend.position="bottom",
        axis.text.x = element_blank(),
        axis.ticks = element_blank())

avdt2 <- ggplot(data = general.info) + 
  geom_violin(aes(x=factor("1"),y=AvDT.grid),fill="gray97")+
  geom_jitter(aes(x=factor("1"),y=AvDT.grid,colour=Ende.WD),size=2)+
  #geom_point(aes(x=1,y=mean(DT.area)),size=5)+
  ylab("")+ xlab("")+
  theme(legend.title=element_blank(),
        legend.position="bottom",
        axis.text.x = element_blank(),
        axis.ticks = element_blank())


grid.arrange(dt1, dt2, pd1, pd2,avdt1,avdt2, ncol=2, nrow =3)

png(filename = "~/Documentos/Omar/Tesis/Taxa/Results/Abril1/Endemic/grid.png",
    width = 2001,
    height=2000,res = 200)
grid.arrange(dt1, dt2, pd1, pd2,avdt1,avdt2, ncol=2, nrow =3)
dev.off()

general.info$

bl.plot <- ggplot()+geom_point(aes(x=1:(length(general.info$BL)-1),
                     y=general.info$BL[which(general.info$BL<max(general.info$BL))],
                     colour=factor(general.info$Ende.WD[-1])),size=5)+
  ylab("Branch Length")+ xlab("Species")+
  theme(legend.position=c("bottom"),
        legend.title=element_text(size = 15,face = "bold"),
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.5, "cm"),
        legend.key = element_rect(fill= "white",colour="black"),
        legend.background = element_rect(fill = "white"))+
  theme(axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.title.y=element_text(size=20),
        axis.title.x=element_text(size=20),
        panel.background = element_rect(fill = "gray97"),
        panel.grid.major = element_line(colour = "white"))+
  scale_colour_discrete(name="Specie",
                       labels=c("Endemic","Widespread"))+
  geom_text(aes(x=90,y=0.4),label="Endemic = 12 spp",size=8)+
  geom_text(aes(x=90,y=0.35),label="Total = 105 spp",size=8)



png(filename = "~/Documentos/Omar/Tesis/Taxa/Results/Abril1/Endemic/BL.png",
    width = 5001,
    height=2000,res = 300)
bl.plot
dev.off()
