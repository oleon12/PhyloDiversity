library(ggplot2)

setwd("~/Documentos/Omar/Tesis/Taxa/Richness/")

GI <- read.csv("General2.info")

head(GI, 5L)

pos <- c("Mammalia","Amphibia","Aves","Insecta","Squamata","Testudines",
         "Magnoliopsida","Liliopsida","Gimnospermae","Bryophyta","Marchantiophyta")

p <- ggplot(data = GI)
p <- p + geom_bar(aes(MajorGroup, fill=MajorPhyllum))+xlab(NULL)+ylab("Number of Species")+
  scale_x_discrete(limits=pos)+scale_y_continuous(expand = c(.01,.01))+
  theme(legend.position=c(.9,.9),
        legend.title=element_blank(),
        legend.text = element_text(size = 60),
        legend.key.size = unit(3.5, "cm"),
        legend.background = element_rect(fill="transparent"),
        legend.key = element_rect(fill="transparent"))+
  theme(axis.text = element_text(size = 37, face="bold"),
        axis.title.y = element_text(size = 55))
p

setwd("/media/omar/MAR/")

png("SpeciesGroup.png",width =3000 ,height =1500)
p
dev.off()
