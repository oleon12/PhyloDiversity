library(ggplot2)
library(classInt)

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Abril1/Endemic/")

end <- read.csv("DT_area.endemic",header = T)

head(end,5L)

png(filename = "DT_area_end.png",width = 1070,height = 699,res = 120)

p <- ggplot(data = end, aes(x=factor(Percentage),y=Values))
p + #geom_violin() +
  geom_jitter(aes(colour=Area),width=0.7)+
  #geom_point(aes(x=factor(Percentage),y=mean(Values)),colour="red",shape=17,size=4)+
  ylab("Taxonomic Distincness")+ xlab("Remove %")+
  ggtitle("Endemic species for areas of endemism (TD)")+
  scale_color_brewer(palette = "Paired")+
  theme(panel.background = element_rect(fill = "gray97"))

dev.off()

#################################################################################

end <- read.csv("PD_area.endemic",header = T)

head(end,5L)

png(filename = "PD_area_end.png",width = 1070,height = 699,res = 120)

p <- ggplot(data = end, aes(x=factor(Percentage),y=Values))
p + geom_violin() +
  geom_jitter(aes(colour=Area))+
  geom_point(aes(x=factor(Percentage),y=mean(Values)),colour="red",shape=17,size=4)+
  ylab("Phylogenetic Diversity")+ xlab("Remove %")+
  ggtitle("Endemic species for areas of endemism (PD)")+
  scale_color_brewer(palette = "Paired")+
  theme(panel.background = element_rect(fill = "gray97"))

dev.off()

##################################################################################

end <- read.csv("AvDT_area.endemic",header = T)

head(end,5L)

png(filename = "AvDT_area_end.png",width = 1070,height = 699,res = 120)

p <- ggplot(data = end, aes(x=factor(Percentage),y=Values))
p + geom_violin() +
  geom_jitter(aes(colour=Area))+
  geom_point(aes(x=factor(Percentage),y=mean(Values)),colour="red",shape=17,size=4)+
  ylab("Average Taxonomic Distincness")+ xlab("Remove %")+
  ggtitle("Endemic species for areas of endemism (AvTD)")+
  scale_color_brewer(palette = "Paired")+
  theme(panel.background = element_rect(fill = "gray97"))

dev.off()

#################################################################################
library(classInt)

end <- read.csv("DT_grid.endemic",header = T)


brks <- classIntervals(end$Values,n=5,style="quantile")
brks <- brks$brks
brks

?classIntervals

q5 <- which(findInterval(end$Values,brks,all.inside = T)==5)
q4 <- which(findInterval(end$Values,brks,all.inside = T)==4)

png(filename = "DT_grid25Q5_end.png",width = 1070,height = 699,res = 120)

p <- ggplot()
p +geom_jitter(aes(x=factor(end$Percentage[q5]),
                  y=end$Values[q5],
                  colour=end$Area[q5]))+
  theme(legend.position = "none")+ 
  geom_violin(aes(x=factor(end$Percentage[q5]),y=end$Values[q5]))+
  ylab("Taxonomic Distincness")+ xlab("Remove %")+
  ggtitle("Endemic species for 25° grid  (TD-Q5)")+
  theme(panel.background = element_rect(fill = "gray97"))

dev.off()

png(filename = "DT_grid25Q4_end.png",width = 1070,height = 699,res = 120)

p <- ggplot()
p +geom_jitter(aes(x=factor(end$Percentage[q4]),
                   y=end$Values[q4],
                   colour=end$Area[q4]))+
  theme(legend.position = "none")+ 
  geom_violin(aes(x=factor(end$Percentage[q4]),y=end$Values[q4]))+
  ylab("Taxonomic Distincness")+ xlab("Remove %")+
  ggtitle("Endemic species for 25° grid  (TD-Q4)")+
  theme(panel.background = element_rect(fill = "gray97"))

dev.off()

##################################################################################

end <- read.csv("PD_grid.endemic",header = T)


brks <- classIntervals(end$Values,n=5,style="quantile")
brks <- brks$brks
brks

q5 <- which(findInterval(end$Values,brks,all.inside = T)==5)
q4 <- which(findInterval(end$Values,brks,all.inside = T)==4)

png(filename = "PD_grid25Q5_end.png",width = 1070,height = 699,res = 120)

p <- ggplot()
p +geom_jitter(aes(x=factor(end$Percentage[q5]),
                   y=end$Values[q5],
                   colour=end$Area[q5]))+
  theme(legend.position = "none")+ 
  geom_violin(aes(x=factor(end$Percentage[q5]),y=end$Values[q5]))+
  ylab("Phylogenetic Diversity")+ xlab("Remove %")+
  ggtitle("Endemic species for 25° grid  (PD-Q5)")+
  theme(panel.background = element_rect(fill = "gray97"))

dev.off()

png(filename = "PD_grid25Q4_end.png",width = 1070,height = 699,res = 120)

p <- ggplot()
p +geom_jitter(aes(x=factor(end$Percentage[q4]),
                   y=end$Values[q4],
                   colour=end$Area[q4]))+
  theme(legend.position = "none")+ 
  geom_violin(aes(x=factor(end$Percentage[q4]),y=end$Values[q4]))+
  ylab("Phylogenetic Diversity")+ xlab("Remove %")+
  ggtitle("Endemic species for 25° grid  (PD-Q4)")+
  theme(panel.background = element_rect(fill = "gray97"))

dev.off()

##################################################################################

end <- read.csv("AvDT_grid.endemic",header = T)

brks <- classIntervals(end$Values,n=5,style="quantile")
brks <- brks$brks
brks

q5 <- which(findInterval(end$Values,brks,all.inside = T)==5)
q4 <- which(findInterval(end$Values,brks,all.inside = T)==4)

png(filename = "AvDT_grid25Q5_end.png",width = 1070,height = 699,res = 120)

p <- ggplot()
p +geom_jitter(aes(x=factor(end$Percentage[q5]),
                   y=end$Values[q5],
                   colour=end$Area[q5]))+
  theme(legend.position = "none")+ 
  geom_violin(aes(x=factor(end$Percentage[q5]),y=end$Values[q5]))+
  ylab("Average Taxonomic Disticness")+ xlab("Remove %")+
  ggtitle("Endemic species for 25° grid  (AvDT-Q5)")+
  theme(panel.background = element_rect(fill = "gray97"))

dev.off()

png(filename = "AvDT_grid25Q4_end.png",width = 1070,height = 699,res = 120)

p <- ggplot()
p +geom_jitter(aes(x=factor(end$Percentage[q4]),
                   y=end$Values[q4],
                   colour=end$Area[q4]))+
  theme(legend.position = "none")+ 
  geom_violin(aes(x=factor(end$Percentage[q4]),y=end$Values[q4]))+
  ylab("Average Taxonomic Distincness")+ xlab("Remove %")+
  ggtitle("Endemic species for 25° grid  (AvDT-Q4)")+
  theme(panel.background = element_rect(fill = "gray97"))

dev.off()


p <- ggplot(data = end, aes(x=factor(Percentage),y=Values))
p + geom_bar(aes(fill=Area),stat = "identity",position = "fill")


p <- ggplot(data = end, aes(x=factor(Area),y=Values))
p + geom_bar(aes(fill=factor(Percentage)),stat="identity",position = "dodge")
