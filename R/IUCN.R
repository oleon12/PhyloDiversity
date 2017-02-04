library(letsR)
library(ggplot2)

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/")

End <- read.csv("End.csv")
WD <- read.csv("WD.csv")

colnames(End) <- c("Sp", "BL")
colnames(WD) <- c("Sp", "BL")

head(End)
head(WD)

##########################################################3

Status <- rep(NA, length(End$Sp))
Criteria <- rep(NA, length(End$Sp))
Population <- rep(NA, length(End$Sp))

for(i in 1:length(End$Sp)){
  
  print(paste(paste("Species",i, sep = " "),paste("of",length(End$Sp), sep = " "), sep = " "))
  
  IUCN <- lets.iucn(End$Sp[i])
  
  Status[i] <- as.character(IUCN$Status)
  
  Criteria[i] <- as.character(IUCN$Criteria)
  
  Population[i] <- as.character(IUCN$Population)
  
}

End.End <- cbind(End,Status,Criteria,Population)

ggplot(End.End)+geom_bar(aes(Status))

length(grep("NT", End.End$Status))

##########################################################3

Status <- rep(NA, length(WD$Sp))
Criteria <- rep(NA, length(WD$Sp))
Population <- rep(NA, length(WD$Sp))

for(i in 913:length(WD$Sp)){
  
  print(paste(paste("Species",i, sep = " "),paste("of",length(WD$Sp), sep = " "), sep = " "))
  
  IUCN <- lets.iucn(WD$Sp[i])
  
  Status[i] <- as.character(IUCN$Status)
  
  Criteria[i] <- as.character(IUCN$Criteria)
  
  Population[i] <- as.character(IUCN$Population)
  
}


WD.WD <- cbind(WD,Status,Criteria,Population)

ggplot(WD.WD)+geom_bar(aes(Status))

length(grep("NE", WD.WD$Status))


End <- rep("Endemic", length(End.End$Sp))
NonEnd <- rep("Non-Endemic", length(WD.WD$Sp))

SpType <- as.matrix(c(End,NonEnd))

AllSp <- rbind(End.End, WD.WD)

AllSp <- cbind(AllSp, SpType)

head(AllSp)

AllSp <- AllSp[-grep("LR", AllSp$Status),]

p <- ggplot(AllSp)+geom_bar(aes(Status, fill=SpType), position="dodge")+
  xlab("IUCN Category")+ylab("Species")+
  scale_x_discrete(limit=c("CR","EN","VU","LC","NE","NT","DD"))+
  scale_fill_discrete(name="Species")+
  theme(legend.position = c(.2,.7),
        legend.title=element_text(size = 30,face = "bold"),
        legend.text = element_text(size = 30),
        legend.key.size = unit(1, "cm"),
        legend.key = element_rect(fill= "white",colour="black"),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15),
        axis.title.y=element_text(size=30),
        axis.title.x=element_text(size=30))+  
  guides(colour = guide_legend(override.aes = list(size=2)))

p

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/Endemic/")


png("IUCN.png", width =2000 , height =1300 ,res = 200)
p
dev.off()
