library(ggplot2)

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final/Raw_IndexR/")

PNN <- data.frame(Area=read.csv("DT.PNN")$area,
                  TD=read.csv("DT.PNN")$W,
                  PD=read.csv("PD.PNN")$PD,
                  AvTD=read.csv("AvTD.PNN")$Dplus)
PNN

## Tres calculos
##             /PNNout+PNNin  /PNNout+PNN+ /NABin
## PNNin
## PNN+
## PNNout-NAB


## For TD

PNNin <- PNN$TD[grep("PNNin",PNN$Area)]
PNN1 <- sum(PNN$TD[1:117])
PNNNAB <- PNN$TD[grep("NABin",PNN$Area)] - PNN$TD[grep("PNNout",PNN$Area)] 
NAB <- PNN$TD[grep("NABin",PNN$Area)]

PNN.td <- matrix(0, ncol = 3, nrow = 3)

PNN.td[1,1] <- PNNin/(PNN$TD[grep("PNNout",PNN$Area)]+PNNin)
PNN.td[2,1] <- PNN1/(PNN$TD[grep("PNNout",PNN$Area)]+PNNin)
PNN.td[3,1] <- PNNNAB/(PNN$TD[grep("PNNout",PNN$Area)]+PNNin)

PNN.td[1,2] <- PNNin/(PNN$TD[grep("PNNout",PNN$Area)]+PNN1)
PNN.td[2,2] <- PNN1/(PNN$TD[grep("PNNout",PNN$Area)]+PNN1)
PNN.td[3,2] <- PNNNAB/(PNN$TD[grep("PNNout",PNN$Area)]+PNN1)

PNN.td[1,3] <- PNNin/NAB
PNN.td[2,3] <- PNN1/NAB
PNN.td[3,3] <- PNNNAB/NAB

colnames(PNN.td) <- c("/PNNout+PNNin","/PNNout+PNN+","/NAB")
rownames(PNN.td) <- c("PNNin","PNN+","PNNout-NAB")

PNN.td

## For PD

PNNin <- PNN$PD[grep("PNNin",PNN$Area)]
PNN1 <- sum(PNN$PD[1:117])
PNNNAB <- PNN$PD[grep("NABin",PNN$Area)] - PNN$PD[grep("PNNout",PNN$Area)] 
NAB <- PNN$PD[grep("NABin",PNN$Area)]

PNN.pd <- matrix(0, ncol = 3, nrow = 3)

PNN.pd[1,1] <- PNNin/(PNN$PD[grep("PNNout",PNN$Area)]+PNNin)
PNN.pd[2,1] <- PNN1/(PNN$PD[grep("PNNout",PNN$Area)]+PNNin)
PNN.pd[3,1] <- PNNNAB/(PNN$PD[grep("PNNout",PNN$Area)]+PNNin)

PNN.pd[1,2] <- PNNin/(PNN$PD[grep("PNNout",PNN$Area)]+PNN1)
PNN.pd[2,2] <- PNN1/(PNN$PD[grep("PNNout",PNN$Area)]+PNN1)
PNN.pd[3,2] <- PNNNAB/(PNN$PD[grep("PNNout",PNN$Area)]+PNN1)

PNN.pd[1,3] <- PNNin/NAB
PNN.pd[2,3] <- PNN1/NAB
PNN.pd[3,3] <- PNNNAB/NAB

colnames(PNN.pd) <- c("/PNNout+PNNin","/PNNout+PNN+","/NAB")
rownames(PNN.pd) <- c("PNNin","PNN+","PNNout-NAB")

PNN.pd


## For AvTD

PNNin <- PNN$AvTD[grep("PNNin",PNN$Area)]
PNN1 <- sum(PNN$AvTD[1:117])
PNNNAB <- (PNN$AvTD[grep("NABin",PNN$Area)] - PNN$AvTD[grep("PNNout",PNN$Area)]) 
NAB <- PNN$AvTD[grep("NABin",PNN$Area)]

PNN.avtd <- matrix(0, ncol = 3, nrow = 3)

PNN.avtd[1,1] <- PNNin/(PNN$AvTD[grep("PNNout",PNN$Area)]+PNNin)
PNN.avtd[2,1] <- PNN1/(PNN$AvTD[grep("PNNout",PNN$Area)]+PNNin)
PNN.avtd[3,1] <- PNNNAB/(PNN$AvTD[grep("PNNout",PNN$Area)]+PNNin)

PNN.avtd[1,2] <- PNNin/(PNN$AvTD[grep("PNNout",PNN$Area)]+PNN1)
PNN.avtd[2,2] <- PNN1/(PNN$AvTD[grep("PNNout",PNN$Area)]+PNN1)
PNN.avtd[3,2] <- PNNNAB/(PNN$AvTD[grep("PNNout",PNN$Area)]+PNN1)

PNN.avtd[1,3] <- PNNin/NAB
PNN.avtd[2,3] <- PNN1/NAB
PNN.avtd[3,3] <- PNNNAB/NAB

colnames(PNN.avtd) <- c("/PNNout+PNNin","/PNNout+PNN+","/NAB")
rownames(PNN.avtd) <- c("PNNin","PNN+","PNNout-NAB")

PNN.avtd



setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final/Raw_IndexR/")

write.table(PNN.td,
            file = "TD_PNN.r",
            col.names = T, row.names = T, quote = F, sep = ",")

write.table(PNN.pd,
            file = "PD_PNN.r",
            col.names = T, row.names = T, quote = F, sep = ",")

write.table(PNN.avtd,
            file = "AvTD_PNN.r",
            col.names = T, row.names = T, quote = F, sep = ",")

####################################################################################
####################################################################################

PNN2 <- PNN[-c(118,119,120), ]

PNN2.zeroTD <- length(which(PNN2$TD==0))
PNN2.zeroPD <- length(which(PNN2$PD==0))
PNN2.zeroAvTD <- length(which(PNN2$AvTD==0))

(PNN2.zeroTD/length(PNN2$Area))*100
(PNN2.zeroPD/length(PNN2$Area))*100
(PNN2.zeroAvTD/length(PNN2$Area))*100

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final/")

PNN.tr0 <- read.csv("PNN.dist.matrix",header = T)

PNN.tr1 <- PNN.tr0[,-1]

PNN.sumTR <- as.matrix(apply(PNN.tr1, 2, sum))

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final/Raw_IndexR/")

PNNTD <- read.csv("DT.PNN", header = T)
PNNTD <- as.data.frame(cbind(PNNTD,PNN.sumTR))
PNNTD <- PNNTD[-c(118,119,120),]
PNNTD <- PNNTD[-which(PNNTD$W==0),]
PNNTD <- PNNTD[order(PNNTD$W, decreasing = T),]

h <- sqrt(length(PNNTD$W))

k <- diff(range(PNNTD$W))/h

ggplot()+geom_histogram(aes(PNNTD$W), bins = k, binwidth = h)

ggplot()+geom_point(aes(x=PNNTD$rich,y=PNNTD$W))+xlab("Phylo-Richness")+ylab("W index")

cor.test(PNNTD$rich,PNNTD$W, method = 'spearman')$estimate

ggplot()+geom_point(aes(x=PNNTD$PNN.sumTR,y=PNNTD$W))+xlab("Total-Richness")+xlab("W index")

cor.test(PNNTD$PNN.sumTR,PNNTD$W, method = 'spearman')$estimate

png("PNN_DT2.png",width = 3000,height = 1240)

ggplot() + geom_point(aes(x=1:length(PNNTD$area),y=PNNTD$W))+
  xlab("NNP")+ylab("W index")+
  geom_point(aes(x=1:length(PNNTD$area),y=PNNTD$rich),colour="red")
#annotate("text", label = PNNTD$area, y = (PNNTD$W)+0.001, 
#           x = 1:length(PNNTD$area), size = 12, colour = "black",angle=0)

dev.off()
#####################################################################

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final/Raw_IndexR/")

PNNPD <- read.csv("PD.PNN", header = T) ## Read the PD values
PNNPD <- as.data.frame(cbind(PNNPD,PNN.sumTR)) # Combine PD matrix with Total Richness
PNNPD <- PNNPD[-c(118,119,120),] # Remove the NABin, PNNin and PNNout
PNNPD <- PNNPD[-which(PNNPD$PD==0),] # Remove PNN with PD value == 0
PNNPD <- PNNPD[order(PNNPD$PD, decreasing = T),] # Order the matrix in decreasing order

h <- sqrt(length(PNNPD$PD))

k <- diff(range(PNNPD$PD))/h

ggplot()+geom_histogram(aes(PNNPD$PD), bins = k, binwidth = h)

ggplot()+geom_point(aes(x=PNNPD$SR,y=PNNPD$PD))+xlab("Phylo-Richness")+xlab("PD index")

cor.test(PNNPD$SR,PNNPD$PD, method = 'spearman')$estimate

ggplot()+geom_point(aes(x=PNNPD$PNN.sumTR,y=PNNPD$PD))+xlab("Total-Richness")+xlab("PD index")

cor.test(PNNPD$PNN.sumTR,PNNPD$PD, method = 'spearman')$estimate

png("PNN_PD2.png",width = 3000,height = 1240)

ggplot() + geom_point(aes(x=1:length(PNNPD$PD),y=PNNPD$PD))+
  xlab("NNP")+ylab("PD index")+
  geom_point(aes(x=1:length(PNNPD$PD),y=PNNPD$SR),colour="red")
#annotate("text", label = PNNTD$area, y = (PNNTD$W)+0.001, 
#           x = 1:length(PNNTD$area), size = 12, colour = "black",angle=0)

dev.off()
#####################################################################


setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final/Raw_IndexR/")

PNNAvTD <- read.csv("AvTD.PNN", header = T)
PNNAvTD <- as.data.frame(cbind(PNNAvTD,PNN.sumTR))
PNNAvTD <- PNNAvTD[-c(118,119,120),]
PNNAvTD <- PNNAvTD[-which(PNNAvTD$Dplus==0),]
PNNAvTD <- PNNAvTD[order(PNNAvTD$Dplus, decreasing = T),]

h <- sqrt(length(PNNAvTD$Dplus))

k <- diff(range(PNNAvTD$Dplus))/h

ggplot()+geom_histogram(aes(PNNAvTD$Dplus), bins = k, binwidth = h)

ggplot()+geom_point(aes(x=PNNAvTD$Species,y=PNNAvTD$Dplus))+xlab("Phylo-Richness")+ylab("AvTD index")

cor.test(PNNAvTD$Species,PNNAvTD$Dplus, method = 'spearman')$estimate

ggplot()+geom_point(aes(x=PNNAvTD$PNN.sumTR,y=PNNAvTD$Dplus))+xlab("Total-Richness")+ylab("AvTD index")

cor.test(PNNAvTD$PNN.sumTR,PNNAvTD$Dplus, method = 'spearman')$estimate

png("PNN_AvDT2.png",width = 3000,height = 1240)

ggplot() + geom_point(aes(x=1:length(PNNAvTD$Species),y=log(PNNAvTD$Dplus)))+
  xlab("NNP")+ylab("AvTD index")+
  geom_point(aes(x=1:length(PNNAvTD$Species),y=PNNAvTD$Species),colour="red")
#annotate("text", label = PNNTD$area, y = (PNNTD$W)+0.001, 
#           x = 1:length(PNNTD$area), size = 12, colour = "black",angle=0)

dev.off()
#####################################################################

Areas <- cbind(as.character(PNNTD$area),rownames(PNNPD), rownames(PNNAvTD))
Areas <- Areas[1:length(rownames(PNNAvTD)), ]
colnames(Areas) <- c("TD", "PD", "AvTD")
Areas <- as.data.frame(Areas)

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final/")

GI <- read.csv("General.info", header = T)
match.sp <- read.csv("Match.sp")

PNN.m <- PNN.tr0[which(PNN.tr0$especie%in%GI$Sp), ]

for(i in 2:length(PNN.m$especie)){
  
  bl <- GI$BL[grep(PNN.m$especie[i],GI$Sp)]
  
  if(length(bl)>1){
    bl <- sum(bl)/length(bl)
  }
  
  PNN.m[i,grep(1,PNN.m[i,])] <- bl
  
}

PNN.m <- PNN.m[,-1]

BL.total <- as.matrix(apply(PNN.m, 2, sum))

BL <- matrix(NA, nrow = length(Areas$TD), ncol = 3)

for(i in 1:length(colnames(Areas))){
  
  for (j in 1:length(Areas$TD)) {
    
    bls <- BL.total[grep(Areas[j,i],rownames(BL.total))]
    if(length(bls)>1){ bls <- bls[1]}
    BL[j,i] <- bls
  }
  
}  

colnames(BL) <- c("BL_TD","BL_PD","BL_AvTD")
BL <- as.data.frame(BL)

Areas.f <- as.data.frame(cbind(as.character(Areas$TD),as.character(BL$BL_TD),
                               as.character(Areas$PD),as.character(BL$BL_PD),
                               as.character(Areas$AvTD),as.character(BL$BL_AvTD)))

colnames(Areas.f) <- c("TD","BL","PD","BL","AvTD","BL")

Areas.f

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Julio1/RawIndex/")

write.table(Areas.f, "PNN_Class.table",quote = F, col.names = T, row.names = F)
