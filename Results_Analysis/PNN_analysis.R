library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)

# Set the working directory
setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

# Read the data, and create a unique matrix
PNN <- data.frame(Area=read.csv("TD.PNN")$area, #Read the PA names
                  TD=read.csv("TD.PNN")$W, # Read TD index values
                  PD=read.csv("PD.PNN")$PD, # Read PD index values
                  AvTD=read.csv("AvTD.PNN")$Dplus) # Read AvTD index values

# Check the resulting matrix
PNN 

## Three different calculates
##             /PNNout+PNNin  /PNNout+PNN+ /NABin
## PNNin
## PNN+
## PNNout-NAB

#################
##   For TD    ##
#################

PNNin <- PNN$TD[grep("PNNin",PNN$Area)] # Find the value for the PA
PNN1 <- sum(PNN$TD[1:117]) # Sum the value of each PA to make a total value
PNNNAB <- PNN$TD[grep("NABin",PNN$Area)] - PNN$TD[grep("PNNout",PNN$Area)] # Difference between NAB value and Outside PA value
NAB <- PNN$TD[grep("NABin",PNN$Area)] # Find the value for the NAB

#Create a outcome matrix
PNN.td <- matrix(0, ncol = 3, nrow = 3)

# Make the calculations

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

#################
##   For PD    ##
#################

PNNin <- PNN$PD[grep("PNNin",PNN$Area)]# Find the value for the PA
PNN1 <- sum(PNN$PD[1:117])# Sum the value of each PA to make a total value
PNNNAB <- PNN$PD[grep("NABin",PNN$Area)] - PNN$PD[grep("PNNout",PNN$Area)] # Difference between NAB value and Outside PA value
NAB <- PNN$PD[grep("NABin",PNN$Area)]# Find the value for the NAB

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

#################
##  For AvTD   ##
#################

PNNin <- PNN$AvTD[grep("PNNin",PNN$Area)]# Find the value for the PA
PNN1 <- sum(PNN$AvTD[1:117])# Sum the value of each PA to make a total value
PNNNAB <- (PNN$AvTD[grep("NABin",PNN$Area)] - PNN$AvTD[grep("PNNout",PNN$Area)])# Difference between NAB value and Outside PA value 
NAB <- PNN$AvTD[grep("NABin",PNN$Area)]# Find the value for the NAB

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

#################################################################################
#################################################################################

PNN.td
PNN.pd
PNN.avtd

#################################################################################
#################################################################################


setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

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
# Now, another question to answer: How many PA have index values equal zero for all three indices ??
# Does the same PA for all three indices ??


# Firsts, remove the rows of PNNin, PNNout and NABin, these rows do not interest us
PNN2 <- PNN[-c(116,117,118), ]

# Find the number of PA with index value equal 0
PNN2.zeroTD <- length(which(PNN2$TD==0))
PNN2.zeroPD <- length(which(PNN2$PD==0))
PNN2.zeroAvTD <- length(which(PNN2$AvTD==0))

# Now, calculate the percentage
(PNN2.zeroTD/length(PNN2$Area))*100 # For TD
(PNN2.zeroPD/length(PNN2$Area))*100 # For PD
(PNN2.zeroAvTD/length(PNN2$Area))*100 # For AvTD

#####################################################################################

############################################################################
#                   PNN indices values vs Richness                         #
############################################################################

# Set a new working directory
setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/")

# Read the absence/presence matrix for the PA
PNN.tr0 <- read.csv("PNN.dist.matrix",header = T)

# Remove the species column
PNN.tr1 <- PNN.tr0[,-1]

#Make the summatory of all species who inhabit each PA
PNN.sumTR <- as.matrix(apply(PNN.tr1, 2, sum))

# Set a new wokking directory
setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R//")

# Read the matrix for the TD index
PNNTD <- read.csv("TD.PNN", header = T)
# join the TD and the total richness matrices
PNNTD <- as.data.frame(cbind(PNNTD,PNN.sumTR))
# Remove the PNNin, PNNout and NABin rows
PNNTD <- PNNTD[-c(118,117,116),]
# Remove the PA with index values equal zero
PNNTD <- PNNTD[-which(PNNTD$W==0),]
# Ordeing the matrix in decreasing order given the TD index
PNNTD <- PNNTD[order(PNNTD$W, decreasing = T),]


## Here the h and k variables will be calculate for the histogram of the TD index
h <- sqrt(length(PNNTD$W))

k <- diff(range(PNNTD$W))/h

# Plot the histogram
ggplot()+geom_histogram(aes(PNNTD$W, fill=..count..), bins = k, binwidth = h)+xlab("TD Index")+ylab("Count")

# Spearman Correlation TD index vs PhyloRichness
cor.test(PNNTD$rich,PNNTD$W, method = 'spearman', exact = F)$estimate

est <- paste("R=",round(cor.test(PNNTD$rich,PNNTD$W, method = 'spearman', exact = F)$estimate, digits = 3),sep="")

#Plot TD index vs PhyloRichnes
cor1 <- ggplot()+geom_point(aes(x=PNNTD$rich,y=PNNTD$W))+xlab("Phylo-Richness")+ylab("W index")+
  geom_label(aes(x=25,y=300), label=est, size=10)+
  geom_smooth(aes(x=PNNTD$rich,y=PNNTD$W), method = "lm")

#Spearman correlation TD index vs TotalRichness
cor.test(PNNTD$PNN.sumTR,PNNTD$W, method = 'spearman', exact = F)$estimate

est2 <- paste("R=",round(cor.test(PNNTD$PNN.sumTR,PNNTD$W, method = 'spearman', exact = F)$estimate, digits = 3),sep="")

#Plot TD index vs TotalRichness
cor2 <- ggplot()+geom_point(aes(x=PNNTD$PNN.sumTR,y=PNNTD$W))+xlab("Total-Richness")+ylab("W index")+
  geom_label(aes(x=50,y=300), label=est2, size=10)+
  geom_smooth(aes(x=PNNTD$PNN.sumTR,y=PNNTD$W), method = "lm")



png("PNN2RichVsIndexTD.png")
grid.arrange(cor1,cor2,heights = c(2/4, 2/4))
dev.off()


p1 <- ggplot() + geom_point(aes(x=1:length(PNNTD$area),y=log(PNNTD$W)),colour="black")+
  xlab("Protected Areas")+ylab("log TD Index")+
  theme(legend.title=element_blank(),
        legend.text=element_text(face="bold"),
        legend.position="bottom",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0,5,1,1),units="points"))
p1  

p2 <- ggplot()+geom_point(aes(x=1:length(PNNTD$area),y=log(PNNTD$rich)),colour="red")+
  xlab(NULL)+ylab("log Phylogenetic-Richness")+
  theme(legend.title=element_blank(),
        legend.text=element_text(face="bold"),
        legend.position="bottom",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0,5,1,1),units="points"))
p2

png("PNN_DT_RichVsIndex.png")
grid.arrange(p2,p1,heights = c(2/4, 2/4))
dev.off()

#####################################################################

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R//")

PNNPD <- read.csv("PD.PNN", header = T) ## Read the PD values
PNNPD <- as.data.frame(cbind(PNNPD,PNN.sumTR)) # Combine PD matrix with Total Richness
PNNPD <- PNNPD[-c(118,117,116),] # Remove the NABin, PNNin and PNNout
PNNPD <- PNNPD[-which(PNNPD$PD==0),] # Remove PNN with PD value == 0
PNNPD <- PNNPD[order(PNNPD$PD, decreasing = T),] # Order the matrix in decreasing order

# Calcualte the h and k variable for the histogram 
h <- sqrt(length(PNNPD$PD))

k <- diff(range(PNNPD$PD))/h
#Plot the hhistogram
ggplot()+geom_histogram(aes(PNNPD$PD, fill=..count..), bins = k, binwidth = h)+xlab("TD Index")+ylab("Count")

cor.test(PNNPD$SR,PNNPD$PD, method = 'spearman', exact = F)$estimate

est <- paste("R=",round(cor.test(PNNPD$SR,PNNPD$PD, method = 'spearman', exact = F)$estimate, digits = 3),sep="")

#Plot the PD index vs Phylogenetic Richness
cor1 <- ggplot()+geom_point(aes(x=PNNPD$SR,y=PNNPD$PD))+xlab("Phylo-Richness")+ylab("PD index")+
  geom_label(aes(x=35,y=90), label=est, size=10)+
  geom_smooth(aes(x=PNNPD$SR,y=PNNPD$PD), method = "lm")


cor.test(PNNPD$PNN.sumTR,PNNPD$PD, method = 'spearman', exact = F)$estimate

est2 <- paste("R=",round(cor.test(PNNPD$PNN.sumTR,PNNPD$PD, method = 'spearman', exact = F)$estimate, digits = 3),sep="")


# Spearman corraltion for PD and Phylogenetic richness
cor2 <- ggplot()+geom_point(aes(x=PNNPD$PNN.sumTR,y=PNNPD$PD))+xlab("Total-Richness")+ylab("PD index")+
  geom_label(aes(x=65,y=90), label=est2, size=10)+
  geom_smooth(aes(x=PNNPD$PNN.sumTR,y=PNNPD$PD), method = "lm")


png("PNN2RichVsIndexPD.png")
grid.arrange(cor1,cor2,heights = c(2/4, 2/4))
dev.off()


p1 <- ggplot() + geom_point(aes(x=1:length(PNNTD$area),y=log(PNNPD$PD)),colour="black")+
  xlab("Protected Areas")+ylab("log PD Index")+
  theme(legend.title=element_blank(),
        legend.text=element_text(face="bold"),
        legend.position="bottom",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0,5,1,1),units="points"))
p1  

p2 <- ggplot()+geom_point(aes(x=1:length(PNNTD$area),y=log(PNNPD$SR)),colour="red")+
  xlab(NULL)+ylab("log Phylogenetic-Richness")+
  theme(legend.title=element_blank(),
        legend.text=element_text(face="bold"),
        legend.position="bottom",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0,5,1,1),units="points"))
p2

png("PNN_PD_RichVsIndex.png")
grid.arrange(p2,p1,heights = c(2/4, 2/4))
dev.off()
#####################################################################


setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

PNNAvTD <- read.csv("AvTD.PNN", header = T) ## Read the AvTD values
PNNAvTD <- as.data.frame(cbind(PNNAvTD,PNN.sumTR))# Combine AvTD matrix with Total Richness
PNNAvTD <- PNNAvTD[-c(118,117,116),]# Remove the NABin, PNNin and PNNout
PNNAvTD <- PNNAvTD[-which(PNNAvTD$Dplus==0),]# Remove PNN with PD value == 0
PNNAvTD <- PNNAvTD[order(PNNAvTD$Dplus, decreasing = T),] # Order the matrix in decreasing order

#Calculate the h and k variables for the histogram
h <- sqrt(length(PNNAvTD$Dplus))

k <- diff(range(PNNAvTD$Dplus))/h

#Plot the histogram
ggplot()+geom_histogram(aes(PNNAvTD$Dplus, fill=..count..), bins = k, binwidth = h)+xlab("AvTD Index")+ylab("Count")

#Spearman corraltion for AvTD and Phylogenetic Richness
cor.test(PNNAvTD$Species,PNNAvTD$Dplus, method = 'spearman', exact = F)$estimate

est <- paste("R=",round(cor.test(PNNAvTD$Species,PNNAvTD$Dplus, method = 'spearman', exact = F)$estimate, digits = 3),sep="")

#Plot AvTD vs Phylogenetic Richness
cor1 <- ggplot()+geom_point(aes(x=PNNAvTD$Species,y=PNNAvTD$Dplus))+xlab("Phylo-Richness")+ylab("AvTD index")+
  geom_label(aes(x=50,y=50), label=est, size=10)+
  geom_smooth(aes(x=PNNAvTD$Species,y=PNNAvTD$Dplus), method = "lm")

#Spearman correlation for AvTD and Total richness
cor.test(PNNAvTD$PNN.sumTR,PNNAvTD$Dplus, method = 'spearman', exact = F)$estimate

est2 <- paste("R=",round(cor.test(PNNAvTD$PNN.sumTR,PNNAvTD$Dplus, method = 'spearman', exact = F)$estimate, digits = 3),sep="")

#Plot AvTD vs Total Richness
cor2 <- ggplot()+geom_point(aes(x=PNNAvTD$PNN.sumTR,y=PNNAvTD$Dplus))+xlab("Phylo-Richness")+ylab("AvTD index")+
  geom_label(aes(x=70,y=50), label=est2, size=10)+
  geom_smooth(aes(x=PNNAvTD$PNN.sumTR,y=PNNAvTD$Dplus), method = "lm")


png("PNN2RichVsIndexAvTD.png")
grid.arrange(cor1,cor2,heights = c(2/4, 2/4))
dev.off()


p1 <- ggplot() + geom_point(aes(x=1:length(PNNAvTD$Species),y=log(PNNAvTD$Dplus)),colour="black")+
  xlab("Protected Areas")+ylab("log AvTD Index")+
  theme(legend.title=element_blank(),
        legend.text=element_text(face="bold"),
        legend.position="bottom",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0,5,1,1),units="points"))
p1  

p2 <- ggplot()+geom_point(aes(x=1:length(PNNAvTD$Species),y=log(PNNAvTD$Species)),colour="red")+
  xlab(NULL)+ylab("log Phylogenetic-Richness")+
  theme(legend.title=element_blank(),
        legend.text=element_text(face="bold"),
        legend.position="bottom",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0,5,1,1),units="points"))
p2

png("PNN_AvTD_RichVsIndex.png")
grid.arrange(p2,p1,heights = c(2/4, 2/4))
dev.off()

#####################################################################
#####################################################################
# Now, lets make a matrix with some util information from the Protected Areas

#bind the name areas from the three different indices
Areas <- cbind(as.character(PNNTD$area),rownames(PNNPD), rownames(PNNAvTD))
#Cause, AvTD have less protected areas with values > zero, the number of PA for AvTD
#will be the row length of the mmatrix
Areas <- Areas[1:length(rownames(PNNAvTD)), ]
#Assing the column names and treat as data frame
colnames(Areas) <- c("TD", "PD", "AvTD")
Areas <- as.data.frame(Areas)

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/")

# Read the table ith the general information
GI <- read.csv("General.info", header = T)
# And the list of species that are shared between phylogenies and absence/presence matrices
match.sp <- read.csv("Match.sp")

#Extract from the original absence/presence matrix the species from the General Information matrix
PNN.m <- PNN.tr0[which(PNN.tr0$especie%in%GI$Sp), ]


# Here for each species, their Branch Length value will be assignef in the places where the species have presence
# There, is in a place the Specie 1 have a presence (==1), there, that 1 will be replace for the Species's branch length


for(i in 2:length(PNN.m$especie)){
  
  bl <- GI$BL[grep(PNN.m$especie[i],GI$Sp)]
  
  if(length(bl)>1){ print(i)
    bl <- sum(bl)/length(bl)
  }
  
  PNN.m[i,grep(1,PNN.m[i,])] <- bl
  
}

warnings()

#Remove the first column, this contain the species name, and is not for our interest
PNN.m <- PNN.m[,-1]

# Make the total sumatory for each PA
BL.total <- as.matrix(apply(PNN.m, 2, sum))

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

write.csv(BL.total, "PNN_BL_Total", quote = F, row.names = T)

# Create and empty matrix
BL <- matrix(NA, nrow = length(Areas$TD), ncol = 3)
# And do the same make it below but for each index
for(i in 1:length(colnames(Areas))){
  
  for (j in 1:length(Areas$TD)) {
    
    bls <- BL.total[grep(Areas[j,i],rownames(BL.total))]
    if(length(bls)>1){ bls <- bls[1]}
    BL[j,i] <- bls
  }
  
}  

#Assing the column names
colnames(BL) <- c("BL_TD","BL_PD","BL_AvTD")
BL <- as.data.frame(BL)

#Create a new matrix with the BL values and the name of the PA
Areas.f <- as.data.frame(cbind(as.character(Areas$TD),as.character(BL$BL_TD),
                               as.character(Areas$PD),as.character(BL$BL_PD),
                               as.character(Areas$AvTD),as.character(BL$BL_AvTD)))

#Column naes
colnames(Areas.f) <- c("TD","BL","PD","BL","AvTD","BL")
#Check it
Areas.f

#Save it
setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/RawIndex_R/")

write.table(Areas.f, "PNN_Class.table",quote = F, col.names = T, row.names = F)
