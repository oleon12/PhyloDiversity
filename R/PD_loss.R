#Load libraries

library(ape)
library(phytools)
library(picante)
library(phangorn)
library(jrich)
library(rgeos)
library(rgdal)
library(factoextra)
library(cluster)
library(ggplot2)
library(rgeos)
library(maps)
library(maptools)
library(rgdal)
library(ggplot2)
library(gridExtra)


source("...phylo.id.table.R")
source("...pd2 .b.R")
source("...toNum.R")
source("...match.phy.dist.R")
source("...phylo.beta.R")
source("...phylo.id.table.R")
source("...get.cluster.n.R")


#### Cluster Btotal
b.total <- readOGR("../Data/shp/Grid25_Cluster.t.shp")@data

#### PD Data
pd <- readOGR("../Data/shp/Grid25_PD.shp")@data
Cells <- as.character(1:length(pd$SP_ID))
pd <- data.frame(Id=Cells[which(pd$Index=="Q5")],
                 Class=pd$Index[which(pd$Index=="Q5")],
                 PD=pd$IndexVals[which(pd$Index=="Q5")],
                 Clust =b.total$Cluster[which(pd$Index=="Q5")])
head(pd)

#################################################################
#                   Remotion for all Q5 Cells                   #
#################################################################

#Remotion %
per <- c(0.125,0.25,0.50)
#Number of Iterations
it <- 1:100

#Output matrices
pd.it <- matrix(NA,nrow = length(it), ncol = length(per) )
cell.12 <- matrix(NA, ncol = round(length(pd$Id)*per[1], digits = 0), nrow = it)
cell.25 <- matrix(NA, ncol = round(length(pd$Id)*per[2], digits = 0), nrow = it)
cell.50 <- matrix(NA, ncol = round(length(pd$Id)*per[3], digits = 0), nrow = it)

cell.id <- list(cell.12, cell.25, cell.50)

ind <- pd

for(i in 1:length(per)){
  
  for(j in 1:length(it)){
    
    #Sample i percentage of Q5 cells
    samp <-sample(1:length(pd$Id),round(length(pd$Id)*per[i], digits = 0))
    #Remove their PD values
    val <- pd$PD[-samp]
    #Save the id of those cells removed
    cells <- pd$Id[samp]
    
    #Calculate total PD for the remain Q5 cells
    sum.pd <- sum(val)
    #Save in the respective position
    pd.it[j,i] <- sum.pd
    
    cell.id[[i]] <- rbind(cell.id[[i]], t(as.data.frame(cells)))
  }
}

pd.it <- as.data.frame(pd.it)
colnames(pd.it) <- c("p12","p25","p50")
head(pd.it)

#Here, total PD for each iteration removing X% of Q5 cells
summary(pd.it)
summary(pd$PD)

###################################################
#Plot1
p12 <- ggplot()+
  geom_point(aes(x=it,y=pd.it$p12))+
  geom_line(aes(x=it,y=pd.it$p12))+
  xlab(NULL)+ylab("12.5%")
p25 <- ggplot()+  
  geom_point(aes(x=it,y=pd.it$p25))+
  geom_line(aes(x=it,y=pd.it$p25))+
  xlab(NULL)+ylab("25%")
p50 <- ggplot()+  
  geom_point(aes(x=it,y=pd.it$p50))+
  geom_line(aes(x=it,y=pd.it$p50))+
  xlab("Iterations")+ylab("50%")


grid.arrange(p12,p25,p50, ncol=1, nrow=3)

#####################################################
#Plot 2
ggplot()+
  geom_boxplot(aes(x="12.5%", y=pd.it$p12))+
  geom_boxplot(aes(x="25%", y=pd.it$p25))+
  geom_boxplot(aes(x="50%", y=pd.it$p50))+
  xlab("Remotion percentage")+ylab("PD values per iteration")

#####################################################
#Plot3

pds <- c(sum(pd$PD), mean(pd.it$p12), mean(pd.it$p25), mean(pd.it$p50))
col <- factor(c("Total PD", "-12.5% (mean)", "-25% (mean)", "-50% (mean)"))

ggplot()+
  geom_bar(aes(y=pds, x=col), stat="identity")+
  scale_x_discrete(limits=col)+
  xlab(NULL)+ylab("PD")

####################################################
# 

minpdloss <-c(0,(sum(pd$PD)-min(pd.it$p12))/sum(pd$PD),
              (sum(pd$PD)-min(pd.it$p25))/sum(pd$PD),
              (sum(pd$PD)-min(pd.it$p50))/sum(pd$PD))
maxpdloss <-c(0,(sum(pd$PD)-max(pd.it$p12))/sum(pd$PD),
              (sum(pd$PD)-max(pd.it$p25))/sum(pd$PD),
              (sum(pd$PD)-max(pd.it$p50))/sum(pd$PD))
meanpdloss <-c(0,(sum(pd$PD)-mean(pd.it$p12))/sum(pd$PD),
               (sum(pd$PD)-mean(pd.it$p25))/sum(pd$PD),
               (sum(pd$PD)-mean(pd.it$p50))/sum(pd$PD))

ggplot()+
  geom_bar(aes(x=c("12.5%","25%","50%"),y=meanpdloss[-1]), stat="identity")+
  #geom_errorbar(aes(ymin=minpdloss[-1], ymax=maxpdloss[-1], x=c("12.5%","25%","50%")), width=.1)+
  xlab("% of Remotion")+ylab("% of PD loss")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size=15))

ggsave(filename = "PDloss_All.pdf", device = "pdf", width = 20, height = 10, dpi = 300, units = "cm")

########################################################
# Same as above, but for each cluster area

d1 <- pd[which(pd$Clust==1),]
d2 <- pd[which(pd$Clust==2),]
d3 <- pd[which(pd$Clust==3),]

########################################################

#Remotion %
per <- c(0.125,0.25,0.50)
#Number of Iterations
it <- 1:100

#Output matrices
pd.it1 <- matrix(NA,nrow = length(it), ncol = length(per) )
cell.12 <- matrix(NA, ncol = round(length(d1$Id)*per[1], digits = 0), nrow = it)
cell.25 <- matrix(NA, ncol = round(length(d1$Id)*per[2], digits = 0), nrow = it)
cell.50 <- matrix(NA, ncol = round(length(d1$Id)*per[3], digits = 0), nrow = it)

cell.id <- list(cell.12, cell.25, cell.50)

ind <- d1

for(i in 1:length(per)){
  
  for(j in 1:length(it)){
    
    #Sample i percentage of Q5 cells
    samp <-sample(1:length(d1$Id),round(length(d1$Id)*per[i], digits = 0))
    #Remove their PD values
    val <- d1$PD[-samp]
    #Save the id of those cells removed
    cells <- d1$Id[samp]
    
    #Calculate total PD for the remain Q5 cells
    sum.pd <- sum(val)
    #Save in the respective position
    pd.it1[j,i] <- sum.pd
    
    cell.id[[i]] <- rbind(cell.id[[i]], t(as.data.frame(cells)))
  }
}

pd.it1 <- as.data.frame(pd.it1)
colnames(pd.it1) <- c("p12","p25","p50")
head(pd.it1)

#Here, total PD for each iteration removing X% of Q5 cells
summary(pd.it1)
summary(pd$PD)

pds1 <- c(sum(d1$PD), mean(pd.it1$p12), mean(pd.it1$p25), mean(pd.it1$p50))
col1 <- factor(c("Total PD", "-12.5% (mean)", "-25% (mean)", "-50% (mean)"))

ggplot()+
  geom_bar(aes(y=pds1, x=col1), stat="identity")+
  scale_x_discrete(limits=col1)+
  xlab(NULL)+ylab("PD")

########################################################

#Remotion %
per <- c(0.125,0.25,0.50)
#Number of Iterations
it <- 1:100

#Output matrices
pd.it2 <- matrix(NA,nrow = length(it), ncol = length(per) )
cell.12 <- matrix(NA, ncol = round(length(d2$Id)*per[1], digits = 0), nrow = it)
cell.25 <- matrix(NA, ncol = round(length(d2$Id)*per[2], digits = 0), nrow = it)
cell.50 <- matrix(NA, ncol = round(length(d2$Id)*per[3], digits = 0), nrow = it)

cell.id <- list(cell.12, cell.25, cell.50)

ind <- d2

for(i in 1:length(per)){
  
  for(j in 1:length(it)){
    
    #Sample i percentage of Q5 cells
    samp <-sample(1:length(d2$Id),round(length(d2$Id)*per[i], digits = 0))
    #Remove their PD values
    val <- d2$PD[-samp]
    #Save the id of those cells removed
    cells <- d2$Id[samp]
    
    #Calculate total PD for the remain Q5 cells
    sum.pd <- sum(val)
    #Save in the respective position
    pd.it2[j,i] <- sum.pd
    
    cell.id[[i]] <- rbind(cell.id[[i]], t(as.data.frame(cells)))
  }
}

pd.it2 <- as.data.frame(pd.it2)
colnames(pd.it2) <- c("p12","p25","p50")
head(pd.it2)

#Here, total PD for each iteration removing X% of Q5 cells
summary(pd.it2)
summary(pd$PD)

pds2 <- c(sum(d2$PD), mean(pd.it2$p12), mean(pd.it2$p25), mean(pd.it2$p50))
col2 <- factor(c("Total PD", "-12.5% (mean)", "-25% (mean)", "-50% (mean)"))

ggplot()+
  geom_bar(aes(y=pds2, x=col2), stat="identity")+
  scale_x_discrete(limits=col2)+
  xlab(NULL)+ylab("PD")

########################################################

#Remotion %
per <- c(0.125,0.25,0.50)
#Number of Iterations
it <- 1:100

#Output matrices
pd.it3 <- matrix(NA,nrow = length(it), ncol = length(per) )
cell.12 <- matrix(NA, ncol = round(length(d3$Id)*per[1], digits = 0), nrow = it)
cell.25 <- matrix(NA, ncol = round(length(d3$Id)*per[2], digits = 0), nrow = it)
cell.50 <- matrix(NA, ncol = round(length(d3$Id)*per[3], digits = 0), nrow = it)

cell.id <- list(cell.12, cell.25, cell.50)

ind <- d3

for(i in 1:length(per)){
  
  for(j in 1:length(it)){
    
    #Sample i percentage of Q5 cells
    samp <-sample(1:length(d3$Id),round(length(d3$Id)*per[i], digits = 0))
    #Remove their PD values
    val <- d3$PD[-samp]
    #Save the id of those cells removed
    cells <- d3$Id[samp]
    
    #Calculate total PD for the remain Q5 cells
    sum.pd <- sum(val)
    #Save in the respective position
    pd.it3[j,i] <- sum.pd
    
    cell.id[[i]] <- rbind(cell.id[[i]], t(as.data.frame(cells)))
  }
}

pd.it3 <- as.data.frame(pd.it3)
colnames(pd.it3) <- c("p12","p25","p50")
head(pd.it3)

#Here, total PD for each iteration removing X% of Q5 cells
summary(pd.it3)
summary(pd$PD)

pds3 <- c(sum(d3$PD), mean(pd.it3$p12), mean(pd.it3$p25), mean(pd.it3$p50))
col3 <- factor(c("Total PD", "-12.5% (mean)", "-25% (mean)", "-50% (mean)"))

ggplot()+
  geom_bar(aes(y=pds3, x=col3), stat="identity")+
  scale_x_discrete(limits=col3)+
  xlab(NULL)+ylab("PD")

##########################################################################

#Results
pds1
pds2
pds3

pdss <- as.data.frame(rbind(pds1, pds2, pds3))
colnames(pdss) <- c("Total PD", "-12.5% (mean)", "-25% (mean)", "-50% (mean)")

pdss

out.all <- data.frame(pd=c(pdss$`Total PD`,pdss$`-12.5% (mean)`,pdss$`-25% (mean)`,pdss$`-50% (mean)`),
                      Clust= rep(c("d1","d2","d3"), 4),
                      ID=rep(colnames(pdss), each=3))

ggplot()+
  geom_bar(aes(x=out.all$ID, y=out.all$pd, fill=out.all$Clust), stat="identity",position=position_dodge())+
  scale_x_discrete(limits=col3)+
  xlab(NULL)+ylab("PD")

#################################################################

pdloss1 <- (pds1[1]-pds1[-1])/pds1[1]
pdloss2 <- (pds2[1]-pds2[-1])/pds2[1]
pdloss3 <- (pds3[1]-pds3[-1])/pds3[1]


pdss2 <- as.data.frame(rbind(pdloss1, pdloss2, pdloss3))
colnames(pdss2) <- c("-12.5% (mean)", "-25% (mean)", "-50% (mean)")

pdss2

out.all2 <- data.frame(pd=c(pdss2$`-12.5% (mean)`,pdss2$`-25% (mean)`,pdss2$`-50% (mean)`),
                      Clust= rep(c("C1","C2","C3"), 3),
                      ID=rep(colnames(pdss2), each=3))
out.all2

ggplot()+
  geom_bar(aes(x=out.all2$ID, y=out.all2$pd, fill=out.all2$Clust), stat="identity",position=position_dodge())+
  scale_x_discrete(limits=col3[-1])+
  scale_fill_manual(name="Groups", values = c("#1b9e77","#1f78b4","#d95f02"))+
  xlab("% of Remotion")+ylab("% PDloss")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size=15),
        legend.position = c(.15,.65),
        legend.key.size = unit(1,"cm"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=15))

ggsave(filename = "PDloss_Areas.pdf", device = "pdf", width = 20, height = 10, dpi = 300, units = "cm")
