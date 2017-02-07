library(ape)
library(phytools)
library(phangorn)
library(ggplot2)


setwd("~/Documentos/Omar/Tesis/Taxa/Trees/")

# Read the directory.
dir.tree <- dir()[grep("_tree",dir())]
#dir.tree <- dir.tree[-grep("~",dir.tree)]

# Creat a empty list where we will put the trees.
multi.phylo <- list()

for (i in 1:length(dir.tree)){
  tree <- read.tree(dir.tree[i]) # Read each tree and...
  #plot.phylo(tree)
  multi.phylo[[i]] <- tree # Put inside the list, at the enda we create a multiphylo object.
}

names(multi.phylo) <- dir.tree

#########################################################################################
##                                                                                     ##
##                   First Approach NodeDist = NnodeSP / NnodeAll                      ##
##                                                                                     ##
#########################################################################################

out <- matrix(0, nrow = 1, ncol = 3)
head(out)


for(i in 1:length(multi.phylo)){
  
  count1 <- matrix(NA, nrow = length(multi.phylo[[i]]$tip.label), ncol = 3)
  
  count1[ , 1] <- multi.phylo[[i]]$tip.label
  
  for(j in 1:length(multi.phylo[[i]]$tip.label)){
    
    nNodes <- multi.phylo[[i]]$Nnode
    
    dNodes <- Ancestors(multi.phylo[[i]], j)
    
    proportion <- length(dNodes)/nNodes
    
    count1[j,2] <- proportion
    
    count1[j,3] <- nNodes
    
  }
  
  out <- rbind(out, count1)
  
}

colnames(out) <- c("Species", "Proportion", "Nnodes")

out <- as.data.frame(out[-1,])

head(out)

###########################################################

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/")

Gen.Info <- read.csv("General.info")

head(Gen.Info)

out <- out[which(out$Species%in%Gen.Info$Sp),]

Gen.Info <- cbind(Gen.Info,out)

head(Gen.Info)

Gen.Info

write.csv(Gen.Info, "Gen.Info3.csv",quote = F, row.names = F, col.names = T)

##########################################################
##########################################################

GI <- read.csv("Gen.Info3.csv")

str(GI)

##########################################################
##########################################################

WD <- GI[grep("WD",GI$Ende.WD), ]
WD <- WD[order(WD$BL),]

End <- GI[grep("End",GI$Ende.WD), ]
End <- End[order(End$BL),]

ggplot()+
  geom_histogram(aes(log(GI$BL), fill=GI$Ende.WD,position="dodge"))+
  xlab("log Branch Lengths")+
  scale_fill_discrete(name="Species", labels=c("Endemic", "Widespread"))+
  theme(legend.position = c(.2,.8),
        legend.background = element_blank(),
        legend.key.size = unit(1, "cm"))
  
  
  
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
  



##########################################################
##########################################################


smmr <- summary(GI$Proportion)
smmr

ggplot()+
  geom_histogram(aes(GI$Proportion))+
  geom_vline(aes(xintercept=smmr[2]))+ #1stQu
  geom_vline(aes(xintercept=smmr[5]))+ #3rdQu
  geom_vline(aes(xintercept=smmr[3], colour="Median"))+ #Median
  geom_vline(aes(xintercept=smmr[4], colour="Mean")) #Mean

ggplot()+
  geom_point(aes(x=GI$Nnodes, y=GI$Proportion))+
  geom_hline(aes(yintercept=smmr[2]), linetype = "longdash")+ #1stQu
  geom_hline(aes(yintercept=smmr[5]), linetype = "longdash")+ #3rdQu
  geom_hline(aes(yintercept=smmr[3], colour="Median"))+ #Median
  geom_hline(aes(yintercept=smmr[4], colour="Mean")) #Mean


GI2 <- GI[grep("End", GI$Ende.WD),]

smmr <- summary(GI2$Proportion)
smmr

ggplot()+
  geom_histogram(aes(GI2$Proportion))+
  geom_vline(aes(xintercept=smmr[2]))+ #1stQu
  geom_vline(aes(xintercept=smmr[5]))+ #3rdQu
  geom_vline(aes(xintercept=smmr[3], colour="Median"))+ #Median
  geom_vline(aes(xintercept=smmr[4], colour="Mean")) #Mean

ggplot()+
  geom_point(aes(x=GI2$Nnodes, y=GI2$Proportion))+
  geom_hline(aes(yintercept=smmr[2]), linetype = "longdash")+ #1stQu
  geom_hline(aes(yintercept=smmr[5]), linetype = "longdash")+ #3rdQu
  geom_hline(aes(yintercept=smmr[3], colour="Median"))+ #Median
  geom_hline(aes(yintercept=smmr[4], colour="Mean")) #Mean


#############################################################

GI2 <- GI[grep("WD", GI$Ende.WD),]

smmr <- summary(GI2$Proportion)
smmr

ggplot()+
  geom_histogram(aes(GI2$Proportion))+
  geom_vline(aes(xintercept=smmr[2]))+ #1stQu
  geom_vline(aes(xintercept=smmr[5]))+ #3rdQu
  geom_vline(aes(xintercept=smmr[3], colour="Median"))+ #Median
  geom_vline(aes(xintercept=smmr[4], colour="Mean")) #Mean

ggplot()+
  geom_point(aes(x=GI2$Nnodes, y=GI2$Proportion))+
  geom_hline(aes(yintercept=smmr[2]), linetype = "longdash")+ #1stQu
  geom_hline(aes(yintercept=smmr[5]), linetype = "longdash")+ #3rdQu
  geom_hline(aes(yintercept=smmr[3], colour="Median"))+ #Median
  geom_hline(aes(yintercept=smmr[4], colour="Mean")) #Mean


##########################################################################
##########################################################################
