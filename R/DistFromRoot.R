library(ape)
library(phytools)
library(phangorn)
library(ggplot2)

find.root <- function(Phylo){
  
  library(ape, verbose = F,quietly = T)
  library(phangorn,verbose = F, quietly = T)
  
  ##########################################
  # Filter
  if(class(Phylo)!="phylo"){
    stop("Your input must be a class phylo")
  }
  if(is.rooted(Phylo)==F){
    stop("Your phylogeny must be rooted")
  }
  ###########################################
  
  RootNode <- length(Phylo$tip.label)+1
  
  Option <- c(1,length(Phylo$tip.label))
  
  RootTip <- c()
  
  for(i in 1:2){
    
    tip <- Option[i]
    
    Anc <- Ancestors(Phylo,tip)
    
    if(length(Anc)==1){
      
      RootTip <- Phylo$tip.label[Option[i]]
    }
  }
  
  if(is.null(RootTip)){
    
    warning("Your root terminal are not a unique specie, thus, just one will choose")
    
    for(i in 1:2){
      
      tip <- Option[i]
      
      Anc <- Ancestors(Phylo,tip)
      
      if(length(Anc)==2){
        
        RootTip <- Phylo$tip.label[Option[i]]
      }
    }
  }
  return(RootTip)
}

DistNodes1 <- function(BasalN, NiNode, Nnodes){
  
  Dist <- (BasalN * NiNode)/Nnodes
  
  return(Dist)
  
}

DistNodes2 <- function(BasalN, NiNode, Nnodes){
  
  Dist <- (BasalN - NiNode)/Nnodes
  
  return(Dist)
  
}

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

#######################################################################
#######################################################################

out <- matrix(0, nrow = 1, ncol = 5) #Species, Proportion, Nnodes, Dist1, Dist2
head(out)

for(i in 1:length(multi.phylo)){
  
  count1 <- matrix(NA, nrow = length(multi.phylo[[i]]$tip.label), ncol = 5)
  
  count1[ , 1] <- multi.phylo[[i]]$tip.label
  
  for(j in 1:length(multi.phylo[[i]]$tip.label)){
    
    nNodes <- multi.phylo[[i]]$Nnode
    
    dNodes <- length(Ancestors(multi.phylo[[i]], j))
    
    RootSp <- find.root(multi.phylo[[i]])
    
    PosRoot <- grep(RootSp, multi.phylo[[i]]$tip.label)
    
    RootNodes <- length(Ancestors(multi.phylo[[i]], PosRoot))
    
    Dist1 <- DistNodes1(RootNodes,dNodes,nNodes)
    
    Dist2 <- DistNodes2(RootNodes,dNodes,nNodes)
    
    proportion <- dNodes-RootNodes
    
    #proportion <- proportion*nNodes
    
    count1[j,2] <- proportion
    
    count1[j,3] <- nNodes
    
    count1[j,4] <- Dist1
    
    count1[j,5] <- Dist2
    
  }
  
  out <- rbind(out, count1)
  
}

colnames(out) <- c("Species", "Proportion", "Nnodes", "Dist1", "Dist2")

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

write.csv(Gen.Info, "Gen.Info4.csv",quote = F, row.names = F, col.names = T)

##########################################################
##########################################################

GI <- read.csv("Gen.Info4.csv")

str(GI)

##########################################################
##########################################################


smmr <- summary(GI$Proportion)
smmr

p <- ggplot()+
  geom_density(aes(GI$Proportion, fill=GI$Ende.WD), alpha=0.5)+
  xlab("Distance from the root (Proportion)")+ylab(NULL)+
  scale_fill_discrete(name="Species", labels=c("Endemic","non-Endemic"))+
  theme(legend.position = c(.8,.8),
        legend.background = element_blank(),
        legend.title = element_text(size=35),
        legend.text = element_text(size = 30),
        legend.key.size = unit(1, "cm"))+
  theme(axis.text.x=element_text(size=30),
        axis.text.y=element_text(size=15),
        axis.title.y=element_text(size=35),
        axis.title.x=element_text(size=35))

p

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/Endemic/")

png("DistanceFromRoot.png",width = 700, height = 400)
p
dev.off()

p <- ggplot()+
  geom_point(aes(x=GI$Nnodes, y=GI$Proportion, colour=GI$Ende.WD))+
  geom_hline(aes(yintercept=smmr[2]), linetype = "longdash")+ #1stQu
  geom_hline(aes(yintercept=smmr[5]), linetype = "longdash")+ #3rdQu
  geom_hline(aes(yintercept=smmr[3]), colour="gray")+ #Median
  geom_hline(aes(yintercept=smmr[4]), colour="black")+ #Mean
  ylab("Distance from root")+
  xlab("N° of internal nodes in the phylogeny")+
  scale_color_discrete(name="Species", labels=c("Endemic","non-Endemic"))+
  theme(legend.position = c(.2,.8),
        legend.background = element_blank(),
        legend.title = element_text(size=35),
        legend.text = element_text(size = 30),
        legend.key.size = unit(1, "cm"))+
  theme(axis.text.x=element_text(size=30),
        axis.text.y=element_text(size=30),
        axis.title.y=element_text(size=35),
        axis.title.x=element_text(size=35))+
  guides(colour = guide_legend(override.aes = list(size=10)))

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final2/Endemic/")

png("DistanceVSInNodes.png",width = 700, height = 400)
p
dev.off()


##################################################3

GI2 <- GI[grep("End", GI$Ende.WD),]

smmr <- summary(GI2$Proportion)
smmr

ggplot()+
  geom_density(aes(GI2$Proportion), fill="black", alpha=.5)+
  geom_vline(aes(xintercept=smmr[2]))+ #1stQu
  geom_vline(aes(xintercept=smmr[5]))+ #3rdQu
  geom_vline(aes(xintercept=smmr[3], colour="Median"))+ #Median
  geom_vline(aes(xintercept=smmr[4], colour="Mean"))+ #Mean


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
  geom_density(aes(GI2$Proportion), fill="black", alpha=.5)+
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
