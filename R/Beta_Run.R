library(ape)
library(phytools)
library(picante)
library(phangorn)
library(jrich)
library(fpc)
library(cluster)
library(maptools)
library(rgdal)
library(ggplot2)
library(gridExtra)

source("...phylo.id.table.R")
source("...pd2.b.R")
source("...toNum.R")
source("...match.phy.dist.R")
source("...cstats.table.R")
source("...cluster.stats.2.R")
source("...phylo.beta.R")
source("...get.cluster.n.R")

########################################################################
#                               Tree List                              #
########################################################################

setwd("../Data/Trees")

dir.tree <- as.list(dir()[grep("_tree",dir())])
pam(out.beta$pd.beta.rich,5)

multi.phylo <- lapply(dir.tree, function(x){read.tree(x)})
  names(multi.phylo) <- dir.tree

########################################################################
#                              Areas List                              #
########################################################################

area <- read.csv("../Data/Distributions/grid_25g.dist.matrix")
  
out.beta <- phylo.beta(multi.phylo = multi.phylo, dist = area)
#save(out.beta, file = "out.beta.rda")
load("out.beta.rda")
str(out.beta)

beta <- out.beta$pd.beta.rich

#########################################################################

stat.clust <- get.cluster.n(dist = beta, k = 9)

plot(stat.clust$N.groups, stat.clust$tot.withinss)
lines(stat.clust$N.groups, stat.clust$tot.withinss)

dend <- hclust(beta,method = "ward.D2")

#Beta total 3 clusters
#Beta rep 3 clusters
#Beta rich 4 clusters
plot(dend)
rect.hclust(dend, 4)

class.cluster <- cutree(dend, k=4)

#########################################################################

grid.25 <- readOGR("../Data/shp/grid_25g.shp")

head(grid.25$ID, 10L)

cells.n <- rep("NO", length(area[1,])-1)
names(cells.n) <- colnames(area)[-1]


cells.n[which(names(cells.n)%in%names(class.cluster))] <- class.cluster


grid.25$Cluster <- cells.n

writeOGR(obj = grid.25,driver = "ESRI Shapefile",dsn = "..Data/shp", layer = "Grid25_Cluster.rich")


###########################################################
#                      Plot Results                       #
###########################################################


clu.t <- readOGR("..Data/shp/Grid25_Cluster.t.shp")
clu.rich <- readOGR("..Data/shp/Grid25_Cluster.rich.shp")
clu.rep <- readOGR("..Data/shp/Grid25_Cluster.rep.shp")
nab <- fortify(readOGR("../Data/shp/NAB.shp"))
andes <- fortify(readOGR("../Data/shp/Andes.shp"))


## Total

cols.t <- rep(clu.t$Cluster, each=5)
cols.rich <- rep(clu.rich$Cluster, each=5)
cols.rep <- rep(clu.rep$Cluster, each=5)


p.t <- fortify(clu.t)
p.rich <- fortify(clu.rich)
p.rep <- fortify(clu.rep)


## Total

bt <- ggplot()+
  geom_polygon(aes(p.t$long, p.t$lat, group=p.t$group, fill=cols.t), colour=NA)+
  geom_polygon(aes(nab$long, nab$lat, group=nab$group), fill=NA, colour="black")+
  geom_polygon(aes(andes$long, andes$lat, group=andes$group), fill=NA, colour="black", size=.25)+
  coord_fixed()+
  scale_y_continuous(limits = c(-3.1,13))+
  scale_fill_manual(values=c("NO"=NA, "1"="#d69b3c","2"="#87e073","3"="#32c8ca"))+
  xlab("Long")+ylab("Lat")+
  theme(legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_rect(fill = "transparent", colour = "black"),
        legend.key.size = unit(0.8, "cm"),
        legend.text = element_text(size = 12))+
  xlab(NULL)+ylab(NULL)+ggtitle("β total")

bt
ggsave(filename = "../Indices/Beta.total_map.png", device = "png", width = 10, height = 15, units = "cm")

## Rep

brep <- ggplot()+
  geom_polygon(aes(p.rep$long, p.rep$lat, group=p.rep$group, fill=cols.rep), colour=NA)+
  geom_polygon(aes(nab$long, nab$lat, group=nab$group), fill=NA, colour="black")+
  geom_polygon(aes(andes$long, andes$lat, group=andes$group), fill=NA, colour="black", size=.25)+
  coord_fixed()+
  scale_y_continuous(limits = c(-3.1,13))+
  scale_fill_manual(values=c("NO"=NA, "1"="#d69b3c","2"="#87e073","3"="#32c8ca"))+
  xlab("Long")+ylab("Lat")+
  theme(legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_rect(fill = "transparent", colour = "black"),
        legend.key.size = unit(0.8, "cm"),
        legend.text = element_text(size = 12))+
  xlab(NULL)+ylab(NULL)+ggtitle("β rep")

brep
ggsave(filename = "../Indices/Beta.rep_map.png", device = "png", width = 10, height = 15, units = "cm")

## Rich

brich <- ggplot()+
  geom_polygon(aes(p.rich$long, p.rich$lat, group=p.rich$group, fill=cols.rich), colour=NA)+
  geom_polygon(aes(nab$long, nab$lat, group=nab$group), fill=NA, colour="black")+
  geom_polygon(aes(andes$long, andes$lat, group=andes$group), fill=NA, colour="black", size=.25)+
  coord_fixed()+
  scale_y_continuous(limits = c(-3.1,13))+
  scale_fill_manual(values=c("NO"=NA, "1"="#d69b3c","2"="#87e073","3"="#32c8ca","4"="#22a884"))+
  xlab("Long")+ylab("Lat")+
  theme(legend.position = "none",
        legend.background = element_blank(),
        legend.key = element_rect(fill = "transparent", colour = "black"),
        legend.key.size = unit(0.8, "cm"),
        legend.text = element_text(size = 12))+
  xlab(NULL)+ylab(NULL)+ggtitle("β rich")

brich
ggsave(filename = "../Indices/Beta.rich_map.png", device = "png", width = 10, height = 15, units = "cm")

#################################################
