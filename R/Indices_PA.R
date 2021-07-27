########################################################################
#                               Tree List                              #
########################################################################

dir.tree <- as.list(dir()[grep("_tree",dir())])

# Creat a empty list where we will put the trees.
multi.phylo <- lapply(dir.tree, function(x){read.tree(x)})

names(multi.phylo) <- dir.tree

########################################################################
#                              Areas List                              #
########################################################################

setwd("~/Documentos/Omar/Projects/Tesis2/Richness/")

multi.data <- list(read.csv("PNN.inn_25g.matrix"),read.csv("PNN.out_25g.matrix"))
names(multi.data) <- c("PNN.in","PNN.out")

########################################################################
#                           IUCN Categories                            #
########################################################################

iucn <- read.csv("../Richness/IUCN_2.csv")

iucn <- iucn[ , 1:2]

########################################################################
#                   Phylogenetic Diversity Indices                     #
########################################################################

##PD
pd.out <- pd2(Phylogeny = multi.phylo, Distribution = multi.data, pdRoot = T)

names(pd.out) <- names(multi.data)

##AvTD
avtd.out <- avtd2(Phylogeny = multi.phylo, Distribution = multi.data, avtdRoot = T)

names(avtd.out) <- names(multi.data)

##EDGE
edge.r2 <- edge2(Phylo = multi.phylo, IUCN.cat = iucn, na.Cat = c("LC","NT"), freq.Val = c(0.5,0.5), n.iter = 30)

edge.r2V <- edge.r2$Resampled.EDGE
edge.r2V <- data.frame(Species = edge.r2V$Species, 
                       EDGE = as.vector(apply(edge.r2V[,2:31], 1, mean)))

edge.out <- edge.area(edge.val = edge.r2V, matrix = multi.data)

names(edge.out) <- names(multi.data)


save(pd.out, file = "../Indices/pd.out_PNN.rda")
save(avtd.out, file = "../Indices/avtd.out_PNN.rda")
save(edge.out, file = "../Indices/edge.out_PNN.rda")

###########################################################

PNN.in <- data.frame(Ricness = pd.out$PNN.in$SR,
                     PD = pd.out$PNN.in$PD,
                     AvTD = avtd.out$PNN.in$Dplus,
                     EDGE = apply(t(edge.out$PNN.in),1,sum))

write.csv(PNN.in, "../Indices/Indices_PNN.in_Table.csv", quote = F, row.names = F)

PNN.out <- data.frame(Ricness = pd.out$PNN.out$SR,
                      PD = pd.out$PNN.out$PD,
                      AvTD = avtd.out$PNN.out$Dplus,
                      EDGE = apply(t(edge.out$PNN.out),1,sum))

write.csv(PNN.out, "../Indices/Indices_PNN.out_Table.csv", quote = F, row.names = F)

###########################################################

inPNN <- read.csv("../Indices/Indices_PNN.in_Table.csv")
outPNN <- read.csv("../Indices/Indices_PNN.out_Table.csv")
Total <- read.csv("../Indices/Indices_Table.csv")

outTotal <- apply(outPNN, 2, sum)
inTotal <- apply(inPNN, 2, sum)
Total <- apply(Total, 2, sum)

Total
outTotal
inTotal

Total; (outTotal+inTotal)

OutTable <- rbind(round(outTotal/Total, digits = 2),round(inTotal/Total, digits = 2))
rownames(OutTable) <- c("Outside_PNN","Within_PNN")
OutTable <- as.data.frame(OutTable)

OutTable

write.csv(OutTable, "../Indices/PNN_Results.csv", quote = F, row.names = T)
