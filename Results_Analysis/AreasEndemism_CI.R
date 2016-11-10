library(vegan)

setwd("~/Documentos/Omar/Tesis/Taxa/Results/Final/")

areas <- read.csv("Area.dist.matrix")

name <- areas$especie

areas <- areas[,-1]

rownames(areas) <- name

areas <- t(areas)

areas

d <- as.matrix(vegdist(areas, method = "jaccard"))

magdalenaCI <- d[grep("Magdalena",rownames(areas)),]

magdalenaCI <- as.matrix(magdalenaCI[order(magdalenaCI,decreasing = T)])

colnames(magdalenaCI) <- "CI"

magdalenaCI

write.csv(magdalenaCI,"~/Documentos/Omar/Tesis/Taxa/Results/Final/Raw_IndexR/AreaEndemis.CI",
          quote = F, row.names = T)
