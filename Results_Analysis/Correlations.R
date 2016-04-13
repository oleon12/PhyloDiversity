library(ggplot2)

## Set working directory
setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Richness/")

## Read the absence/presence tables for the total richness
area.rich <- read.csv("Area.rich")
grid.rich <- read.csv("grid_25g.rich")
## Create two vectors where the total richness will are
area.t.rich <- c()
grid.t.rich <- c()
## Make the summatory of the richness (This could be made by a apply function, check it later)
for(i in 1:length(colnames(area.rich))){
  total <- sum(area.rich[,i])
  area.t.rich[i] <- total
}

for(i in 1:length(colnames(grid.rich))){
  total <- sum(grid.rich[,i])
  grid.t.rich[i] <- total
}

# Create the two data frames with the index values and richness values
# f.Rich are the the filogenetic richness
# t.Rich are the whole richness

setwd("/home/omar/Documentos/Omar/Tesis/Taxa/Results/Abril1/")

index.area <- data.frame(Area=read.csv("DTenar.total")$area,
                         DT=read.csv("DTenar.total")$W,
                         PD=read.csv("PDenar.total")$PD,
                         AvDT=read.csv("AVDTenar.total")$Dplus,
                         f.Rich=read.csv("DTenar.total")$rich,
                         t.Rich=area.t.rich)

index.grid <- data.frame(Cell.Grid=read.csv("DTg25.total")$area,
                         DT=read.csv("DTg25.total")$W,
                         PD=read.csv("PDg25.total")$PD,
                         AvDT=read.csv("AVDTg25.total")$Dplus,
                         f.Rich=read.csv("DTg25.total")$rich,
                         t.Rich=grid.t.rich)

## Check first f.Rich vs t.Rich

f.t  <- ggplot()
f.t <- f.t + geom_point(aes(x=index.grid$t.Rich,y=index.grid$f.Rich))
f.t

## f.Rich vs Index

fi <- ggplot()
fi <- fi + geom_point(aes(x=index.grid$f.Rich,y = log(index.grid$DT)),colour="red",shape=1)
fi <- fi + geom_point(aes(x=index.grid$f.Rich,y = log(index.grid$PD)),colour="green",shape=1)
fi <- fi + geom_point(aes(x=index.grid$f.Rich,y = index.grid$AvDT),colour="orange",shape=1)
fi

## t.Rich vs Index

fi <- ggplot()
fi <- fi + geom_point(aes(x=index.grid$t.Rich,y = log(index.grid$DT)),colour="red",shape=2)
fi <- fi + geom_point(aes(x=index.grid$t.Rich,y = log(index.grid$PD)),colour="green",shape=2)
fi <- fi + geom_point(aes(x=index.grid$t.Rich,y = index.grid$AvDT),colour="orange",shape=2)
fi

# t.Rich+f.rich vs Index

fi <- ggplot()
fi <- fi + geom_point(aes(x=index.grid$f.Rich,y = log(index.grid$DT)),colour="red",shape=1)
fi <- fi + geom_point(aes(x=index.grid$f.Rich,y = log(index.grid$PD)),colour="green",shape=1)
fi <- fi + geom_point(aes(x=index.grid$f.Rich,y = index.grid$AvDT),colour="orange",shape=1)
fi <- fi + geom_point(aes(x=index.grid$t.Rich,y = log(index.grid$DT)),colour="red",shape=2)
fi <- fi + geom_point(aes(x=index.grid$t.Rich,y = log(index.grid$PD)),colour="green",shape=2)
fi <- fi + geom_point(aes(x=index.grid$t.Rich,y = index.grid$AvDT),colour="orange",shape=2)
fi

## Create outcome matrix

correlations <- matrix(NA, nrow=2,ncol=3)
rownames(correlations) <- c("Total.Rich","Filo.Rich")
colnames(correlations) <- c("DT", "PD", "AVDT")

## Correlations t.rich
##DT vs Total Richness
dt <- cor.test(index.grid$t.Rich,index.grid$DT,method = "spearman")
correlations[1,1] <- dt$estimate
## PD vs Total Richness
pd <- cor.test(index.grid$t.Rich,index.grid$PD,method = "spearman")
correlations[1,2] <- pd$estimate
## AVDT vs Total Richness
avdt <- cor.test(index.grid$t.Rich,index.grid$AvDT,method = "spearman")
correlations[1,3] <- avdt$estimate

## Correlations f.rich
## DT vs Filogenetic Richness
dt <- cor.test(index.grid$f.Rich,index.grid$DT,method = "spearman")
correlations[2,1] <- dt$estimate
## PD vs Filogenetic Richness
pd <- cor.test(index.grid$f.Rich,index.grid$PD,method = "spearman")
correlations[2,2] <- pd$estimate
## AVDT vs Filogenetic Richness
avdt <- cor.test(index.grid$f.Rich,index.grid$AvDT,method = "spearman")
correlations[2,3] <- avdt$estimate

correlations

write.table(correlations,file="Corr.RichvIndex",quote = F, row.names = T,col.names = T)

ind.cor <- cor(index.grid[,2:4],method = "spearman")

write.table(ind.cor,file = "Corr.IndvInd", quote = F, row.names = T,col.names = T)
